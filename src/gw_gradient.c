#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <volume_io.h>
#include <bicpl.h>


// TODO:
// cleanup: remove dist array (or use it in a second pass to remove
//          3mm distance for estimation of cortical thickness)
//          remove debug stuff
// 

#define DBG 0

#define EVAL_NEAREST -1
#define EVAL_LINEAR 0
#define EVAL_CUBIC 2

#define CSF 1
#define GM 2
#define WM 3

static void smooth_1d( double *, int, int, double );
static void deriv_1d( double *, double *, int );
static void assign_missing_values( int, int, int [], int *[], Real *, short * );
static void smooth_surface_field( int, int, int [], int *[], int, Real, 
                                  Real *, Real * );

// Determine the magnitude of the distance between the vertices
// of the white surface to the positions of maximum t1-gradient.

public void adjust_gm_wm_gradient( Volume t1, Volume classified, 
                                   Real search_dist, Real search_inc,
                                   int start_point, int end_point,
                                   int n_ngh[], int * ngh[],  
                                   Real * coords, Real * t1grad ) {

  int i, j, k, sizes[N_DIMENSIONS];
  Vector normal;
  int n_points = end_point - start_point;

  short * mask = malloc( n_points * sizeof( short ) );
  if( !mask ) {
    printf( "Error allocating memory for mask array.\n" );
    exit( 1 );
  }

  Real * dist = malloc( n_points * sizeof( Real ) );
  if( !dist ) {
    printf( "Error allocating memory for dist.\n" );
    exit( 1 );
  }

  int max_neighbours = 0;
  for( i = start_point; i < end_point; i++ ) {
    if( n_ngh[i] > max_neighbours ) max_neighbours = n_ngh[i];
  }
  Point * neigh_points = (Point *)malloc( max_neighbours * sizeof( Point ) );

  int ndist = (int)( 2.0 * search_dist / search_inc ) + 1;
  Real * vals = (Real*)malloc( ndist * sizeof( Real ) );
  Real * deriv = (Real*)malloc( ndist * sizeof( Real ) );
  int * cls = (int*)malloc( ndist * sizeof( int ) );

  int count_gm = 0;
  Real avg_gm = 0.0;
  get_volume_sizes( classified, sizes );
  for( k = 0;  k < sizes[0];  k++ ) {
    for( j = 0;  j < sizes[1];  j++ ) {
      for( i = 0;  i < sizes[2];  i++ ) {
        int val = rint( get_volume_real_value(classified,k,j,i,0,0) );
        if( val == GM ) {
          int t1val = get_volume_real_value(t1,k,j,i,0,0);
          avg_gm += t1val;
          count_gm++;
        }
      }
    }
  }
  avg_gm /= (Real)count_gm;

#if DBG
  FILE * fp = fopen( "debug1.txt", "w+t" );
#endif

  for( i = start_point; i < end_point; i++ ) {
    mask[i-start_point] = 0;
    dist[i-start_point] = 0.0;
    t1grad[i-start_point] = 0.0;

    // get normal vector pointing outwards
    for( j = 0; j < n_ngh[i]; j++ ) {
      int neigh = ngh[i][j];
      fill_Point( neigh_points[j], coords[3*neigh+0],
                  coords[3*neigh+1], coords[3*neigh+2] );
    }
    find_polygon_normal( n_ngh[i], neigh_points, &normal );

    Real nx = normal.coords[0];
    Real ny = normal.coords[1];
    Real nz = normal.coords[2];

    // scan in +/- normal vector direction for intensity
    Real d = -search_dist;
    for( j = 0; j < ndist; j++ ) {
      Real x = coords[3*i+0] + d * nx;
      Real y = coords[3*i+1] + d * ny;
      Real z = coords[3*i+2] + d * nz;
      evaluate_volume_in_world( t1, x, y, z, EVAL_CUBIC, FALSE,
                                0.0, &vals[j], NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL );
      Real tmp;
      evaluate_volume_in_world( classified, x, y, z, EVAL_NEAREST, FALSE,
                                0.0, &tmp, NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL );
      cls[j] = rint( tmp );
      // ignore possible labelled WM that could be CSF after filling the ventricles
      if( cls[j] == WM && vals[j] < avg_gm ) cls[j] = CSF;
      d += search_inc;
    }
    t1grad[i-start_point] = vals[(ndist-1)/2];

    // compute intensity gradient (derivative) along normal line
    smooth_1d( vals, ndist, 5, 0.25 );
    deriv_1d( vals, deriv, ndist );

#if DBG
    fprintf( fp, "# i=%d\n", i-start_point );
    fprintf( fp, "#   j    x    v    dv\n" );
    for( j = 0; j < ndist; j++ ) {
      fprintf( fp, "%d %g %d %g %g\n", j, -search_dist + j * search_inc,
               cls[j], vals[j], deriv[j] );
    }
#endif

    // gradient threshold in pure WM.
    int count = 0;
    Real dwm_thresh = 0.0;
    for( j = 1; j < (int)(ndist/2); j++ ) {
      if( cls[j] == WM && deriv[j] < 0.0 ) {
        count++;
        dwm_thresh += deriv[j];
      }
    }
    if( count > 0 ) dwm_thresh /= (Real)count;

    // detect gm-wm maximum gradient 

    // classified gm-csf border in the outside direction
    int j_csf = ndist;

    // remove trailing low-intensity voxels where CSF should be
    for( j_csf = (int)(ndist/2); j_csf < ndist; j_csf++ ) {
      if( vals[j_csf] < avg_gm ) break;
    }
    // if no csf found, could be in a narrow sulcus, so find 
    // smallest local min
    if( j_csf == ndist ) {
      int global_min = -1;
      for( j_csf = (int)(ndist/2); j_csf < ndist-1; j_csf++ ) {
        if( vals[j_csf] < vals[j_csf-1] && vals[j_csf] < vals[j_csf+1] ) {
          if( global_min == -1 ) {
            global_min = j_csf;
          } else {
            if( vals[j_csf] < vals[global_min] ) global_min = j_csf;
          }
        }
      }
      if( global_min > 0 ) j_csf = global_min;
    }

    // skip leading non-WM and positive derivative in WM.
    int j_white = 1;
    for( j_white = 1; j_white < j_csf; j_white++ ) {
      if( cls[j_white] == WM && deriv[j_white] < 0.0 ) break;
    }

    // guess gm-wm border 3.0mm inside csf border.
    int j_wm = j_csf - (int)(3.0/search_inc);
    if( j_wm < j_white ) j_wm = j_white+1;
 
    // classified gm-wm border
    int j_gmwm = -1;
    for( j = j_white; j <= j_csf; j++ ) {
      if( cls[j-1] == WM && cls[j] == GM ) j_gmwm = j;
    }
    if( j_gmwm != -1 ) {
      t1grad[i-start_point] = vals[j_gmwm];
    }

    // max negative first gradient closest to j_gmwm
    int j_max_grad = -1;
    for( j = j_white; j < j_csf-1; j++ ) {
      if( deriv[j-1] > deriv[j] && deriv[j+1] > deriv[j] &&
          vals[j] > avg_gm && deriv[j] < dwm_thresh ) {
        if( j_max_grad < -1 ) {
          j_max_grad = j;
        } else {
          if( abs( j - j_wm ) < abs( j_max_grad - j_wm ) ) {
            j_max_grad = j;
          }
        }
      }
    }
#if DBG
    fprintf( fp, "# j_csf = %d j_gmwm = %d j_wm = %d j_max=%d j_* = %d th = %g\n\n", 
             j_csf, j_gmwm, j_white, j_max_grad, j_wm, dwm_thresh );
#endif

    if( j_max_grad > -1 ) {
      mask[i-start_point] = 1;
      // dist[i-start_point] = -(-search_dist + j_max_grad * search_inc);
      dist[i-start_point] = (j_csf - j_max_grad + 1) * search_inc;
      t1grad[i-start_point] = vals[j_max_grad];
    }


  }
  free( vals );
  free( deriv );
  free( neigh_points );
#if DBG
  fclose(fp);
#endif

  // Fill in the values that were not assigned (cannot find gradient position).

  assign_missing_values( start_point, end_point, n_ngh, ngh, dist, mask );
  free( mask );

  // Smooth the displacements to the gradient position.
#if DBG
  FILE * fp0 = fopen( "debug4_pre.txt", "w+t" );
  for( i = 0; i < n_points; i++ ) {
    fprintf( fp0, "%g\n", dist[i] );
  }
  fclose(fp0);
#endif

  int n_iters = 3;
  Real relax = 0.25;
  smooth_surface_field( start_point, end_point, n_ngh, ngh, n_iters, relax, 
                        coords, dist );
  free( dist );

#if DBG
  FILE * fp1 = fopen( "debug4_post.txt", "w+t" );
  for( i = 0; i < n_points; i++ ) {
    fprintf( fp1, "%g\n", dist[i] );
  }
  fclose(fp1);
  FILE * fp2 = fopen( "debug3_pre.txt", "w+t" );
  for( i = 0; i < n_points; i++ ) {
    fprintf( fp2, "%g\n", t1grad[i] );
  }
  fclose(fp2);
#endif

  n_iters = 3;
  smooth_surface_field( start_point, end_point, n_ngh, ngh, n_iters, relax, 
                        coords, t1grad );

// debug stuff

#if DBG
  fp2 = fopen( "debug3_post.txt", "w+t" );
  for( i = 0; i < n_points; i++ ) {
    fprintf( fp2, "%g\n", t1grad[i] );
  }
  fclose(fp2);
  exit(1);
#endif

}

static void smooth_1d( double val[], int n, int iter, double alpha ) {

  int      it, i;

  if( iter < 1 ) return;
  if( n < 5 ) return;

  double * tmp = malloc( n * sizeof( double ) );
  alpha /= 12.0;

  for( it = 1; it <= iter; it++ ) {
    tmp[0] = -( -35.0 * val[0] + 104.0 * val[1] - 114.0 * val[2] + 
               56.0 * val[3] - 11.0 * val[4] );  // sign should be +35
    tmp[1] = ( 11.0 * val[0] - 20.0 * val[1] + 6.0 * val[2] + 
               4.0 * val[3] - val[4] );
    for( i = 2; i < n-2; i++ ) {
      tmp[i] = ( -val[i-2] + 16.0 * val[i-1] - 30.0 * val[i] + 
                 16.0 * val[i+1] - val[i+2] );
    }
    tmp[n-2] = ( 11.0 * val[n-1] - 20.0 * val[n-2] + 6.0 * val[n-3] + 
                 4.0 * val[n-4] - val[n-5] );
    tmp[n-1] = -( -35.0 * val[n-1] + 104.0 * val[n-2] - 114.0 * val[n-3] + 
                 56.0 * val[n-4] - 11.0 * val[n-5] );  // sign should be +35
    for( i = 0; i < n; i++ ) {
      val[i] += alpha * tmp[i];
    }
  }
  free( tmp );
}

// Compute first derivatives of a vector.

static void deriv_1d( double val[], double deriv[], int n ) {

  int      it, i;

  if( n < 5 ) return;

  deriv[0] = ( -25.0 * val[0] + 48.0 * val[1] - 36.0 * val[2] + 
               16.0 * val[3] - 3.0 * val[4] ) / 12.0;
  deriv[1] = ( -3.0 * val[0] - 10.0 * val[1] + 18.0 * val[2] - 
               6.0 * val[3] + val[4] ) / 12.0;
  for( i = 2; i < n-2; i++ ) {
    deriv[i] = ( val[i-2] - 8.0 * val[i-1] +
               8.0 * val[i+1] - val[i+2] ) / 12.0;
  }
  deriv[n-2] = -( -3.0 * val[n-1] - 10.0 * val[n-2] + 18.0 * val[n-3] - 
                  6.0 * val[n-4] + val[n-5] ) / 12.0;
  deriv[n-1] = -( -25.0 * val[n-1] + 48.0 * val[n-2] - 36.0 * val[n-3] + 
               16.0 * val[n-4] - 3.0 * val[n-5] ) / 12.0;
}

// 
// The mask defines vertices with assigned values. Fill in
// the unassigned values by averaging the assigned neighbours.

static void assign_missing_values( int start_point, int end_point, 
                                   int n_ngh[], int * ngh[], 
                                   Real * scalar, short * mask ) {


  int i, j, n, kk, count;

  do {
    count = 0;
    for( i = start_point; i < end_point; i++ ) {
      if( mask[i-start_point] ) continue;

      kk = 0;
      scalar[i-start_point] = 0.0;
      for( n = 0; n < n_ngh[i]; n++ ) {
        if( mask[ngh[i][n]-start_point] ) {
          kk++;
          scalar[i-start_point] += scalar[ngh[i][n]-start_point];
        }
      }
      if( kk >= 3 ) {
        mask[i-start_point] = 1;
        scalar[i-start_point] /= (Real)kk;
      } else {
        count++;
      }
    }
  } while( count > 0 );
}


// Do smoothing of a scalar field on a surface (simple averaging).
// The mask defines vertices with initialized values.

static void smooth_surface_field( int start_point, int end_point,
                                  int n_ngh[], int * ngh[], int n_iters, 
                                  Real relax, Real * coords,
                                  Real * scalar ) {

  int i, j, k, kk;
  int n_points = end_point - start_point;

  Real * new_scalar = (Real *)malloc( n_points * sizeof( Real ) );
  if( !new_scalar ) {
    printf( "Error allocating memory for new_scalar.\n" );
    exit( 1 );
  }

  Real * weight = (Real *)malloc( n_points * sizeof( Real ) );
  if( !weight ) {
    printf( "Error allocating memory for weight.\n" );
    exit( 1 );
  }

  for( i = 0; i < n_points; i++ ) {
    weight[i] = 0.0;
  }

  for( i = start_point; i < end_point; i++ ) {
    for( j = 0; j < n_ngh[i]; j++ ) {
      int k1 = ngh[i][j];
      int k2 = ngh[i][(j+1)%n_ngh[i]];
      if( i < k1 && i < k2 ) {
        Real e01 = sqrt( ( coords[3*i+0] - coords[3*k1+0] ) *
                         ( coords[3*i+0] - coords[3*k1+0] ) +
                         ( coords[3*i+1] - coords[3*k1+1] ) *
                         ( coords[3*i+1] - coords[3*k1+1] ) +
                         ( coords[3*i+2] - coords[3*k1+2] ) *
                         ( coords[3*i+2] - coords[3*k1+2] ) );
        Real e12 = sqrt( ( coords[3*k1+0] - coords[3*k2+0] ) *
                         ( coords[3*k1+0] - coords[3*k2+0] ) +
                         ( coords[3*k1+1] - coords[3*k2+1] ) *
                         ( coords[3*k1+1] - coords[3*k2+1] ) +
                         ( coords[3*k1+2] - coords[3*k2+2] ) *
                         ( coords[3*k1+2] - coords[3*k2+2] ) );
        Real e20 = sqrt( ( coords[3*k2+0] - coords[3*i+0] ) *
                         ( coords[3*k2+0] - coords[3*i+0] ) +
                         ( coords[3*k2+1] - coords[3*i+1] ) *
                         ( coords[3*k2+1] - coords[3*i+1] ) +
                         ( coords[3*k2+2] - coords[3*i+2] ) *
                         ( coords[3*k2+2] - coords[3*i+2] ) );
        Real s = 0.5 * ( e01 + e12 + e20 );
        Real area = sqrt( fabs( s * ( s - e01 ) * ( s - e12 ) *
                              ( s - e20 ) ) + 1.0e-10 );
        weight[i-start_point] += area;
        weight[k1-start_point] += area;
        weight[k2-start_point] += area;
      }
    }
  }

  for( kk = 0; kk < n_iters; kk++ ) {
    for( i = start_point; i < end_point; i++ ) {
      new_scalar[i-start_point] = 0.0;
      if( n_ngh[i] > 0 ) {
        Real sumw = 0.0;
        for( j = 0; j < n_ngh[i]; j++ ) {
          int jj = ngh[i][j]-start_point;
          sumw += weight[jj];
          new_scalar[i-start_point] += weight[jj] * scalar[jj];
        }
        new_scalar[i-start_point] /= sumw;
      }
    }
    for( i = 0; i < n_points; i++ ) {
      scalar[i] = relax * scalar[i] + ( 1.0 - relax ) * new_scalar[i];
    }
  }
  free( new_scalar );
  free( weight );
}

