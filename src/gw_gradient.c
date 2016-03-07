#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <volume_io.h>
#include <bicpl.h>

#define DBG 0

#define EVAL_NEAREST -1
#define EVAL_LINEAR 0
#define EVAL_CUBIC 2

#define BG 0
#define CSF 1
#define GM 2
#define WM 3
#define L1 5

enum { UNKNOWN, T1, HISTO, T2 };

static void smooth_1d( double *, int, int, double );
static void deriv_1d( double *, double *, int );
static void assign_missing_values( int, int, int [], int *[], Real *, short * );
static void smooth_surface_field( int, int, int [], int *[], int, Real, 
                                  Real *, Real * );

// Determine the magnitude of the distance between the vertices
// of the white surface to the positions of maximum t1-gradient.
// This logic will work for t1 in-vivo MRI and histology image.
//   t1: bg < csf < gm < wm
//   histo: bg/csf < wm < gm (need black background/csf)
// This logic will not work for t2/pd with bright csf (bg < wm < gm < csf ).

public void adjust_gm_wm_gradient( Volume t1, Volume classified, 
                                   Real search_dist, Real search_inc,
                                   int start_point, int end_point,
                                   int n_ngh[], int * ngh[],  
                                   Real * coords, int * mask, 
                                   Real * t1grad,
                                   int * image_type ) {

  int i, j, k, sizes[N_DIMENSIONS];
  Vector normal;
  int n_points = end_point - start_point;

  short * defined = malloc( n_points * sizeof( short ) );
  if( !defined ) {
    printf( "Error allocating memory for defined array.\n" );
    exit( 1 );
  }

  Real * dist = malloc( n_points * sizeof( Real ) );
  if( !dist ) {
    printf( "Error allocating memory for dist array.\n" );
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

  int count_gm = 0, count_wm = 0, count_csf = 0;
  Real avg_gm = 0.0, avg_wm = 0.0, avg_csf = 0.0;
  get_volume_sizes( classified, sizes );
  for( k = 0;  k < sizes[0];  k++ ) {
    for( j = 0;  j < sizes[1];  j++ ) {
      for( i = 0;  i < sizes[2];  i++ ) {
        int val = rint( get_volume_real_value(classified,k,j,i,0,0) );
        int t1val = get_volume_real_value(t1,k,j,i,0,0);
        if( val == CSF ) {
          avg_csf += t1val; count_csf++;
        }
        if( val == GM ) {
          avg_gm += t1val; count_gm++;
        }
        if( val == WM ) {
          avg_wm += t1val; count_wm++;
        }
      }
    }
  }
  if( count_csf > 0 ) avg_csf /= (Real)count_csf;
  if( count_gm > 0 ) avg_gm /= (Real)count_gm;
  if( count_wm > 0 ) avg_wm /= (Real)count_wm;
#if DBG
  printf( "avg_csf = %f\n", avg_csf );
  printf( "avg_gm = %f\n", avg_gm );
  printf( "avg_wm = %f\n", avg_wm );
#endif

  *image_type = UNKNOWN;
  if( ( avg_csf < avg_gm ) && ( avg_csf < avg_wm ) && 
      ( avg_gm < avg_wm ) ) *image_type = T1;

  if( ( avg_csf < avg_gm ) && ( avg_csf < avg_wm ) && 
      ( avg_gm > avg_wm ) ) *image_type = HISTO;

  if( ( avg_csf > avg_gm ) && ( avg_csf > avg_wm ) &&
      ( avg_gm > avg_wm ) ) *image_type = T2;
  if( *image_type == UNKNOWN ) {
    printf( "Cannot determine image type for gradient correction.\n" );
    printf( "gm = %f wm = %f bg = %f csf = %f\n", avg_gm, avg_wm,
            avg_csf );
    exit( 1 );
  }
  Real ventricle_thresh = avg_csf;
  if( *image_type == T1 ) {
    // this is in-vivo T1: csf < gm < wm
    ventricle_thresh = 0.5 * ( avg_csf + avg_gm );
  } else if( *image_type == HISTO ) {
    // this is histology: csf < wm < gm
    ventricle_thresh = 0.5 * ( avg_csf + avg_wm );
  } else if( *image_type == T2 ) {
    // this is T2: wm < gm < csf
    ventricle_thresh = 0.5 * ( avg_csf + avg_gm );
  }

#if DBG
  FILE * fp = fopen( "debug1.txt", "w+t" );
#endif

  for( i = start_point; i < end_point; i++ ) {
    defined[i-start_point] = 0;
    dist[i] = 0.0;
    t1grad[i-start_point] = 0.0;

    // is the vertex on the medial cut plane or under the brainstem?
    if( mask[i-start_point] ) {
      defined[i-start_point] = 1;
      dist[i-start_point] = 0.0;
      evaluate_volume_in_world( t1, coords[3*i+0], coords[3*i+1], coords[3*i+2],
                                EVAL_CUBIC, FALSE, 0.0, &t1grad[i-start_point],
                                NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL );
#if DBG
      fprintf( fp, "# i=%d (masked) \n", i-start_point );
#endif
      continue;
    }

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
      if( cls[j] == L1 ) cls[j] = CSF;
      // ignore possibly labelled WM that could be CSF after filling 
      // the ventricles

      if( *image_type == T1 ) {
        // this is in-vivo T1: csf < gm < wm
        if( cls[j] == WM && vals[j] < ventricle_thresh ) cls[j] = CSF;
      } else if( *image_type == HISTO ) {
        // this is histology: csf < wm < gm
        if( cls[j] == WM && vals[j] < ventricle_thresh ) cls[j] = CSF;
      } else if( *image_type == T2 ) {
        // this is T2: wm < gm < csf
        if( cls[j] == WM && vals[j] > ventricle_thresh ) cls[j] = CSF;
      }
      d += search_inc;
    }
    t1grad[i-start_point] = vals[(ndist-1)/2];

    // compute intensity gradient (derivative) along normal line
    smooth_1d( vals, ndist, 5, 0.25 );
#if 0
    for( j = 0; j < 100; j++ ) {
      smooth_1d( vals, ndist, 1, 0.20 );
      smooth_1d( vals, ndist, 1, -0.204 );
    }

#endif
    deriv_1d( vals, deriv, ndist );

#if DBG
    fprintf( fp, "# i=%d\n", i-start_point );
    fprintf( fp, "#   j    x    v    dv\n" );
    for( j = 0; j < ndist; j++ ) {
      fprintf( fp, "%d %g %d %g %g\n", j, -search_dist + j * search_inc,
               cls[j], vals[j], deriv[j] );
    }
#endif

    // gradient threshold in pure WM inside current white surface.
    int count = 0;
    Real wm_thresh = 0.0;
    Real dwm_thresh = 0.0;
    for( j = 1; j < (int)(ndist/2); j++ ) {
      if( cls[j] == WM ) {
        if( *image_type == T1 ) {
          if( deriv[j] < 0.0 ) {
            count++;
            wm_thresh += vals[j];
            dwm_thresh += deriv[j];
          }
        } else if( *image_type == HISTO || *image_type == T2 ) {
          if( deriv[j] > 0.0 ) {
            count++;
            wm_thresh += vals[j];
            dwm_thresh += deriv[j];
          }
        }
      }
    }
    if( count > 0 ) wm_thresh /= (Real)count;
    if( count > 0 ) dwm_thresh /= (Real)count;

    // GM threshold outside current white surface.
    count = 0;
    Real gm_thresh = 0.0;
    for( j = (int)(ndist/2); j < ndist; j++ ) {
      if( cls[j] == GM ) {
        count++;
        gm_thresh += vals[j];
      }
    }
    if( count > 0 ) gm_thresh /= (Real)count;

    // detect gm-wm maximum gradient 

    // classified gm-csf border in the outside direction
    int j_csf = ndist;

    // remove trailing low-intensity voxels where CSF should be
    for( j_csf = (int)(ndist/2); j_csf < ndist; j_csf++ ) {
      if( *image_type == T1 ) {
        // if( vals[j_csf] < avg_gm ) break;
        if( vals[j_csf] < ventricle_thresh ) break;
      } else if( *image_type == HISTO ) {
        if( vals[j_csf] < ventricle_thresh || cls[j_csf] == CSF ) break;
      } else if( *image_type == T2 ) {
        // if( vals[j_csf] > avg_gm ) break;  // this one is tricky
        if( vals[j_csf] > ventricle_thresh ) break;  // this one is tricky
      }
    }
    // if no csf found, could be in a narrow sulcus, so find 
    // smallest local min
    if( j_csf == ndist ) {
      if( *image_type == T1 || *image_type == HISTO ) {
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
      } else if( *image_type == T2 ) {
        // don't really know what to do. if no csf is visible, 
        // we probably have a local max. do nothing for now.
      }
    }

    int j_max_grad = -1;
    if( *image_type == HISTO ) {
      for( j = 1; j < j_csf; j++ ) {
        if( deriv[j] >= dwm_thresh && vals[j] >= wm_thresh &&
            vals[j] <= gm_thresh ) {
          if( j_max_grad == -1 ) {
            j_max_grad = j;
          } else {
            if( deriv[j] > deriv[j_max_grad] ) {
              j_max_grad = j;
            }
          }
        }
      }
    }

    if( j_max_grad == -1 ) {
      // skip leading non-WM and positive derivative in WM (for t1).
      int j_white = 1;
      for( j_white = 1; j_white < j_csf; j_white++ ) {
        if( cls[j_white] == WM ) {
          if( *image_type == T1 ) {
            if( deriv[j_white] < 0.0 ) break;
          } else if( *image_type == HISTO || *image_type == T2 ) {
            if( deriv[j_white] > 0.0 ) break;
          }
        }
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

// CLAUDE: in below, should avg_gm be gm_thresh?????

      // max negative first gradient closest to j_gmwm (for t1)
      for( j = j_white; j < j_csf-1; j++ ) {
        if( ( *image_type == T1 && deriv[j-1] > deriv[j] && 
              deriv[j+1] > deriv[j] && vals[j] > avg_gm && 
              deriv[j] < dwm_thresh ) ||
            ( ( *image_type == HISTO || *image_type == T2 ) && 
              deriv[j-1] < deriv[j] && deriv[j+1] < deriv[j] && 
              vals[j] > avg_wm && deriv[j] > dwm_thresh ) ) {
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
      fprintf( fp, "# j_csf = %d j_gmwm = %d j_wm = %d j_max=%d j_* = %d wm_th = %g dwm_th = %g gm_th = %g\n\n", 
               j_csf, j_gmwm, j_white, j_max_grad, j_wm, wm_thresh, 
               dwm_thresh, gm_thresh  );
#endif
    } else {
#if DBG
      fprintf( fp, "# j_csf = %d j_max=%d wm_th = %g dwm_th = %g gm_th = %g\n\n", 
               j_csf, j_max_grad, wm_thresh, dwm_thresh, gm_thresh );
#endif
    }

    if( j_max_grad > -1 ) {
      defined[i-start_point] = 1;
      dist[i-start_point] = -(-search_dist + j_max_grad * search_inc);
      // dist[i-start_point] = (j_csf - j_max_grad + 1) * search_inc;
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
  // Note: don't fill in the missing values for t1grad, since missing
  //       values use the image intensity value as a default.

  assign_missing_values( start_point, end_point, n_ngh, ngh, dist, defined );
  free( defined );
 
  // assign_missing_values( start_point, end_point, n_ngh, ngh, t1grad, mask );

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
  free( dist );

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
  // exit(1);
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

