/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
/*
   inflate_to_sphere.c

   Inflate a .obj as a sphere.

   inflate_to_sphere in.obj out.obj

   Values: in.obj = input object file
           out.obj = output object file

   Author: Claude Lepage, May 2011.

*/

#include <math.h>
#include <stdio.h>
#include <volume_io.h>
#include <bicpl.h>
#include "octree.h"

// Prototypes of functions in this file.

static void usage( char * );
static Status read_surface_obj( STRING, int *, Point *[],
                                Vector *[], int *, int *[], int *[], int **[] );
static Status get_surface_neighbours( polygons_struct *, int *[],
                                      int ** [] );
static void save_surface_obj( STRING, int, Point *, Vector *, int, int []);

static void compute_triangle_angles( int, int *, Point * );
static void compute_edge_stats( int, int *, Point * );
static void compute_surface_normals( int, int, int *, Point *, Vector * );
static void compute_triangle_normal( Point, Point, Point, Real[3] );
void smooth( int, int, int, Real, int *, Point *, float );

// Main program.

int main( int argc, char * argv[] ) {

  int      i, j, k, jj, kk, pp, v1, v2, opp1, opp2, t1, t2, found3;
  int      target_nodes, changed, num_swapped, histo[101];
  Real     thresholdlen, minlen, maxlen, factor;

  int      n_points;           // number of grid points per object
  int      n_elems;            // number of triangles per object
  Point  * coords;             // coordinates
  Vector * normals;            // normal vectors
  int    * connec;             // connectivity

  if( argc != 3 ) {
    usage( argv[0] );
    return( 1 );
  }

  // Read the surface file.
  if( read_surface_obj( argv[1], &n_points, &coords, &normals,
                        &n_elems, &connec, NULL, NULL ) != OK ) {
    return 1;
  }

  // Re-center the surface.
  for( j = 0; j < 3; j++ ) {
    Real  sum = 0.0;
    for( i = 1; i < n_points; i++ ) {
      sum += coords[i].coords[j];
    }
    sum /= (Real)n_points;
    for( i = 0; i < n_points; i++ ) {
      coords[i].coords[j] -= sum;
    }
  }

  // Preliminary smooth.
  Real relax = 0.10;
  printf( "preliminary smoothing 1000 iterations...\n" );
  smooth( n_points, n_elems, 1000, relax, connec, coords, 1000.0 );

  // Estimation of average radius.

  Real Rad = 0.0;
  for( i = 0; i < n_points; i++ ) {
    Rad += sqrt( coords[i].coords[0] * coords[i].coords[0] +
                 coords[i].coords[1] * coords[i].coords[1] +
                 coords[i].coords[2] * coords[i].coords[2] );
  }
  Rad /= (Real)n_points;
  printf( "Estimated radius is %f.\n", Rad );

  printf( "inflating to unit sphere...\n" );
  printf( "# iter   mean    stddev     min      max\n" );
  Real alpha0 = 0.05;

  for( kk = 0; kk < 100; kk++ ) {

    // Do some smoothing on the coordinates (simple averaging).
    smooth( n_points, n_elems, 100, relax, connec, coords, 1000.0 );
    printf( "  %3d ", kk );
    compute_edge_stats( n_elems, connec, coords );

    // Normalize the coordinates towards a sphere centered at zero.
    Real alpha = alpha0 + ( 0.5 - alpha0 ) * kk / 25.0;
    if( kk > 25 ) alpha = 1.0;
    for( i = 0; i < n_points; i++ ) {
      Real mag = sqrt( coords[i].coords[0] * coords[i].coords[0] +
                       coords[i].coords[1] * coords[i].coords[1] +
                       coords[i].coords[2] * coords[i].coords[2] );
      mag = 1.0 + alpha * ( Rad / mag - 1.0 );
      coords[i].coords[0] *= mag;
      coords[i].coords[1] *= mag;
      coords[i].coords[2] *= mag;
    }
  }

  // Normalize the coordinates onto the unit sphere centered at zero.
  for( i = 0; i < n_points; i++ ) {
    Real mag = sqrt( coords[i].coords[0] * coords[i].coords[0] +
                     coords[i].coords[1] * coords[i].coords[1] +
                     coords[i].coords[2] * coords[i].coords[2] );
    coords[i].coords[0] /= mag;
    coords[i].coords[1] /= mag;
    coords[i].coords[2] /= mag;
  } 
  compute_triangle_angles( n_elems, connec, coords );

  // Save the inflated sphere.
  compute_surface_normals( n_elems, n_points, connec, coords, normals );
  save_surface_obj( argv[2], n_points, coords, normals, n_elems, connec );

  FREE( connec );
  FREE( coords );
  FREE( normals );

  return 0;
}

// Do smoothing on the coordinates (simple averaging).

void smooth( int n_points, int n_elems, int n_iters, Real relax, 
             int * connec, Point * coords, float minAR ) {


  int i, j, k, kk;
  Real len[3];
  Real * new_coords = (Real *)malloc( 3 * n_points * sizeof( Real ) );
  if( !new_coords ) {
    printf( "Error allocating memory for new_coords.\n" );
    exit( 1 );
  }
  int * countNgh = (int *)malloc( n_points * sizeof( int ) );
  if( !countNgh ) {
    printf( "Error allocating memory for countNgh.\n" );
    exit( 1 );
  }
  Real * weight = (Real *)malloc( n_points * sizeof( Real ) );
  if( !weight ) {
    printf( "Error allocating memory for weight.\n" );
    exit( 1 );
  }
  Real * sumweight = (Real *)malloc( n_points * sizeof( Real ) );
  if( !sumweight ) {
    printf( "Error allocating memory for sumweight.\n" );
    exit( 1 );
  }

  for( i = 0; i < n_points; i++ ) {
    countNgh[i] = 0;
  }
  for( i = 0; i < n_elems; i++ ) {
    countNgh[connec[3*i]]++;
    countNgh[connec[3*i+1]]++;
    countNgh[connec[3*i+2]]++;
  }

  int num_active = n_points;
  if( minAR < 1.0 ) {
    Real invminAR = 1.0 / minAR;
    char * flag = (char *)malloc( n_points * sizeof( char ) );
    if( !flag ) {
      printf( "Error allocating memory for flag.\n" );
      exit( 1 );
    }
    for( i = 0; i < n_points; i++ ) {
      flag[i] = 0;
    }

    for( i = 0; i < n_elems; i++ ) {
      for( j = 0; j < 3; j++ ) {
        Real dx = coords[connec[3*i+j]].coords[0] - coords[connec[3*i+(j+1)%3]].coords[0];
        Real dy = coords[connec[3*i+j]].coords[1] - coords[connec[3*i+(j+1)%3]].coords[1];
        Real dz = coords[connec[3*i+j]].coords[2] - coords[connec[3*i+(j+1)%3]].coords[2];
        len[j] = sqrt( dx * dx + dy * dy + dz * dz );
      }
      Real s = 0.5 * ( len[0] + len[1] + len[2] );
      Real invAR = 0.125 * ( len[0] * len[1] * len[2] ) /
                   ( ( s - len[0] ) * ( s - len[1] ) * ( s - len[2] ) );
      if( invAR > invminAR ) {
        for( j = 0; j < 3; j++ ) {
          flag[connec[3*i+j]] = 1;
        }
      }
    }

    num_active = 0;
    for( i = 0; i < n_points; i++ ) {
      countNgh[i] *= flag[i];
      num_active += flag[i];
    }
    free( flag );
    printf( "Smoothing for AR on %d vertices.\n", num_active );
  }

  for( kk = 0; kk < n_iters && num_active > 0; kk++ ) {
    for( i = 0; i < n_points; i++ ) {
      weight[i] = 0.0;
      sumweight[i] = 0.0;
    }
    for( i = 0; i < 3*n_points; i++ ) {
      new_coords[i] = 0.0;
    }
    for( i = 0; i < n_elems; i++ ) {
      for( j = 0; j < 3; j++ ) {
        int k0 = connec[3*i+j];
        int k1 = connec[3*i+(j+1)%3];
        len[j] = sqrt( ( coords[k0].coords[0] - coords[k1].coords[0] ) *
                       ( coords[k0].coords[0] - coords[k1].coords[0] ) +
                       ( coords[k0].coords[1] - coords[k1].coords[1] ) *
                       ( coords[k0].coords[1] - coords[k1].coords[1] ) +
                       ( coords[k0].coords[2] - coords[k1].coords[2] ) *
                       ( coords[k0].coords[2] - coords[k1].coords[2] ) );
      }
      Real s = 0.5 * ( len[0] + len[1] + len[2] );
      Real area = sqrt( fabs( s * ( s - len[0] ) * ( s - len[1] ) * ( s - len[2] ) ) + 1.0e-10 );
      for( j = 0; j < 3; j++ ) {
        weight[connec[3*i+j]] += area;
      }
    }

    for( i = 0; i < n_elems; i++ ) {
      for( k = 0; k < 3; k++ ) {       // k is 3 verts of triangle
        int k0 = connec[3*i+k];
        int k1 = connec[3*i+(k+1)%3];
        int k2 = connec[3*i+(k+2)%3];
        for( j = 0; j < 3; j++ ) {     // j is x,y,z
          new_coords[3*k0+j] += weight[k1] * coords[k1].coords[j] +
                                weight[k2] * coords[k2].coords[j];
        }
        sumweight[k0] += weight[k1] + weight[k2];
      }
    }

    for( i = 0; i < n_points; i++ ) {
      if( countNgh[i] > 0 ) {
        for( j = 0; j < 3; j++ ) {
          new_coords[3*i+j] = new_coords[3*i+j] / (Real)(sumweight[i]);
          coords[i].coords[j] = relax * coords[i].coords[j] + ( 1.0 - relax ) * new_coords[3*i+j];
        }
      }
    }
  }
  free( countNgh );
  free( new_coords );
  free( weight );
  free( sumweight );
}


// Recompute the surface normals at the nodes. Simply average the normal
// vector of the neighbouring faces at each node.
//
static void compute_surface_normals( int n_elems, int n_points, int * connec, 
                                     Point * coords, Vector * normals ) {

    int  i, j, v0, v1, v2;
    Real norm[3], mag;

    for( i = 0; i < n_points; i++ ) {
      normals[i].coords[0] = 0.0;
      normals[i].coords[1] = 0.0;
      normals[i].coords[2] = 0.0;
    }

    for( i = 0; i < n_elems; i++ ) {
      v0 = connec[3*i];
      v1 = connec[3*i+1];
      v2 = connec[3*i+2];
      compute_triangle_normal( coords[v0], coords[v1], coords[v2], norm );
      for( j = 0; j < 3; j++ ) {
        normals[v0].coords[j] += norm[j];
        normals[v1].coords[j] += norm[j];
        normals[v2].coords[j] += norm[j];
      }
    }

    for( i = 0; i < n_points; i++ ) {
      mag = sqrt( normals[i].coords[0] * normals[i].coords[0] +
                  normals[i].coords[1] * normals[i].coords[1] +
                  normals[i].coords[2] * normals[i].coords[2] );
      normals[i].coords[0] /= mag;
      normals[i].coords[1] /= mag;
      normals[i].coords[2] /= mag;
    }
}

// Compute the normal vector to a triangular face.

static void compute_triangle_normal( Point v1, Point v2, Point v3, Real norm[3] ) {

  Real a1x, a1y, a1z, a2x, a2y, a2z, mag;

  a1x = v2.coords[0] - v1.coords[0];
  a1y = v2.coords[1] - v1.coords[1];
  a1z = v2.coords[2] - v1.coords[2];

  a2x = v3.coords[0] - v1.coords[0];
  a2y = v3.coords[1] - v1.coords[1];
  a2z = v3.coords[2] - v1.coords[2];

  norm[0] = a1y * a2z - a1z * a2y;
  norm[1] = a1z * a2x - a1x * a2z;
  norm[2] = a1x * a2y - a1y * a2x;
  mag = sqrt( norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2] );
  norm[0] /= mag;
  norm[1] /= mag;
  norm[2] /= mag;
}

// Evaluate the quality of the triangles.
//
static void compute_triangle_angles( int n_elems, int * connec, Point * coords ) {

    int  i, j, pp, v1, v2;
    Real a1x, a1y, a1z, a2x, a2y, a2z, mag, cc;
    int  obtuse_flag;

    // Compute current min/max angles of triangles.
    int  num_obtuse = 0;
    Real min_angle = -1.0;
    Real max_angle = 1.0;

    for( i = 0; i < n_elems; i++ ) {
      obtuse_flag = 0;
      for( j = 0; j < 3; j++ ) {
        pp = connec[3*i+j];
        v1 = connec[3*i+(j+1)%3];
        v2 = connec[3*i+(j+2)%3];
        a1x = coords[v1].coords[0] - coords[pp].coords[0];
        a1y = coords[v1].coords[1] - coords[pp].coords[1];
        a1z = coords[v1].coords[2] - coords[pp].coords[2];
        a2x = coords[v2].coords[0] - coords[pp].coords[0];
        a2y = coords[v2].coords[1] - coords[pp].coords[1];
        a2z = coords[v2].coords[2] - coords[pp].coords[2];
        mag = sqrt( ( a1x * a1x + a1y * a1y + a1z * a1z ) * 
                    ( a2x * a2x + a2y * a2y + a2z * a2z ) );
        cc = ( a1x * a2x + a1y * a2y + a1z * a2z ) / mag;
        if( cc > min_angle ) min_angle = cc;
        if( cc < max_angle ) max_angle = cc;
        if( cc <= 0.0 ) obtuse_flag = 1;
      }
      num_obtuse += obtuse_flag;
    }
    min_angle = 57.29577951 * acos( min_angle );
    max_angle = 57.29577951 * acos( max_angle );
    printf( "Min angle = %6.2f deg  Max angle = %6.2f deg\n",
            min_angle, max_angle );
    printf( "%d triangles with angles > 90 degrees\n", num_obtuse );
}

// Evaluate min, max, mean and stddev of edges on the mesh.
//
static void compute_edge_stats( int n_elems, int * connec, Point * coords ) {

  int  i, j;

  Real min_len = 1.0e20;
  Real max_len = 0.0;
  Real sum = 0.0;
  Real sum_sq = 0.0;

  for( i = 0; i < n_elems; i++ ) {
    for( j = 0; j < 3; j++ ) {
      int v1 = connec[3*i+j%3];
      int v2 = connec[3*i+(j+1)%3];
      if( v1 < v2 ) {
        Real dx = coords[v1].coords[0] - coords[v2].coords[0];
        Real dy = coords[v1].coords[1] - coords[v2].coords[1];
        Real dz = coords[v1].coords[2] - coords[v2].coords[2];
        Real len_sq = dx * dx + dy * dy + dz * dz;
        sum_sq += len_sq;
        sum += sqrt( len_sq );
        if( len_sq < min_len ) min_len = len_sq;
        if( len_sq > max_len ) max_len = len_sq;
      }
    }
  }
  int n_edges = ( n_elems * 3 ) / 2;
  Real mean = sum / n_edges;
  Real stddev = sqrt( sum_sq / n_edges - mean * mean );
  printf( " %9.6f %9.6f %9.6f %9.6f\n", (float)mean, 
          (float)stddev, (float)sqrt(min_len), (float)sqrt(max_len) );
}

// -------------------------------------------------------------------
// Help message on how to use this module.
//
static void usage( char * executable_name ) {

  STRING  usage_format = "\
Usage: %s in.obj out.obj\n\
Values: in.obj = input object file\n\
        out.obj = output object file\n\n\
Copyright Alan C. Evans\n\
Professor of Neurology\n\
McGill University\n\n";

  print_error( usage_format, executable_name );
}


// -------------------------------------------------------------------
// Load the cortical surface.
//
// filename: name of the .obj file
// n_points: the number of the vertices
// points: (x,y,z) coordinates
// normals: normal vectors
// n_elem: number of triangles
// connec: connectivity of triangles
// n_neighbours: number of vertices around each node
// neighbours: the set of ordered triangle consisting of the vertices
//
static Status read_surface_obj( STRING filename,
                                 int * n_points,
                                 Point * points[],
                                 Vector * normals[],
                                 int * n_elem,
                                 int * connec[],
                                 int * n_neighbours[],
                                 int ** neighbours[] ) {

  int               i, n_objects;
  object_struct  ** object_list;
  polygons_struct * surface;
  File_formats      format;
  STRING            expanded;

  expanded = expand_filename( filename );   // why?????

  int err = input_graphics_file( expanded, &format, &n_objects,
                                 &object_list );

  if( err != OK ) {
    print_error( "Error reading file %s\n", expanded );
    return( ERROR );
  }

  if( n_objects != 1 || 
      ( n_objects == 1 && get_object_type(object_list[0]) != POLYGONS ) ) {
    print_error( "Error in contents of file %s\n", expanded );
    return( ERROR );
  }

  delete_string( expanded );

  surface = get_polygons_ptr( object_list[0] );

  int ntri = 0, nquad = 0, unknown = 0;
  int start_ind = 0;
  for( i = 0; i < surface->n_items; i++ ) {
    int nn = surface->end_indices[i] - start_ind;
    start_ind = surface->end_indices[i];
    if( nn == 3 ) {
      ntri++;
    } else {
     if( nn == 4 ) {
       nquad++;
     } else {
       unknown++;
       printf( "face with %d nodes\n", nn );
     }
   }
  }
  printf( "%d triangles, %d quads, %d unknown faces in mesh\n", ntri, nquad, unknown );

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Error: Surface must contain only triangular polygons.\n" );
    delete_object_list( n_objects, object_list );
    return ERROR;
  }

  // Make a copy of the coordinates, the normals, and the
  // connectivity since delete_object_list will destroy them.

  *n_points = surface->n_points;
  *n_elem = surface->n_items;
  ALLOC( *points, surface->n_points );
  if( !(*points) ) {
    printf( "Error allocating memory for points.\n" );
    exit( 1 );
  }
  ALLOC( *normals, surface->n_points );
  if( !(*normals) ) {
    printf( "Error allocating memory for normals.\n" );
    exit( 1 );
  }

  for( i = 0; i < *n_points; i++ ) {
    (*points)[i].coords[0] = surface->points[i].coords[0];
    (*points)[i].coords[1] = surface->points[i].coords[1];
    (*points)[i].coords[2] = surface->points[i].coords[2];
    (*normals)[i].coords[0] = surface->normals[i].coords[0];
    (*normals)[i].coords[1] = surface->normals[i].coords[1];
    (*normals)[i].coords[2] = surface->normals[i].coords[2];
  }

  if( connec ) {
    ALLOC( *connec, 3*surface->n_items );
    if( !(*connec) ) {
      printf( "Error allocating memory for connec.\n" );
      exit( 1 );
    }
    for( i = 0; i < 3*surface->n_items; i++ ) {
      (*connec)[i] = surface->indices[i];
    }
  }

  if( n_neighbours && neighbours ) {
    get_surface_neighbours( surface, n_neighbours, neighbours );
  }

  delete_object_list( n_objects, object_list );

  return( OK );
}


// -------------------------------------------------------------------
// Construct the edges around each node. The edges are sorted to
// make an ordered closed loop.
//
private Status get_surface_neighbours( polygons_struct * surface,
                                       int * n_neighbours_return[],
                                       int ** neighbours_return[] ) {

  int    i, j, k, jj;
  int  * tri;
  int  * n_ngh;
  int ** ngh;
  int  * ngh_array;

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Surface must contain only triangular polygons.\n" );
    return ERROR;
  }

  // Check if the node numbering starts at 0 or 1.

  int min_idx, max_idx;

  min_idx = 100*surface->n_points;  // anything big
  max_idx = 0;                      // anything small

  for( i = 0; i < 3*surface->n_items; i++ ) {
    if( surface->indices[i] < min_idx ) min_idx = surface->indices[i];
    if( surface->indices[i] > max_idx ) max_idx = surface->indices[i];
  }

  // Shift numbering to start at zero, for array indexing. Note
  // that we don't care if surface->indices array is modified.

  if( min_idx != 0 ) {
    for( i = 0; i < 3*surface->n_items; i++ ) {
      surface->indices[i] -= min_idx;
    }
  }

  // Count number of triangles attached to each node.

  ALLOC( n_ngh, surface->n_points );
  if( !n_ngh ) {
    printf( "Error allocating memory for n_ngh.\n" );
    exit( 1 );
  }
  ALLOC( ngh, surface->n_points );
  if( !ngh ) {
    printf( "Error allocating memory for ngh.\n" );
    exit( 1 );
  }
  ALLOC( ngh_array, 3*surface->n_items );
  if( !ngh_array ) {
    printf( "Error allocating memory for ngh_array.\n" );
    exit( 1 );
  }

  for( i = 0; i < surface->n_points; i++ ) {
    n_ngh[i] = 0;
  }

  for( i = 0; i < 3*surface->n_items; i++ ) {
    n_ngh[surface->indices[i]]++;
    ngh_array[i] = -1;
  }

  int max_ngh = 0;
  int sum_ngh = 0;
  for( i = 0; i < surface->n_points; i++ ) {
    ngh[i] = &(ngh_array[sum_ngh]);
    sum_ngh += n_ngh[i];
    max_ngh = MAX( max_ngh, n_ngh[i] );
  }

  // At first, store the indices of the triangles in the neighbours.
  for( i = 0; i < surface->n_items; i++ ) {
    for( j = 0; j < 3; j++ ) {
      jj = surface->indices[3*i+j];
      for( k = 0; k < n_ngh[jj]; k++ ) {
        if( ngh[jj][k] == -1 ) {
          ngh[jj][k] = i;
          break;
        }
      }
    }
  }

  // Now create a sort closed loop of the node neighbours.
  // This is needed by the parametric=0 FEM algorithm.
  //
  //         1 ----- 2
  //          /\   /\
  //         /  \ /  \
  //       0 ----P---- 3
  //         \  / \  /
  //          \/   \/
  //         5 ----- 4
  //

  int * tmp;
  ALLOC( tmp, 2*max_ngh );
  if( !tmp ) {
    printf( "Error allocating memory for tmp.\n" );
    exit( 1 );
  }

  for( i = 0; i < surface->n_points; i++ ) {
    for( k = 0; k < n_ngh[i]; k++ ) {
      tri = &(surface->indices[3*ngh[i][k]]);
      for( j = 0; j < 3; j++ ) {
        if( tri[j] == i ) break;
      }
      tmp[2*k+0] = tri[(j+1)%3];
      tmp[2*k+1] = tri[(j+2)%3];
    }

    ngh[i][0] = tmp[0];
    ngh[i][1] = tmp[1];
    for( k = 2; k < n_ngh[i]; k++ ) {
      for( j = 1; j < n_ngh[i]; j++ ) {
        if( tmp[2*j] == ngh[i][k-1] || tmp[2*j+1] == ngh[i][k-1] ) {
          if( tmp[2*j] == ngh[i][k-1] ) {
            ngh[i][k] = tmp[2*j+1];
          } else {
            ngh[i][k] = tmp[2*j];
          }
          tmp[2*j] = -1;
          tmp[2*j+1] = -1;
          break;
        }
      }
    }
  }

  *n_neighbours_return = n_ngh;
  *neighbours_return = ngh;

  FREE( tmp );

  return OK;

}

 
// Save an .obj file.
 
static void save_surface_obj( STRING filename,
                              int n_points,
                              Point coords[],
                              Vector normals[],
                              int n_elems,
                              int connec[] ) {

  int               i, j;

  FILE * fp = fopen( filename, "w" );
  fprintf( fp, "P 0.3 0.3 0.4 10 1 %d\n", n_points );

  // print the coords
  for( j = 0; j < n_points; j++ ) {
    fprintf( fp, "%g %g %g\n", coords[j].coords[0], 
             coords[j].coords[1], coords[j].coords[2] );
  }
  fprintf( fp, "\n" );

  // print the normals
  for( j = 0; j < n_points; j++ ) {
    fprintf( fp, "%g %g %g\n", normals[j].coords[0], 
             normals[j].coords[1], normals[j].coords[2] );
  }

  // The connectivity - part 1.
  fprintf( fp, "\n" );
  fprintf( fp, "%d\n", n_elems );
  fprintf( fp, "0 1 1 1 1\n\n" );

  for( i = 1; i <= n_elems; i++ ) {
    fprintf( fp, "%d ", 3*i );
    if( i%8 == 0 ) fprintf( fp, "\n" );
  }

  // The connectivity - part 2.

  int count = 0;
  for( j = 0; j < 3*n_elems; j++ ) {
    if( count%8 == 0 ) fprintf( fp, "\n" );
    fprintf( fp, "%d ", connec[j] );
    count++;
  }
  fclose( fp );
}
