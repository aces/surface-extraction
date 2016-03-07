/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
/*
   surf_surf_interpolate.c

   Surface-to-surface interpolation to find the projected
   distance map between two surfaces of the same subject or 
   to interpolate a field from one of these surfaces to the 
   other. The output is in the space of the first surface.

   interpolate_sphere surf1.obj surf2.obj [map.txt] out.txt 

   Values: surf1.obj = input object file
           surf2.obj = input object file
           dist.txt = projected distance map (based on surf1.obj)

   or

   Values: surf1.obj = input object file
           surf2.obj = input object file
           field.txt = input field in the space of surf2.obj
           out.txt = interpolated output field in the space of surf1.obj

   Author: Claude Lepage, May 2015
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
static void read_field( STRING, int, Real * field );
static void save_field( STRING, int, Real * field );

// Main program.

int main( int argc, char * argv[] ) {

  int      i, j, k, jj, kk, pp, v1, v2, opp1, opp2, t1, t2, found3;
  int      target_nodes, changed, num_swapped, histo[101];
  Real     thresholdlen, minlen, maxlen, factor;

  int      n_points, n_points_1;  // number of grid points per object
  int      n_elems, n_elems_1;    // number of triangles per object
  Point  * original_coords;       // coordinates
  Point  * coords_1;              // coordinates
  Point  * coords;                // coordinates
  Vector * normals;               // normal vectors
  int    * connec;                // connectivity
  int    * n_ngh = NULL;          // node neighbours (inverse connectivity)
  int   ** ngh = NULL;
  Real   * field_2 = NULL;        // field associated to surf2.obj
  Real   * new_field = NULL;      // interpolated field on surf1.obj

  if( !( argc == 4 || argc == 5 ) ) {
    usage( argv[0] );
    return( 1 );
  }

  // Read the second surface.
  if( read_surface_obj( argv[2], &n_points_1, &coords_1, &normals,
                        &n_elems_1, &connec, &n_ngh, &ngh ) != OK ) {
    return 1;
  }
  FREE( normals );
  FREE( connec );

  // Build an octree for the second surface.
  SurfaceOctree surf_tree;
  Real * newxyz = (Real*)malloc( 3 * n_points_1 * sizeof( Real ) );
  if( !newxyz ) {
    printf( "Error allocating memory for newxyz.\n" );
    exit( 1 );
  }
  for( i = 0; i < n_points_1; i++ ) {
    for( j = 0; j < 3; j++ ) {
      newxyz[3*i+j] = coords_1[i].coords[j];
    }
  }
  initialize_surface_octree( n_points_1, n_ngh, ngh, newxyz, &surf_tree );
  FREE( n_ngh );
  FREE( ngh[0] );
  FREE( ngh );

  // Read the first surface.
  if( read_surface_obj( argv[1], &n_points, &coords, &normals,
                        &n_elems, &connec, NULL, NULL ) != OK ) {
    return 1;
  }
//  if( n_points != n_points_1 || n_elems != n_elems_1 ) {
//    fprintf( stderr, "Mis-matched sizes of object files.\n" );
//    return 1;
//  }

  // Read the surface field for the second surface, if any.
  new_field = (Real *)malloc( n_points * sizeof( Real ) );
  if( argc == 5 ) {
    printf( "Mapping %s onto %s...\n", argv[3], argv[1] );
    field_2 = (Real *)malloc( n_points_1 * sizeof( Real ) );
    read_field( argv[3], n_points_1, field_2 );
  }

  // Surface-to-surface interpolation.

  printf( "surface-to-surface interpolation...\n" );
  for( i = 0; i < n_points; i++ ) {

    int n0, n1, n2;
    Real w0, w1, w2;
    Real dx, dy, dz, mag;
    Real fit = search_point_octree( coords[i].coords[0],
                                    coords[i].coords[1],
                                    coords[i].coords[2],
                                    &n0, &n1, &n2, &w0, &w1, &w2,
                                    &surf_tree );
    if( field_2 ) {
      new_field[i] = w0 * field_2[n0] + w1 * field_2[n1] + w2 * field_2[n2];
    } else {
      dx = w0 * coords_1[n0].coords[0] + w1 * coords_1[n1].coords[0] +
           w2 * coords_1[n2].coords[0] - coords[i].coords[0];
      dy = w0 * coords_1[n0].coords[1] + w1 * coords_1[n1].coords[1] +
           w2 * coords_1[n2].coords[1] - coords[i].coords[1];
      dz = w0 * coords_1[n0].coords[2] + w1 * coords_1[n1].coords[2] +
           w2 * coords_1[n2].coords[2] - coords[i].coords[2];
      mag = sqrt( dx * dx + dy * dy + dz * dz );
      if( dx * normals[i].coords[0] + dy * normals[i].coords[1] +
          dz * normals[i].coords[2] > 0.0 ) {
        mag = -mag;
      }
      new_field[i] = mag;
    }

  }

  free_surface_octree( &surf_tree );
  free( newxyz );
  FREE( coords_1 );

  if( field_2 ) {
    save_field( argv[4], n_points, new_field );
    free( field_2 );
  } else {
    save_field( argv[3], n_points, new_field );
  }

  free( new_field );
  FREE( coords );
  FREE( normals );
  FREE( connec );

  return 0;
}


// -------------------------------------------------------------------
// Help message on how to use this module.
//
static void usage( char * executable_name ) {

  STRING  usage_format = "\
   Surface-to-surface interpolation to find the projected\n\
   distance map between two surfaces of the same subject or \n\
   to interpolate a field from one of these surfaces to the \n\
   other. The output is in the space of the first surface.\n\n\
Usage: %s surf1.obj surf2.obj [map.txt] out.txt\n\n\
Values: surf1.obj = input object file\n\
        surf2.obj = input object file\n\
        dist.txt = projected distance map (based on surf1.obj)\n\n\
or\n\n\
Values: surf1.obj = input object file\n\
        surf2.obj = input object file\n\
        field.txt = input field in the space of surf2.obj\n\
        out.txt = interpolated output field in the space of surf1.obj\n\n\
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
static Status get_surface_neighbours( polygons_struct * surface,
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
  if( fp ) {
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
  } else {
    printf( "Cannot open obj file %s for writing\n", filename );
    exit(1);
  }
}

// Read a .txt file.

static void read_field( STRING filename, 
                        int n_points,
                        Real * field ) {

  int  j;
  double val;

  FILE * fp = fopen( filename, "r" );
  if( fp ) {
    for( j = 0; j < n_points; j++ ) {
      fscanf( fp, "%lg", &val );
      field[j] = (Real)val;
    }
    fclose( fp );
  } else {
    printf( "Cannot open txt file %s for reading\n", filename );
    exit(1);
  }
}

// Save a .txt file.

static void save_field( STRING filename, 
                        int n_points,
                        Real * field ) {
  int  j;

  FILE * fp = fopen( filename, "w" );
  if( fp ) {
    for( j = 0; j < n_points; j++ ) {
      fprintf( fp, "%g\n", field[j] );
    }
    fclose( fp );
  } else {
    printf( "Cannot open txt file %s for writing\n", filename );
    exit(1);
  }

}
