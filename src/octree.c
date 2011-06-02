/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
/*
   Functions for surface-to-surface interpolation using an octree.

   Author: Claude Lepage, May 2011
*/

#include "octree.h"

#define MAX_ELEM_BRANCH  200

// Prototypes for private functions

private void initialize_branch( SurfaceOctree *, OctreeBranch *,
                        int *, int );
private void free_octree_branch( OctreeBranch * );
private Real search_edge_branch( int, int, OctreeBranch *, Real,
                                 SurfaceOctree * );
private int edge_triangle_intersect( Real [3], Real [3],
                                     Real [3], Real [3], Real [3],
                                     Real [4], Real [6], Real * );
private Real edge_edge_distance( Real [3], Real [3], Real [3], Real [3],
                                 Real );

private Real point_triangle_distance( Real [3], Real [3], Real [3],
                                      Real [3], Real [3], Real [3] );

private Real single_point_triangle_distance( Real, Real, Real,
                                             Real *, Real *, Real *,
                                             Real [3], Real [3], Real [3] );
private Real search_point_branch( Real, Real, Real, int *, int *, int *,
                                  Real *, Real *, Real *,
                                  OctreeBranch *, Real, SurfaceOctree * );


// Static variables

static int num_evals[5];
static int max_elem_branch;
static int min_elem_branch;

// -------------------------------------------------------------------
// Initialize the octree.
//
public void initialize_surface_octree( int n_points,
                                       int n_neighbours[],
                                       int * neighbours[],
                                       Real parameters[],
                                       SurfaceOctree * tree ) {

  int    i, j, k;

  num_evals[0] = 0;
  num_evals[1] = 0;
  num_evals[2] = 0;
  num_evals[3] = 0;
  num_evals[4] = 0;

  max_elem_branch = MAX_ELEM_BRANCH;
  min_elem_branch = MAX_ELEM_BRANCH / 2;
  char * s = getenv( "MAX_ELEM_BRANCH" );
  if( s ) {
    if( strlen( s ) > 0 ) {
      if( sscanf( s, "%d", &max_elem_branch ) != 1 ) {
        max_elem_branch = MAX_ELEM_BRANCH;
      }
      if( max_elem_branch < 0 ) max_elem_branch = MAX_ELEM_BRANCH;
      max_elem_branch = MAX( 10, max_elem_branch );  // no less than 10
      min_elem_branch = MAX( 10, max_elem_branch / 2 );
    }
  }

  printf( "max_elem_branch = %d\n", max_elem_branch );
  printf( "min_elem_branch = %d\n", min_elem_branch );

  tree->n_points = n_points;
  tree->xyz = parameters;
 
  // Create the list of polygons.
  int n_polygons = 0;
  for( i = 0; i < n_points; i++ ) {
    n_polygons += n_neighbours[i];
  }
  tree->n_polygons = n_polygons / 3;

  ALLOC( tree->polys, 3*tree->n_polygons );
  n_polygons = 0;
  int n_edges = 0;
  for( i = 0; i < n_points; i++ ) {
    for( k = 0; k < n_neighbours[i]; k++ ) {
      int p1 = neighbours[i][k];
      int p2 = neighbours[i][(k+1)%n_neighbours[i]];
      if( i < p1 && i < p2 ) {
        tree->polys[3*n_polygons] = i;
        tree->polys[3*n_polygons+1] = p1;
        tree->polys[3*n_polygons+2] = p2;
        n_polygons++;
      }
      if( i < p1 ) n_edges++;
    }
  }

#if 0
  // Check that all triangles are oriented in the same way.
  for( i = 0; i < n_polygons; i++ ) {
    int p0 = tree->polys[3*i];
    int p1 = tree->polys[3*i+1];
    int p2 = tree->polys[3*i+2];
    for( k = 0; k < n_neighbours[p1]; k++ ) {
      if( neighbours[p1][k] == p2 ) {
        if( neighbours[p1][(k+1)%n_neighbours[p1]] == p0 ) {
          // orientation is ok
        } else {
          // orientation is flipped
          printf( "flipped orientation at node %d :" );
          for( j = 0; j < n_neighbours[p1]; j++ ) {
            printf( " %d", neighbours[p1][j] );
          }
          printf( "for tri %d %d %d\n", p0, p1, p2 );
          break;
        }
      }
    }
  }
#endif

  // Compute the equation of the plane for each face.
  //     A*x + B*y + C*z = D
  // (A,B,C) is the normal vector, D = const.

  Real  v1[3], v2[3], vcross[3], mag;
  ALLOC( tree->coeffs, 4*n_polygons );

  for( i = 0; i < n_polygons; i++ ) {

    for( j = 0; j < 3; j++ ) {
      v1[j] = tree->xyz[3*tree->polys[3*i+1]+j] - tree->xyz[3*tree->polys[3*i]+j];
      v2[j] = tree->xyz[3*tree->polys[3*i+2]+j] - tree->xyz[3*tree->polys[3*i]+j];
    }
    vcross[0] = v1[1]*v2[2] - v1[2]*v2[1];
    vcross[1] = v1[2]*v2[0] - v1[0]*v2[2];
    vcross[2] = v1[0]*v2[1] - v1[1]*v2[0];

    tree->coeffs[4*i] = vcross[0];
    tree->coeffs[4*i+1] = vcross[1];
    tree->coeffs[4*i+2] = vcross[2];
    tree->coeffs[4*i+3] = vcross[0]*tree->xyz[3*tree->polys[3*i]] +
                          vcross[1]*tree->xyz[3*tree->polys[3*i]+1] +
                          vcross[2]*tree->xyz[3*tree->polys[3*i]+2];
  }

  // Compute the bounding box for each triangle.

  ALLOC( tree->polybbox, 6*n_polygons );
  for( i = 0; i < n_polygons; i++ ) {

    int j0 = 3*tree->polys[3*i];
    int j1 = 3*tree->polys[3*i+1];
    int j2 = 3*tree->polys[3*i+2];

    Real bmin , bmax;

    for( j = 0; j < 3; j++ ) {
      if( tree->xyz[j0+j] < tree->xyz[j1+j] ) {
        bmin = tree->xyz[j0+j];
        bmax = tree->xyz[j1+j];
      } else {
        bmin = tree->xyz[j1+j];
        bmax = tree->xyz[j0+j];
      }
      if( tree->xyz[j2+j] < bmin ) {
        bmin = tree->xyz[j2+j];
      } else if( tree->xyz[j2+j] > bmax ) {
        bmax = tree->xyz[j2+j];
      }
      tree->polybbox[6*i+j] = bmin;
      tree->polybbox[6*i+3+j] = bmax;
    }
  }

  // Make a lookup table for the edges - to know in which terminal
  // octree leaf an edge is contained.
  ALLOC( tree->edge_index, n_points+1 );
  ALLOC( tree->edge_table, n_edges );
  ALLOC( tree->edge_ptr, n_edges );
  n_edges = 0;
  tree->edge_index[0] = 0;
  for( i = 0; i < n_points; i++ ) {
    for( k = 0; k < n_neighbours[i]; k++ ) {
      if( i < neighbours[i][k] ) {
        tree->edge_table[n_edges] = neighbours[i][k];
        tree->edge_ptr[n_edges] = NULL;
        n_edges++;
      }
    }
    tree->edge_index[i+1] = n_edges;
  }

  // Start creating the octree. List of available triangles (all of them).

  int * list_polys;
  ALLOC( list_polys, tree->n_polygons );
  for( i = 0; i < tree->n_polygons; i++ ) {
    list_polys[i] = i;
  }

  tree->n_id = 0;
  ALLOC( tree->root, 1 );
  initialize_branch( tree, tree->root, list_polys, tree->n_polygons );

  FREE( list_polys );
}

// -------------------------------------------------------------------
// Initialize a branch of the octree (recursive).
// A triangle can only be contained in on branch (no overlap). 
// This means any edge can be in at most 2 boxes.
//
private void initialize_branch( SurfaceOctree * tree,
                                OctreeBranch * branch,
                                int * list_polys,
                                int N ) {

  int    i, k, j, p0, p1;

  branch->id = tree->n_id;
  tree->n_id++;

  // Find min/max bounding box of current domain.

  // First box of first triangle.

  int pp = 6*list_polys[0];
  for( j = 0; j < 3; j++ ) {
    branch->minbox[j] = tree->polybbox[pp+j];
    branch->maxbox[j] = tree->polybbox[pp+3+j];
  }

  for( i = 1; i < N; i++ ) {
    int pp = 6*list_polys[i];
    for( j = 0; j < 3; j++ ) {
      if( tree->polybbox[pp+j] < branch->minbox[j] ) {
        branch->minbox[j] = tree->polybbox[pp+j];
      }
      if( tree->polybbox[pp+3+j] > branch->maxbox[j] ) {
        branch->maxbox[j] = tree->polybbox[pp+3+j];
      }
    }
  }

  int last_right;
  int terminal_branch = 1;

  if( N > max_elem_branch ) {

    Real mid_plane[3];
    mid_plane[0] = 0.0;
    mid_plane[1] = 0.0;
    mid_plane[2] = 0.0;
    for( i = 0; i < N; i++ ) {
      for( j = 0; j < 3; j++ ) {
        p0 = 3*tree->polys[3*list_polys[i]+j];  // node j of triangle i
        for( k = 0; k < 3; k++ ) {
          mid_plane[k] += tree->xyz[p0+k];
        }
      }
    }
    mid_plane[0] /= ( 3.0 * (Real)N );
    mid_plane[1] /= ( 3.0 * (Real)N );
    mid_plane[2] /= ( 3.0 * (Real)N );

    // Find the dimension with the best split in terms of 
    // number of full elements on each side.

    int count[3];
    count[0] = 0;
    count[1] = 0;
    count[2] = 0;

    for( i = 0; i < N; i++ ) {
      for( k = 0; k < 3; k++ ) {
        int check_lo = 0;
        for( j = 0; j < 3; j++ ) {
          p0 = tree->polys[3*list_polys[i]+j];  // node j of triangle i
          if( tree->xyz[3*p0+k] < mid_plane[k] ) check_lo++;
        }
        if( check_lo == 1 || check_lo == 2 ) count[k]++;
      }
    }
    int split_dir = 0;
    if( count[1] < count[split_dir] ) split_dir = 1;
    if( count[2] < count[split_dir] ) split_dir = 2;
    mid_plane[split_dir] *= 3.0;

    last_right = N;

    for( i = 0; i < last_right; i++ ) {
      Real split_avg = 0.0;
      for( j = 0; j < 3; j++ ) {
        p0 = tree->polys[3*list_polys[i]+j];  // node j of triangle i
        split_avg += tree->xyz[3*p0+split_dir];
      }
      if( split_avg > mid_plane[split_dir] ) {
        last_right--;
        int temp = list_polys[i];
        list_polys[i] = list_polys[last_right];
        list_polys[last_right] = temp;
        i--;   // repeat this i
      }
    }

    // Make sure the terminal box is not too small.

    if( last_right < min_elem_branch || N - last_right < min_elem_branch ) {
      terminal_branch = 1;
    } else {
      terminal_branch = 0;
    }

#if 0
// double check if list is ok

    for( i = 0; i < N; i++ ) {
      Real split_avg = 0.0;
      for( j = 0; j < 3; j++ ) {
        p0 = tree->polys[3*list_polys[i]+j];  // node j of triangle i
        split_avg += tree->xyz[3*p0+split_dir];
      }
      if( i < last_right ) {
        if( split_avg > xyz_avg[split_dir] ) {
          printf( "err: i = %d last_right = %d avg = %f split = %f\n",
                  i, last_right, split_avg/3.0, xyz_avg[split_dir]/3.0 );
        }
      } else {
        if( split_avg <= xyz_avg[split_dir] ) {
          printf( "err: i = %d last_right = %d avg = %f split = %f\n",
                  i, last_right, split_avg/3.0, xyz_avg[split_dir]/3.0 );
        }
      }
    }
#endif

  }

  // Branch size is small enough. No more sub-division.
  if( terminal_branch ) {
    // This branch becomes a terminal leaf.
    branch->n_size = N;
    branch->left = NULL;
    branch->right = NULL;
    ALLOC( branch->list_polys, N );
    for( i = 0; i < N; i++ ) {
      branch->list_polys[i] = list_polys[i];
      for( j = 0; j < 3; j++ ) {
        p0 = tree->polys[3*list_polys[i]+j];        // node j of triangle i
        p1 = tree->polys[3*list_polys[i]+(j+1)%3];  // node j+1 of triangle i
        if( p0 < p1 ) {
          for( k = tree->edge_index[p0]; k < tree->edge_index[p0+1]; k++ ) {
            if( tree->edge_table[k] == p1 ) {
              tree->edge_ptr[k] = branch;
              break;
            }
          }
        } else {
          for( k = tree->edge_index[p1]; k < tree->edge_index[p1+1]; k++ ) {
            if( tree->edge_table[k] == p0 ) {
              tree->edge_ptr[k] = branch;
              break;
            }
          }
        }
      }
    }
  } else {
    branch->n_size = 0;
    branch->list_polys = NULL;

    ALLOC( branch->left, 1 );
    initialize_branch( tree, branch->left, &(list_polys[0]),
                       last_right );

    ALLOC( branch->right, 1 );
    initialize_branch( tree, branch->right, &(list_polys[last_right]),
                       N-last_right );

  }
}

// -------------------------------------------------------------------
// Free the memory allocated by the octree.
//
public void free_surface_octree( SurfaceOctree * tree ) {

  printf( "Number of quick edge-tri intersection evaluations = %d\n",
          num_evals[0] );
  printf( "Number of full edge-tri intersection evaluations = %d\n",
          num_evals[1] );
  printf( "Number of edge-edge intersection evaluations = %d\n",
          num_evals[2] );
  num_evals[0] = 0;
  num_evals[1] = 0;
  num_evals[2] = 0;
  num_evals[3] = 0;
  num_evals[4] = 0;

  tree->n_points = 0;
  tree->n_polygons = 0;
  tree->xyz = NULL;         // coords are simply a pointer, not a copy
  if( tree->polys ) {
    FREE( tree->polys );
    tree->polys = NULL;
  }

  if( tree->root ) {
    free_octree_branch( tree->root );
    FREE( tree->root );
    tree->root = NULL;
  }

  if( tree->edge_index ) {
    FREE( tree->edge_index );
    tree->edge_index = NULL;
  }

  if( tree->edge_table ) {
    FREE( tree->edge_table );
    tree->edge_table = NULL;
  }

  if( tree->edge_ptr ) {
    FREE( tree->edge_ptr );
    tree->edge_ptr = NULL;
  }

  if( tree->coeffs ) {
    FREE( tree->coeffs );
    tree->coeffs = NULL;
  }

  if( tree->polybbox ) {
    FREE( tree->polybbox );
    tree->polybbox = NULL;
  }
}

// -------------------------------------------------------------------
// Free the memory allocated by a branch of the octree (recursive).
//
private void free_octree_branch( OctreeBranch * branch ) {

  if( branch == NULL ) return;

  branch->n_size = 0;
  if( branch->list_polys ) {
    FREE( branch->list_polys );
    branch->list_polys = NULL;
  }

  if( branch->left ) {
    free_octree_branch( branch->left );
    FREE( branch->left );
    branch->left = NULL;
  }

  if( branch->right ) {
    free_octree_branch( branch->right );
    FREE( branch->right );
    branch->right = NULL;
  }

}


// -------------------------------------------------------------------
// Find closest distance of the edge (p0:p1) in the octree.
//
// Return: minimal distance found from this edge to any triangle.
//         A negative value (-1.0) means intersection.
//
public Real search_edge_octree( int p0, int p1,
                                SurfaceOctree * tree ) {

  int   i;
  Real  min_dist;

  // For convenience, set p0 < p1.

  if( p0 > p1 ) {
    int tmp = p0;
    p0 = p1;
    p1 = tmp;
  }

  // As a starting point, find the box that contains both p0 and p1.
  // This should give a pretty good initial minimum distance.

  OctreeBranch * branch_ptr = NULL;
  for( i = tree->edge_index[p0]; i < tree->edge_index[p0+1]; i++ ) {
    if( tree->edge_table[i] == p1 ) {
      branch_ptr = tree->edge_ptr[i];
      break;
    }
  }

  min_dist = 1.0e10;
  tree->branch_cache = NULL;
  if( branch_ptr ) {
    min_dist = search_edge_branch( 3*p0, 3*p1, branch_ptr, min_dist, tree );
    if( min_dist < 0.0 ) {
      return( min_dist );
    }
  } else {
    printf( "no starting branch for edge %d %d\n", p0, p1 );
  }

  // Global search in the rest of the tree for the minimum distance.

  tree->branch_cache = branch_ptr;

  min_dist = search_edge_branch( 3*p0, 3*p1, tree->root, min_dist, tree );

//  printf( "edge %d %d  d = %g\n", p0, p1, min_dist );

  return( min_dist );

}

// 

private Real search_edge_branch( int p0, int p1, 
                                 OctreeBranch * branch,
                                 Real curr_dist,
                                 SurfaceOctree * tree ) {

  int  i;

  // Check if we have already visited this branch for 
  // the initial search.

  if( tree->branch_cache == branch ) return( curr_dist );

  // Make sure that the edge can intersect this box within curr_dist.
  // Does the bounding box of the edge closer than curr_dist from 
  // the bounding box of the branch?

  for( i = 0; i < 3; i++ ) {
    Real  bmin, bmax;
    if( tree->xyz[p0+i] < tree->xyz[p1+i] ) {
      bmin = tree->xyz[p0+i];
      bmax = tree->xyz[p1+i];
    } else {
      bmin = tree->xyz[p1+i];
      bmax = tree->xyz[p0+i];
    }
    if( bmax + curr_dist < branch->minbox[i] ) break;
    if( bmin - curr_dist > branch->maxbox[i] ) break;
  }

  int check_box = ( i == 3 );

  printf( "  edge %d %d  d = %g\n", p0/3, p1/3, curr_dist );

  if( check_box ) {

    Real dist_sq, curr_dist_sq;
    curr_dist_sq = curr_dist * curr_dist;

    if( branch->n_size > 0 ) {
      for( i = 0; i < branch->n_size; i++ ) {
        int pp = branch->list_polys[i];
        int t0 = 3*tree->polys[3*pp];
        int t1 = 3*tree->polys[3*pp+1];
        int t2 = 3*tree->polys[3*pp+2];

        // Make sure that the edge does not share a node with the triangle.
        int touching = 0;

        if( p0 == t0 || p0 == t1 || p0 == t2 ) touching = 1;
        if( p1 == t0 || p1 == t1 || p1 == t2 ) touching = 1;

        printf( "check tri %d %d %d ", t0/3, t1/3, t2/3 );

        if( !touching ) {
          if( edge_triangle_intersect( &(tree->xyz[p0]), &(tree->xyz[p1]),
                                       &(tree->xyz[t0]), &(tree->xyz[t1]),
                                       &(tree->xyz[t2]),
                                       &(tree->coeffs[4*pp]),
                                       &(tree->polybbox[6*pp]),
                                       &curr_dist ) ) {
            return( -1.0 );   // intersection
          }
        }
        printf( "  d = %g\n", curr_dist );

        // If no intersection, find the closest distance from edge to
        // any other edge of the triangle. To process the oriented
        // edges of the triangle only once, we make sure that:
        //    t0 < t1 : always true (no need to check)
        //    t1 < t2 : true half the time
        //    t2 < t0 : always false (no need to check)
        // By the construction of the triangles, t0 < t1 and t0 < t2.

        touching = 0;
        if( p0 == t0 || p0 == t1 || p1 == t0 || p1 == t1 ) touching = 1;

        if( !touching ) {
          dist_sq = edge_edge_distance( &(tree->xyz[p0]), &(tree->xyz[p1]),
                                        &(tree->xyz[t0]), &(tree->xyz[t1]),
                                        curr_dist );
          if( dist_sq < curr_dist_sq ) {
            curr_dist = sqrt( dist_sq );
            curr_dist_sq = dist_sq;
          }
        }

        if( t1 < t2 ) {
          touching = 0;
          if( p0 == t1 || p0 == t2 || p1 == t1 || p1 == t2 ) touching = 1;
          if( !touching ) {
            dist_sq = edge_edge_distance( &(tree->xyz[p0]), &(tree->xyz[p1]),
                                          &(tree->xyz[t1]), &(tree->xyz[t2]),
                                          curr_dist );
            if( dist_sq < curr_dist_sq ) {
              curr_dist = sqrt( dist_sq );
              curr_dist_sq = dist_sq;
            }
          }
        }
      }
    } else {
      if( branch->left ) {
        curr_dist = search_edge_branch( p0, p1, branch->left,
                                        curr_dist, tree );
        if( curr_dist < 0.0 ) return( curr_dist );
      }
      if( branch->right ) {
        curr_dist = search_edge_branch( p0, p1, branch->right,
                                        curr_dist, tree );
        if( curr_dist < 0.0 ) return( curr_dist );
      }
    }
  }

  return( curr_dist );
}


public Real tri_tri_distance( Real p0[3], Real p1[3], Real p2[3],
                              Real t0[3], Real t1[3], Real t2[3] ) {

  Real dist_sq, min_dist_sq;

  // Project each point onto the other triangle.

  min_dist_sq = point_triangle_distance( p0, p1, p2, t0, t1, t2 );
  dist_sq = point_triangle_distance( t0, t1, t2, p0, p1, p2 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;

  // Project each edge onto the edges of the other triangle.

  dist_sq = edge_edge_distance( p0, p1, t0, t1, 100.0 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;
  dist_sq = edge_edge_distance( p0, p1, t0, t2, 100.0 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;
  dist_sq = edge_edge_distance( p0, p1, t1, t2, 100.0 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;
  //printf( "p0:p1 on t0:t1 = %g\n", edge_edge_distance( p0, p1, t0, t1, 100.0 ) );
  //printf( "p0:p1 on t0:t2 = %g\n", edge_edge_distance( p0, p1, t0, t2, 100.0 ) );
  //printf( "p0:p1 on t1:t2 = %g\n", edge_edge_distance( p0, p1, t1, t2, 100.0 ) );

  dist_sq = edge_edge_distance( p0, p2, t0, t1, 100.0 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;
  dist_sq = edge_edge_distance( p0, p2, t0, t2, 100.0 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;
  dist_sq = edge_edge_distance( p0, p2, t1, t2, 100.0 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;
  //printf( "p0:p2 on t0:t1 = %g\n", edge_edge_distance( p0, p2, t0, t1, 100.0 ) );
  //printf( "p0:p2 on t0:t2 = %g\n", edge_edge_distance( p0, p2, t0, t2, 100.0 ) );
  //printf( "p0:p2 on t1:t2 = %g\n", edge_edge_distance( p0, p2, t1, t2, 100.0 ) );

  dist_sq = edge_edge_distance( p1, p2, t0, t1, 100.0 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;
  dist_sq = edge_edge_distance( p1, p2, t0, t2, 100.0 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;
  dist_sq = edge_edge_distance( p1, p2, t1, t2, 100.0 );
  if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;
  //printf( "p1:p2 on t0:t1 = %g\n", edge_edge_distance( p1, p2, t0, t1, 100.0 ) );
  //printf( "p1:p2 on t0:t2 = %g\n", edge_edge_distance( p1, p2, t0, t2, 100.0 ) );
  //printf( "p1:p2 on t1:t2 = %g\n", edge_edge_distance( p1, p2, t1, t2, 100.0 ) );

  return( min_dist_sq );
}

private Real point_triangle_distance( Real p0[3], Real p1[3], Real p2[3],
                                      Real t0[3], Real t1[3], Real t2[3] ) {

  // Points p0, p1, p2 on triangle (t0,t1,t2).

  Real A_tt = ( t1[0] - t0[0] ) * ( t1[0] - t0[0] ) +
              ( t1[1] - t0[1] ) * ( t1[1] - t0[1] ) +
              ( t1[2] - t0[2] ) * ( t1[2] - t0[2] );
  Real A_st = ( t1[0] - t0[0] ) * ( t2[0] - t0[0] ) +
              ( t1[1] - t0[1] ) * ( t2[1] - t0[1] ) +
              ( t1[2] - t0[2] ) * ( t2[2] - t0[2] );
  Real A_ss = ( t2[0] - t0[0] ) * ( t2[0] - t0[0] ) +
              ( t2[1] - t0[1] ) * ( t2[1] - t0[1] ) +
              ( t2[2] - t0[2] ) * ( t2[2] - t0[2] );
  Real det = A_tt * A_ss - A_st * A_st;

  if( det == 0.0 ) {
    // ????? something is wrong???
    // printf( "fatal error in intersection\n" );
    // exit(1);
    // looks like a flat triangle, so project point to
    // a different closest triangle instead.
    return( 1.0e20 );
  } else {
    Real rhs_t, rhs_s, xi, eta, dx, dy, dz, dist_sq, min_dist_sq;
    // Point p0.
    rhs_t = ( t1[0] - t0[0] ) * ( p0[0] - t0[0] ) +
            ( t1[1] - t0[1] ) * ( p0[1] - t0[1] ) +
            ( t1[2] - t0[2] ) * ( p0[2] - t0[2] );
    rhs_s = ( t2[0] - t0[0] ) * ( p0[0] - t0[0] ) +
            ( t2[1] - t0[1] ) * ( p0[1] - t0[1] ) +
            ( t2[2] - t0[2] ) * ( p0[2] - t0[2] );
    xi = ( A_ss * rhs_t - A_st * rhs_s ) / det;
    eta = ( A_tt * rhs_s - A_st * rhs_t ) / det;
    // clip to nearest side of triangle
    if( xi + eta > 1.0 ) {    // clip on line xi+eta=1
      xi = xi / ( xi + eta );
      if( xi < 0.0 ) {
        xi = 0.0;
      } else if( xi > 1.0 ) {
        xi = 1.0;
      }
      eta = 1.0 - xi;
    } else if( xi < 0.0 ) {      // clip on line xi=0
      eta = eta / ( 1.0 - xi );
      xi = 0.0;
      if( eta < 0.0 ) {
        eta = 0.0;
      } else if( eta > 1.0 ) {
        eta = 1.0;
      }
    } else if( eta < 0.0 ) {      // clip on line eta=0
      xi = xi / ( 1.0 - eta );
      eta = 0.0;
      if( xi < 0.0 ) {
        xi = 0.0;
      } else if( xi > 1.0 ) {
        xi = 1.0;
      }
    }
    dx = ( 1.0 - xi - eta ) * t0[0] + xi * t1[0] + eta * t2[0] - p0[0];
    dy = ( 1.0 - xi - eta ) * t0[1] + xi * t1[1] + eta * t2[1] - p0[1];
    dz = ( 1.0 - xi - eta ) * t0[2] + xi * t1[2] + eta * t2[2] - p0[2];
    min_dist_sq = dx * dx + dy * dy + dz * dz;

    // Point p1.
    rhs_t = ( t1[0] - t0[0] ) * ( p1[0] - t0[0] ) +
            ( t1[1] - t0[1] ) * ( p1[1] - t0[1] ) +
            ( t1[2] - t0[2] ) * ( p1[2] - t0[2] );
    rhs_s = ( t2[0] - t0[0] ) * ( p1[0] - t0[0] ) +
            ( t2[1] - t0[1] ) * ( p1[1] - t0[1] ) +
            ( t2[2] - t0[2] ) * ( p1[2] - t0[2] );
    xi = ( A_ss * rhs_t - A_st * rhs_s ) / det;
    eta = ( A_tt * rhs_s - A_st * rhs_t ) / det;
    // clip to nearest side of triangle
    if( xi + eta > 1.0 ) {    // clip on line xi+eta=1
      xi = xi / ( xi + eta );
      if( xi < 0.0 ) {
        xi = 0.0;
      } else if( xi > 1.0 ) {
        xi = 1.0;
      }
      eta = 1.0 - xi;
    } else if( xi < 0.0 ) {      // clip on line xi=0
      eta = eta / ( 1.0 - xi );
      xi = 0.0;
      if( eta < 0.0 ) {
        eta = 0.0;
      } else if( eta > 1.0 ) {
        eta = 1.0;
      }
    } else if( eta < 0.0 ) {      // clip on line eta=0
      xi = xi / ( 1.0 - eta );
      eta = 0.0;
      if( xi < 0.0 ) {
        xi = 0.0;
      } else if( xi > 1.0 ) {
        xi = 1.0;
      }
    }
    dx = ( 1.0 - xi - eta ) * t0[0] + xi * t1[0] + eta * t2[0] - p1[0];
    dy = ( 1.0 - xi - eta ) * t0[1] + xi * t1[1] + eta * t2[1] - p1[1];
    dz = ( 1.0 - xi - eta ) * t0[2] + xi * t1[2] + eta * t2[2] - p1[2];
    dist_sq = dx * dx + dy * dy + dz * dz;
    if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;

    // Point p2.
    rhs_t = ( t1[0] - t0[0] ) * ( p2[0] - t0[0] ) +
            ( t1[1] - t0[1] ) * ( p2[1] - t0[1] ) +
            ( t1[2] - t0[2] ) * ( p2[2] - t0[2] );
    rhs_s = ( t2[0] - t0[0] ) * ( p2[0] - t0[0] ) +
            ( t2[1] - t0[1] ) * ( p2[1] - t0[1] ) +
            ( t2[2] - t0[2] ) * ( p2[2] - t0[2] );
    xi = ( A_ss * rhs_t - A_st * rhs_s ) / det;
    eta = ( A_tt * rhs_s - A_st * rhs_t ) / det;
    // clip to nearest side of triangle
    if( xi + eta > 1.0 ) {    // clip on line xi+eta=1
      xi = xi / ( xi + eta );
      if( xi < 0.0 ) {
        xi = 0.0;
      } else if( xi > 1.0 ) {
        xi = 1.0;
      }
      eta = 1.0 - xi;
    } else if( xi < 0.0 ) {      // clip on line xi=0
      eta = eta / ( 1.0 - xi );
      xi = 0.0;
      if( eta < 0.0 ) {
        eta = 0.0;
      } else if( eta > 1.0 ) {
        eta = 1.0;
      }
    } else if( eta < 0.0 ) {      // clip on line eta=0
      xi = xi / ( 1.0 - eta );
      eta = 0.0;
      if( xi < 0.0 ) {
        xi = 0.0;
      } else if( xi > 1.0 ) {
        xi = 1.0;
      }
    }
    dx = ( 1.0 - xi - eta ) * t0[0] + xi * t1[0] + eta * t2[0] - p2[0];
    dy = ( 1.0 - xi - eta ) * t0[1] + xi * t1[1] + eta * t2[1] - p2[1];
    dz = ( 1.0 - xi - eta ) * t0[2] + xi * t1[2] + eta * t2[2] - p2[2];
    dist_sq = dx * dx + dy * dy + dz * dz;
    if( dist_sq < min_dist_sq ) min_dist_sq = dist_sq;

    return( min_dist_sq );
  }
  return( 0.0 );
}



// Return: 0 if no intersection
//         1 if intersection

private int edge_triangle_intersect( Real p0[3], Real p1[3],
                                     Real t0[3], Real t1[3], Real t2[3],
                                     Real coeffs[4], Real tribbox[6],
                                     Real * curr_dist ) {

  int  j;

  // Check the bounding box first.

  for( j = 0; j < 3; j++ ) {
    if( p0[j] < p1[j] ) {
      if( p1[j] + (*curr_dist) < tribbox[j] ||
          p0[j] - (*curr_dist) > tribbox[3+j] ) return( 0 );
    } else {
      if( p0[j] + (*curr_dist) < tribbox[j] ||
          p1[j] - (*curr_dist) > tribbox[3+j] ) return( 0 );
    }
  }

  // Check first if the edge intersects the plane of the triangle.
  Real denom;
  denom = coeffs[0] * ( p1[0] - p0[0] ) +
          coeffs[1] * ( p1[1] - p0[1] ) +
          coeffs[2] * ( p1[2] - p0[2] );
  if( ABS(denom) < 1.0e-20 ) {
    printf( "denom = %g\n", denom );
    // edge is essentially parallel to triangle - assume no intersection
  } else {
    Real t = ( coeffs[3] - coeffs[0] * p0[0] - coeffs[1] * p0[1] -
               coeffs[2] * p0[2] ) / denom;
    num_evals[0]++;

    Real inter[3];

    int inside = 0;
    if( t >= 0.0 && t <= 1.0 ) {
      // The plane intersects inside the edge. Must check that
      // this intersection point projects inside the triangle.

      for( j = 0; j < 3; j++ ) {
        inter[j] = p0[j] + t * ( p1[j] - p0[j] );
      }
      inside = 1;
    } else {
      // Project the closest end point onto the triangle.
      if( t <= 0.0 ) {
        inter[0] = p0[0];
        inter[1] = p0[1];
        inter[2] = p0[2];
      } else {
        inter[0] = p1[0];
        inter[1] = p1[1];
        inter[2] = p1[2];
      }
    }

    num_evals[1]++;

    Real A_tt = ( t1[0] - t0[0] ) * ( t1[0] - t0[0] ) +
                ( t1[1] - t0[1] ) * ( t1[1] - t0[1] ) +
                ( t1[2] - t0[2] ) * ( t1[2] - t0[2] );
    Real A_st = ( t1[0] - t0[0] ) * ( t2[0] - t0[0] ) +
                ( t1[1] - t0[1] ) * ( t2[1] - t0[1] ) +
                ( t1[2] - t0[2] ) * ( t2[2] - t0[2] );
    Real A_ss = ( t2[0] - t0[0] ) * ( t2[0] - t0[0] ) +
                ( t2[1] - t0[1] ) * ( t2[1] - t0[1] ) +
                ( t2[2] - t0[2] ) * ( t2[2] - t0[2] );
    Real rhs_t = ( t1[0] - t0[0] ) * ( inter[0] - t0[0] ) +
                 ( t1[1] - t0[1] ) * ( inter[1] - t0[1] ) +
                 ( t1[2] - t0[2] ) * ( inter[2] - t0[2] );
    Real rhs_s = ( t2[0] - t0[0] ) * ( inter[0] - t0[0] ) +
                 ( t2[1] - t0[1] ) * ( inter[1] - t0[1] ) +
                 ( t2[2] - t0[2] ) * ( inter[2] - t0[2] );
    Real det = A_tt * A_ss - A_st * A_st;

    if( det == 0.0 ) {
      // ????? something is wrong???
      // printf( "fatal error in intersection\n" );
      // exit(1);
      // looks like a flat triangle, so project point to
      // a different closest triangle instead.
      return( 0 );
    } else {
      Real xi = ( A_ss * rhs_t - A_st * rhs_s ) / det;
      Real eta = ( A_tt * rhs_s - A_st * rhs_t ) / det;
      if( ( xi >= 0.0 ) && ( eta >= 0.0 ) && ( xi + eta <= 1.0 ) ) {
        if( inside ) return( 1 );    // bad news: intersection found
      } else {
        // clip to nearest side of triangle
        if( xi + eta > 1.0 ) {    // clip on line xi+eta=1
          xi = xi / ( xi + eta );
          if( xi < 0.0 ) {
            xi = 0.0;
          } else if( xi > 1.0 ) {
            xi = 1.0;
          }
          eta = 1.0 - xi;
        } else if( xi < 0.0 ) {      // clip on line xi=0
          eta = eta / ( 1.0 - xi );
          xi = 0.0;
          if( eta < 0.0 ) {
            eta = 0.0;
          } else if( eta > 1.0 ) {
            eta = 1.0;
          }
        } else if( eta < 0.0 ) {      // clip on line eta=0
          xi = xi / ( 1.0 - eta );
          eta = 0.0;
          if( xi < 0.0 ) {
            xi = 0.0;
          } else if( xi > 1.0 ) {
            xi = 1.0;
          }
        }
      }
      Real dx = ( 1.0 - xi - eta ) * t0[0] + xi * t1[0] + eta * t2[0] - inter[0];
      Real dy = ( 1.0 - xi - eta ) * t0[1] + xi * t1[1] + eta * t2[1] - inter[1];
      Real dz = ( 1.0 - xi - eta ) * t0[2] + xi * t1[2] + eta * t2[2] - inter[2];
      Real dist_sq = dx * dx + dy * dy + dz * dz;
      if( dist_sq < (*curr_dist)*(*curr_dist) ) {
        *curr_dist = sqrt( dist_sq );
      }
    }
  }
  return( 0 );
}


// -------------------------------------------------------------------
// Compute the minimum distance between two line segments.
// Return: dist**2
//
private Real edge_edge_distance( Real p0[3], Real p1[3],
                                 Real q0[3], Real q1[3],
                                 Real curr_dist ) {

  // Check bounding box.

  int  j;
  Real pmin, pmax, qmin, qmax;
  for( j = 0; j < 3; j++ ) {
    if( p0[j] < p1[j] ) {
      pmin = p0[j];
      pmax = p1[j];
    } else {
      pmin = p1[j];
      pmax = p0[j];
    }
    if( q0[j] < q1[j] ) {
      qmin = q0[j];
      qmax = q1[j];
    } else {
      qmin = q1[j];
      qmax = q0[j];
    }
    if( ( qmin - pmax > curr_dist ) || ( pmin - qmax > curr_dist ) ) break;
  }

  if( j < 3 ) return( curr_dist * curr_dist );

  num_evals[2]++;

  Real A_tt = ( p1[0] - p0[0] ) * ( p1[0] - p0[0] ) +
              ( p1[1] - p0[1] ) * ( p1[1] - p0[1] ) +
              ( p1[2] - p0[2] ) * ( p1[2] - p0[2] );
  Real A_st = - ( ( p1[0] - p0[0] ) * ( q1[0] - q0[0] ) +
                  ( p1[1] - p0[1] ) * ( q1[1] - q0[1] ) +
                  ( p1[2] - p0[2] ) * ( q1[2] - q0[2] ) );
  Real A_ss = ( q1[0] - q0[0] ) * ( q1[0] - q0[0] ) +
              ( q1[1] - q0[1] ) * ( q1[1] - q0[1] ) +
              ( q1[2] - q0[2] ) * ( q1[2] - q0[2] );
  Real rhs_t = -( ( p1[0] - p0[0] ) * ( p0[0] - q0[0] ) +
                  ( p1[1] - p0[1] ) * ( p0[1] - q0[1] ) +
                  ( p1[2] - p0[2] ) * ( p0[2] - q0[2] ) );
  Real rhs_s = ( q1[0] - q0[0] ) * ( p0[0] - q0[0] ) +
               ( q1[1] - q0[1] ) * ( p0[1] - q0[1] ) +
               ( q1[2] - q0[2] ) * ( p0[2] - q0[2] );

  Real det = A_tt * A_ss - A_st * A_st;
  Real dx, dy, dz, dist_sq;

  if( det == 0.0 ) {
    // edges are parallel - there are an infinity of solution (s,t).
    // Take t = 0 and project p0 on line q0:q1.

    Real s = rhs_s / A_ss;
    if( s < 0.0 ) s = 0.0;
    if( s > 1.0 ) s = 1.0;

    dx = q0[0] + s * ( q1[0] - q0[0] ) - p0[0];
    dy = q0[1] + s * ( q1[1] - q0[1] ) - p0[1];
    dz = q0[2] + s * ( q1[2] - q0[2] ) - p0[2];
    dist_sq = dx * dx + dy * dy + dz * dz;
    return( dist_sq );

  } else {

    Real t = ( A_ss * rhs_t - A_st * rhs_s ) / det;
    Real s = ( A_tt * rhs_s - A_st * rhs_t ) / det;

    Real dist_s, dist_t;
    int  clip = 0;
    if( t < 0.0 || t > 1.0 ) {
      Real * pp = p0;
      if( t >= 1.0 ) {
        rhs_s = ( q1[0] - q0[0] ) * ( p1[0] - q0[0] ) +
                ( q1[1] - q0[1] ) * ( p1[1] - q0[1] ) +
                ( q1[2] - q0[2] ) * ( p1[2] - q0[2] );
        pp = p1;
      }
      Real ss = rhs_s / A_ss;
      if( ss < 0.0 ) ss = 0.0;
      if( ss > 1.0 ) ss = 1.0;

      dx = q0[0] + ss * ( q1[0] - q0[0] ) - pp[0];
      dy = q0[1] + ss * ( q1[1] - q0[1] ) - pp[1];
      dz = q0[2] + ss * ( q1[2] - q0[2] ) - pp[2];
      dist_s = dx * dx + dy * dy + dz * dz;
      clip = 1;
    }
    if( s < 0.0 || s > 1.0 ) {
      Real * qq = q0;
      if( s >= 1.0 ) {
        rhs_t = -( ( p1[0] - p0[0] ) * ( p0[0] - q1[0] ) +
                   ( p1[1] - p0[1] ) * ( p0[1] - q1[1] ) +
                   ( p1[2] - p0[2] ) * ( p0[2] - q1[2] ) );
        qq = q1;
      }
      Real tt = rhs_t / A_tt;
      if( tt < 0.0 ) tt = 0.0;
      if( tt > 1.0 ) tt = 1.0;

      dx = p0[0] + tt * ( p1[0] - p0[0] ) - qq[0];
      dy = p0[1] + tt * ( p1[1] - p0[1] ) - qq[1];
      dz = p0[2] + tt * ( p1[2] - p0[2] ) - qq[2];
      dist_t = dx * dx + dy * dy + dz * dz;
      clip += 2;
    }

    if( !clip ) {
      dx = q0[0] + s * ( q1[0] - q0[0] ) - ( p0[0] + t * ( p1[0] - p0[0] ) );
      dy = q0[1] + s * ( q1[1] - q0[1] ) - ( p0[1] + t * ( p1[1] - p0[1] ) );
      dz = q0[2] + s * ( q1[2] - q0[2] ) - ( p0[2] + t * ( p1[2] - p0[2] ) );
      dist_sq = dx * dx + dy * dy + dz * dz;
      return( dist_sq );
    } else {
      if( clip == 1 ) {
        return( dist_s );
      } else if( clip == 2 ) {
        return( dist_t );
      } else {
        return( MIN( dist_s, dist_t ) );
      }
    }

  }

}

// -------------------------------------------------------------------
// Find closest distance of the point (x,y,z) in the octree.
// Return the vertices of this closest triangle and its local
// coordinates at (x,y,z).
//
// Return: minimal distance found from this point to any triangle.
//
public Real search_point_octree( Real x, Real y, Real z,
                                 int * t0, int * t1, int * t2,
                                 Real * w0, Real * w1, Real * w2,
                                 SurfaceOctree * tree ) {

  Real min_dist = 1.0e10;

  // Global search in the tree for the minimum distance, starting at the root.

  min_dist = search_point_branch( x, y, z, t0, t1, t2, w0, w1, w2,
                                  tree->root, min_dist, tree );

  return( min_dist );

}

// -------------------------------------------------------------------
// 
//
private Real search_point_branch( Real x, Real y, Real z,
                                  int * t0, int * t1, int * t2,
                                  Real * w0, Real * w1, Real * w2,
                                  OctreeBranch * branch, Real curr_dist,
                                  SurfaceOctree * tree ) {

  int  i, outside;

  // Make sure that the point can close to this box within curr_dist.
  // Does the bounding box of the point closer than curr_dist from 
  // the bounding box of the branch?

  outside = 0;
  if( x + curr_dist < branch->minbox[0] ) outside = 1;
  if( x - curr_dist > branch->maxbox[0] ) outside = 1;
  if( y + curr_dist < branch->minbox[1] ) outside = 1;
  if( y - curr_dist > branch->maxbox[1] ) outside = 1;
  if( z + curr_dist < branch->minbox[2] ) outside = 1;
  if( z - curr_dist > branch->maxbox[2] ) outside = 1;

  // Check minimum distance between this point and all
  // triangles in this box.

  if( !outside ) {

    Real dist_sq, curr_dist_sq, tmp_w0, tmp_w1, tmp_w2;

    curr_dist_sq = curr_dist * curr_dist;

    if( branch->n_size > 0 ) {

      // it's a terminal leaf - search it
      for( i = 0; i < branch->n_size; i++ ) {
        int pp = branch->list_polys[i];
        int n0 = 3*tree->polys[3*pp];
        int n1 = 3*tree->polys[3*pp+1];
        int n2 = 3*tree->polys[3*pp+2];

        // Find the closest distance from point to the triangle.
        dist_sq = single_point_triangle_distance( x, y, z, &tmp_w0, &tmp_w1,
                                                  &tmp_w2, &(tree->xyz[n0]),
                                                  &(tree->xyz[n1]), &(tree->xyz[n2]) );
        if( dist_sq < curr_dist_sq ) {
          *t0 = n0/3;
          *t1 = n1/3;
          *t2 = n2/3;
          *w0 = tmp_w0;
          *w1 = tmp_w1;
          *w2 = tmp_w2;
          curr_dist_sq = dist_sq;
        }
      }
      curr_dist = sqrt( curr_dist_sq );
    } else {
      // it's a branch - look at its sub-branches
      if( branch->left ) {
        curr_dist = search_point_branch( x, y, z, t0, t1, t2, w0, w1, w2,
                                         branch->left, curr_dist, tree );
      }
      if( branch->right ) {
        curr_dist = search_point_branch( x, y, z, t0, t1, t2, w0, w1, w2,
                                         branch->right, curr_dist, tree );
      }
    }
  }

  return( curr_dist );
}


//
// Find the minimum distance from the point (x,y,z) to the 
// triangle (t0,t1,t2). Also compute the location of the point
// in terms of the local coordinates (w0,w1,w2) of the triangle.
//
private Real single_point_triangle_distance( Real x, Real y, Real z,
                                             Real * w0, Real * w1, Real * w2,
                                             Real t0[3], Real t1[3], Real t2[3] ) {

  // Point (x,y,z) on triangle (t0,t1,t2).

  Real A_tt = ( t1[0] - t0[0] ) * ( t1[0] - t0[0] ) +
              ( t1[1] - t0[1] ) * ( t1[1] - t0[1] ) +
              ( t1[2] - t0[2] ) * ( t1[2] - t0[2] );
  Real A_st = ( t1[0] - t0[0] ) * ( t2[0] - t0[0] ) +
              ( t1[1] - t0[1] ) * ( t2[1] - t0[1] ) +
              ( t1[2] - t0[2] ) * ( t2[2] - t0[2] );
  Real A_ss = ( t2[0] - t0[0] ) * ( t2[0] - t0[0] ) +
              ( t2[1] - t0[1] ) * ( t2[1] - t0[1] ) +
              ( t2[2] - t0[2] ) * ( t2[2] - t0[2] );
  Real det = A_tt * A_ss - A_st * A_st;

  if( det == 0.0 ) {
    // ????? something is wrong???
    // printf( "fatal error in intersection\n" );
    // exit(1);
    // looks like a flat triangle, so project point to
    // a different closest triangle instead.
    return( 1.0e20 );
  } else {
    Real rhs_t, rhs_s, xi, eta, dx, dy, dz, dist_sq, min_dist_sq;
    // Point p0.
    rhs_t = ( t1[0] - t0[0] ) * ( x - t0[0] ) +
            ( t1[1] - t0[1] ) * ( y - t0[1] ) +
            ( t1[2] - t0[2] ) * ( z - t0[2] );
    rhs_s = ( t2[0] - t0[0] ) * ( x - t0[0] ) +
            ( t2[1] - t0[1] ) * ( y - t0[1] ) +
            ( t2[2] - t0[2] ) * ( z - t0[2] );
    xi = ( A_ss * rhs_t - A_st * rhs_s ) / det;
    eta = ( A_tt * rhs_s - A_st * rhs_t ) / det;
    // clip to nearest side of triangle
    if( xi + eta > 1.0 ) {    // clip on line xi+eta=1
      xi = xi / ( xi + eta );
      if( xi < 0.0 ) {
        xi = 0.0;
      } else if( xi > 1.0 ) {
        xi = 1.0;
      }
      eta = 1.0 - xi;
    } else if( xi < 0.0 ) {      // clip on line xi=0
      eta = eta / ( 1.0 - xi );
      xi = 0.0;
      if( eta < 0.0 ) {
        eta = 0.0;
      } else if( eta > 1.0 ) {
        eta = 1.0;
      }
    } else if( eta < 0.0 ) {      // clip on line eta=0
      xi = xi / ( 1.0 - eta );
      eta = 0.0;
      if( xi < 0.0 ) {
        xi = 0.0;
      } else if( xi > 1.0 ) {
        xi = 1.0;
      }
    }

    *w0 = ( 1.0 - xi - eta );
    *w1 = xi;
    *w2 = eta;

    dx = ( 1.0 - xi - eta ) * t0[0] + xi * t1[0] + eta * t2[0] - x;
    dy = ( 1.0 - xi - eta ) * t0[1] + xi * t1[1] + eta * t2[1] - y;
    dz = ( 1.0 - xi - eta ) * t0[2] + xi * t1[2] + eta * t2[2] - z;
    min_dist_sq = dx * dx + dy * dy + dz * dz;

    return( min_dist_sq );
  }
  return( 0.0 );   // will never get here
}

