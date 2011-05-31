/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
/*
   Functions for surface-to-surface interpolation using an octree.

   Author: Claude Lepage, May 2011
*/

#ifndef  _DEF_OCTREE_H
#define  _DEF_OCTREE_H

#include  <math.h>
#include  <volume_io.h>
#include  <bicpl.h>

typedef struct OctreeBranch {
  int              id;           // for debugging purposes
  Real             minbox[3];    // bounding box
  Real             maxbox[3];
  int              n_size;       // number of triangles in box
  int            * list_polys;   // list of triangles in box
  struct OctreeBranch  * left;
  struct OctreeBranch  * right;
} OctreeBranch;


typedef struct {
  int             n_id;
  int             n_points;
  int             n_polygons;
  Real          * xyz;
  int           * polys;
  Real          * coeffs;
  Real          * polybbox;
  int           * edge_index;
  int           * edge_table;
  OctreeBranch ** edge_ptr;
  OctreeBranch  * root;
  OctreeBranch  * branch_cache;  // box of first search
} SurfaceOctree;

// Prototypes for public functions

void initialize_surface_octree( int, int [], int * [], Real [],
                                SurfaceOctree * );
void free_surface_octree( SurfaceOctree * );

Real search_edge_octree( int, int, SurfaceOctree * );
Real search_point_octree( Real, Real, Real, int *, int *, int *,
                          Real *, Real *, Real *, SurfaceOctree * );

#endif
