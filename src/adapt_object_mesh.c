/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
/*
   adapt_object_mesh.c

   Perform mesh smoothing, coarsening, swapping to improve
   the quality of a surface mesh with spherical topology.

   adapt_object_mesh input.obj output.obj target_nodes [n_iters]

   Values: input.obj = input object file
           output.obj = output object file after adaptation
           target_nodes = target number of points
           n_iterations = number of iterations

   Author: Claude Lepage, May 2011
*/

#include <math.h>
#include <stdio.h>
#include <volume_io.h>
#include <bicpl.h>

// Prototypes of functions in this file.

static void usage( char * );
static Status read_surface_obj( STRING, int *, Point *[],
                                Vector *[], int *, int *[] );
static void compute_triangle_angles( int, int *, Point * );
static void compute_surface_normals( int, int, int *, Point *, Vector * );
static int compute_triangle_normal( Point, Point, Point, Real [3] );
void smooth( int, int, int, Real, int *, Point *, float );

// Main program.

int main( int argc, char * argv[] ) {

  int      i, j, k, jj, kk, pp, v1, v2, opp1, opp2, t1, t2, count, found3;
  int      target_nodes, max_num_iters, max_sm_iters, changed, 
           num_swapped, histo[101];
  Real     thresholdlen, minlen, maxlen, factor;

  int      n_points;           // number of grid points per object
  int      n_elems;            // number of triangles per object
  Point  * coords;             // coordinates
  Vector * normals;            // normal vectors
  int    * connec;             // connectivity

  int      n_edges;            // total number of edges
  int    * edges;              // list of edges
  int    * countNgh;           // number of neighbours per point
  int    * cumulNgh;           // cumulative number of neighbours per point
  int    * triNgh;             // triangle neighbours per point
  int    * flags;              // active points
  int    * renum;              // new node numbers
  Real   * edgelen;            // length of edges

  FILE   * fp;

  if( argc < 4 ) {
    usage( argv[0] );
    return( 1 );
  }

  target_nodes = atoi( argv[3] );
  if( argc == 5 ) {
    max_sm_iters = atoi( argv[4] );
  } else {
    max_sm_iters = 10;
  }

  // Read the surface file.
  if( read_surface_obj( argv[1], &n_points, &coords, &normals,
                        &n_elems, &connec ) != OK ) {
    return 1;
  }
  printf( "Initial number of nodes = %d\n", n_points );
  printf( "Initial number of triangles = %d\n", n_elems );

  factor = 1.0;
  max_num_iters = 100;

  if( n_points > target_nodes ) {

  do {

    printf( "Iteration %d...\n", max_num_iters );

    changed = 0;
    num_swapped = 0;

    // Count number of neighbours per point.
    countNgh = (int *)malloc( n_points * sizeof( int ) );
    cumulNgh = (int *)malloc( ( n_points+1 ) * sizeof( int ) );
    triNgh = (int *)malloc( 3 * n_elems * sizeof( int ) );

    for( i = 0; i < n_points; i++ ) {
      countNgh[i] = 0;
    }
    for( i = 0; i < n_elems; i++ ) {
      triNgh[3*i] = -1;
      triNgh[3*i+1] = -1;
      triNgh[3*i+2] = -1;
      countNgh[connec[3*i]]++;
      countNgh[connec[3*i+1]]++;
      countNgh[connec[3*i+2]]++;
    }
    cumulNgh[0] = 0;
    for( i = 0; i < n_points; i++ ) {
      cumulNgh[i+1] = cumulNgh[i] + countNgh[i];
    }
    // Make a histogram for number of neighbours.
    j = 1;
    count = 0;
    found3 = 0;
    do {
      k = 0;
      for( i = 0; i < n_points; i++ ) {
        if( countNgh[i] == j ) k++;
      }
      if( k > 0 ) printf( "%d points with %d neighbours\n", k, j );
      if( j == 3 && k > 0 ) found3 = 1;
      j++;
      count += k;
    } while( count < n_points && j < 50 );

    // Find triangle neighbours at each point.
    for( i = 0; i < n_elems; i++ ) {
      for( j = 0; j < 3; j++ ) {
        pp = connec[3*i+j];
        for( k = cumulNgh[pp]; k < cumulNgh[pp+1]; k++ ) {
          if( triNgh[k] == i ) break;
          if( triNgh[k] == -1 ) {
            triNgh[k] = i;
            break;
          }
        }
      }
    }

    // Check for pinched vertices on the surface. 
    int bad_connect = 0;
    int ngh_list[50];
    int ngh_done[50];
    for( i = 0; i < n_points; i++ ) {
      for( j = 0; j < cumulNgh[i+1]-cumulNgh[i]; j++ ) {
        ngh_done[j] = 0;
      }
      for( jj = 0; jj < 3; jj++ ) {
        if( connec[3*triNgh[cumulNgh[i]]+jj] == i ) break;
      }
      ngh_list[0] = connec[3*triNgh[cumulNgh[i]]+(jj+1)%3];
      ngh_list[1] = connec[3*triNgh[cumulNgh[i]]+(jj+2)%3];
      int ngh_count = 1;
      ngh_done[0] = 1;

      do {
        int changed = 0;
        for( k = cumulNgh[i]; k < cumulNgh[i+1]; k++ ) {
          if( ngh_done[k-cumulNgh[i]] ) continue;
          for( jj = 0; jj < 3; jj++ ) {
            if( connec[3*triNgh[k]+jj] == ngh_list[ngh_count] ) {
              ngh_count++;
              ngh_list[ngh_count] = connec[3*triNgh[k]+(jj+1)%3];
              ngh_done[k-cumulNgh[i]] = 1;
              changed = 1;
              break;
            }
          }
        }
        if( !changed ) break;
      } while( 1  );

      if( ngh_list[0] == ngh_list[ngh_count] ) {
        if( ngh_count < cumulNgh[i+1] - cumulNgh[i] ) {
          printf( "A:Wrong connectivity for vertex %d at %f,%f,%f\n", i,
                  coords[i].coords[0], coords[i].coords[1],
                  coords[i].coords[2] );
          for( k = cumulNgh[i]; k < cumulNgh[i+1]; k++ ) {
            printf( "  %d %d %d\n", connec[3*triNgh[k]],
                    connec[3*triNgh[k]+1], connec[3*triNgh[k]+2] );
          }
          bad_connect++;
        }
      } else {
        printf( "B:Wrong connectivity for vertex %d at %f,%f,%f\n", i,
                coords[i].coords[0], coords[i].coords[1],
                coords[i].coords[2] );
        for( k = cumulNgh[i]; k < cumulNgh[i+1]; k++ ) {
          printf( "  %d %d %d\n", connec[3*triNgh[k]],
                  connec[3*triNgh[k]+1], connec[3*triNgh[k]+2] );
        }
        printf( "ngh_count = %d\n", ngh_count );
        printf( "list: " );
        for( k = 0; k <= ngh_count; k++ ) {
          printf( "%d ", ngh_list[k] );
        }
        printf( "\n\n" );
        bad_connect++;
      }
    }
    if( bad_connect ) {
      printf( "%d vertices with bad connectivity.\n", bad_connect );
      exit(1);
    }

    // Create the edges from the triangles.

    n_edges = ( 3 * n_elems ) / 2;
    edges = (int *)malloc( 4 * n_edges * sizeof( int ) );
    edgelen = (Real *)malloc( n_edges * sizeof( Real ) );

    count = 0;
    minlen = 1.0e20;
    maxlen = -1.0e20;
    for( i = 0; i < n_elems; i++ ) {
      for( j = 0; j < 3; j++ ) {
        v1 = connec[3*i+j];
        v2 = connec[3*i+(j+1)%3];

        if( v1 < v2 ) {
          edges[4*count] = v1;
          edges[4*count+1] = v2;
          edges[4*count+2] = i;
          edges[4*count+3] = -1;
          for( k = cumulNgh[v1]; k < cumulNgh[v1+1]; k++ ) {
            if( triNgh[k] != i ) {
              if( connec[3*triNgh[k]] == v2 ||
                  connec[3*triNgh[k]+1] == v2 ||
                  connec[3*triNgh[k]+2] == v2 ) {
                edges[4*count+3] = triNgh[k];
                break;
              }
            }
          }
          if( edges[4*count+3] == -1 ) {
            printf( "Edge %d:%d has only one triangle %d attached\n",
                    v1, v2, i );
            exit(1);
          }

          edgelen[count] = sqrt( ( coords[v1].coords[0] - coords[v2].coords[0] ) *
                                 ( coords[v1].coords[0] - coords[v2].coords[0] ) +
                                 ( coords[v1].coords[1] - coords[v2].coords[1] ) *
                                 ( coords[v1].coords[1] - coords[v2].coords[1] ) +
                                 ( coords[v1].coords[2] - coords[v2].coords[2] ) *
                                 ( coords[v1].coords[2] - coords[v2].coords[2] ) );
          if( edgelen[count] < minlen ) minlen = edgelen[count];
          if( edgelen[count] > maxlen ) maxlen = edgelen[count];
          count++;
        }
      }
    }

    if( count != n_edges ) {
      printf( "count = %d  n_edges = %d\n", count, n_edges );
      return 1;
    }

    // Make histogram of edge lengths to find coarsening threshold.
    Real delta = ( maxlen - minlen ) / 100.0;
    for( i = 0; i < 101; i++ ) {
      histo[i] = 0;
    }
    for( i = 0; i < n_edges; i++ ) {
      j = (int)( ( edgelen[i] - minlen ) / delta );
      histo[j]++;
    }
    int target_edges = (int)( (float)n_edges * (float)target_nodes / (float)n_points );
    target_edges = ( n_edges - target_edges ) / 2;

    for( i = 0; i < 101; i++ ) {
      target_edges -= histo[i];
      if( target_edges < 0 ) break;
    }
    thresholdlen = minlen + ((Real)i+0.5) * delta * factor;   // close enough

    printf( "Found %d internal edges\n", n_edges );
    printf( "Edge length threshold = %g  factor = %g\n", thresholdlen, factor );
    printf( "Minimum edge length   = %g\n", minlen );
    printf( "Maximum edge length   = %g\n", maxlen );

    // Loop over edges to coarsen.

    flags = (int *)malloc( n_points * sizeof( int ) );
    renum = (int *)malloc( n_points * sizeof( int ) );

    for( i = 0; i < n_points; i++ ) {
      flags[i] = 0;       // 0 active, 1 off
      renum[i] = i;
    }

    for( i = 0; i < n_edges; i++ ) {
      v1 = edges[4*i];
      v2 = edges[4*i+1];     // we know v1 < v2.
      t1 = edges[4*i+2];
      t2 = edges[4*i+3];
      if( flags[v1] || flags[v2] ) continue;
      if( connec[3*t1] != v1 && connec[3*t1] != v2 ) opp1 = connec[3*t1];
      if( connec[3*t1+1] != v1 && connec[3*t1+1] != v2 ) opp1 = connec[3*t1+1];
      if( connec[3*t1+2] != v1 && connec[3*t1+2] != v2 ) opp1 = connec[3*t1+2];
      if( connec[3*t2] != v1 && connec[3*t2] != v2 ) opp2 = connec[3*t2];
      if( connec[3*t2+1] != v1 && connec[3*t2+1] != v2 ) opp2 = connec[3*t2+1];
      if( connec[3*t2+2] != v1 && connec[3*t2+2] != v2 ) opp2 = connec[3*t2+2];
      if( opp1 == opp2 ) {
        printf( "edge %d : v1 = %d  v2 = %d\n", i, v1, v2 );
        printf( "v1 at %g %g %g\n", coords[v1].coords[0], coords[v1].coords[1], 
                coords[v1].coords[2] );
        printf( "v2 at %g %g %g\n", coords[v2].coords[0], coords[v2].coords[1], 
                coords[v2].coords[2] );
        printf( "opp1 = opp2 = %d\n", opp1 );
        printf( "t1 = %d %d %d\n", connec[3*t1], connec[3*t1+1], connec[3*t1+2] );
        printf( "t2 = %d %d %d\n", connec[3*t2], connec[3*t2+1], connec[3*t2+2] );
        printf( "ngh of t1:\n" );
        for( j = cumulNgh[v1]; j < cumulNgh[v1+1]; j++ ) {
          int tt = triNgh[j];
          printf( "tri %d with nodes %d %d %d\n", j, connec[3*tt], 
                  connec[3*tt+1], connec[3*tt+2] );
        }
        printf( "ngh of t2:\n" );
        for( j = cumulNgh[v2]; j < cumulNgh[v2+1]; j++ ) {
          int tt = triNgh[j];
          printf( "tri %d with nodes %d %d %d\n", j, connec[3*tt], 
                  connec[3*tt+1], connec[3*tt+2] );
        }
        return 1;
      }
      if( flags[opp1] || flags[opp2] ) continue;
      if( countNgh[v1] + countNgh[v2] < 7 ) {
        printf( "Strange edge %d:%d with nghrs %d:%d\n", v1, v2, countNgh[v1], countNgh[v2] );
        printf( "node %d at %g %g %g\n", v1, coords[v1].coords[0], coords[v1].coords[1],
                coords[v1].coords[2] );
        printf( "node %d at %g %g %g\n", v2, coords[v2].coords[0], coords[v2].coords[1],
                coords[v2].coords[2] );
        continue;
      }

      // Special configuration to eliminate: a triangle with exactly 4 neighbours
      // for each of its 3 vertices, all enclosed in a larger triangle.

      if( countNgh[v1] == 4 && countNgh[v2] == 4 && 
          ( countNgh[opp1] == 4 || countNgh[opp2] == 4 ) ) {

        if( countNgh[opp1] == 4 && countNgh[opp2] == 4 ) {
          continue;
          printf( "Problem with connectivity: 4 neighbours for verts %d %d %d %d\n",
                  v1, v2, opp1, opp2 );

          printf( "v1 = %d (%d) v2 = %d (%d) opp1 = %d (%d) opp2 = %d (%d)\n", 
                  v1, countNgh[v1], v2, countNgh[v2], opp1, countNgh[opp1], 
                  opp2, countNgh[opp2] );
          printf( "Triangles of v1:\n" );
          for( j = cumulNgh[v1]; j < cumulNgh[v1+1]; j++ ) {
            printf( "  tri:%d %d(%d) %d(%d) %d(%d)\n", triNgh[j], connec[3*triNgh[j]],
                    flags[connec[3*triNgh[j]]], connec[3*triNgh[j]+1], 
                    flags[connec[3*triNgh[j]+1]], connec[3*triNgh[j]+2],
                    flags[connec[3*triNgh[j]+2]] );
          }
          printf( "Triangles of v2:\n" );
          for( j = cumulNgh[v2]; j < cumulNgh[v2+1]; j++ ) {
            printf( "  tri:%d %d(%d) %d(%d) %d(%d)\n", triNgh[j], connec[3*triNgh[j]],
                    flags[connec[3*triNgh[j]]], connec[3*triNgh[j]+1], 
                    flags[connec[3*triNgh[j]+1]], connec[3*triNgh[j]+2],
                    flags[connec[3*triNgh[j]+2]] );
          }
          printf( "Triangles of opp1:\n" );
          for( j = cumulNgh[opp1]; j < cumulNgh[opp1+1]; j++ ) {
            printf( "  tri:%d %d(%d) %d(%d) %d(%d)\n", triNgh[j], connec[3*triNgh[j]],
                    flags[connec[3*triNgh[j]]], connec[3*triNgh[j]+1], 
                    flags[connec[3*triNgh[j]+1]], connec[3*triNgh[j]+2],
                    flags[connec[3*triNgh[j]+2]] );
          }
          printf( "Triangles of opp2:\n" );
          for( j = cumulNgh[opp2]; j < cumulNgh[opp2+1]; j++ ) {
            printf( "  tri:%d %d(%d) %d(%d) %d(%d)\n", triNgh[j], connec[3*triNgh[j]],
                    flags[connec[3*triNgh[j]]], connec[3*triNgh[j]+1], 
                    flags[connec[3*triNgh[j]+1]], connec[3*triNgh[j]+2],
                    flags[connec[3*triNgh[j]+2]] );
          }

          exit(1);
        }

        // Find the 7 triangles in the cluster inside the larger triangle.

        int vt1, vt2, vt3, other;
        if( countNgh[opp1] == 4 ) {
          other = opp1;
          vt1 = opp2;
        } else {
          other = opp2;
          vt1 = opp1;
        }

        int   cluster[7];
        int   ncluster = 0;
        for( j = cumulNgh[v1]; j < cumulNgh[v1+1]; j++ ) {
          cluster[ncluster] = triNgh[j];
          ncluster++;
        }
        for( j = cumulNgh[v2]; j < cumulNgh[v2+1]; j++ ) {
          for( k = 0; k < ncluster; k++ ) {
            if( cluster[k] == triNgh[j] ) break;
          }
          if( k == ncluster ) {
            cluster[ncluster] = triNgh[j];
            ncluster++;
          }
          if( ncluster > 7 ) {
            printf( "Problem with cluster for v2: 4 neighbours for verts %d %d %d %d\n",
                    v1, v2, opp1, opp2 );
            exit(1);
          }
        }
        for( j = cumulNgh[other]; j < cumulNgh[other+1]; j++ ) {
          for( k = 0; k < ncluster; k++ ) {
            if( cluster[k] == triNgh[j] ) break;
          }
          if( k == ncluster ) {
            cluster[ncluster] = triNgh[j];
            ncluster++;
          }
          if( ncluster > 7 ) {
            printf( "Problem with cluster for v2: 4 neighbours for verts %d %d %d %d\n",
                    v1, v2, opp1, opp2 );
            exit(1);
          }
        }

        if( ncluster == 7 ) {
          int global_flag = 0;
          for( k = 0; k < ncluster; k++ ) {
            global_flag += flags[connec[3*cluster[k]]] + 
                           flags[connec[3*cluster[k]+1]] + 
                           flags[connec[3*cluster[k]+2]];
          }
          if( global_flag == 0 ) {

//          printf( "Remove this 7-cluster of triangles...\n" );
//          printf( "Cluster: " );
//          for( k = 0; k < ncluster; k++ ) {
//            printf( "%d ", cluster[k] );
//          }
//          printf( "\n" );

            int largerTri[3], nlT = 0;
            for( k = 0; k < ncluster; k++ ) {
              for( j = 0; j < 3; j++ ) {
                if( connec[3*cluster[k]+j] != v1 &&
                    connec[3*cluster[k]+j] != v2 &&
                    connec[3*cluster[k]+j] != other ) {
                  for( jj = 0; jj < nlT; jj++ ) {
                    if( largerTri[jj] == connec[3*cluster[k]+j] ) break;
                  }
                  if( jj == nlT ) {
                    largerTri[nlT] = connec[3*cluster[k]+j];
                    nlT++;
                  }
                }
              }
            }
            if( countNgh[largerTri[0]] > 5 && countNgh[largerTri[1]] > 5 &&
                countNgh[largerTri[2]] > 5 ) {

//            printf( "Larger triangle (%d): ", nlT );
//            for( k = 0; k < nlT; k++ ) {
//              printf( "%d (%d)", largerTri[k], countNgh[largerTri[k]] );
//            }
//            printf( "\n" );
          
              renum[v1] = -1;
              renum[v2] = -1;
              renum[other] = -1;
              countNgh[v1] = 0;
              countNgh[v2] = 0;
              countNgh[other] = 0;
              for( k = 0; k < ncluster; k++ ) {
                int newTri = 0;
                for( j = 0; j < 3; j++ ) {
                  flags[connec[3*cluster[k]+j]] = 1;
                  if( connec[3*cluster[k]+j] == largerTri[0] ||
                      connec[3*cluster[k]+j] == largerTri[1] ) newTri++;
                }
                if( newTri == 2 ) {
                  // re-use this triangle for the larger triangle
                  for( j = 0; j < 3; j++ ) {
                    if( connec[3*cluster[k]+j] != largerTri[0] &&
                        connec[3*cluster[k]+j] != largerTri[1] ) {
                      connec[3*cluster[k]+j] = largerTri[2];
                      break;
                    }
                  }
                } else {
                  // delete this inner triangle
                  for( j = 0; j < 3; j++ ) {
                    connec[3*cluster[k]+j] = -1;
                  }
                }
              }
              for( j = 0; j < 3; j++ ) {
                countNgh[largerTri[j]] -= 2;
              }

              changed++;
              continue;
            }
          }
        }

      }

      int coarsen = 0;
      if( found3 ) {
        if( countNgh[v1] == 3 || countNgh[v2] == 3 ) coarsen = 1;
      } else {
        if( countNgh[v1] <= 4 || countNgh[v2] <= 4 ) coarsen = 1;
        // Coarsen if surface is flat enough and edgelen < threshold.
        if( !coarsen ) {
          if( edgelen[i] <= thresholdlen ) {
            Real n1x = normals[v1].coords[0];
            Real n1y = normals[v1].coords[1];
            Real n1z = normals[v1].coords[2];
            Real n2x = normals[v2].coords[0];
            Real n2y = normals[v2].coords[1];
            Real n2z = normals[v2].coords[2];
            Real mag = n1x * n2x + n1y * n2y + n1z * n2z;
            if( mag > 0.99 ) {   // cosine(angle)
              coarsen = 1;
            }
          }
        }
        // This would create 3 neighbours on opp1 or opp2 after collapsed.
        if( countNgh[opp1] <= 4 || countNgh[opp2] <= 4 ) coarsen = 0;
      }

      // make sure remaining connectivity is ok: any outer triangle of v1
      // must not touch any outer triangle of v2.
      if( coarsen && countNgh[v1] > 3 && countNgh[v2] > 3 ) {
        for( j = cumulNgh[v1]; j < cumulNgh[v1+1]; j++ ) {
          int tt = triNgh[j];
          if( tt != t1 && tt != t2 ) {
            for( k = 0; k < 3; k++ ) {
              pp = renum[connec[3*tt+k]];
              if( pp != v1 && pp != opp1 && pp != opp2 ) {
                for( jj = cumulNgh[v2]; jj < cumulNgh[v2+1]; jj++ ) {
                  if( triNgh[jj] != t1 && triNgh[jj] != t2 ) {
                    for( kk = 0; kk < 3; kk++ ) {
                      if( renum[connec[3*triNgh[jj]+kk]] == pp ) coarsen = 0;
                    }
                  }
                  if( !coarsen ) break;
                }
              }
              if( !coarsen ) break;
            }
          }
          if( !coarsen ) break;
        }
      }

      if( coarsen ) {
        changed++;

        renum[v2] = v1;      // we know v1 < v2.
        flags[v1] = 1;
        flags[v2] = 1;
        flags[opp1] = 1;
        flags[opp2] = 1;
        connec[3*t1] = -1;
        connec[3*t1+1] = -1;
        connec[3*t1+2] = -1;
        connec[3*t2] = -1;
        connec[3*t2+1] = -1;
        connec[3*t2+2] = -1;

        // Make the new point at the geometric center of the neighbours.

        count = 0;
        coords[v1].coords[0] = 0.0;
        coords[v1].coords[1] = 0.0;
        coords[v1].coords[2] = 0.0;
        normals[v1].coords[0] = 0.0;
        normals[v1].coords[1] = 0.0;
        normals[v1].coords[2] = 0.0;
        for( j = cumulNgh[v1]; j < cumulNgh[v1+1]; j++ ) {
          int tt = triNgh[j];
          if( tt != t1 && tt != t2 ) {
            for( k = 0; k < 3; k++ ) {
              pp = connec[3*tt+k];
              if( pp != v1 && pp != v2 ) {
                coords[v1].coords[0] += coords[pp].coords[0];
                coords[v1].coords[1] += coords[pp].coords[1];
                coords[v1].coords[2] += coords[pp].coords[2];
                normals[v1].coords[0] += normals[pp].coords[0];
                normals[v1].coords[1] += normals[pp].coords[1];
                normals[v1].coords[2] += normals[pp].coords[2];
                count++;
              }
            }
          }
        }
        for( j = cumulNgh[v2]; j < cumulNgh[v2+1]; j++ ) {
          int tt = triNgh[j];
          if( tt != t1 && tt != t2 ) {
            for( k = 0; k < 3; k++ ) {
              pp = connec[3*tt+k];
              if( pp != v1 && pp != v2 ) {
                coords[v1].coords[0] += coords[pp].coords[0];
                coords[v1].coords[1] += coords[pp].coords[1];
                coords[v1].coords[2] += coords[pp].coords[2];
                normals[v1].coords[0] += normals[pp].coords[0];
                normals[v1].coords[1] += normals[pp].coords[1];
                normals[v1].coords[2] += normals[pp].coords[2];
                count++;
              }
            }
          }
        }
        coords[v1].coords[0] /= (Real)count;
        coords[v1].coords[1] /= (Real)count;
        coords[v1].coords[2] /= (Real)count;
        Real mag = sqrt( normals[v1].coords[0] * normals[v1].coords[0] +
                         normals[v1].coords[1] * normals[v1].coords[1] + 
                         normals[v1].coords[2] * normals[v1].coords[2] );
        normals[v1].coords[0] /= mag;
        normals[v1].coords[1] /= mag;
        normals[v1].coords[2] /= mag;

        countNgh[opp1]--;
        countNgh[opp2]--;
        countNgh[v1] = countNgh[v1] + countNgh[v2] - 4;
        countNgh[v2] = countNgh[v1];
      }

      // Try to swap is no coarsening is available. We can swap to improve:
      //    - connectivity (get closer to 6 neighbours per node)
      //    - eliminate obtuse angles (> 90 deg)
      if( !coarsen ) {
      // if( !coarsen && !found3 ) {
        int swap = 0;
        Real a1x, a1y, a1z, a2x, a2y, a2z, mag, cc1, cc2, cc3, cc4, nn1, nn2;
        Real n1x, n1y, n1z, n2x, n2y, n2z;
        Real norm1[3], norm2[3];

        // check if connectivity would benefit from a swap.
        if( 4 + countNgh[opp1] + countNgh[opp2] <= countNgh[v1] + countNgh[v2] ) swap = 1;
        if( ( countNgh[opp1] <= 4 || countNgh[opp2] <= 4 ) &&
            ( countNgh[v1] > 4 && countNgh[v2] > 4 ) ) swap = 1;

        // Check if triangles have obtuse angles (> 120 deg) before swap.
        if( !swap ) {
          a1x = coords[v1].coords[0] - coords[opp1].coords[0];
          a1y = coords[v1].coords[1] - coords[opp1].coords[1];
          a1z = coords[v1].coords[2] - coords[opp1].coords[2];
          a2x = coords[v2].coords[0] - coords[opp1].coords[0];
          a2y = coords[v2].coords[1] - coords[opp1].coords[1];
          a2z = coords[v2].coords[2] - coords[opp1].coords[2];
          mag = sqrt( ( a1x * a1x + a1y * a1y + a1z * a1z ) * 
                      ( a2x * a2x + a2y * a2y + a2z * a2z ) );
          cc1 = ( a1x * a2x + a1y * a2y + a1z * a2z ) / mag;
  
          a1x = coords[v1].coords[0] - coords[opp2].coords[0];
          a1y = coords[v1].coords[1] - coords[opp2].coords[1];
          a1z = coords[v1].coords[2] - coords[opp2].coords[2];
          a2x = coords[v2].coords[0] - coords[opp2].coords[0];
          a2y = coords[v2].coords[1] - coords[opp2].coords[1];
          a2z = coords[v2].coords[2] - coords[opp2].coords[2];
          mag = sqrt( ( a1x * a1x + a1y * a1y + a1z * a1z ) * 
                      ( a2x * a2x + a2y * a2y + a2z * a2z ) );
          cc2 = ( a1x * a2x + a1y * a2y + a1z * a2z ) / mag;
          if( cc1 <= -0.5 || cc2 <= -0.5 ) swap = 1;
        }

        // Make sure than new edge (opp1,opp2) does not exist.
        if( swap ) {
          for( j = cumulNgh[opp1]; j < cumulNgh[opp1+1]; j++ ) {
            int tt = triNgh[j];
            for( k = 0; k < 3; k++ ) {
              pp = connec[3*tt+k];
              if( renum[pp] == renum[opp2] ) swap = 0;
            }
          }
        }

        // Make sure that swap does not introduce a vertex with
        // only 3 vertices.
        if( swap ) {
          if( countNgh[v1] <= 4 || countNgh[v2] <= 4 ) swap = 0;
        }

        // Check triangle angles after. Are there worse obtuse angles than before? 
        // Ignore improvement in angles if an opp node has only 4 neighbours, in
        // which case, swap it.
        if( swap && countNgh[opp1] > 4 && countNgh[opp2] > 4 ) {
          a1x = coords[opp1].coords[0] - coords[v1].coords[0];
          a1y = coords[opp1].coords[1] - coords[v1].coords[1];
          a1z = coords[opp1].coords[2] - coords[v1].coords[2];
          a2x = coords[opp2].coords[0] - coords[v1].coords[0];
          a2y = coords[opp2].coords[1] - coords[v1].coords[1];
          a2z = coords[opp2].coords[2] - coords[v1].coords[2];
          mag = sqrt( ( a1x * a1x + a1y * a1y + a1z * a1z ) * 
                      ( a2x * a2x + a2y * a2y + a2z * a2z ) );
          cc3 = ( a1x * a2x + a1y * a2y + a1z * a2z ) / mag;

          a1x = coords[opp1].coords[0] - coords[v2].coords[0];
          a1y = coords[opp1].coords[1] - coords[v2].coords[1];
          a1z = coords[opp1].coords[2] - coords[v2].coords[2];
          a2x = coords[opp2].coords[0] - coords[v2].coords[0];
          a2y = coords[opp2].coords[1] - coords[v2].coords[1];
          a2z = coords[opp2].coords[2] - coords[v2].coords[2];
          mag = sqrt( ( a1x * a1x + a1y * a1y + a1z * a1z ) * 
                      ( a2x * a2x + a2y * a2y + a2z * a2z ) );
          cc4 = ( a1x * a2x + a1y * a2y + a1z * a2z ) / mag;

          if( cc2 < cc1 ) cc1 = cc2;  // worse angle before
          if( cc4 < cc3 ) cc3 = cc4;  // worse angle after
          if( cc3 <= -0.5 && cc3 < cc1 ) swap = 0;
        }

        // Make sure swap does not create a bad configuration. Check
        // the normals of the swapped triangles. Both normals must be
        // in the same direction and no worse than 60 deg. We cannot
        // have flipped triangles which will lead to self-intersections.
        if( swap ) {
          // surface normals before and after. is it worse?
          if( compute_triangle_normal( coords[v1], coords[v2], coords[opp1], norm1 ) &&
              compute_triangle_normal( coords[v1], coords[opp2], coords[v2], norm2 ) ) {
            nn1 = norm1[0] * norm2[0] + norm1[1] * norm2[1] + norm1[2] * norm2[2];
            if( compute_triangle_normal( coords[v1], coords[opp2], coords[opp1], norm1 ) &&
                compute_triangle_normal( coords[v2], coords[opp1], coords[opp2], norm2 ) ) {
              nn2 = norm1[0] * norm2[0] + norm1[1] * norm2[1] + norm1[2] * norm2[2];
              // faces are inverted after, but allow if it makes the surface less worse.
              // (60 deg).
              if( nn2 <= 0.5 && nn2 < nn1 ) swap = 0;
            } else {
              swap = 0;
            }
          } else {
            swap = 0;
          }
        }

        if( swap ) {
          flags[v1] = 1;
          flags[v2] = 1;
          flags[opp1] = 1;
          flags[opp2] = 1;
          if( connec[3*t1] == v2 ) connec[3*t1] = opp2;
          if( connec[3*t1+1] == v2 ) connec[3*t1+1] = opp2;
          if( connec[3*t1+2] == v2 ) connec[3*t1+2] = opp2;
          if( connec[3*t2] == v1 ) connec[3*t2] = opp1;
          if( connec[3*t2+1] == v1 ) connec[3*t2+1] = opp1;
          if( connec[3*t2+2] == v1 ) connec[3*t2+2] = opp1;
          countNgh[v1]--;
          countNgh[v2]--;
          countNgh[opp1]++;
          countNgh[opp2]++;
          num_swapped++;
        }
      }
    }

    if( changed ) {
      printf( "%d edges collapsed\n", changed );
    }
    if( num_swapped ) {
      printf( "%d edges swapped\n", num_swapped );
    }

    // Adjust global data structures after removal of nodes.

    if( changed || num_swapped ) {
      count = 0;
      for( i = 0; i < n_points; i++ ) {
        if( renum[i] == i ) {
          coords[count].coords[0] = coords[i].coords[0];
          coords[count].coords[1] = coords[i].coords[1];
          coords[count].coords[2] = coords[i].coords[2];
          normals[count].coords[0] = normals[i].coords[0];
          normals[count].coords[1] = normals[i].coords[1];
          normals[count].coords[2] = normals[i].coords[2];
          renum[i] = count;
          count++;
        } else {
          renum[i] = renum[renum[i]];   // good since v1 < v2.
        }
      }
      n_points = count;
      printf( "New number of nodes = %d\n", n_points );

      count = 0;
      for( i = 0; i < n_elems; i++ ) {
        if( connec[3*i] >= 0 ) {
          connec[3*count] = renum[connec[3*i]];
          connec[3*count+1] = renum[connec[3*i+1]];
          connec[3*count+2] = renum[connec[3*i+2]];
          count++;
        }
      }
      n_elems = count;
      printf( "New number of triangles = %d\n", n_elems );
      printf( "\n" );

    } else {
      if( found3 ) {
        printf( "No change while removing nodes with 3 neighbours.\n" );
        // return 1;
      }
      factor *= 1.05;
    }

    // Smoothing on triangles with aspect ratio less than 0.10.
    smooth( n_points, n_elems, 25, 0.25, connec, coords, 0.20 );

    compute_surface_normals( n_elems, n_points, connec, coords, normals );
    compute_triangle_angles( n_elems, connec, coords );

    free( countNgh );
    free( cumulNgh );
    free( triNgh );
    free( edges );
    free( flags );
    free( renum );
    free( edgelen );

    max_num_iters--;
    if( max_num_iters <= 0 ) break;

    if( (float)changed / (float)target_nodes < 0.05 ) {
      if( n_points > target_nodes ) {
        factor *= ( 1.0 + 0.05 * (float)( n_points - target_nodes ) / target_nodes );
      }
    }

  } while( ( n_points > target_nodes ) /* || ( n_points <= target_nodes && changed ) */ ); }

  // Do a little bit of smoothing on the coordinates (simple averaging).
  Real   relax = 0.75;
  smooth( n_points, n_elems, max_sm_iters, relax, connec, coords, 1000.0 );

  compute_surface_normals( n_elems, n_points, connec, coords, normals );
  compute_triangle_angles( n_elems, connec, coords );

  // Quick and dirty output into .obj format

  fp = fopen( argv[2], "w" );
  fprintf( fp, "P 0.3 0.3 0.4 10 1 %d\n", n_points );

  // print the coords
  for( j = 0; j < n_points; j++ ) {
    fprintf( fp, "%g %g %g\n", coords[j].coords[0], 
             coords[j].coords[1], coords[j].coords[2] );
  }
  fprintf( fp, "\n" );
  FREE( coords );

  // print the normals
  for( j = 0; j < n_points; j++ ) {
    // fprintf( fp, "%g %g %g\n", -normals[j].coords[0], 
    //          -normals[j].coords[1], -normals[j].coords[2] );
    fprintf( fp, "%g %g %g\n", normals[j].coords[0], 
             normals[j].coords[1], normals[j].coords[2] );
  }
  FREE( normals );

  // The connectivity - part 1.
  fprintf( fp, "\n" );
  fprintf( fp, "%d\n", n_elems );
  fprintf( fp, "0 1 1 1 1\n\n" );

  for( i = 1; i <= n_elems; i++ ) {
    fprintf( fp, "%d ", 3*i );
    if( i%8 == 0 ) fprintf( fp, "\n" );
  }

  // The connectivity - part 2.

  count = 0;
  for( j = 0; j < 3*n_elems; j++ ) {
    if( count%8 == 0 ) fprintf( fp, "\n" );
    fprintf( fp, "%d ", connec[j] );
    count++;
  }
  FREE( connec );

  fclose( fp );

  return 0;
}

// Do smoothing on the coordinates (simple averaging).

void smooth( int n_points, int n_elems, int n_iters, Real relax, 
             int * connec, Point * coords, float minAR ) {


  int i, j, k, kk;
  Real edgeLen[3];
  Real * new_coords = (Real *)malloc( 3 * n_points * sizeof( Real ) );
  if( !new_coords ) {
    printf( "Error allocating memory for new_coords.\n" );
    exit( 1 );
  }

  char * countNgh = (char *)malloc( n_points * sizeof( char ) );
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
        edgeLen[j] = sqrt( dx * dx + dy * dy + dz * dz );
      }
      Real s = 0.5 * ( edgeLen[0] + edgeLen[1] + edgeLen[2] );
      Real invAR = 0.125 * ( edgeLen[0] * edgeLen[1] * edgeLen[2] ) /
                   ( ( s - edgeLen[0] ) * ( s - edgeLen[1] ) * ( s - edgeLen[2] ) );
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
        edgeLen[j] = sqrt( ( coords[k0].coords[0] - coords[k1].coords[0] ) *
                           ( coords[k0].coords[0] - coords[k1].coords[0] ) +
                           ( coords[k0].coords[1] - coords[k1].coords[1] ) *
                           ( coords[k0].coords[1] - coords[k1].coords[1] ) +
                           ( coords[k0].coords[2] - coords[k1].coords[2] ) *
                           ( coords[k0].coords[2] - coords[k1].coords[2] ) );
      }
      Real s = 0.5 * ( edgeLen[0] + edgeLen[1] + edgeLen[2] );
      Real area = sqrt( fabs( s * ( s - edgeLen[0] ) * ( s - edgeLen[1] ) * 
                              ( s - edgeLen[2] ) ) + 1.0e-10 );
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
      if( compute_triangle_normal( coords[v0], coords[v1], coords[v2], norm ) ) {
        for( j = 0; j < 3; j++ ) {
          normals[v0].coords[j] += norm[j];
          normals[v1].coords[j] += norm[j];
          normals[v2].coords[j] += norm[j];
        }
      }
    }

    for( i = 0; i < n_points; i++ ) {
      mag = sqrt( normals[i].coords[0] * normals[i].coords[0] +
                  normals[i].coords[1] * normals[i].coords[1] +
                  normals[i].coords[2] * normals[i].coords[2] );
      if( mag > 1.0e-10 ) {
        normals[i].coords[0] /= mag;
        normals[i].coords[1] /= mag;
        normals[i].coords[2] /= mag;
      } else {
        normals[i].coords[0] = 0.0;
        normals[i].coords[1] = 0.0;
        normals[i].coords[2] = 0.0;
      }
    }
}

// Compute the normal vector to a triangular face.

static int compute_triangle_normal( Point v1, Point v2, Point v3, Real norm[3] ) {

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
  if( mag > 1.0e-10 ) {
    norm[0] /= mag;
    norm[1] /= mag;
    norm[2] /= mag;
    return 1;
  } else {
    return 0;
  }
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

// -------------------------------------------------------------------
// Help message on how to use this module.
//
static void usage( char * executable_name ) {

  STRING  usage_format = "\
Usage: %s in.obj out.obj target_points [n_iterations]\n\
Values: in.obj = input object file\n\
        out.obj = output object file\n\
        target_points = target number of points\n\
        n_iterations = number of iterations\n\n\
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
//
static Status read_surface_obj( STRING filename,
                                 int * n_points,
                                 Point * points[],
                                 Vector * normals[],
                                 int * n_elem,
                                 int * connec[] ) {

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
  ALLOC( *normals, surface->n_points );
  ALLOC( *connec, 3*surface->n_items );

  for( i = 0; i < *n_points; i++ ) {
    (*points)[i].coords[0] = surface->points[i].coords[0];
    (*points)[i].coords[1] = surface->points[i].coords[1];
    (*points)[i].coords[2] = surface->points[i].coords[2];
    (*normals)[i].coords[0] = surface->normals[i].coords[0];
    (*normals)[i].coords[1] = surface->normals[i].coords[1];
    (*normals)[i].coords[2] = surface->normals[i].coords[2];
  }

  for( i = 0; i < 3*surface->n_items; i++ ) {
    (*connec)[i] = surface->indices[i];
  }

  delete_object_list( n_objects, object_list );

  return( OK );
}

