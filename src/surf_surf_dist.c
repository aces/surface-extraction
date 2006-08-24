#include  <fit_3d.h>

#define  MAX_THRESHOLD_N_POLYGONS   200

#define  THRESHOLD_N_POLYGONS   65

#ifdef PRINT_DIST
#define PRINT_DIST
static Real sum_dist;
#endif


private  void    recursive_find_close_pairs( 
    Real               min_distance,
    Real               search_distance,
    Real               *max_step_size,
    Real               parameters1[],
    poly_info_struct   poly_info1[],
    int                **list_of_polys1,
    unsigned int       **poly_classes1,
    int                *n_alloced1,
    int                offset_index1,
    int                n_polygons1,
    int                next_offset_index1,
    Real               line_dir1[],
    Real               parameters2[],
    poly_info_struct   poly_info2[],
    int                **list_of_polys2,
    unsigned int       **poly_classes2,
    int                *n_alloced2,
    int                offset_index2,
    int                n_polygons2,
    int                next_offset_index2,
    Real               line_dir2[],
    int                *n_pairs_ptr,
    int                *n_pairs_alloc,
    int                *p1s_ptr[],
    int                *p2s_ptr[],
    unsigned char      *cases_ptr[],
    float              *min_line_dists_ptr[],
    Real               *closest_distance )
{
    int               parms1[3], parms2[3];
    unsigned char     which_case;
    Real              mid_plane, min_distance_sq, search_distance_sq;
    Real              movement, movement_sq, dx, dy, dz, d, t_min, dist;
    int               which1, which2, n1, n2;
    float             mid_plane_plus, mid_plane_minus;
    Real              min_range[N_DIMENSIONS], max_range[N_DIMENSIONS];
    float             fl_search_distance;
    float             x1_low, x1_high, y1_low, y1_high, z1_low, z1_high;
    int               split_axis, dim, *list1, *list2, poly1;
    int               p1, p2, poly;
    unsigned int      *classes1, *src_classes1, class_bit;
    unsigned int      *new_classes11, *new_classes12;
    unsigned int      *classes2, *src_classes2;
    unsigned int      *new_classes21, *new_classes22;
    int               i, i1, i2, diff;
    int               *src_list1, *new_list11, *new_list12;
    int               new_n_polygons11, new_n_polygons12;
    int               *src_list2, *new_list21, *new_list22;
    int               new_n_polygons21, new_n_polygons22;
    int               n11, n12, n21, n22;
    int               parm_p1, parm_p2, parm_n11, parm_n12, parm_n21, parm_n22;
    BOOLEAN           on_left, on_right;
    unsigned int      cl;
    int               count1[17], count2[17], c1, c2;
    int               static_list1[MAX_THRESHOLD_N_POLYGONS];
    int               *sorted_list1;
    int               static_list2[MAX_THRESHOLD_N_POLYGONS];
    int               *sorted_list2;
    int               n_to_compare2;
    poly_info_struct  *static_compare2[MAX_THRESHOLD_N_POLYGONS];
    int               static_compare2_index[MAX_THRESHOLD_N_POLYGONS];
    int               *compare2_index;
    poly_info_struct  **compare2, *info_ptr;
    static int        threshold_n_polygons = THRESHOLD_N_POLYGONS;
    static BOOLEAN    first = TRUE;

    if( first ) {
        first = FALSE;
        if( getenv( "THRESHOLD_N_POLYGONS" ) == NULL ||
            sscanf( getenv("THRESHOLD_N_POLYGONS" ), "%d",
                    &threshold_n_polygons ) != 1  ||
            threshold_n_polygons > MAX_THRESHOLD_N_POLYGONS )
            threshold_n_polygons = THRESHOLD_N_POLYGONS;
    }

    if( n_polygons1 == 0 || n_polygons2 == 0 )
        return;

    if( n_polygons1 >= threshold_n_polygons ||
        n_polygons2 >= threshold_n_polygons ) {

      // Compute the true min/max bounding box of this list of polygons
      // (both surfaces together).
      list1 = &(*list_of_polys1)[offset_index1];
      for( dim = 0; dim < N_DIMENSIONS; dim++ ) {
        min_range[dim] = poly_info1[list1[0]].low_limits[dim];
        max_range[dim] = poly_info1[list1[0]].high_limits[dim];
      }

      for_less( i, 1, n_polygons1 ) {
        poly = list1[i];
        for( dim = 0; dim < N_DIMENSIONS; dim++ ) {
          if( poly_info1[poly].low_limits[dim] < min_range[dim] ) {
            min_range[dim] = poly_info1[poly].low_limits[dim];
          } else {
            if( poly_info1[poly].high_limits[dim] > max_range[dim] ) {
              max_range[dim] = poly_info1[poly].high_limits[dim];
            }
          }
        }
      }

      list2 = &(*list_of_polys1)[offset_index1];
      for_less( i, 0, n_polygons2 ) {
        poly = list2[i];
        for( dim = 0; dim < N_DIMENSIONS; dim++ ) {
          if( poly_info2[poly].low_limits[dim] < min_range[dim] ) {
            min_range[dim] = poly_info2[poly].low_limits[dim];
          } else {
            if( poly_info2[poly].high_limits[dim] > max_range[dim] ) {
              max_range[dim] = poly_info2[poly].high_limits[dim];
            }
          }
        }
      }

      split_axis = 0;

      for_less( dim, 1, N_DIMENSIONS ) {
        if( max_range[dim] - min_range[dim] >
            max_range[split_axis] - min_range[split_axis] ) {
          split_axis = dim;
        }
      }

      if( max_range[split_axis] - min_range[split_axis] >= 2.0*search_distance ) {
        if( next_offset_index1 + 2 * n_polygons1 > *n_alloced1 ) {
            *n_alloced1 = next_offset_index1 + 2 * n_polygons1;
            REALLOC( *list_of_polys1, *n_alloced1 );
            REALLOC( *poly_classes1, *n_alloced1 );
        }

        if( next_offset_index2 + 2 * n_polygons2 > *n_alloced2 ) {
            *n_alloced2 = next_offset_index2 + 2 * n_polygons2;
            REALLOC( *list_of_polys2, *n_alloced2 );
            REALLOC( *poly_classes2, *n_alloced2 );
        }

        // Compute the mid_plane by averaging the centroids
        // of the bounding boxes of the triangles. This gives
        // a much better balance for the number of polygons
        // per branch.
        mid_plane = 0.0;
        for_less( i, 0, n_polygons1 ) {
          poly = list1[i];
          mid_plane += poly_info1[poly].low_limits[split_axis] +
                       poly_info1[poly].high_limits[split_axis];
        }
        for_less( i, 0, n_polygons2 ) {
          poly = list2[i];
          mid_plane += poly_info2[poly].low_limits[split_axis] +
                       poly_info2[poly].high_limits[split_axis];
        }

        mid_plane = 0.5 * mid_plane / ( (float)(n_polygons1+n_polygons2) );
        mid_plane_plus = (float) (mid_plane + search_distance / 2.0);
        mid_plane_minus = (float) (mid_plane - search_distance / 2.0);

        /*--- do everything to left */

        src_list1 = &(*list_of_polys1)[offset_index1];
        src_classes1 = &(*poly_classes1)[offset_index1];

        new_list11 = &(*list_of_polys1)[next_offset_index1];
        new_classes11 = &(*poly_classes1)[next_offset_index1];
        new_n_polygons11 = 0;

        new_list12 = &(*list_of_polys1)[next_offset_index1+n_polygons1];
        new_classes12 = &(*poly_classes1)[next_offset_index1+n_polygons1];
        new_n_polygons12 = 0;

        src_list2 = &(*list_of_polys2)[offset_index2];
        src_classes2 = &(*poly_classes2)[offset_index2];

        new_list21 = &(*list_of_polys2)[next_offset_index2];
        new_classes21 = &(*poly_classes2)[next_offset_index2];
        new_n_polygons21 = 0;

        new_list22 = &(*list_of_polys2)[next_offset_index2+n_polygons2];
        new_classes22 = &(*poly_classes2)[next_offset_index2+n_polygons2];
        new_n_polygons22 = 0;

        class_bit = (1u << split_axis);

        for_less( i, 0, n_polygons1 ) {
            poly = src_list1[i];
            on_left = poly_info1[poly].low_limits[split_axis] <= mid_plane_plus;
            on_right = poly_info1[poly].high_limits[split_axis] >=
                                              mid_plane_minus;
            if( on_left ) {
                new_list11[new_n_polygons11] = poly;
                new_classes11[new_n_polygons11] = src_classes1[i];
                ++new_n_polygons11;
            }
            if( on_right ) {
                new_list12[new_n_polygons12] = poly;
                if( on_left )
                    new_classes12[new_n_polygons12] = (src_classes1[i] | class_bit);
                else
                    new_classes12[new_n_polygons12] = src_classes1[i];
                ++new_n_polygons12;
            }
        }

        diff = n_polygons1 - new_n_polygons11;
        for_less( i, 0, new_n_polygons12 ) {
            new_list12[i-diff] = new_list12[i];
            new_classes12[i-diff] = new_classes12[i];
        }

        for_less( i, 0, n_polygons2 ) {
            poly = src_list2[i];
            on_left = poly_info2[poly].low_limits[split_axis] <= mid_plane_plus;
            on_right = poly_info2[poly].high_limits[split_axis] >=
                                 mid_plane_minus;
            if( on_left ) {
                new_list21[new_n_polygons21] = poly;
                new_classes21[new_n_polygons21] = src_classes2[i];
                ++new_n_polygons21;
            }
            if( on_right ) {
                new_list22[new_n_polygons22] = poly;
                if( on_left )
                    new_classes22[new_n_polygons22] = (src_classes2[i] |
                                                       class_bit);
                else
                    new_classes22[new_n_polygons22] = src_classes2[i];
                ++new_n_polygons22;
            }
        }

        diff = n_polygons2 - new_n_polygons21;
        for_less( i, 0, new_n_polygons22 ) {
            new_list22[i-diff] = new_list22[i];
            new_classes22[i-diff] = new_classes22[i];
        }

        // Avoid creating a branch with too few triangles or with
        // nearly the same number of triangles as its parent. We
        // split if in between 15% and 85% of n_polygons (that's
        // the corresponding value for 0.1275).

        if( ( new_n_polygons11*(n_polygons1-new_n_polygons11) >
              0.1275*n_polygons1*n_polygons1 ) &&
            ( new_n_polygons12*(n_polygons1-new_n_polygons12) >
              0.1275*n_polygons1*n_polygons1 ) &&
            ( new_n_polygons21*(n_polygons2-new_n_polygons21) >
              0.1275*n_polygons2*n_polygons2 ) &&
            ( new_n_polygons22*(n_polygons2-new_n_polygons22) >
              0.1275*n_polygons2*n_polygons2 ) ) {

            recursive_find_close_pairs( min_distance, search_distance, max_step_size,
                        parameters1, poly_info1, list_of_polys1, poly_classes1,
                        n_alloced1, next_offset_index1, new_n_polygons11,
                        next_offset_index1 + new_n_polygons11 +new_n_polygons12,
                        line_dir1, parameters2,
                        poly_info2, list_of_polys2, poly_classes2,
                        n_alloced2, next_offset_index2, new_n_polygons21,
                        next_offset_index2 + new_n_polygons21 +new_n_polygons22,
                        line_dir2, n_pairs_ptr, n_pairs_alloc, 
                        p1s_ptr, p2s_ptr, cases_ptr, min_line_dists_ptr,
                        closest_distance );

           recursive_find_close_pairs( min_distance, search_distance,
                        max_step_size, parameters1,
                        poly_info1, list_of_polys1, poly_classes1,
                        n_alloced1, next_offset_index1 + new_n_polygons11,
                        new_n_polygons12,
                        next_offset_index1 + new_n_polygons11 +new_n_polygons12,
                        line_dir1,
                        parameters2,
                        poly_info2, list_of_polys2, poly_classes2,
                        n_alloced2, next_offset_index2 + new_n_polygons21,
                        new_n_polygons22,
                        next_offset_index2 + new_n_polygons21 +new_n_polygons22,
                        line_dir2, n_pairs_ptr, n_pairs_alloc,
                        p1s_ptr, p2s_ptr, cases_ptr, min_line_dists_ptr,
                        closest_distance );

            return;
        }
      }
    }

    min_distance_sq = min_distance * min_distance;
    search_distance_sq = search_distance * search_distance;
    fl_search_distance = (float) search_distance;

    list1 = &(*list_of_polys1)[offset_index1];
    classes1 = &(*poly_classes1)[offset_index1];

    if( n_polygons1 > MAX_THRESHOLD_N_POLYGONS ) {
        ALLOC( sorted_list1, n_polygons1 );
    } else {
        sorted_list1 = static_list1;
    }
        
    for_less( i1, 0, 16 )
        count1[i1] = 0;

    for_less( i1, 0, n_polygons1 )
        ++count1[classes1[i1]];

    for_less( i1, 1, 16 )
        count1[i1] += count1[i1-1];

    for_less( i1, 0, n_polygons1 ) {
        cl = classes1[i1];
        --count1[cl];
        sorted_list1[count1[cl]] = list1[i1];
    }

    count1[16] = n_polygons1;

    list2 = &(*list_of_polys2)[offset_index2];
    classes2 = &(*poly_classes2)[offset_index2];

    if( n_polygons2 > MAX_THRESHOLD_N_POLYGONS ) {
        ALLOC( sorted_list2, n_polygons2 );
        ALLOC( compare2, n_polygons2 );
        ALLOC( compare2_index, n_polygons2 );
    } else {
        sorted_list2 = static_list2;
        compare2 = static_compare2;
        compare2_index = static_compare2_index;
    }
        
    for_less( i1, 0, 16 )
        count2[i1] = 0;

    for_less( i1, 0, n_polygons2 )
        ++count2[classes2[i1]];

    for_less( i1, 1, 16 )
        count2[i1] += count2[i1-1];

    for_less( i1, 0, n_polygons2 ) {
        cl = classes2[i1];
        --count2[cl];
        sorted_list2[count2[cl]] = list2[i1];
    }

    count2[16] = n_polygons2;

    for_less( c1, 0, 16 ) {
      if( count1[c1] == count1[c1+1] ){
        continue;
      }

        n_to_compare2 = 0;

        for_less( c2, 0, 16 ) {
            if( (c1 & c2) != 0 || count2[c2] == count2[c2+1] )
                continue;

            for_less( i2, count2[c2], count2[c2+1] ) {
                compare2[n_to_compare2] = &poly_info2[sorted_list2[i2]];
                compare2_index[n_to_compare2] = sorted_list2[i2];
                ++n_to_compare2;
            }
        }

        for_less( i1, count1[c1], count1[c1+1] ) {
            poly1 = sorted_list1[i1];
            p1 = poly_info1[poly1].p1;
            n11 = poly_info1[poly1].n11;
            n12 = poly_info1[poly1].n12;

            parm_p1 = IJ( p1, 0, 3 );
            parm_n11 = IJ( n11, 0, 3 );
            parm_n12 = IJ( n12, 0, 3 );

            x1_low = poly_info1[poly1].low_limits[X] - fl_search_distance;
            x1_high = poly_info1[poly1].high_limits[X] + fl_search_distance;
            y1_low = poly_info1[poly1].low_limits[Y] - fl_search_distance;
            y1_high = poly_info1[poly1].high_limits[Y] + fl_search_distance;
            z1_low = poly_info1[poly1].low_limits[Z] - fl_search_distance;
            z1_high = poly_info1[poly1].high_limits[Z] + fl_search_distance;

            for_less( i2, 0, n_to_compare2 ) {
                info_ptr = compare2[i2];

                if( info_ptr->low_limits[X] > x1_high ||
                    info_ptr->high_limits[X] < x1_low ||
                    info_ptr->low_limits[Y] > y1_high ||
                    info_ptr->high_limits[Y] < y1_low ||
                    info_ptr->low_limits[Z] > z1_high ||
                    info_ptr->high_limits[Z] < z1_low ) {
                    continue;
                }

                p2 = info_ptr->p1;
                n21 = info_ptr->n11;
                n22 = info_ptr->n12;
                parm_p2 = 3 * p2;
                parm_n21 = 3 * n21;
                parm_n22 = 3 * n22;

                if( sq_triangle_triangle_dist_estimate(
                                                &parameters1[parm_p1],
                                                &parameters1[parm_n11],
                                                &parameters1[parm_n12],
                                                &parameters2[parm_p2],
                                                &parameters2[parm_n21],
                                                &parameters2[parm_n22],
                                                search_distance_sq ) ) {
                    continue;
                }

                which_case = 0;
                dist = sq_triangle_triangle_dist( &parameters1[parm_p1],
                                             &parameters1[parm_n11],
                                             &parameters1[parm_n12],
                                             &parameters2[parm_p2],
                                             &parameters2[parm_n21],
                                             &parameters2[parm_n22],
                                             &which_case );

                if( dist >= search_distance_sq ) {
                    continue;
                }

#ifdef PRINT_DIST
sum_dist += dist;
#endif

                if( *closest_distance < 0.0 || dist < *closest_distance )
                    *closest_distance = dist;

                if( dist <= min_distance_sq || line_dir1 == NULL ) {
                    t_min = 0.0;
                } else {
                    parms1[0] = parm_p1;
                    parms1[1] = parm_n11;
                    parms1[2] = parm_n12;
                    parms2[0] = parm_p2;
                    parms2[1] = parm_n21;
                    parms2[2] = parm_n22;

                    movement_sq = 0.0;
                    for_less( which1, 0, 3 ) {
                      for_less( which2, 0, 3 ) {
                        dx = line_dir1[parms1[which1]+0] -
                             line_dir2[parms2[which2]+0];
                        dy = line_dir1[parms1[which1]+1] -
                             line_dir2[parms2[which2]+1];
                        dz = line_dir1[parms1[which1]+2] -
                             line_dir2[parms2[which2]+2];
                        d = dx * dx + dy * dy + dz * dz;

                        if( which1 == 0 && which2 == 0 ||
                            d > movement_sq )
                            movement_sq = d;
                      }
                    }

                    if( movement_sq == 0.0 ) {
                        t_min = 1.0e30;
                    } else {
                        movement = sqrt( movement_sq );
                        dist = sqrt( dist );
                        t_min = (dist - min_distance) / movement;
                        *max_step_size = MIN( *max_step_size, dist/movement );
                    }
                }

                if( (*n_pairs_ptr) >= (*n_pairs_alloc) ) {
                  // Re-allocate more space and copy. Add 10% more.
                  *n_pairs_alloc = (int)( 1.1 * (float)(*n_pairs_alloc) );
                  REALLOC( *p1s_ptr, *n_pairs_alloc );
                  REALLOC( *p2s_ptr, *n_pairs_alloc );
                  REALLOC( *cases_ptr, *n_pairs_alloc );
                  REALLOC( *min_line_dists_ptr, *n_pairs_alloc );
                }

                (*min_line_dists_ptr)[*n_pairs_ptr] = (float)
                                                 (t_min * 0.999);
                (*p1s_ptr)[*n_pairs_ptr] = poly1;
                (*p2s_ptr)[*n_pairs_ptr] = compare2_index[i2];
                (*cases_ptr)[*n_pairs_ptr] = which_case;
                ++(*n_pairs_ptr);
            }
        }
    }

    if( n_polygons1 > MAX_THRESHOLD_N_POLYGONS ) {
        FREE( sorted_list1 );
    }

    if( n_polygons2 > MAX_THRESHOLD_N_POLYGONS ) {
        FREE( sorted_list2 );
        FREE( compare2 );
        FREE( compare2_index );
    }
}


private  void   find_surf_intersect_candidates(
    Real               min_distance,
    int                *n_pairs_ptr,
    int                *n_pairs_alloc,
    int                *p1s_ptr[],
    int                *p2s_ptr[],
    unsigned char      *cases_ptr[],
    float              *min_line_dists_ptr[],
    int                n_points1,
    Real               parameters1[],
    Real               line_dir1[],
    int                n_triangles1,
    int                triangles1[],
    int                n_points2,
    Real               parameters2[],
    Real               line_dir2[],
    int                n_triangles2,
    int                triangles2[],
    Real               max_movement,
    Real               *max_step_size,
    Real               *closest_distance )
{
    poly_info_struct  *poly_info1, *poly_info2;
    int               poly1, poly2, n_polygons1;
    int               n_polygons2;
    int               p1, dim, p, n;
    int               n11, n12;
    int               *list_of_polys1, n_alloced1;
    unsigned int      *poly_classes1;
    int               *list_of_polys2, n_alloced2;
    unsigned int      *poly_classes2;
    Real              search_distance;

#ifdef PRINT_DIST
sum_dist = 0.0;
#endif

    search_distance = min_distance + max_movement;

    /*--- prepare surface 1 */

    n_polygons1 = n_triangles1;

    // Store the triangles in poly_info in the same order
    // as they were read. Also, find the min/max bounding
    // box for each triangle.

    ALLOC( poly_info1, n_polygons1 );

    for( p = 0; p < n_polygons1; p++ ) {
      n = 0;
      if( triangles1[3*p+1] < triangles1[3*p+n] ) n = 1;
      if( triangles1[3*p+2] < triangles1[3*p+n] ) n = 2;
      p1 = triangles1[3*p+n];
      n11 = triangles1[3*p+(n+1)%3];
      n12 = triangles1[3*p+(n+2)%3];
      poly_info1[p].p1 = p1;           // want p1 < n11 && p1 < n12
      poly_info1[p].n11 = n11;
      poly_info1[p].n12 = n12;
      for_less( dim, 0, N_DIMENSIONS ) {
        poly_info1[p].low_limits[dim] = (float) MIN3( parameters1[IJ(p1,dim,3)],
                                                      parameters1[IJ(n11,dim,3)],
                                                      parameters1[IJ(n12,dim,3)] );
        poly_info1[p].high_limits[dim] = (float) MAX3( parameters1[IJ(p1,dim,3)],
                                                       parameters1[IJ(n11,dim,3)],
                                                       parameters1[IJ(n12,dim,3)] );
      }
    }

    // Set-up arrays for the partitionning of the domain.

    n_alloced1 = 3 * n_polygons1;
    ALLOC( list_of_polys1, n_alloced1 );
    ALLOC( poly_classes1, n_alloced1 );

    for_less( poly1, 0, n_polygons1 ) {
        list_of_polys1[poly1] = poly1;
        p1 = poly_info1[poly1].p1;
        n11 = poly_info1[poly1].n11;
        n12 = poly_info1[poly1].n12;

        if( line_dir1[IJ(p1,0,3)] == 0.0 &&
            line_dir1[IJ(p1,1,3)] == 0.0 &&
            line_dir1[IJ(p1,2,3)] == 0.0 &&
            line_dir1[IJ(n11,0,3)] == 0.0 &&
            line_dir1[IJ(n11,1,3)] == 0.0 &&
            line_dir1[IJ(n11,2,3)] == 0.0 &&
            line_dir1[IJ(n12,0,3)] == 0.0 &&
            line_dir1[IJ(n12,1,3)] == 0.0 &&
            line_dir1[IJ(n12,2,3)] == 0.0 ) {
            poly_classes1[poly1] = 8;
        } else {
            poly_classes1[poly1] = 0;
        }
    }

    /*--- prepare surface 2 */

    n_polygons2 = n_triangles2;

    // Store the triangles in poly_info in the same order
    // as they were read. Also, find the min/max bounding
    // box for each triangle.

    ALLOC( poly_info2, n_polygons2 );

    for( p = 0; p < n_polygons2; p++ ) {
      n = 0;
      if( triangles2[3*p+1] < triangles2[3*p+n] ) n = 1;
      if( triangles2[3*p+2] < triangles2[3*p+n] ) n = 2;
      p1 = triangles2[3*p+n];
      n11 = triangles2[3*p+(n+1)%3];
      n12 = triangles2[3*p+(n+2)%3];
      poly_info2[p].p1 = p1;           // want p1 < n11 && p1 < n12
      poly_info2[p].n11 = n11;
      poly_info2[p].n12 = n12;
      for_less( dim, 0, N_DIMENSIONS ) {
        poly_info2[p].low_limits[dim] = (float) MIN3( parameters2[IJ(p1,dim,3)],
                                                      parameters2[IJ(n11,dim,3)],
                                                      parameters2[IJ(n12,dim,3)] );
        poly_info2[p].high_limits[dim] = (float) MAX3( parameters2[IJ(p1,dim,3)],
                                                       parameters2[IJ(n11,dim,3)],
                                                       parameters2[IJ(n12,dim,3)] );
      }
    }

    n_alloced2 = 3 * n_polygons2;
    ALLOC( list_of_polys2, n_alloced2 );
    ALLOC( poly_classes2, n_alloced2 );

    for_less( poly2, 0, n_polygons2 ) {
        list_of_polys2[poly2] = poly2;
        p1 = poly_info2[poly2].p1;
        n11 = poly_info2[poly2].n11;
        n12 = poly_info2[poly2].n12;

        if( line_dir2[IJ(p1,0,3)] == 0.0 &&
            line_dir2[IJ(p1,1,3)] == 0.0 &&
            line_dir2[IJ(p1,2,3)] == 0.0 &&
            line_dir2[IJ(n11,0,3)] == 0.0 &&
            line_dir2[IJ(n11,1,3)] == 0.0 &&
            line_dir2[IJ(n11,2,3)] == 0.0 &&
            line_dir2[IJ(n12,0,3)] == 0.0 &&
            line_dir2[IJ(n12,1,3)] == 0.0 &&
            line_dir2[IJ(n12,2,3)] == 0.0 ) {
            poly_classes2[poly2] = 8;
        } else {
            poly_classes2[poly2] = 0;
        }
    }

    /*--- now set up list */

    *closest_distance = -1.0;
    *n_pairs_ptr = 0;

    recursive_find_close_pairs( min_distance, search_distance, max_step_size,
                        parameters1, poly_info1,
                        &list_of_polys1, &poly_classes1,
                        &n_alloced1, 0, n_polygons1, n_polygons1,
                        line_dir1,
                        parameters2, poly_info2,
                        &list_of_polys2, &poly_classes2,
                        &n_alloced2, 0, n_polygons2, n_polygons2,
                        line_dir2, n_pairs_ptr, n_pairs_alloc, 
                        p1s_ptr, p2s_ptr, cases_ptr, min_line_dists_ptr,
                        closest_distance );

    if( n_alloced1 > 0 ) {
        FREE( list_of_polys1 );
        FREE( poly_classes1 );
    }

    if( n_alloced2 > 0 ) {
        FREE( list_of_polys2 );
        FREE( poly_classes2 );
    }

    FREE( poly_info1 );
    FREE( poly_info2 );

    if( *closest_distance > 0.0 )
        *closest_distance = sqrt( *closest_distance );

#ifdef PRINT_DIST
print( "Sum dist: %.15g\n", sum_dist );
#endif
}

public  Real  recompute_surf_surfs(
    Deform_struct                 *deform,
    int                           start_parameter[],
    Real                          parameters[],
    Real                          max_movement,
    Real                          line_dir[],
    Real                          *max_step_size,
    surf_surf_lookup_struct       *ss_lookup )
{
    int                    i, w, s1, s2;
    Real                   *this_parms1, *this_parms2, *this_line1, *this_line2;
    Real                   closest_distance;
    surf_surf_struct       *surf;
    Real                   min_distance, close;

    closest_distance = -1.0;
    *max_step_size = 1.0e30;

    for_less( i, 0, deform->n_surf_surfs ) {
        surf = &deform->surf_surfs[i];
        s1 = surf->surface_index1;
        s2 = surf->surface_index2;

        this_parms1 = &parameters[start_parameter[s1]];
        this_parms2 = &parameters[start_parameter[s2]];

        if( line_dir == NULL ) {
            this_line1 = NULL;
            this_line2 = NULL;
        } else {
            this_line1 = &line_dir[start_parameter[s1]];
            this_line2 = &line_dir[start_parameter[s2]];
        }

        min_distance = 0.0;
        for_less( w, 0, surf->n_weights ) {
            min_distance = MAX( min_distance, surf->min_distances[w] );
        }

        find_surf_intersect_candidates(
                       min_distance,
                       &ss_lookup[i].n_pairs,
                       &ss_lookup[i].n_pairs_alloc,
                       &ss_lookup[i].p1s,
                       &ss_lookup[i].p2s,
                       &ss_lookup[i].cases,
                       &ss_lookup[i].min_line_dists,
                       deform->surfaces[s1].surface.n_points,
                       this_parms1, this_line1,
                       deform->surfaces[s1].surface.n_polygons,
                       deform->surfaces[s1].surface.triangles,
                       deform->surfaces[s2].surface.n_points,
                       this_parms2, this_line2,
                       deform->surfaces[s2].surface.n_polygons,
                       deform->surfaces[s2].surface.triangles,
                       max_movement, max_step_size, &close );

        if( close >= 0.0 &&
            (closest_distance < 0.0 || close < closest_distance) )
            closest_distance = close;
    }

    return( closest_distance );
}

public  void  delete_surf_surf_lookup(
    surf_surf_lookup_struct  *ss_lookup ) {

    if( ss_lookup->n_pairs > 0 ) {
        FREE( ss_lookup->p1s );
        FREE( ss_lookup->p2s );
        FREE( ss_lookup->cases );
        FREE( ss_lookup->min_line_dists );
        ss_lookup->n_pairs = 0;
        ss_lookup->n_pairs_alloc = 0;
    }
}

public  int  get_n_surf_surf_candidate(
    surf_surf_lookup_struct  *ss_lookup )
{
    return( ss_lookup->n_pairs );
}

public  Real   get_surf_surf_deriv_factor( Real min_distance,
                                           Real weight,
                                           Real dist ) {

    Real  diff, factor;

    if( dist == 0.0 ) {
        factor = 0.0;
    } else {
        diff = min_distance - dist;
#ifdef USE_CUBE
#define USE_CUBE
        factor = -1.5 * diff * diff * weight / ( dist * min_distance );
#else
        factor = -diff * weight / dist;
#endif
    }
    return( factor );
}

public  void   get_surf_surf_deriv(
    int                      p1,
    int                      n11,
    int                      n12,
    int                      p2,
    int                      n21,
    int                      n22,
    Real                     factor,
    Real                     deriv1[],
    Real                     deriv2[],
    Real                     tri_tri_deriv[6][3] ) {

    deriv1[IJ(p1,0,3)] += factor * tri_tri_deriv[0][0];
    deriv1[IJ(p1,1,3)] += factor * tri_tri_deriv[0][1];
    deriv1[IJ(p1,2,3)] += factor * tri_tri_deriv[0][2];

    deriv1[IJ(n11,0,3)] += factor * tri_tri_deriv[1][0];
    deriv1[IJ(n11,1,3)] += factor * tri_tri_deriv[1][1];
    deriv1[IJ(n11,2,3)] += factor * tri_tri_deriv[1][2];

    deriv1[IJ(n12,0,3)] += factor * tri_tri_deriv[2][0];
    deriv1[IJ(n12,1,3)] += factor * tri_tri_deriv[2][1];
    deriv1[IJ(n12,2,3)] += factor * tri_tri_deriv[2][2];

    deriv2[IJ(p2,0,3)] += factor * tri_tri_deriv[3][0];
    deriv2[IJ(p2,1,3)] += factor * tri_tri_deriv[3][1];
    deriv2[IJ(p2,2,3)] += factor * tri_tri_deriv[3][2];

    deriv2[IJ(n21,0,3)] += factor * tri_tri_deriv[4][0];
    deriv2[IJ(n21,1,3)] += factor * tri_tri_deriv[4][1];
    deriv2[IJ(n21,2,3)] += factor * tri_tri_deriv[4][2];

    deriv2[IJ(n22,0,3)] += factor * tri_tri_deriv[5][0];
    deriv2[IJ(n22,1,3)] += factor * tri_tri_deriv[5][1];
    deriv2[IJ(n22,2,3)] += factor * tri_tri_deriv[5][2];
}
