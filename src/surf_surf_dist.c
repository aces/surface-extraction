#include  <fit_3d.h>

#define  MAX_THRESHOLD_N_POLYGONS   200

#define  THRESHOLD_N_POLYGONS   65

#ifdef PRINT_DIST
#define PRINT_DIST
static Real sum_dist;
#endif

typedef  struct
{
    float   low_limits[N_DIMENSIONS];
    float   high_limits[N_DIMENSIONS];
    int     p1;
    int     n11;
    int     n12;
} poly_info_struct;

private  void    recursive_find_close_pairs( 
    Real               min_range[N_DIMENSIONS],
    Real               max_range[N_DIMENSIONS],
    Real               min_distance,
    Real               search_distance,
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
    int                *p1s_ptr[],
    int                *p2s_ptr[],
    unsigned char      *cases_ptr[],
    float              *min_line_dists_ptr[],
    int                n_points1,
    int                n_neighbours1[],
    int                *neighbours1[],
    int                n_points2,
    int                n_neighbours2[],
    int                *neighbours2[],
    Real               *closest_distance )
{
    int               parms1[3], parms2[3], which_case;
    Real              mid_plane, min_distance_sq, search_distance_sq;
    Real              movement, movement_sq, dx, dy, dz, d, t_min, dist;
    int               which1, which2, n1, n2;
    float             mid_plane_plus, mid_plane_minus;
    Real              new_min[N_DIMENSIONS], new_max[N_DIMENSIONS];
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
    poly_info_struct  **compare2, *info_ptr;
    static int        threshold_n_polygons = THRESHOLD_N_POLYGONS;
    static BOOLEAN    first = TRUE;

    if( first )
    {
        first = FALSE;
        if( getenv( "THRESHOLD_N_POLYGONS" ) == NULL ||
            sscanf( getenv("THRESHOLD_N_POLYGONS" ), "%d",
                    &threshold_n_polygons ) != 1  ||
            threshold_n_polygons > MAX_THRESHOLD_N_POLYGONS )
            threshold_n_polygons = THRESHOLD_N_POLYGONS;
    }

    if( n_polygons1 == 0 || n_polygons2 == 0 )
        return;

    split_axis = -1;

    for_less( dim, 0, N_DIMENSIONS )
    {
        if( split_axis < 0 ||
            max_range[dim] - min_range[dim] >
            max_range[split_axis] - min_range[split_axis] )
            split_axis = dim;
    }

    if( split_axis >= 0 &&
        max_range[split_axis] - min_range[split_axis] >= 2.0*search_distance &&
        (n_polygons1 >= threshold_n_polygons ||
         n_polygons2 >= threshold_n_polygons) )
    {
        if( next_offset_index1 + 2 * n_polygons1 > *n_alloced1 )
        {
            *n_alloced1 = next_offset_index1 + 2 * n_polygons1;
            REALLOC( *list_of_polys1, *n_alloced1 );
            REALLOC( *poly_classes1, *n_alloced1 );
        }

        if( next_offset_index2 + 2 * n_polygons2 > *n_alloced2 )
        {
            *n_alloced2 = next_offset_index2 + 2 * n_polygons2;
            REALLOC( *list_of_polys2, *n_alloced2 );
            REALLOC( *poly_classes2, *n_alloced2 );
        }

        mid_plane = (max_range[split_axis] + min_range[split_axis]) / 2.0;
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

        for_less( i, 0, n_polygons1 )
        {
            poly = src_list1[i];
            on_left = poly_info1[poly].low_limits[split_axis] <= mid_plane_plus;
            on_right = poly_info1[poly].high_limits[split_axis] >=
                                              mid_plane_minus;
            if( on_left )
            {
                new_list11[new_n_polygons11] = poly;
                new_classes11[new_n_polygons11] = src_classes1[i];
                ++new_n_polygons11;
            }
            if( on_right )
            {
                new_list12[new_n_polygons12] = poly;
                if( on_left )
                    new_classes12[new_n_polygons12] = (src_classes1[i] | class_bit);
                else
                    new_classes12[new_n_polygons12] = src_classes1[i];
                ++new_n_polygons12;
            }
        }

        diff = n_polygons1 - new_n_polygons11;
        for_less( i, 0, new_n_polygons12 )
        {
            new_list12[i-diff] = new_list12[i];
            new_classes12[i-diff] = new_classes12[i];
        }

        for_less( i, 0, n_polygons2 )
        {
            poly = src_list2[i];
            on_left = poly_info2[poly].low_limits[split_axis] <= mid_plane_plus;
            on_right = poly_info2[poly].high_limits[split_axis] >=
                                 mid_plane_minus;
            if( on_left )
            {
                new_list21[new_n_polygons21] = poly;
                new_classes21[new_n_polygons21] = src_classes2[i];
                ++new_n_polygons21;
            }
            if( on_right )
            {
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
        for_less( i, 0, new_n_polygons22 )
        {
            new_list22[i-diff] = new_list22[i];
            new_classes22[i-diff] = new_classes22[i];
        }

        if( (new_n_polygons11 < n_polygons1 || new_n_polygons21 < n_polygons2)&&
            (new_n_polygons12 < n_polygons1 || new_n_polygons22 < n_polygons2) )
        {
            new_max[0] = max_range[0];
            new_max[1] = max_range[1];
            new_max[2] = max_range[2];
            new_max[split_axis] = (Real) mid_plane_plus;

            recursive_find_close_pairs( min_range, new_max,
                        min_distance, search_distance,
                        parameters1,
                        poly_info1, list_of_polys1, poly_classes1,
                        n_alloced1, next_offset_index1, new_n_polygons11,
                        next_offset_index1 + new_n_polygons11 +new_n_polygons12,
                        line_dir1,
                        parameters2,
                        poly_info2, list_of_polys2, poly_classes2,
                        n_alloced2, next_offset_index2, new_n_polygons21,
                        next_offset_index2 + new_n_polygons21 +new_n_polygons22,
                        line_dir2,
                        n_pairs_ptr, p1s_ptr, p2s_ptr, cases_ptr,
                        min_line_dists_ptr,
                        n_points1, n_neighbours1, neighbours1,
                        n_points2, n_neighbours2, neighbours2,
                        closest_distance );

           new_min[0] = min_range[0];
           new_min[1] = min_range[1];
           new_min[2] = min_range[2];
           new_min[split_axis] = (Real) mid_plane_minus;

           recursive_find_close_pairs( new_min, max_range,
                        min_distance, search_distance,
                        parameters1,
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
                        line_dir2,
                        n_pairs_ptr, p1s_ptr, p2s_ptr, cases_ptr,
                        min_line_dists_ptr,
                        n_points1, n_neighbours1, neighbours1,
                        n_points2, n_neighbours2, neighbours2,
                        closest_distance );

            return;
        }
    }

    min_distance_sq = min_distance * min_distance;
    search_distance_sq = search_distance * search_distance;
    fl_search_distance = (float) search_distance;

    list1 = &(*list_of_polys1)[offset_index1];
    classes1 = &(*poly_classes1)[offset_index1];

    if( n_polygons1 > threshold_n_polygons )
    {
        ALLOC( sorted_list1, n_polygons1 );
    }
    else
    {
        sorted_list1 = static_list1;
    }
        
    for_less( i1, 0, 16 )
        count1[i1] = 0;

    for_less( i1, 0, n_polygons1 )
        ++count1[classes1[i1]];

    for_less( i1, 1, 16 )
        count1[i1] += count1[i1-1];

    for_less( i1, 0, n_polygons1 )
    {
        cl = classes1[i1];
        --count1[cl];
        sorted_list1[count1[cl]] = list1[i1];
    }

    count1[16] = n_polygons1;

    list2 = &(*list_of_polys2)[offset_index2];
    classes2 = &(*poly_classes2)[offset_index2];

    if( n_polygons2 > threshold_n_polygons )
    {
        ALLOC( sorted_list2, n_polygons2 );
        ALLOC( compare2, n_polygons2 );
    }
    else
    {
        sorted_list2 = static_list2;
        compare2 = static_compare2;
    }
        
    for_less( i1, 0, 16 )
        count2[i1] = 0;

    for_less( i1, 0, n_polygons2 )
        ++count2[classes2[i1]];

    for_less( i1, 1, 16 )
        count2[i1] += count2[i1-1];

    for_less( i1, 0, n_polygons2 )
    {
        cl = classes2[i1];
        --count2[cl];
        sorted_list2[count2[cl]] = list2[i1];
    }

    count2[16] = n_polygons2;

    for_less( c1, 0, 16 )
    {
      if( count1[c1] == count1[c1+1] ){
        continue;
      }

        n_to_compare2 = 0;

        for_less( c2, 0, 16 )
        {
            if( (c1 & c2) != 0 || count2[c2] == count2[c2+1] )
                continue;

            for_less( i2, count2[c2], count2[c2+1] )
            {
                compare2[n_to_compare2] = &poly_info2[sorted_list2[i2]];
                ++n_to_compare2;
            }
        }

        for_less( i1, count1[c1], count1[c1+1] )
        {
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

            for_less( i2, 0, n_to_compare2 )
            {
                info_ptr = compare2[i2];

                if( info_ptr->low_limits[X] > x1_high ||
                    info_ptr->high_limits[X] < x1_low ||
                    info_ptr->low_limits[Y] > y1_high ||
                    info_ptr->high_limits[Y] < y1_low ||
                    info_ptr->low_limits[Z] > z1_high ||
                    info_ptr->high_limits[Z] < z1_low )
                {
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
                                                search_distance_sq ) >=
                                                search_distance_sq )
                {
                    continue;
                }

                which_case = 0;
                dist = sq_triangle_triangle_dist( &parameters1[parm_p1],
                                             &parameters1[parm_n11],
                                             &parameters1[parm_n12],
                                             &parameters2[parm_p2],
                                             &parameters2[parm_n21],
                                             &parameters2[parm_n22],
                                             &which_case,0 );

                if( dist >= search_distance_sq )
                {
                    continue;
                }

#ifdef PRINT_DIST
sum_dist += dist;
#endif

                if( *closest_distance < 0.0 || dist < *closest_distance )
                    *closest_distance = dist;

                if( dist <= min_distance_sq || line_dir1 == NULL )
                    t_min = 0.0;
                else
                {
                    parms1[0] = parm_p1;
                    parms1[1] = parm_n11;
                    parms1[2] = parm_n12;
                    parms2[0] = parm_p2;
                    parms2[1] = parm_n21;
                    parms2[2] = parm_n22;

                    movement_sq = 0.0;
                    for_less( which1, 0, 3 )
                    for_less( which2, 0, 3 )
                    {
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

                    if( movement_sq == 0.0 )
                        t_min = 1.0e30;
                    else
                    {
                        movement = sqrt( movement_sq );
                        dist = sqrt( dist );
                        t_min = (dist - min_distance) / movement;
                    }
                }

                for_less( n1, 0, n_neighbours1[p1] )
                {
                    if( neighbours1[p1][n1] == n11 &&
                        neighbours1[p1][(n1+1)%n_neighbours1[p1]] == n12 )
                        break;
                }
                if( n1 >= n_neighbours1[p1] )
                    handle_internal_error( "n1 >= n_neighbours1[p1]" );

                for_less( n2, 0, n_neighbours2[p2] )
                {
                    if( neighbours2[p2][n2] == n21 &&
                        neighbours2[p2][(n2+1)%n_neighbours2[p2]] == n22 )
                        break;
                }
                if( n2 >= n_neighbours2[p2] )
                    handle_internal_error( "n2 >= n_neighbours2[p2]" );

                SET_ARRAY_SIZE( *p1s_ptr, (*n_pairs_ptr),
                                (*n_pairs_ptr) + 1,
                                DEFAULT_CHUNK_SIZE );
                SET_ARRAY_SIZE( *p2s_ptr, (*n_pairs_ptr),
                                (*n_pairs_ptr) + 1,
                                DEFAULT_CHUNK_SIZE );
                SET_ARRAY_SIZE( *cases_ptr, (*n_pairs_ptr),
                                (*n_pairs_ptr) + 1,
                                DEFAULT_CHUNK_SIZE );
                SET_ARRAY_SIZE( *min_line_dists_ptr, (*n_pairs_ptr),
                                (*n_pairs_ptr) + 1, DEFAULT_CHUNK_SIZE);

                (*min_line_dists_ptr)[*n_pairs_ptr] = (float)
                                                 (t_min * 0.999);
                (*p1s_ptr)[*n_pairs_ptr] = IJ( n1, p1, n_points1 );
                (*p2s_ptr)[*n_pairs_ptr] = IJ( n2, p2, n_points2 );
                (*cases_ptr)[*n_pairs_ptr] = (unsigned char) which_case;
                ++(*n_pairs_ptr);
            }
        }
    }

    if( n_polygons1 > threshold_n_polygons )
    {
        FREE( sorted_list1 );
    }

    if( n_polygons2 > threshold_n_polygons )
    {
        FREE( sorted_list2 );
        FREE( compare2 );
    }
}

private  void   find_surf_intersect_candidates(
    Real               min_distance,
    int                *n_pairs_ptr,
    int                *p1s_ptr[],
    int                *p2s_ptr[],
    unsigned char      *cases_ptr[],
    float              *min_line_dists_ptr[],
    int                n_points1,
    Real               parameters1[],
    Real               line_dir1[],
    int                n_neighbours1[],
    int                *neighbours1[],
    int                n_points2,
    Real               parameters2[],
    Real               line_dir2[],
    int                n_neighbours2[],
    int                *neighbours2[],
    Real               max_movement,
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
    Real              min_range[N_DIMENSIONS], max_range[N_DIMENSIONS];
    Real              search_distance;

#ifdef PRINT_DIST
sum_dist = 0.0;
#endif

    search_distance = min_distance + max_movement;

    /*--- prepare surface 1 */

    n_polygons1 = 0;
    for_less( p, 0, n_points1 )
    {
        for_less( n, 0, n_neighbours1[p] )
        {
            n11 = neighbours1[p][n];
            n12 = neighbours1[p][(n+1) % n_neighbours1[p]];
            if( n11 > p && n12 > p )
                ++n_polygons1;
        }
    }

    ALLOC( poly_info1, n_polygons1 );

    poly1 = 0;
    for_less( p, 0, n_points1 )
    {
        for_less( n, 0, n_neighbours1[p] )
        {
            n11 = neighbours1[p][n];
            n12 = neighbours1[p][(n+1) % n_neighbours1[p]];
            if( n11 > p && n12 > p )
            {
                poly_info1[poly1].p1 = p;
                poly_info1[poly1].n11 = n11;
                poly_info1[poly1].n12 = n12;
                for_less( dim, 0, N_DIMENSIONS )
                {
                    poly_info1[poly1].low_limits[dim] = (float) MIN3(
                                                   parameters1[IJ(p,dim,3)],
                                                   parameters1[IJ(n11,dim,3)],
                                                   parameters1[IJ(n12,dim,3)] );
                    poly_info1[poly1].high_limits[dim] = (float) MAX3(
                                                   parameters1[IJ(p,dim,3)],
                                                   parameters1[IJ(n11,dim,3)],
                                                   parameters1[IJ(n12,dim,3)] );
                }

                ++poly1;
            }
        }
    }

    n_alloced1 = 3 * n_polygons1;
    ALLOC( list_of_polys1, n_alloced1 );
    ALLOC( poly_classes1, n_alloced1 );

    for_less( poly1, 0, n_polygons1 )
    {
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
            line_dir1[IJ(n12,2,3)] == 0.0 )
        {
            poly_classes1[poly1] = 8;
        }
        else
        {
            poly_classes1[poly1] = 0;
        }
    }

    /*--- prepare surface 2 */

    n_polygons2 = 0;
    for_less( p, 0, n_points2 )
    {
        for_less( n, 0, n_neighbours2[p] )
        {
            n11 = neighbours2[p][n];
            n12 = neighbours2[p][(n+1) % n_neighbours2[p]];
            if( n11 > p && n12 > p )
                ++n_polygons2;
        }
    }

    ALLOC( poly_info2, n_polygons2 );
    poly2 = 0;
    for_less( p, 0, n_points2 )
    {
        for_less( n, 0, n_neighbours2[p] )
        {
            n11 = neighbours2[p][n];
            n12 = neighbours2[p][(n+1) % n_neighbours2[p]];
            if( n11 > p && n12 > p )
            {
                poly_info2[poly2].p1 = p;
                poly_info2[poly2].n11 = n11;
                poly_info2[poly2].n12 = n12;

                for_less( dim, 0, N_DIMENSIONS )
                {
                    poly_info2[poly2].low_limits[dim] = (float) MIN3(
                                                   parameters2[IJ(p,dim,3)],
                                                   parameters2[IJ(n11,dim,3)],
                                                   parameters2[IJ(n12,dim,3)] );
                    poly_info2[poly2].high_limits[dim] = (float) MAX3(
                                                   parameters2[IJ(p,dim,3)],
                                                   parameters2[IJ(n11,dim,3)],
                                                   parameters2[IJ(n12,dim,3)] );
                }

                ++poly2;
            }
        }
    }

    n_alloced2 = 3 * n_polygons2;
    ALLOC( list_of_polys2, n_alloced2 );
    ALLOC( poly_classes2, n_alloced2 );

    for_less( poly2, 0, n_polygons2 )
    {
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
            line_dir2[IJ(n12,2,3)] == 0.0 )
        {
            poly_classes2[poly2] = 8;
        }
        else
        {
            poly_classes2[poly2] = 0;
        }
    }

    /*--- now set up list */

    for_less( dim, 0, N_DIMENSIONS )
    {
        min_range[dim] = parameters1[IJ(0,dim,3)];
        max_range[dim] = parameters1[IJ(0,dim,3)];
    }

    for_less( dim, 0, N_DIMENSIONS )
    {
        for_less( p, 0, n_points1 )
        {
            if( parameters1[IJ(p,dim,3)] < min_range[dim] )
                min_range[dim] = parameters1[IJ(p,dim,3)];
            if( parameters1[IJ(p,dim,3)] > max_range[dim] )
                max_range[dim] = parameters1[IJ(p,dim,3)];
        }
        for_less( p, 0, n_points2 )
        {
            if( parameters2[IJ(p,dim,3)] < min_range[dim] )
                min_range[dim] = parameters2[IJ(p,dim,3)];
            if( parameters2[IJ(p,dim,3)] > max_range[dim] )
                max_range[dim] = parameters2[IJ(p,dim,3)];
        }
    }

    *closest_distance = -1.0;

    *n_pairs_ptr = 0;
    *p1s_ptr = NULL;
    *p2s_ptr = NULL;
    *cases_ptr = NULL;
    *min_line_dists_ptr = NULL;
    recursive_find_close_pairs( min_range, max_range,
                        min_distance, search_distance,
                        parameters1,
                        poly_info1,
                        &list_of_polys1, &poly_classes1,
                        &n_alloced1, 0, n_polygons1, n_polygons1,
                        line_dir1,
                        parameters2,
                        poly_info2,
                        &list_of_polys2, &poly_classes2,
                        &n_alloced2, 0, n_polygons2, n_polygons2,
                        line_dir2,
                        n_pairs_ptr, p1s_ptr, p2s_ptr, cases_ptr,
                        min_line_dists_ptr,
                        n_points1, n_neighbours1, neighbours1,
                        n_points2, n_neighbours2, neighbours2,
                        closest_distance );

    if( n_alloced1 > 0 )
    {
        FREE( list_of_polys1 );
        FREE( poly_classes1 );
    }

    if( n_alloced2 > 0 )
    {
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
    surf_surf_lookup_struct       *ss_lookup )
{
    int                    i, w, s1, s2;
    Real                   *this_parms1, *this_parms2, *this_line1, *this_line2;
    Real                   closest_distance;
    surf_surf_struct       *surf;
    Real                   min_distance, close;

    closest_distance = -1.0;

    for_less( i, 0, deform->n_surf_surfs )
    {
        surf = &deform->surf_surfs[i];
        s1 = surf->surface_index1;
        s2 = surf->surface_index2;

        this_parms1 = &parameters[start_parameter[s1]];
        this_parms2 = &parameters[start_parameter[s2]];

        if( line_dir == NULL )
        {
            this_line1 = NULL;
            this_line2 = NULL;
        }
        else
        {
            this_line1 = &line_dir[start_parameter[s1]];
            this_line2 = &line_dir[start_parameter[s2]];
        }

        min_distance = 0.0;
        for_less( w, 0, surf->n_weights )
            min_distance = MAX( min_distance, surf->min_distances[w] );

        find_surf_intersect_candidates(
                       min_distance,
                       &ss_lookup[i].n_pairs,
                       &ss_lookup[i].p1s,
                       &ss_lookup[i].p2s,
                       &ss_lookup[i].cases,
                       &ss_lookup[i].min_line_dists,
                       deform->surfaces[s1].surface.n_points,
                       this_parms1, this_line1,
                       deform->surfaces[s1].surface.n_neighbours,
                       deform->surfaces[s1].surface.neighbours,
                       deform->surfaces[s2].surface.n_points,
                       this_parms2, this_line2,
                       deform->surfaces[s2].surface.n_neighbours,
                       deform->surfaces[s2].surface.neighbours,
                       max_movement, &close );

        if( close >= 0.0 &&
            (closest_distance < 0.0 || close < closest_distance) )
            closest_distance = close;
    }

    return( closest_distance );
}

public  void  initialize_surf_surf_lookup(
    surf_surf_lookup_struct  *ss_lookup )
{
    ss_lookup->n_pairs = 0;
    ss_lookup->p1s = NULL;
    ss_lookup->p2s = NULL;
    ss_lookup->cases = NULL;
    ss_lookup->min_line_dists = NULL;
}

public  void  delete_surf_surf_lookup(
    surf_surf_lookup_struct  *ss_lookup )
{
    if( ss_lookup->n_pairs > 0 )
    {
        FREE( ss_lookup->p1s );
        FREE( ss_lookup->p2s );
        FREE( ss_lookup->cases );
        FREE( ss_lookup->min_line_dists );
        ss_lookup->n_pairs = 0;
    }
}

public  int  get_n_surf_surf_candidate(
    surf_surf_lookup_struct  *ss_lookup )
{
    return( ss_lookup->n_pairs );
}

public  BOOLEAN   test_surf_surf_candidate(
    surf_surf_lookup_struct  *ss_lookup,
    int                      which,
    Real                     dist_from_computed_self_intersect,
    Real                     max_distance_sq,
    int                      n_points1,
    Real                     parameters1[],
    Smallest_int             active_flags1[],
    int                      n_neighbours1[],
    int                      *neighbours1[],
    int                      n_points2,
    Real                     parameters2[],
    Smallest_int             active_flags2[],
    int                      n_neighbours2[],
    int                      *neighbours2[],
    Real                     *dist_sq )
{
    int   p1, n1, p2, n2, n11, n12, n21, n22, prev_case, which_case;

    if( dist_from_computed_self_intersect> 0.0 &&
        (Real) ss_lookup->min_line_dists[which] >
        dist_from_computed_self_intersect)
    {
        return( FALSE );
    }

    p1 = ss_lookup->p1s[which] % n_points1;
    n1 = ss_lookup->p1s[which] / n_points1;
    p2 = ss_lookup->p2s[which] % n_points2;
    n2 = ss_lookup->p2s[which] / n_points2;

    prev_case = ss_lookup->cases[which];

    n11 = neighbours1[p1][n1];
    n12 = neighbours1[p1][(n1+1)%n_neighbours1[p1]];

    n21 = neighbours2[p2][n2];
    n22 = neighbours2[p2][(n2+1)%n_neighbours2[p2]];

    if( active_flags1 != NULL &&
        !active_flags1[p1] && !active_flags1[n11] && !active_flags1[n12] &&
        !active_flags2[p2] && !active_flags2[n21] && !active_flags2[n22] )
    {
        return( FALSE );
    }

    which_case = prev_case;
    *dist_sq = sq_triangle_triangle_dist( &parameters1[p1*3],
                                          &parameters1[n11*3],
                                          &parameters1[n12*3],
                                          &parameters2[p2*3],
                                          &parameters2[n21*3],
                                          &parameters2[n22*3], &which_case,0 );

    if( which_case != prev_case )
        ss_lookup->cases[which] = (unsigned char) which_case;

    return( TRUE );
}

public  void   create_surf_surf_deriv_info(
    surf_surf_lookup_struct       *ss_lookup,
    int                           which,
    Real                          dist_sq,
    int                           n_points1,
    Real                          parameters1[],
    int                           n_neighbours1[],
    int                           *neighbours1[],
    int                           n_points2,
    Real                          parameters2[],
    int                           n_neighbours2[],
    int                           *neighbours2[],
    surf_surf_deriv_struct        *deriv )
{
    int   p1, n1, p2, n2, n11, n12, n21, n22;

    p1 = ss_lookup->p1s[which] % n_points1;
    n1 = ss_lookup->p1s[which] / n_points1;
    p2 = ss_lookup->p2s[which] % n_points2;
    n2 = ss_lookup->p2s[which] / n_points2;

    n11 = neighbours1[p1][n1];
    n12 = neighbours1[p1][(n1+1)%n_neighbours1[p1]];
    n21 = neighbours2[p2][n2];
    n22 = neighbours2[p2][(n2+1)%n_neighbours2[p2]];

    if( dist_sq > 0.0 )
        deriv->dist = sqrt( dist_sq );
    else
      //deriv->dist = sqrt( dist_sq );
      // Added by June
      deriv->dist = 1;

    sq_triangle_triangle_dist_deriv( &parameters1[IJ(p1,0,3)],
                                     &parameters1[IJ(n11,0,3)],
                                     &parameters1[IJ(n12,0,3)],
                                     &parameters2[IJ(p2,0,3)],
                                     &parameters2[IJ(n21,0,3)],
                                     &parameters2[IJ(n22,0,3)],
                                     deriv->tri_tri_deriv[0],
                                     deriv->tri_tri_deriv[1],
                                     deriv->tri_tri_deriv[2],
                                     deriv->tri_tri_deriv[3],
                                     deriv->tri_tri_deriv[4],
                                     deriv->tri_tri_deriv[5] );
}

public  void   get_surf_surf_deriv(
    surf_surf_lookup_struct       *ss_lookup,
    int                           which,
    Real                          min_distance,
    Real                          weight,
    int                           n_points1,
    Real                          deriv1[],
    int                           n_neighbours1[],
    int                           *neighbours1[],
    int                           n_points2,
    Real                          deriv2[],
    int                           n_neighbours2[],
    int                           *neighbours2[],
    surf_surf_deriv_struct        *deriv_info )
{
    int   p1, n1, p2, n2, n11, n12, n21, n22;
    Real  diff, factor;

    p1 = ss_lookup->p1s[which] % n_points1;
    n1 = ss_lookup->p1s[which] / n_points1;
    p2 = ss_lookup->p2s[which] % n_points2;
    n2 = ss_lookup->p2s[which] / n_points2;

    n11 = neighbours1[p1][n1];
    n12 = neighbours1[p1][(n1+1)%n_neighbours1[p1]];
    n21 = neighbours2[p2][n2];
    n22 = neighbours2[p2][(n2+1)%n_neighbours2[p2]];

    diff = min_distance - deriv_info->dist;

    if( deriv_info->dist == 0.0 )
        factor = 0.0;
    else
    {
#ifdef USE_CUBE
#define USE_CUBE
        factor = -3.0 * diff * diff * weight / 2.0 / deriv_info->dist /
                 min_distance;
#else
        factor = -2.0 * diff * weight / 2.0 / deriv_info->dist;
#endif
    }

    deriv1[IJ(p1,0,3)] += factor * deriv_info->tri_tri_deriv[0][0];
    deriv1[IJ(p1,1,3)] += factor * deriv_info->tri_tri_deriv[0][1];
    deriv1[IJ(p1,2,3)] += factor * deriv_info->tri_tri_deriv[0][2];

    deriv1[IJ(n11,0,3)] += factor * deriv_info->tri_tri_deriv[1][0];
    deriv1[IJ(n11,1,3)] += factor * deriv_info->tri_tri_deriv[1][1];
    deriv1[IJ(n11,2,3)] += factor * deriv_info->tri_tri_deriv[1][2];

    deriv1[IJ(n12,0,3)] += factor * deriv_info->tri_tri_deriv[2][0];
    deriv1[IJ(n12,1,3)] += factor * deriv_info->tri_tri_deriv[2][1];
    deriv1[IJ(n12,2,3)] += factor * deriv_info->tri_tri_deriv[2][2];

    deriv2[IJ(p2,0,3)] += factor * deriv_info->tri_tri_deriv[3][0];
    deriv2[IJ(p2,1,3)] += factor * deriv_info->tri_tri_deriv[3][1];
    deriv2[IJ(p2,2,3)] += factor * deriv_info->tri_tri_deriv[3][2];

    deriv2[IJ(n21,0,3)] += factor * deriv_info->tri_tri_deriv[4][0];
    deriv2[IJ(n21,1,3)] += factor * deriv_info->tri_tri_deriv[4][1];
    deriv2[IJ(n21,2,3)] += factor * deriv_info->tri_tri_deriv[4][2];

    deriv2[IJ(n22,0,3)] += factor * deriv_info->tri_tri_deriv[5][0];
    deriv2[IJ(n22,1,3)] += factor * deriv_info->tri_tri_deriv[5][1];
    deriv2[IJ(n22,2,3)] += factor * deriv_info->tri_tri_deriv[5][2];
}
