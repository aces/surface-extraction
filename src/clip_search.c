#include  <internal_volume_io.h>
#include  <fit_3d.h>

#define  MAX_POLYGONS_PER_NODE   5

typedef struct  clip_node_struct
{
    float                      pos;
    int                        n_polygons;
    int                        *polygons;
    struct  clip_node_struct   *children[2];
} clip_node_struct;

struct clip_struct
{
    int    n_nodes;
    int    debug_n_polygons;
    int    *n_neighbours;
    int    **neighbours;
    Real   *parameters;

    int                *all_polygons;

    Real               limits[N_DIMENSIONS][2];
    clip_node_struct   *root;
};

private  clip_node_struct  *create_node( void )
{
    clip_node_struct  *node;

    ALLOC( node, 1 );

    node->n_polygons = 0;
    node->children[0] = NULL;
    node->children[1] = NULL;

    return( node );
}

private  void  recursive_subdivide(
    clip_struct       *clip,
    Real              limits[N_DIMENSIONS][2],
    clip_node_struct  *node,
    int               axis,
    int               n_successive_failures,
    int               max_polygons_per_node )
{
    Real      save_limit, mid_plane, left_limit, right_limit;
    Real      min_pos, max_pos, mid_pos;
    int       start, end, p, n, neigh1, neigh2, swap;

    if( node->n_polygons <= max_polygons_per_node ||
        n_successive_failures >= N_DIMENSIONS )
        return;

    mid_plane = (limits[axis][0] + limits[axis][1]) / 2.0;
    left_limit = limits[axis][1];
    right_limit = limits[axis][0];

    start = 0;
    end = node->n_polygons - 1;
    while( start <= end )
    {
        p = node->polygons[start] % clip->n_nodes;;
        n = node->polygons[start] / clip->n_nodes;;
        neigh1 = clip->neighbours[p][n];
        neigh2 = clip->neighbours[p][(n+1)%clip->n_neighbours[p]];

        min_pos = MIN3( clip->parameters[IJ(p,axis,N_DIMENSIONS)],
                        clip->parameters[IJ(neigh1,axis,N_DIMENSIONS)],
                        clip->parameters[IJ(neigh2,axis,N_DIMENSIONS)] );
        max_pos = MAX3( clip->parameters[IJ(p,axis,N_DIMENSIONS)],
                        clip->parameters[IJ(neigh1,axis,N_DIMENSIONS)],
                        clip->parameters[IJ(neigh2,axis,N_DIMENSIONS)] );

        mid_pos = (min_pos + max_pos) / 2.0;

        if( mid_pos <= mid_plane )
        {
            right_limit = MAX( right_limit, max_pos );
            ++start;
        }
        else
        {
            left_limit = MIN( left_limit, min_pos );
            swap = node->polygons[start];
            node->polygons[start] = node->polygons[end];
            node->polygons[end] = swap;
            --end;
        }
    }

    if( start > 0 )
    {
        node->children[0] = create_node();
        node->children[0]->n_polygons = start;
        node->children[0]->polygons = node->polygons;
        node->children[0]->pos = (float) right_limit;

        if( node->children[0]->n_polygons != node->n_polygons )
            n_successive_failures = 0;
        else
            ++n_successive_failures;

        save_limit = limits[axis][1];
        limits[axis][1] = right_limit;
        recursive_subdivide( clip, limits, node->children[0],
                             (axis + 1) % N_DIMENSIONS,
                             n_successive_failures, max_polygons_per_node );
        limits[axis][1] = save_limit;
    }

    if( end < node->n_polygons-1 )
    {
        node->children[1] = create_node();
        node->children[1]->n_polygons = node->n_polygons - start;
        node->children[1]->polygons = &node->polygons[end+1];
        node->children[1]->pos = (float) left_limit;

        if( node->children[1]->n_polygons != node->n_polygons )
            n_successive_failures = 0;
        else
            ++n_successive_failures;

        save_limit = limits[axis][0];
        limits[axis][0] = left_limit;
        recursive_subdivide( clip, limits, node->children[1],
                             (axis + 1) % N_DIMENSIONS,
                             n_successive_failures, max_polygons_per_node );
        limits[axis][0] = save_limit;
    }

    node->n_polygons = 0;
}

public  clip_struct  *initialize_clip_search(
    int    n_nodes,
    int    n_neighbours[],
    int    *neighbours[],
    Real   parameters[] )
{
    int          node, n, neigh, n_polygons, neigh1, neigh2, dim;
    Real         pos;
    clip_struct  *clip;

    ALLOC( clip, 1 );

    clip->n_nodes = n_nodes;
    clip->n_neighbours = n_neighbours;
    clip->neighbours = neighbours;
    clip->parameters = parameters;

    for_less( node, 0, n_nodes )
    {
        for_less( dim, 0, N_DIMENSIONS )
        {
            pos = parameters[IJ(node,dim,N_DIMENSIONS)];
            if( node == 0 || pos < clip->limits[dim][0] )
                clip->limits[dim][0] = pos;
            if( node == 0 || pos > clip->limits[dim][1] )
                clip->limits[dim][1] = pos;
        }
    }

    n_polygons = 0;
    for_less( node, 0, n_nodes )
    {
        for_less( n, 0, n_neighbours[node] )
        {
            neigh1 = neighbours[node][n];
            neigh2 = neighbours[node][(n+1)%n_neighbours[node]];

            if( node < neigh1 && node < neigh2 )
                ++n_polygons;
        }
    }

    ALLOC( clip->all_polygons, n_polygons );
    n_polygons = 0;

    for_less( node, 0, n_nodes )
    {
        for_less( n, 0, n_neighbours[node] )
        {
            neigh1 = neighbours[node][n];
            neigh2 = neighbours[node][(n+1)%n_neighbours[node]];

            if( node > neigh1 || node > neigh2 )
                continue;

            clip->all_polygons[n_polygons] = IJ(n,node,n_nodes);
            ++n_polygons;
        }
    }

    clip->root = create_node();
    clip->root->n_polygons = n_polygons;
    clip->root->polygons = clip->all_polygons;
    clip->debug_n_polygons = n_polygons;

    recursive_subdivide( clip, clip->limits, clip->root, X,
                         0, MAX_POLYGONS_PER_NODE );

    return( clip );
}


private  BOOLEAN   intersect_ray_triangle(
    Point            *ray_origin,
    Vector           *ray_direction,
    Real             point0[],
    Real             point1[],
    Real             point2[],
    Real             *dist )
{
    Real     n_dot_d, d;
    Real     v01x, v01y, v01z, v02x, v02y, v02z, nx, ny, nz, rx, ry, rz;
    Real     vx, vy, vz;
    Real     tx0, ty0, tz0, tx1, ty1, tz1, tx2, ty2, tz2;
    Real     sign;

    tx0 = point0[X] - RPoint_x( *ray_origin );
    ty0 = point0[Y] - RPoint_y( *ray_origin );
    tz0 = point0[Z] - RPoint_z( *ray_origin );

    tx1 = point1[X] - RPoint_x( *ray_origin );
    ty1 = point1[Y] - RPoint_y( *ray_origin );
    tz1 = point1[Z] - RPoint_z( *ray_origin );

    tx2 = point2[X] - RPoint_x( *ray_origin );
    ty2 = point2[Y] - RPoint_y( *ray_origin );
    tz2 = point2[Z] - RPoint_z( *ray_origin );

    v01x = tx1 - tx0;
    v01y = ty1 - ty0;
    v01z = tz1 - tz0;
    v02x = tx2 - tx0;
    v02y = ty2 - ty0;
    v02z = tz2 - tz0;

    nx = v01y * v02z - v01z * v02y;
    ny = v01z * v02x - v01x * v02z;
    nz = v01x * v02y - v01y * v02x;

    rx = RVector_x(*ray_direction);
    ry = RVector_y(*ray_direction);
    rz = RVector_z(*ray_direction);

    n_dot_d = rx * nx + ry * ny + rz * nz;

    if( n_dot_d == 0.0 )
        return( FALSE );

    vx = ty1 * tz0 - tz1 * ty0;
    vy = tz1 * tx0 - tx1 * tz0;
    vz = tx1 * ty0 - ty1 * tx0;
    d = rx * vx + ry * vy + rz * vz;
    if( n_dot_d * d > 0.0 )
        return( FALSE );

    vx = ty2 * tz1 - tz2 * ty1;
    vy = tz2 * tx1 - tx2 * tz1;
    vz = tx2 * ty1 - ty2 * tx1;
    d = rx * vx + ry * vy + rz * vz;
    if( n_dot_d * d > 0.0 )
        return( FALSE );

    vx = ty0 * tz2 - tz0 * ty2;
    vy = tz0 * tx2 - tx0 * tz2;
    vz = tx0 * ty2 - ty0 * tx2;
    d = rx * vx + ry * vy + rz * vz;
    if( n_dot_d * d > 0.0 )
        return( FALSE );

    *dist = (nx * tx0 + ny * ty0 + nz * tz0)/ n_dot_d;

    return( TRUE );
}

private  BOOLEAN  clip_line_to_polygon(
    Point             *origin,
    Vector            *normal,
    Real              *outer_distance,
    Real              *inner_distance,
    Real              p1[],
    Real              p2[],
    Real              p3[] )
{
    Real       t;

    if( !intersect_ray_triangle( origin, normal, p1, p2, p3, &t ) )
        return( FALSE );

    if( t > 0.0 )
    {
        if( t < *outer_distance )
        {
            *outer_distance = t;
            return( TRUE );
        }
    }
    else if( t < 0.0 )
    {
        t = -t;
        if( t < *inner_distance )
        {
            *inner_distance = t;
            return( TRUE );
        }
    }
    return( FALSE );
}

private  void  recursive_clip(
    clip_struct       *clip,
    clip_node_struct  *node,
    int               axis,
    int               n_nodes_to_ignore,
    int               nodes_to_ignore[],
    Point             *origin,
    Vector            *normal,
    Real              t_min,
    Real              t_max,
    Real              *outer_distance,
    Real              *inner_distance )
{
    int      i, j, p, n, neigh1, neigh2, poly;
    Real     save_limit, t1, t2, delta, t;
    BOOLEAN  search;

    if( node->n_polygons != 0 )
    {
        for_less( i, 0, node->n_polygons )
        {
            poly = node->polygons[i];
            p = poly % clip->n_nodes;
            n = poly / clip->n_nodes;
            neigh1 = clip->neighbours[p][n];
            neigh2 = clip->neighbours[p][(n+1)%clip->n_neighbours[p]];

            for_less( j, 0, n_nodes_to_ignore )
            {
                if( nodes_to_ignore[j] != p &&
                    nodes_to_ignore[j] != neigh1 &&
                    nodes_to_ignore[j] != neigh2 )
                {
                    break;
                }
            }

            if( j >= n_nodes_to_ignore )
                continue;

            (void) clip_line_to_polygon( origin, normal, outer_distance,
                                     inner_distance,
                                     &clip->parameters[p*N_DIMENSIONS],
                                     &clip->parameters[neigh1*N_DIMENSIONS],
                                     &clip->parameters[neigh2*N_DIMENSIONS] );
        }

        return;
    }

    delta = RVector_coord(*normal,axis);

    if( node->children[0] != NULL )
    {
        t1 = t_min;
        t2 = t_max;
        if( delta != 0.0 )
        {
            t = ((Real) node->children[0]->pos - RPoint_coord(*origin,axis)) /
                delta;

            if( delta > 0.0 && t < t2 )
                t2 = t;
            else if( delta < 0.0 && t > t1 )
                t1 = t;

            search = t1 <= t2;
        }
        else
        {
            search = RPoint_coord(*origin,axis) <
                     (Real) node->children[0]->pos;
        }

        if( search )
        {
            recursive_clip( clip, node->children[0],
                            (axis + 1) % N_DIMENSIONS,
                            n_nodes_to_ignore, nodes_to_ignore, origin, normal,
                            t1, t2, outer_distance, inner_distance );
        }
    }

    if( node->children[1] != NULL )
    {
        t1 = t_min;
        t2 = t_max;
        if( delta != 0.0 )
        {
            t = ((Real) node->children[1]->pos - RPoint_coord(*origin,axis)) /
                delta;

            if( delta > 0.0 && t > t1 )
                t1 = t;
            else if( delta < 0.0 && t < t2 )
                t2 = t;

            search = t1 <= t2;
        }
        else
        {
            search = RPoint_coord(*origin,axis) >
                     (Real) node->children[1]->pos;
        }

        if( search )
        {
            recursive_clip( clip, node->children[1],
                            (axis + 1) % N_DIMENSIONS,
                            n_nodes_to_ignore, nodes_to_ignore, origin, normal,
                            t1, t2, outer_distance, inner_distance );
        }
    }
}

private  void  brute_force_clip(
    clip_struct       *clip,
    int               n_polygons,
    int               n_nodes_to_ignore,
    int               nodes_to_ignore[],
    Point             *origin,
    Vector            *normal,
    Real              *outer_distance,
    Real              *inner_distance )
{
    int      i, j, p, n, neigh1, neigh2, poly;
    Real     save_limit;

    for_less( i, 0, n_polygons )
    {
        poly = clip->all_polygons[i];
        p = poly % clip->n_nodes;
        n = poly / clip->n_nodes;
        neigh1 = clip->neighbours[p][n];
        neigh2 = clip->neighbours[p][(n+1)%clip->n_neighbours[p]];

        for_less( j, 0, n_nodes_to_ignore )
        {
            if( nodes_to_ignore[j] != p &&
                nodes_to_ignore[j] != neigh1 &&
                nodes_to_ignore[j] != neigh2 )
            {
                break;
            }
        }

        if( j >= n_nodes_to_ignore )
            continue;

        (void) clip_line_to_polygon( origin, normal, outer_distance,
                                     inner_distance,
                                     &clip->parameters[p*N_DIMENSIONS],
                                     &clip->parameters[neigh1*N_DIMENSIONS],
                                     &clip->parameters[neigh2*N_DIMENSIONS] );
    }
}

public  void  clip_search_line(
    clip_struct  *clip,
    int          n_nodes_to_ignore,
    int          nodes_to_ignore[],
    Point        *origin,
    Vector       *normal,
    Real         outer_distance,
    Real         inner_distance,
    Real         *clipped_outer_distance,
    Real         *clipped_inner_distance )
{
    recursive_clip( clip, clip->root, X,
                    n_nodes_to_ignore, nodes_to_ignore,
                    origin, normal, -inner_distance, outer_distance,
                    clipped_outer_distance, clipped_inner_distance );
}

private  void  recursive_delete_tree(
    clip_node_struct  *node )
{
    if( node->children[0] != NULL )
        recursive_delete_tree( node->children[0] );

    if( node->children[1] != NULL )
        recursive_delete_tree( node->children[1] );

    FREE( node );
}

public  void  delete_clip_search(
    clip_struct  *clip )
{
    recursive_delete_tree( clip->root );
    FREE( clip->all_polygons );
    FREE( clip );
}
