#include  <fit_3d.h>

private  Real  get_max_movement_per_unit_line_method_bb(
    int    n_parameters,
    Real   line_dir[],
    Real   *lower_bound )
{
    int    p1;
    int    i, j, indices[6];
    Real   x_min, x_max, y_min, y_max, z_min, z_max, dx, dy, dz, max_move;
    Real   x, y, z, dist;

    x_min = 0.0;
    x_max = 0.0;
    y_min = 0.0;
    y_max = 0.0;
    z_min = 0.0;
    z_max = 0.0;
    for_less( p1, 0, n_parameters / N_DIMENSIONS )
    {
        if( p1 == 0 || line_dir[IJ(p1,0,3)] < x_min )
        {
            x_min = line_dir[IJ(p1,0,3)];
            indices[0] = p1;
        }
        if( p1 == 0 || line_dir[IJ(p1,0,3)] > x_max )
        {
            x_max = line_dir[IJ(p1,0,3)];
            indices[1] = p1;
        }
        if( p1 == 0 || line_dir[IJ(p1,1,3)] < y_min )
        {
            y_min = line_dir[IJ(p1,1,3)];
            indices[2] = p1;
        }
        if( p1 == 0 || line_dir[IJ(p1,1,3)] > y_max )
        {
            y_max = line_dir[IJ(p1,1,3)];
            indices[3] = p1;
        }
        if( p1 == 0 || line_dir[IJ(p1,2,3)] < z_min )
        {
            z_min = line_dir[IJ(p1,2,3)];
            indices[4] = p1;
        }
        if( p1 == 0 || line_dir[IJ(p1,2,3)] > z_max )
        {
            z_max = line_dir[IJ(p1,2,3)];
            indices[5] = p1;
        }
    }

    dx = x_max - x_min;
    dy = y_max - y_min;
    dz = z_max - z_min;

    max_move = dx * dx + dy * dy + dz * dz;

    if( max_move > 0.0 )
        max_move = sqrt( max_move );

    /*--- find a lower bound */

    *lower_bound = 0.0;
    for_less( i, 0, 5 )
    {
        x = line_dir[IJ(indices[i],0,3)];
        y = line_dir[IJ(indices[i],1,3)];
        z = line_dir[IJ(indices[i],2,3)];
        for_less( j, i+1, 6 )
        {
            dx = x - line_dir[IJ(indices[j],0,3)];
            dy = y - line_dir[IJ(indices[j],1,3)];
            dz = z - line_dir[IJ(indices[j],2,3)];
            dist = dx * dx + dy * dy + dz * dz;
            if( dist > *lower_bound )
                *lower_bound = dist;
        }
    }

    if( *lower_bound > 0.0 )
        *lower_bound = sqrt( *lower_bound );

    return( max_move );
}

private  Real  get_max_movement_per_unit_line_method_bs(
    int    n_parameters,
    Real   line_dir[] )
{
    int    p1;
    Real   dx, dy, dz, max_move, move;

    max_move = 0.0;
    for_less( p1, 0, n_parameters / N_DIMENSIONS )
    {
        dx = line_dir[IJ(p1,0,3)];
        dy = line_dir[IJ(p1,1,3)];
        dz = line_dir[IJ(p1,2,3)];

        move = dx * dx + dy * dy + dz * dz;

        if( p1 == 0 || move > max_move )
            max_move = move;
    }

    if( max_move > 0.0 )
        max_move = sqrt( max_move );

    return( 2.0 * max_move );
}

private  Real  get_max_movement_per_unit_line_method3(
    int    n_parameters,
    Real   line_dir[] )
{
    int   n_points, ****grid, ***n_in_grid, p, p1, p2;
    int   grid_size, *list, x, y, z, i;
    int   x_grid_min, x_grid_max, y_grid_min, y_grid_max;
    int   z_grid_min, z_grid_max;
    Real  x_min, x_max, y_min, y_max, z_min, z_max, dx, dy, dz;
    Real  move, max_move;
    Real  x1, x2, y1, y2, z1, z2, xp, yp, zp;

    n_points = n_parameters / N_DIMENSIONS;

    x_min = 0.0;
    x_max = 0.0;
    y_min = 0.0;
    y_max = 0.0;
    z_min = 0.0;
    z_max = 0.0;
    for_less( p, 0, n_points )
    {
        if( p == 0 || line_dir[IJ(p,0,3)] < x_min )
            x_min = line_dir[IJ(p,0,3)];
        if( p == 0 || line_dir[IJ(p,0,3)] > x_max )
            x_max = line_dir[IJ(p,0,3)];
        if( p == 0 || line_dir[IJ(p,1,3)] < y_min )
            y_min = line_dir[IJ(p,1,3)];
        if( p == 0 || line_dir[IJ(p,1,3)] > y_max )
            y_max = line_dir[IJ(p,1,3)];
        if( p == 0 || line_dir[IJ(p,2,3)] < z_min )
            z_min = line_dir[IJ(p,2,3)];
        if( p == 0 || line_dir[IJ(p,2,3)] > z_max )
            z_max = line_dir[IJ(p,2,3)];
    }

    grid_size = (int) pow( (Real) n_points / 5.0, 0.3333 ) + 1;

    ALLOC3D( grid, grid_size, grid_size, grid_size );
    ALLOC3D( n_in_grid, grid_size, grid_size, grid_size );

    for_less( x, 0, grid_size )
    for_less( y, 0, grid_size )
    for_less( z, 0, grid_size )
        n_in_grid[x][y][z] = 0;

    for_less( p, 0, n_points )
    {
        x = (int) ((Real) grid_size * (line_dir[IJ(p,0,3)] - x_min) /
                                                   (x_max - x_min));
        y = (int) ((Real) grid_size * (line_dir[IJ(p,1,3)] - y_min) /
                                                   (y_max - y_min));
        z = (int) ((Real) grid_size * (line_dir[IJ(p,2,3)] - z_min) /
                                                   (z_max - z_min));

        if( x == grid_size )
            x = grid_size - 1;
        if( y == grid_size )
            y = grid_size - 1;
        if( z == grid_size )
            z = grid_size - 1;

        ++n_in_grid[x][y][z];
    }

    ALLOC( list, n_points );

    for_less( x, 0, grid_size )
    for_less( y, 0, grid_size )
    for_less( z, 0, grid_size )
    {
        grid[x][y][z] = list;
        list = &list[n_in_grid[x][y][z]];
        n_in_grid[x][y][z] = 0;
    }

    for_less( p, 0, n_points )
    {
        x = (int) ((Real) grid_size * (line_dir[IJ(p,0,3)] - x_min) /
                               (x_max - x_min));
        y = (int) ((Real) grid_size * (line_dir[IJ(p,1,3)] - y_min) /
                               (y_max - y_min));
        z = (int) ((Real) grid_size * (line_dir[IJ(p,2,3)] - z_min) /
                               (z_max - z_min));

        if( x == grid_size )
            x = grid_size - 1;
        if( y == grid_size )
            y = grid_size - 1;
        if( z == grid_size )
            z = grid_size - 1;

        grid[x][y][z][n_in_grid[x][y][z]] = p;
        ++n_in_grid[x][y][z];
    }

    max_move = 0.0;

    for_less( p1, 0, n_points-1 )
    {
        xp = line_dir[IJ(p1,0,3)];
        yp = line_dir[IJ(p1,1,3)];
        zp = line_dir[IJ(p1,2,3)];

        x_grid_min = (int) ((Real) grid_size *
                     (xp - max_move - x_min) / (x_max - x_min));
        x_grid_max = (int) ((Real) grid_size *
                     (xp + max_move - x_min) / (x_max - x_min));

        y_grid_min = (int) ((Real) grid_size *
                     (yp - max_move - y_min) / (y_max - y_min));
        y_grid_max = (int) ((Real) grid_size *
                     (yp + max_move - y_min) / (y_max - y_min));

        z_grid_min = (int) ((Real) grid_size *
                     (zp - max_move - z_min) / (z_max - z_min));
        z_grid_max = (int) ((Real) grid_size *
                     (zp + max_move - z_min) / (z_max - z_min));

        x_grid_min = MAX( 0, x_grid_min+1 );
        x_grid_max = MIN( grid_size-1, x_grid_max-1 );
        y_grid_min = MAX( 0, y_grid_min+1 );
        y_grid_max = MIN( grid_size-1, y_grid_max-1 );
        z_grid_min = MAX( 0, z_grid_min+1 );
        z_grid_max = MIN( grid_size-1, z_grid_max-1 );

        for_less( x, 0, grid_size )
        for_less( y, 0, grid_size )
        for_less( z, 0, grid_size )
        {
            if( x >= x_grid_min && x <= x_grid_max &&
                y >= y_grid_min && y <= y_grid_max &&
                z >= z_grid_min && z <= z_grid_max )
                continue;

            x1 = (Real) grid_size * ((Real) x - x_min) / (x_max - x_min);
            x2 = (Real) grid_size * ((Real) (x+1) - x_min) / (x_max - x_min);
            y1 = (Real) grid_size * ((Real) y - y_min) / (y_max - y_min);
            y2 = (Real) grid_size * ((Real) (y+1) - y_min) / (y_max - y_min);
            z1 = (Real) grid_size * ((Real) z - z_min) / (z_max - z_min);
            z2 = (Real) grid_size * ((Real) (z+1) - z_min) / (z_max - z_min);

            dx = MAX3( 0.0, x2 - xp, xp - x1 );
            dy = MAX3( 0.0, y2 - yp, yp - y1 );
            dz = MAX3( 0.0, z2 - zp, zp - z1 );

            move = dx * dx + dy * dy + dz * dz;

            if( move < max_move )
                continue;

            for_less( i, 0, n_in_grid[x][y][z] )
            {
                p2 = grid[x][y][z][i];
                if( p2 <= p1 )
                    continue;

                dx = line_dir[IJ(p1,0,3)] - line_dir[IJ(p2,0,3)];
                dy = line_dir[IJ(p1,1,3)] - line_dir[IJ(p2,1,3)];
                dz = line_dir[IJ(p1,2,3)] - line_dir[IJ(p2,2,3)];

                move = dx * dx + dy * dy + dz * dz;

                if( move > max_move )
                    max_move = move;
            }
        }
    }

    FREE3D( n_in_grid );
    FREE3D( grid );

    if( max_move > 0.0 )
        max_move = sqrt( max_move );

    return( max_move );
}

#ifdef DEBUG
#define DEBUG
private  Real  get_max_movement_per_unit_line_method4(
    int    n_parameters,
    Real   line_dir[] )
{
    int    p1, p2, n_points;
    Real   dx, dy, dz, max_move, move;

    n_points = n_parameters / N_DIMENSIONS;
    max_move = 0.0;
    for_less( p1, 0, n_points-1 )
    for_less( p2, p1+1, n_points )
    {
        dx = line_dir[IJ(p1,0,3)] - line_dir[IJ(p2,0,3)];
        dy = line_dir[IJ(p1,1,3)] - line_dir[IJ(p2,1,3)];
        dz = line_dir[IJ(p1,2,3)] - line_dir[IJ(p2,2,3)];

        move = dx * dx + dy * dy + dz * dz;

        if( p1 == 0 && p2 == 1 || move > max_move )
            max_move = move;
    }

    if( max_move > 0.0 )
        max_move = sqrt( max_move );

    return( max_move );
}
#endif

#define  MAX_POINTS   5000

private  Real  get_max_movement_per_unit_line_method5(
    int    n_parameters,
    Real   line_dir[],
    Real   lower_bound,
    Real   max_dist_from_centre )
{
    int   n_points, p, i, j;
    Real  x, y, z, dx, dy, dz, reject_dist, diameter, dist;
    int   n_outer_points;
    Real  points[MAX_POINTS][N_DIMENSIONS];

    reject_dist = lower_bound - max_dist_from_centre;
    if( reject_dist < 0.0 )
        return( lower_bound );

    reject_dist = reject_dist * reject_dist;

    n_points = n_parameters / N_DIMENSIONS;
    n_outer_points = 0;

    for_less( p, 0, n_points )
    {
        x = line_dir[IJ(p,0,3)];
        y = line_dir[IJ(p,1,3)];
        z = line_dir[IJ(p,2,3)];

        dist = x * x + y * y + z * z;

        if( dist > reject_dist )
        {
            if( n_outer_points > MAX_POINTS )
            {
                return( 2.0 * max_dist_from_centre );
            }
            points[n_outer_points][0] = x;
            points[n_outer_points][1] = y;
            points[n_outer_points][2] = z;
            ++n_outer_points;
        }
    }

/*
    print( "N points after rejection: %d/%d\n", n_outer_points, n_points );
*/

    diameter = 0.0;
    for_less( i, 0, n_outer_points-1 )
    {
        x = points[i][0];
        y = points[i][1];
        z = points[i][2];
        for_less( j, i+1, n_outer_points )
        {
            dx = x - points[j][0];
            dy = y - points[j][1];
            dz = z - points[j][2];
            dist = dx * dx + dy * dy + dz * dz;
            if( dist > diameter )
                diameter = dist;
        }
    }

    if( diameter > 0.0 )
        diameter = sqrt( diameter );

    return( diameter );
}

#ifdef NO
private  void  output_points(
    int    n_parameters,
    Real   line_dir[] )
{
    int            point, n_points;
    object_struct  **objects;

    n_points = n_parameters / 3;

    ALLOC( objects, n_points );
    for_less( point, 0, n_points )
    {
        objects[point] = create_object( MARKER );
        initialize_marker( get_marker_ptr(objects[point]),
                           BOX_MARKER, WHITE );
        fill_Point( get_marker_ptr(objects[point])->position,
                    line_dir[IJ(point,X,3)],
                    line_dir[IJ(point,Y,3)],
                    line_dir[IJ(point,Z,3)] );
    }

    (void) output_graphics_file( "deriv.obj", BINARY_FORMAT, n_points,
                                 objects );

    delete_object_list( n_points, objects );
}
#endif

public  Real  get_max_movement_per_unit_line(
    int    n_parameters,
    Real   line_dir[] )
{
    Real   meth_bb, meth_bs, meth5, largest_axial_dist;
#ifdef DEBUG
    Real   meth3, meth4;
#endif

    if( line_dir == NULL )
        return( 0.0 );

    meth_bb = get_max_movement_per_unit_line_method_bb( n_parameters, line_dir,
                                                        &largest_axial_dist );
    meth_bs = get_max_movement_per_unit_line_method_bs( n_parameters, line_dir);

    meth5 = get_max_movement_per_unit_line_method5( n_parameters, line_dir,
                                                    largest_axial_dist,
                                                    meth_bs/2.0 );

#ifdef DEBUG

    if( getenv( "DEBUG" ) != NULL )
    {
        static  BOOLEAN  first = TRUE;
        if( first )
        {
            first = FALSE;
            output_points( n_parameters, line_dir );
        }
/*
    meth3 = get_max_movement_per_unit_line_method3( n_parameters, line_dir );
*/
        meth4 = get_max_movement_per_unit_line_method4( n_parameters, line_dir );
        print( "Movement 1: %g\n", meth_bb );
        print( "Movement 2: %g\n", meth_bs );
        print( "Movement 4: %g\n", meth4 );
        print( "Movement 5: %g\n", meth5 );
        return( MIN3( meth_bb,meth_bs,meth4 ) );
    }
    else
        return( MIN3( meth_bb, meth_bs, meth5 ) );
#else
{
    return( MIN3( meth_bb, meth_bs, meth5 ) );
}
#endif
}
