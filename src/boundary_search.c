#include  <volume_io/internal_volume_io.h>
#include  <fit_3d.h>
#include  <bicpl/numerical.h>

#define  MAX_IN_VOXEL_COEF_LOOKUP  20000

private  void  find_voxel_line_polynomial(
    Real        coefs[],
    int         x,
    int         y,
    int         z,
    Real        x_origin,
    Real        y_origin,
    Real        z_origin,
    Real        x_dir,
    Real        y_dir,
    Real        z_dir,
    Real        line_poly[] );

private  void   get_trilinear_gradient(
    Real   coefs[],
    Real   u,
    Real   v,
    Real   w,
    Real   derivs[] );

private  BOOLEAN  trilinear_voxel_contains_range(
    voxel_coef_struct   *lookup,
    int                 x,
    int                 y,
    int                 z,
    Real                isovalue,
    Real                coefs[] );

private  void  clip_line_to_volume(
    int        sizes[],
    Real       line_origin[],
    Real       line_direction[],
    Real       *t_min,
    Real       *t_max )
{
    Point  origin;
    Vector direction;

    fill_Point( origin, line_origin[X], line_origin[Y], line_origin[Z] );
    fill_Point( direction,
                line_direction[X], line_direction[Y], line_direction[Z] );

    (void) clip_line_to_box( &origin, &direction,
                             0.0, (Real) sizes[X]-1.0,
                             0.0, (Real) sizes[Y]-1.0,
                             0.0, (Real) sizes[Z]-1.0,
                             t_min, t_max );
}

public  BOOLEAN  find_isosurface_boundary_in_direction(
    Volume                      volume,
    voxel_coef_struct           *lookup,
    bitlist_3d_struct           *done_bits,
    bitlist_3d_struct           *surface_bits,
    Point                       *point_ray_origin,
    Vector                      *vector_unit_dir,
    Real                        max_outwards_search_distance,
    Real                        max_inwards_search_distance,
    int                         degrees_continuity,
    Real                        isovalue,
    Normal_directions           normal_direction,
    Real                        *boundary_distance )
{
    private const Real  TOLERANCE = 0.0001;
    BOOLEAN       found, done0, done1, contains;
    Real          next_closest0, next_closest1, coefs[8];
    Real          rx0, ry0, rz0;
    Real          rx1, ry1, rz1, t_min, t_max;
    int           x0, y0, z0, x1, y1, z1;
    int           dx_voxel, dy_voxel, dz_voxel, dim;
    Real          current_distance0, current_distance1;
    Real          t_max0, t_max1;
    Real          max_dist;
    Real          voxel_dir[N_DIMENSIONS];
    Real          voxel_origin[N_DIMENSIONS];
    Real          ox, oy, oz;
    Real          dx, dy, dz, ndx, ndy, ndz;
    Real          next_x0, next_y0, next_z0, next_x1, next_y1, next_z1;
    Real          delta_dist_x, delta_dist_y, delta_dist_z;
    Real          ray_origin[N_DIMENSIONS], unit_dir[N_DIMENSIONS];
    int           word_index, bit_sub_index;
    int           i;
    Real          dot_prod, vx, vy, vz, pos;
    Real          first_deriv[N_DIMENSIONS];
    BOOLEAN       deriv_dir_correct;
    int           n_boundaries;
    Real          boundary_positions[3];
    Real          line_coefs[4];
    int           x_offset, y_offset;
    bitlist_type  word, bit0, *done_flags0, *surface_flags0;
    bitlist_type  bit1, *done_flags1, *surface_flags1;

    if( degrees_continuity != 0 )
    {
        print_error( "find_isosurface_boundary_in_direction(): " );
        print_error( "ignoring degrees_continuity != 0\n" );
    }

    for_less( dim, 0, N_DIMENSIONS )
    {
        ray_origin[dim] = RPoint_coord( *point_ray_origin, dim );
        unit_dir[dim] = RVector_coord( *vector_unit_dir, dim );
    }

    voxel_origin[0] = lookup->trans00 * ray_origin[0] +
                      lookup->trans01 * ray_origin[1] +
                      lookup->trans02 * ray_origin[2] +
                      lookup->trans03;
    voxel_origin[1] = lookup->trans10 * ray_origin[0] +
                      lookup->trans11 * ray_origin[1] +
                      lookup->trans12 * ray_origin[2] +
                      lookup->trans13;
    voxel_origin[2] = lookup->trans20 * ray_origin[0] +
                      lookup->trans21 * ray_origin[1] +
                      lookup->trans22 * ray_origin[2] +
                      lookup->trans23;

    voxel_dir[0] = lookup->trans00 * unit_dir[0] +
                   lookup->trans01 * unit_dir[1] +
                   lookup->trans02 * unit_dir[2];
    voxel_dir[1] = lookup->trans10 * unit_dir[0] +
                   lookup->trans11 * unit_dir[1] +
                   lookup->trans12 * unit_dir[2];
    voxel_dir[2] = lookup->trans20 * unit_dir[0] +
                   lookup->trans21 * unit_dir[1] +
                   lookup->trans22 * unit_dir[2];

    clip_line_to_volume( lookup->sizes, voxel_origin, voxel_dir,
                         &t_min, &t_max );

    t_min += TOLERANCE;
    t_max -= TOLERANCE;

    current_distance0 = MAX( 0.0, t_min );
    t_max0 = MIN( max_outwards_search_distance, t_max );

    current_distance1 = MAX( 0.0, -t_max );
    t_max1 = MIN( max_inwards_search_distance, -t_min );

    dx = voxel_dir[X];
    dy = voxel_dir[Y];
    dz = voxel_dir[Z];
    ndx = -dx;
    ndy = -dy;
    ndz = -dz;

    ox = voxel_origin[X];
    oy = voxel_origin[Y];
    oz = voxel_origin[Z];

    rx0 = ox + current_distance0 * dx;
    ry0 = oy + current_distance0 * dy;
    rz0 = oz + current_distance0 * dz;
    x0 = FLOOR( rx0 );
    y0 = FLOOR( ry0 );
    z0 = FLOOR( rz0 );

    rx1 = ox + current_distance1 * ndx;
    ry1 = oy + current_distance1 * ndy;
    rz1 = oz + current_distance1 * ndz;
    x1 = FLOOR( rx1 );
    y1 = FLOOR( ry1 );
    z1 = FLOOR( rz1 );

    if( dx == 0.0 ) {
        dx_voxel = 0;
        next_x0 = t_max0 + 1.0;
        next_x1 = t_max1 + 1.0;
    } else if( dx > 0.0 ) {
        dx_voxel = 1;
        next_x0 = current_distance0 + ((Real) x0 + 1.0 - rx0) / dx;
        next_x1 = current_distance1 -((Real) x1 - rx1) / dx;
        delta_dist_x = 1.0 / dx;
    } else {
        dx_voxel = -1;
        next_x0 = current_distance0 + ((Real) x0 - rx0) / dx;
        next_x1 = current_distance1 -((Real) x1 + 1.0 - rx1) / dx;
        delta_dist_x = -1.0 / dx;
    }

    if( dy == 0.0 ) {
        dy_voxel = 0;
        next_y0 = t_max0 + 1.0;
        next_y1 = t_max1 + 1.0;
    } else if( dy > 0.0 ) {
        dy_voxel = 1;
        next_y0 = current_distance0 + ((Real) y0 + 1.0 - ry0) / dy;
        next_y1 = current_distance1 -((Real) y1 - ry1) / dy;
        delta_dist_y = 1.0 / dy;
    } else {
        dy_voxel = -1;
        next_y0 = current_distance0 + ((Real) y0 - ry0) / dy;
        next_y1 = current_distance1 -((Real) y1 + 1.0 - ry1) / dy;
        delta_dist_y = -1.0 / dy;
    }

    if( dz == 0.0 ) {
        dz_voxel = 0;
        next_z0 = t_max0 + 1.0;
        next_z1 = t_max1 + 1.0;
    } else if( dz > 0.0 ) {
        dz_voxel = 1;
        next_z0 = current_distance0 + ((Real) z0 + 1.0 - rz0) / dz;
        next_z1 = current_distance1 -((Real) z1 - rz1) / dz;
        delta_dist_z = 1.0 / dz;
    } else {
        dz_voxel = -1;
        next_z0 = current_distance0 + ((Real) z0 - rz0) / dz;
        next_z1 = current_distance1 -((Real) z1 + 1.0 - rz1) / dz;
        delta_dist_z = -1.0 / dz;
    }

    next_closest0 = MIN3( next_x0, next_y0, next_z0 );
    next_closest1 = MIN3( next_x1, next_y1, next_z1 );

    found = FALSE;
    done0 = current_distance0 >= t_max0;
    done1 = current_distance1 >= t_max1;

    if( current_distance0 >= t_max0 )
        current_distance0 = MAX( t_max0, t_max1 );
    if( current_distance1 >= t_max1 )
        current_distance1 = MAX( t_max0, t_max1 );

    if( !done0 ) {
        bit_sub_index = z0 & BITS_PER_BITLIST_WORD_MINUS_1;
        bit0 = ((bitlist_type) 1 << (bitlist_type) bit_sub_index);
        word_index = z0 >> LOG_BITS_PER_BITLIST_WORD;
        done_flags0 = &done_bits->bits[x0][y0][word_index];
        surface_flags0 = &surface_bits->bits[x0][y0][word_index];
    }

    if( !done1 ) {
        bit_sub_index = z1 & BITS_PER_BITLIST_WORD_MINUS_1;
        bit1 = ((bitlist_type) 1 << (bitlist_type) bit_sub_index);
        word_index = z1 >> LOG_BITS_PER_BITLIST_WORD;
        done_flags1 = &done_bits->bits[x1][y1][word_index];
        surface_flags1 = &surface_bits->bits[x1][y1][word_index];
    }

    if( done0 && done1 )
        return( FALSE );

    y_offset = done_bits->n_z_words;
    x_offset = done_bits->ny * y_offset;
    y_offset *= dy_voxel;
    x_offset *= dx_voxel;

    while( TRUE ) {
        if( current_distance0 <= current_distance1 ) {
            max_dist = next_closest0;
            if( max_dist >= t_max0 ) {
                max_dist = t_max0;
                done0 = TRUE;
            }

            word = *done_flags0;

            if( (word & bit0) == 0 ) {
                contains = trilinear_voxel_contains_range( lookup, x0, y0, z0,
                                                           isovalue, coefs );

                if( contains )
                    *surface_flags0 |= bit0;

                *done_flags0 = (word | bit0);
            } else {
                contains = (*surface_flags0 & bit0) != 0;
                if( contains )
                    lookup_volume_coeficients( lookup, x0, y0, z0, coefs );
            }

            if( contains ) {
                find_voxel_line_polynomial( coefs, x0, y0, z0,
                                            ox, oy, oz, dx, dy, dz,
                                            line_coefs );

                line_coefs[0] -= isovalue;

                n_boundaries = solve_cubic( line_coefs[3], line_coefs[2],
                                            line_coefs[1], line_coefs[0],
                                            boundary_positions );

                for_less( i, 0, n_boundaries ) {
                    pos = boundary_positions[i];
                    if( pos >= current_distance0 && pos <= max_dist &&
                        (!found || pos < FABS(*boundary_distance)) ) {
                        vx = ox + pos * dx;
                        vy = oy + pos * dy;
                        vz = oz + pos * dz;
                        deriv_dir_correct = TRUE;
                        if( normal_direction != ANY_DIRECTION ) {
                            get_trilinear_gradient( coefs,
                                                    vx - (Real) x0,
                                                    vy - (Real) y0,
                                                    vz - (Real) z0,
                                                    first_deriv );

                            dot_prod = first_deriv[X] * dx +
                                       first_deriv[Y] * dy +
                                       first_deriv[Z] * dz;
                            if( normal_direction == TOWARDS_LOWER &&
                                dot_prod > 0.0 ||
                                normal_direction == TOWARDS_HIGHER &&
                                dot_prod < 0.0 ) {
                                deriv_dir_correct = FALSE;
                            }
                        }

                        if( deriv_dir_correct ) {
                            *boundary_distance = pos;
                            found = TRUE;
                        }
                    }
                }
            }

            if( done0 ) {
                current_distance0 = MAX( t_max0, t_max1 );
                if( done1 )
                    break;
            } else {
                current_distance0 = next_closest0;

                if( next_x0 <= current_distance0 ) {
                    next_x0 += delta_dist_x;
                    x0 += dx_voxel;
                    done_flags0 += x_offset;
                    surface_flags0 += x_offset;
                } else if( next_y0 <= current_distance0 ) {
                    next_y0 += delta_dist_y;
                    y0 += dy_voxel;
                    done_flags0 += y_offset;
                    surface_flags0 += y_offset;
                } else if( next_z0 <= current_distance0 ) {
                    next_z0 += delta_dist_z;
                    z0 += dz_voxel;
                    if( dz_voxel > 0 ) {
                        bit0 <<= 1;
                        if( bit0 == 0 ) {
                            bit0 = 1;
                            ++done_flags0;
                            ++surface_flags0;
                        }
                    } else if( dz_voxel < 0 ) {
                        bit0 >>= 1;
                        if( bit0 == 0 ) {
                            bit0 = ((bitlist_type) 1 <<
                                 (bitlist_type) BITS_PER_BITLIST_WORD_MINUS_1);
                            --done_flags0;
                            --surface_flags0;
                        }
                    }
                }

                next_closest0 = MIN3( next_x0, next_y0, next_z0 );
            }
        } else {
            max_dist = next_closest1;
            if( max_dist >= t_max1 ) {
                max_dist = t_max1;
                done1 = TRUE;
            }

            word = *done_flags1;

            if( (word & bit1) == 0 ) {
                contains = trilinear_voxel_contains_range( lookup, x1, y1, z1,
                                                           isovalue, coefs );

                if( contains )
                    *surface_flags1 |= bit1;

                *done_flags1 = (word | bit1);
            } else {
                contains = (*surface_flags1 & bit1) != 0;
                if( contains )
                    lookup_volume_coeficients( lookup, x1, y1, z1, coefs );
            }

            if( contains ) {
                find_voxel_line_polynomial( coefs, x1, y1, z1,
                                            ox, oy, oz, ndx, ndy, ndz,
                                            line_coefs );

                line_coefs[0] -= isovalue;

                n_boundaries = solve_cubic( line_coefs[3], line_coefs[2],
                                            line_coefs[1], line_coefs[0],
                                            boundary_positions );

                for_less( i, 0, n_boundaries ) {
                    pos = boundary_positions[i];
                    if( pos >= current_distance1 && pos <= max_dist &&
                        (!found || pos < FABS(*boundary_distance)) ) {
                        vx = ox + pos * ndx;
                        vy = oy + pos * ndy;
                        vz = oz + pos * ndz;
                        deriv_dir_correct = TRUE;
                        if( normal_direction != ANY_DIRECTION ) {
                            get_trilinear_gradient( coefs,
                                                    vx - (Real) x1,
                                                    vy - (Real) y1,
                                                    vz - (Real) z1,
                                                    first_deriv );

                            dot_prod = first_deriv[X] * dx +
                                       first_deriv[Y] * dy +
                                       first_deriv[Z] * dz;
                            if( normal_direction == TOWARDS_LOWER &&
                                dot_prod > 0.0 ||
                                normal_direction == TOWARDS_HIGHER &&
                                dot_prod < 0.0 ) {
                                deriv_dir_correct = FALSE;
                            }
                        }

                        if( deriv_dir_correct ) {
                            *boundary_distance = -pos;
                            found = TRUE;
                        }
                    }
                }
            }

            if( done1 ) {
                current_distance1 = MAX( t_max0, t_max1 );
                if( done0 )
                    break;
            } else {
                current_distance1 = next_closest1;

                if( next_x1 <= current_distance1 ) {
                    next_x1 += delta_dist_x;
                    x1 -= dx_voxel;
                    done_flags1 -= x_offset;
                    surface_flags1 -= x_offset;
                } else if( next_y1 <= current_distance1 ) {
                    next_y1 += delta_dist_y;
                    y1 -= dy_voxel;
                    done_flags1 -= y_offset;
                    surface_flags1 -= y_offset;
                } else if( next_z1 <= current_distance1 ) {
                    next_z1 += delta_dist_z;
                    z1 -= dz_voxel;

                    if( dz_voxel < 0 ) {
                        bit1 <<= 1;
                        if( bit1 == 0 ) {
                            bit1 = 1;
                            ++done_flags1;
                            ++surface_flags1;
                        }
                    } else if( dz_voxel > 0 ) {
                        bit1 >>= 1;
                        if( bit1 == 0 ) {
                            bit1 = ((bitlist_type) 1 <<
                                 (bitlist_type) BITS_PER_BITLIST_WORD_MINUS_1);
                            --done_flags1;
                            --surface_flags1;
                        }
                    }
                }

                next_closest1 = MIN3( next_x1, next_y1, next_z1 );
            }
        }

        if( found &&
            (done0 || FABS(*boundary_distance) <= current_distance0) &&
            (done1 || FABS(*boundary_distance) <= current_distance1) ) {
            break;
        }
    }

    return( found );
}

private  void   get_trilinear_gradient(
    Real   coefs[],
    Real   u,
    Real   v,
    Real   w,
    Real   derivs[] )
{
    Real  du00, du01, du10, du11, c00, c01, c10, c11;
    Real  dv0, dv1, c0, c1, dw, du0, du1;

    /*--- get the 4 differences in the u direction */

    du00 = coefs[4] - coefs[0];
    du01 = coefs[5] - coefs[1];
    du10 = coefs[6] - coefs[2];
    du11 = coefs[7] - coefs[3];

    /*--- reduce to a 2D problem, by interpolating in the u direction */

    c00 = coefs[0] + u * du00;
    c01 = coefs[1] + u * du01;
    c10 = coefs[2] + u * du10;
    c11 = coefs[3] + u * du11;

    /*--- get the 2 differences in the v direction for the 2D problem */

    dv0 = c10 - c00;
    dv1 = c11 - c01;
    /*--- reduce 2D to a 1D problem, by interpolating in the v direction */

    c0 = c00 + v * dv0;
    c1 = c01 + v * dv1;

    /*--- get the 1 difference in the w direction for the 1D problem */

    dw = c1 - c0;

    /*--- reduce the 2D u derivs to 1D */

    du0 = INTERPOLATE( v, du00, du10 );
    du1 = INTERPOLATE( v, du01, du11 );

    /*--- interpolate the 1D problems in w, or for Z deriv, just use dw */

    derivs[X] = INTERPOLATE( w, du0, du1 );
    derivs[Y] = INTERPOLATE( w, dv0, dv1 );
    derivs[Z] = dw;
}



private  BOOLEAN  trilinear_voxel_contains_range(
    voxel_coef_struct   *lookup,
    int                 x,
    int                 y,
    int                 z,
    Real                isovalue,
    Real                values[] )
{
    int      i;

    lookup_volume_coeficients( lookup, x, y, z, values );

    if( values[0] < isovalue )
    {
        for_less( i, 1, 8 )
        {
            if( values[i] >= isovalue )
                return( TRUE );
        }
        return( FALSE );
    }
    else if( values[0] > isovalue )
    {
        for_less( i, 1, 8 )
        {
            if( values[i] <= isovalue )
                return( TRUE );
        }
        return( FALSE );
    }
    else
        return( TRUE );
}

public  void  initialize_lookup_volume_coeficients(
    voxel_coef_struct  *lookup,
    Volume             volume )
{
    int        voxel, x_stride, y_stride, z_stride;
    Transform  *trans, inverse;

    if( get_volume_data_type(volume) != UNSIGNED_BYTE )
        handle_internal_error( "initialize_lookup_volume_coeficients: invalid volume type" );

    get_volume_sizes( volume, lookup->sizes );

    for_less( voxel, 0, 256 )
    {
        lookup->voxel_to_real_values[voxel] = convert_voxel_to_value( volume,
                                                 (Real) voxel );
    }

    x_stride = lookup->sizes[Y] * lookup->sizes[Z];
    y_stride = lookup->sizes[Z];
    z_stride = 1;

    lookup->offset1 =                       z_stride;
    lookup->offset2 =            y_stride;
    lookup->offset3 =            y_stride + z_stride;
    lookup->offset4 = x_stride;
    lookup->offset5 = x_stride +            z_stride;
    lookup->offset6 = x_stride + y_stride;
    lookup->offset7 = x_stride + y_stride + z_stride;

    trans = get_linear_transform_ptr(get_voxel_to_world_transform(volume));

    compute_transform_inverse( trans, &inverse );

    lookup->trans00 = Transform_elem(inverse,0,0);
    lookup->trans01 = Transform_elem(inverse,0,1);
    lookup->trans02 = Transform_elem(inverse,0,2);
    lookup->trans03 = Transform_elem(inverse,0,3);

    lookup->trans10 = Transform_elem(inverse,1,0);
    lookup->trans11 = Transform_elem(inverse,1,1);
    lookup->trans12 = Transform_elem(inverse,1,2);
    lookup->trans13 = Transform_elem(inverse,1,3);

    lookup->trans20 = Transform_elem(inverse,2,0);
    lookup->trans21 = Transform_elem(inverse,2,1);
    lookup->trans22 = Transform_elem(inverse,2,2);
    lookup->trans23 = Transform_elem(inverse,2,3);

    if( Transform_elem(inverse,3,0) != 0.0 ||
        Transform_elem(inverse,3,1) != 0.0 ||
        Transform_elem(inverse,3,2) != 0.0 ||
        Transform_elem(inverse,3,3) != 1.0 )
        handle_internal_error( "initialize_lookup_volume_coeficients" );

    lookup->voxel_ptr = (unsigned char ***) (volume->array.data);
}

public  void  lookup_volume_coeficients(
    voxel_coef_struct  *lookup,
    int                x,
    int                y,
    int                z,
    Real               coefs[] )
{
    unsigned char          *ptr;

    ptr = &lookup->voxel_ptr[x][y][z];

    coefs[0] = lookup->voxel_to_real_values[*ptr];
    coefs[1] = lookup->voxel_to_real_values[ptr[lookup->offset1]];
    coefs[2] = lookup->voxel_to_real_values[ptr[lookup->offset2]];
    coefs[3] = lookup->voxel_to_real_values[ptr[lookup->offset3]];
    coefs[4] = lookup->voxel_to_real_values[ptr[lookup->offset4]];
    coefs[5] = lookup->voxel_to_real_values[ptr[lookup->offset5]];
    coefs[6] = lookup->voxel_to_real_values[ptr[lookup->offset6]];
    coefs[7] = lookup->voxel_to_real_values[ptr[lookup->offset7]];
}

public  void  delete_lookup_volume_coeficients(
    voxel_coef_struct  *lookup )
{
}

private  void  find_voxel_line_polynomial(
    Real        coefs[],
    int         x,
    int         y,
    int         z,
    Real        x_origin,
    Real        y_origin,
    Real        z_origin,
    Real        x_dir,
    Real        y_dir,
    Real        z_dir,
    Real        line_poly[] )
{
    Real  ou, ov, ow, delta0, delta1, delta0t, delta1t, delta0tt;
    Real  d00, d01, d10, d11;
    Real  c00a, c01a, c10a, c11a;
    Real  c00t, c01t, c10t, c11t;
    Real  c0a, c1a, c0t, c1t, c0tt, c1tt;

    d00 = coefs[1] - coefs[0];
    d01 = coefs[3] - coefs[2];
    d10 = coefs[5] - coefs[4];
    d11 = coefs[7] - coefs[6];

    ou = x_origin - (Real) x;
    ov = y_origin - (Real) y;
    ow = z_origin - (Real) z;

    c00a = coefs[0] + ow * d00;
    c01a = coefs[2] + ow * d01;
    c10a = coefs[4] + ow * d10;
    c11a = coefs[6] + ow * d11;
    c00t = z_dir * d00;
    c01t = z_dir * d01;
    c10t = z_dir * d10;
    c11t = z_dir * d11;

    delta0 = c01a - c00a;
    delta1 = c11a - c10a;
    delta0t = c01t - c00t;
    delta1t = c11t - c10t;

    c0a = c00a + ov * delta0;
    c0t = c00t + ov * delta0t + y_dir * delta0;
    c0tt = y_dir * delta0t;
    c1a = c10a + ov * delta1;
    c1t = c10t + ov * delta1t + y_dir * delta1;
    c1tt = y_dir * delta1t;

    delta0 = c1a - c0a;
    delta0t = c1t - c0t;
    delta0tt = c1tt - c0tt;

    line_poly[0] = c0a + ou * delta0;
    line_poly[1] = c0t + ou * delta0t + x_dir * delta0;
    line_poly[2] = c0tt + ou * delta0tt + x_dir * delta0t;
    line_poly[3] = x_dir * delta0tt;
}
