#include  <fit_3d.h>

#define  TOLERANCE_FACTOR   1.0e-10

#define  SUB( d, a, b ) \
          { GLUE(d,x) = GLUE(a,x) - GLUE(b,x); \
            GLUE(d,y) = GLUE(a,y) - GLUE(b,y); \
            GLUE(d,z) = GLUE(a,z) - GLUE(b,z); }

#define  CROSS( d, a, b ) \
          { GLUE(d,x) = GLUE(a,y) * GLUE(b,z) - GLUE(a,z) * GLUE(b,y); \
            GLUE(d,y) = GLUE(a,z) * GLUE(b,x) - GLUE(a,x) * GLUE(b,z); \
            GLUE(d,z) = GLUE(a,x) * GLUE(b,y) - GLUE(a,y) * GLUE(b,x); }

#define  DOT( a, b ) \
            (GLUE(a,x) * GLUE(b,x) + GLUE(a,y) * GLUE(b,y) + \
             GLUE(a,z) * GLUE(b,z))

#define  A0  (1 << 0)
#define  A1  (1 << 1)
#define  A2  (1 << 2)
#define  B0  (1 << 3)
#define  B1  (1 << 4)
#define  B2  (1 << 5)

#define  VERTEX_VERTEX( a0, a1, a2, b0, b1, b2 ) \
        SUB( delta, b0, a0 ); \
        SUB( da0, a1, a0 ); \
        SUB( da1, a2, a0 ); \
        SUB( db0, b1, b0 ); \
        SUB( db1, b2, b0 ); \
 \
        if( DOT( delta, da0 ) <= 0.0 && \
            DOT( delta, da1 ) <= 0.0 && \
            DOT( delta, db0 ) >= 0.0 && \
            DOT( delta, db1 ) >= 0.0 ) \
        { \
            *dist = DOT( delta, delta ); \
            return( TRUE ); \
        }

#define  VERTEX_EDGE( a0, a1, a2, b0, b1, b2 ) \
        SUB( da0, a1, a0 ); \
        len = DOT( da0, da0 ); \
 \
        if( len == 0.0 ) \
            return( FALSE ); \
 \
        SUB( delta, b0, a0 ); \
        t = DOT( delta, da0 ) / len; \
        if( t < 0.0 || t > 1.0 ) \
            return( FALSE ); \
 \
        deltax -= t * da0x; \
        deltay -= t * da0y; \
        deltaz -= t * da0z; \
 \
        SUB( da0, a2, a0 ); \
        SUB( db0, b1, b0 ); \
        SUB( db1, b2, b0 ); \
 \
        if( DOT( delta, da0 ) <= 0.0 && \
            DOT( delta, db0 ) >= 0.0 && DOT( delta, db1 ) >= 0.0 ) \
        { \
            *dist = DOT( delta, delta ); \
            return( TRUE ); \
        }

#define  EDGE_EDGE( a0, a1, a2, b0, b1, b2 ) \
        SUB( da0, a1, a0 ); \
        SUB( db0, b1, b0 ); \
        SUB( da2, a2, a0 ); \
        SUB( db2, b2, b0 ); \
        CROSS( delta, da0, db0 ); \
        len = DOT( delta, delta ); \
        if( len == 0.0 ) \
            return( FALSE ); \
        SUB( d, b0, a0 ); \
        d = DOT( d, delta ); \
        if( d >= 0.0 && DOT( delta, da2 ) <= 0.0 && DOT( delta, db2 ) >= 0.0 ||\
            d  < 0.0 && DOT( delta, da2 ) >= 0.0 && DOT( delta, db2 ) <= 0.0 ) \
        { \
            *dist = d * d / len; \
            return( TRUE ); \
        }

#define  POINT_TRIANGLE( a0, a1, a2, b0, b1, b2 ) \
        SUB( da0, a1, a0 ); \
        SUB( da2, a0, a2 ); \
        SUB( db0, b1, b0 ); \
        SUB( db2, b0, b2 ); \
        CROSS( an, da2, da0 ); \
        d = DOT( b0, an ) - DOT( a0, an ); \
        if( d * DOT( db0, an ) >= 0.0 && d * DOT( db2, an ) <= 0.0 ) \
        { \
            a0_dot_da0 = DOT( a0, da0 ); \
            x_dot_y = DOT( a2, da0 ) - a0_dot_da0; \
            len_da0 = DOT( da0, da0 ); \
            len_da2 = DOT( da2, da2 ); \
 \
            bottom = len_da0 * len_da2 - x_dot_y * x_dot_y; \
 \
            if( bottom == 0.0 ) \
                return( FALSE ); \
 \
            x_dot_v = DOT( b0, da0 ) - a0_dot_da0; \
            y_dot_v = -(DOT( b0, da2 ) - DOT( a0, da2 )); \
 \
            x_pos = (x_dot_v * len_da2 - x_dot_y * y_dot_v) / bottom; \
            y_pos = (y_dot_v * len_da0 - x_dot_y * x_dot_v) / bottom; \
 \
            if( x_pos >= 0.0 && y_pos >= 0.0 && x_pos + y_pos <= 1.0 ) \
            { \
                len_an = DOT( an, an ); \
                *dist = d * d / len_an; \
                return( TRUE ); \
            } \
        }

private  BOOLEAN  test_distance_with_known_case(
    Real            a0[],
    Real            a1[],
    Real            a2[],
    Real            b0[],
    Real            b1[],
    Real            b2[],
    int             which_case,
    Real            *dist )
{
    Real      a0x, a0y, a0z, a1x, a1y, a1z, a2x, a2y, a2z;
    Real      b0x, b0y, b0z, b1x, b1y, b1z, b2x, b2y, b2z;
    Real      da0x, da0y, da0z, da1x, da1y, da1z, da2x, da2y, da2z;
    Real      db2x, db2y, db2z, bnx, bny, bnz, eb0x, eb0y, eb0z;
    Real      anx, any, anz, ea0x, ea0y, ea0z, dot;
    Real      db0x, db0y, db0z, db1x, db1y, db1z;
    Real      deltax, deltay, deltaz, len, t;
    Real      c, d, f, s, denom, len_da0, len_db0, dx, dy, dz;
    Real      a0_dot_da0, x_dot_y, len_da2, bottom, x_dot_v;
    Real      y_dot_v, x_pos, y_pos, len_an;

    a0x = a0[X];     a0y = a0[Y];     a0z = a0[Z];
    a1x = a1[X];     a1y = a1[Y];     a1z = a1[Z];
    a2x = a2[X];     a2y = a2[Y];     a2z = a2[Z];

    b0x = b0[X];     b0y = b0[Y];     b0z = b0[Z];
    b1x = b1[X];     b1y = b1[Y];     b1z = b1[Z];
    b2x = b2[X];     b2y = b2[Y];     b2z = b2[Z];

    switch( which_case )
    {
    case A0 | B0:
        VERTEX_VERTEX( a0, a1, a2, b0, b1, b2 );
        break;

    case A0 | B1:
        VERTEX_VERTEX( a0, a1, a2, b1, b2, b0 );
        break;

    case A0 | B2:
        VERTEX_VERTEX( a0, a1, a2, b2, b0, b1 );
        break;

    case A1 | B0:
        VERTEX_VERTEX( a1, a2, a0, b0, b1, b2 );
        break;

    case A1 | B1:
        VERTEX_VERTEX( a1, a2, a0, b1, b2, b0 );
        break;

    case A1 | B2:
        VERTEX_VERTEX( a1, a2, a0, b2, b0, b1 );
        break;

    case A2 | B0:
        VERTEX_VERTEX( a2, a0, a1, b0, b1, b2 );
        break;

    case A2 | B1:
        VERTEX_VERTEX( a2, a0, a1, b1, b2, b0 );
        break;

    case A2 | B2:
        VERTEX_VERTEX( a2, a0, a1, b2, b0, b1 );
        break;

    case A0 | A1 | B0:
        VERTEX_EDGE( a0, a1, a2, b0, b1, b2 );
        break;

    case A0 | A1 | B1:
        VERTEX_EDGE( a0, a1, a2, b1, b2, b0);
        break;

    case A0 | A1 | B2:
        VERTEX_EDGE( a0, a1, a2, b2, b0, b1 );
        break;

    case A1 | A2 | B0:
        VERTEX_EDGE( a1, a2, a0, b0, b1, b2 );
        break;

    case A1 | A2 | B1:
        VERTEX_EDGE( a1, a2, a0, b1, b2, b0);
        break;

    case A1 | A2 | B2:
        VERTEX_EDGE( a1, a2, a0, b2, b0, b1 );
        break;

    case A2 | A0 | B0:
        VERTEX_EDGE( a2, a0, a1, b0, b1, b2 );
        break;

    case A2 | A0 | B1:
        VERTEX_EDGE( a2, a0, a1, b1, b2, b0);
        break;

    case A2 | A0 | B2:
        VERTEX_EDGE( a2, a0, a1, b2, b0, b1 );
        break;

    case B0 | B1 | A0:
        VERTEX_EDGE( b0, b1, b2, a0, a1, a2 );
        break;

    case B0 | B1 | A1:
        VERTEX_EDGE( b0, b1, b2, a1, a2, a0);
        break;

    case B0 | B1 | A2:
        VERTEX_EDGE( b0, b1, b2, a2, a0, a1 );
        break;

    case B1 | B2 | A0:
        VERTEX_EDGE( b1, b2, b0, a0, a1, a2 );
        break;

    case B1 | B2 | A1:
        VERTEX_EDGE( b1, b2, b0, a1, a2, a0);
        break;

    case B1 | B2 | A2:
        VERTEX_EDGE( b1, b2, b0, a2, a0, a1 );
        break;

    case B2 | B0 | A0:
        VERTEX_EDGE( b2, b0, b1, a0, a1, a2 );
        break;

    case B2 | B0 | A1:
        VERTEX_EDGE( b2, b0, b1, a1, a2, a0);
        break;

    case B2 | B0 | A2:
        VERTEX_EDGE( b2, b0, b1, a2, a0, a1 );
        break;

    case A0 | A1 | B0 | B1:
        EDGE_EDGE( a0, a1, a2, b0, b1, b2 );
        break;

    case A0 | A1 | B1 | B2:
        EDGE_EDGE( a0, a1, a2, b1, b2, b0 );
        break;

    case A0 | A1 | B2 | B0:
        EDGE_EDGE( a0, a1, a2, b2, b0, b1 );
        break;

    case A1 | A2 | B0 | B1:
        EDGE_EDGE( a1, a2, a0, b0, b1, b2 );
        break;

    case A1 | A2 | B1 | B2:
        EDGE_EDGE( a1, a2, a0, b1, b2, b0 );
        break;

    case A1 | A2 | B2 | B0:
        EDGE_EDGE( a1, a2, a0, b2, b0, b1 );
        break;

    case A2 | A0 | B0 | B1:
        EDGE_EDGE( a2, a0, a1, b0, b1, b2 );
        break;

    case A2 | A0 | B1 | B2:
        EDGE_EDGE( a2, a0, a1, b1, b2, b0 );
        break;

    case A2 | A0 | B2 | B0:
        EDGE_EDGE( a2, a0, a1, b2, b0, b1 );
        break;

    /*------- vertex - triangle */

    case A0 | A1 | A2 | B0:
        POINT_TRIANGLE( a0, a1, a2, b0, b1, b2 );
        break;

    case A0 | A1 | A2 | B1:
        POINT_TRIANGLE( a0, a1, a2, b1, b2, b0 );
        break;

    case A0 | A1 | A2 | B2:
        POINT_TRIANGLE( a0, a1, a2, b2, b0, b1 );
        break;

    case B0 | B1 | B2 | A0:
        POINT_TRIANGLE( b0, b1, b2, a0, a1, a2 );
        break;

    case B0 | B1 | B2 | A1:
        POINT_TRIANGLE( b0, b1, b2, a1, a2, a0 );
        break;

    case B0 | B1 | B2 | A2:
        POINT_TRIANGLE( b0, b1, b2, a2, a0, a1 );
        break;
    }

    return( FALSE );
}

private  Real  sq_triangle_triangle_dist_and_weights(
    Real            a0[],
    Real            a1[],
    Real            a2[],
    Real            b0[],
    Real            b1[],
    Real            b2[],
    Real            weights[],
    int             *tri_case )
{
    Real      a0x, a0y, a0z, a1x, a1y, a1z, a2x, a2y, a2z;
    Real      b0x, b0y, b0z, b1x, b1y, b1z, b2x, b2y, b2z;
    Real      da0x, da0y, da0z, da1x, da1y, da1z, da2x, da2y, da2z;
    Real      db0x, db0y, db0z, db1x, db1y, db1z, db2x, db2y, db2z;
    Real      anx, any, anz, bnx, bny, bnz;
    Real      x_dot_y, x_dot_v, y_dot_v;
    Real      bottom, x_pos, y_pos;
    Real      dx, dy, dz;
    Real      s, t, c, d, f;
    Real      denom;
    Real      dist;
    Real      len_an, len_bn;
    Real      v_da0, v_da1, v_da2, v_ea0, v_ea1, v_ea2;
    Real      v_eb0, v_eb1, v_eb2;
    Real      v_db0, v_db1, v_db2;
    BOOLEAN   outside_ea0, outside_ea1, outside_ea2;
    BOOLEAN   outside_eb0, outside_eb1, outside_eb2;
    Real      a0_dot_db0, a1_dot_db0, a2_dot_db0;
    Real      a0_dot_db1, a1_dot_db1, a2_dot_db1;
    Real      a0_dot_db2, a1_dot_db2, a2_dot_db2;
    Real      b0_dot_da0, b1_dot_da0, b2_dot_da0;
    Real      b0_dot_da1, b1_dot_da1, b2_dot_da1;
    Real      b0_dot_da2, b1_dot_da2, b2_dot_da2;
    Real      ea0x, ea0y, ea0z, ea1x, ea1y, ea1z, ea2x, ea2y, ea2z;
    Real      eb0x, eb0y, eb0z, eb1x, eb1y, eb1z, eb2x, eb2y, eb2z;
    Real      a0_dot_ea0, a1_dot_ea1, a2_dot_ea2;
    Real      b0_dot_eb0, b1_dot_eb1, b2_dot_eb2;
    Real      a0_dot_da0, a1_dot_da1, a2_dot_da2;
    Real      b0_dot_db0, b1_dot_db1, b2_dot_db2;
    Real      a0_dot_an;
    Real      b0_dot_bn;
    Real      len_da0, len_da1, len_da2, len_db0, len_db1, len_db2;
    Real      db0_dot_an, db1_dot_an, db2_dot_an;
    Real      da0_dot_bn, da1_dot_bn, da2_dot_bn;
    Real      px, py, pz;
    Real      characteristic_length, tolerance;

    weights[0] = 0.0;
    weights[1] = 0.0;
    weights[2] = 0.0;
    weights[3] = 0.0;
    weights[4] = 0.0;
    weights[5] = 0.0;

    a0x = a0[X];     a0y = a0[Y];     a0z = a0[Z];
    a1x = a1[X];     a1y = a1[Y];     a1z = a1[Z];
    a2x = a2[X];     a2y = a2[Y];     a2z = a2[Z];

    b0x = b0[X];     b0y = b0[Y];     b0z = b0[Z];
    b1x = b1[X];     b1y = b1[Y];     b1z = b1[Z];
    b2x = b2[X];     b2y = b2[Y];     b2z = b2[Z];

    SUB( da0, a1, a0 );
    SUB( da1, a2, a1 );
    SUB( da2, a0, a2 );

    SUB( db0, b1, b0 );
    SUB( db1, b2, b1 );
    SUB( db2, b0, b2 );

    CROSS( an, da2, da0 );
    CROSS( bn, db2, db0 );

    CROSS( ea0, an, da0 );
    CROSS( ea1, an, da1 );
    CROSS( ea2, an, da2 );

    CROSS( eb0, bn, db0 );
    CROSS( eb1, bn, db1 );
    CROSS( eb2, bn, db2 );

    b0_dot_da0 = DOT( b0, da0 );
    b1_dot_da0 = DOT( b1, da0 );
    b2_dot_da0 = DOT( b2, da0 );
    b0_dot_da1 = DOT( b0, da1 );
    b1_dot_da1 = DOT( b1, da1 );
    b2_dot_da1 = DOT( b2, da1 );
    b0_dot_da2 = DOT( b0, da2 );
    b1_dot_da2 = DOT( b1, da2 );
    b2_dot_da2 = DOT( b2, da2 );

    a0_dot_db0 = DOT( a0, db0 );
    a1_dot_db0 = DOT( a1, db0 );
    a2_dot_db0 = DOT( a2, db0 );
    a0_dot_db1 = DOT( a0, db1 );
    a1_dot_db1 = DOT( a1, db1 );
    a2_dot_db1 = DOT( a2, db1 );
    a0_dot_db2 = DOT( a0, db2 );
    a1_dot_db2 = DOT( a1, db2 );
    a2_dot_db2 = DOT( a2, db2 );

    a0_dot_ea0 = DOT( a0, ea0 );
    a1_dot_ea1 = DOT( a1, ea1 );
    a2_dot_ea2 = DOT( a2, ea2 );

    b0_dot_eb0 = DOT( b0, eb0 );
    b1_dot_eb1 = DOT( b1, eb1 );
    b2_dot_eb2 = DOT( b2, eb2 );

    a0_dot_da0 = DOT( a0, da0 );
    a1_dot_da1 = DOT( a1, da1 );
    a2_dot_da2 = DOT( a2, da2 );

    a0_dot_an = DOT( a0, an );

    b0_dot_db0 = DOT( b0, db0 );
    b1_dot_db1 = DOT( b1, db1 );
    b2_dot_db2 = DOT( b2, db2 );

    b0_dot_bn = DOT( b0, bn );

    len_da0 = DOT( da0, da0 );
    len_da1 = DOT( da1, da1 );
    len_da2 = DOT( da2, da2 );

    len_db0 = DOT( db0, db0 );
    len_db1 = DOT( db1, db1 );
    len_db2 = DOT( db2, db2 );

    characteristic_length = (len_da0 + len_da1 + len_da2 +
                             len_db0 + len_db1 + len_db2) / 6.0;
    tolerance = characteristic_length * TOLERANCE_FACTOR;

    len_an = DOT( an, an );
    len_bn = DOT( bn, bn );

    db0_dot_an = DOT( db0, an );
    db1_dot_an = DOT( db1, an );
    db2_dot_an = DOT( db2, an );

    da0_dot_bn = DOT( da0, bn );
    da1_dot_bn = DOT( da1, bn );
    da2_dot_bn = DOT( da2, bn );

    /*--- test b0 against the triangle a0, a1, a2 */

    v_da0 = b0_dot_da0 - a0_dot_da0;
    v_da1 = b0_dot_da1 - a1_dot_da1;
    v_da2 = b0_dot_da2 - a2_dot_da2;

    v_ea0 = DOT( b0, ea0 ) - a0_dot_ea0;
    v_ea1 = DOT( b0, ea1 ) - a1_dot_ea1;
    v_ea2 = DOT( b0, ea2 ) - a2_dot_ea2;

    outside_ea0 = (v_ea0 <= 0.0);
    outside_ea1 = (v_ea1 <= 0.0);
    outside_ea2 = (v_ea2 <= 0.0);

    if( !outside_ea0 && !outside_ea1 && !outside_ea2 )
    {
        d = DOT( b0, an ) - a0_dot_an;
        if( d * db0_dot_an >= -tolerance && d * db2_dot_an <= tolerance )
        {
            x_dot_y = DOT( a2, da0 ) - a0_dot_da0;

            bottom = len_da0 * len_da2 - x_dot_y * x_dot_y;

            if( bottom != 0.0 )
            {
                x_dot_v = v_da0;
                y_dot_v = -(DOT( b0, da2 ) - DOT( a0, da2 ));

                x_pos = (x_dot_v * len_da2 - x_dot_y * y_dot_v) / bottom;
                y_pos = (y_dot_v * len_da0 - x_dot_y * x_dot_v) / bottom;

                weights[0] = 1.0 - x_pos - y_pos;
                weights[1] = x_pos;
                weights[2] = y_pos;
                weights[3] = 1.0;

                dist = d * d / len_an;
                *tri_case = A0 | A1 | A2 | B0;

                return( dist );
            }
        }
    }
    else if( outside_ea0 && v_da0 > 0.0 && v_da0 < len_da0 )
    {
        d = v_da0 / len_da0;
        px = a0x + d * da0x;
        py = a0y + d * da0y;
        pz = a0z + d * da0z;
        SUB( d, b0, p );

        if( DOT( d, db0 ) >= -tolerance && DOT( d, db2 ) <= tolerance )
        {
            weights[0] = 1.0 - d;
            weights[1] = d;
            weights[3] = 1.0;
            dist = DOT( d, d );
            *tri_case = A0 | A1 | B0;
            return( dist );
        }
    }
    else if( outside_ea1 && v_da1 > 0.0 && v_da1 < len_da1 )
    {
        d = v_da1 / len_da1;
        px = a1x + d * da1x;
        py = a1y + d * da1y;
        pz = a1z + d * da1z;
        SUB( d, b0, p );

        if( DOT( d, db0 ) >= -tolerance && DOT( d, db2 ) <= tolerance )
        {
            weights[1] = 1.0 - d;
            weights[2] = d;
            weights[3] = 1.0;
            dist = DOT( d, d );
            *tri_case = A1 | A2 | B0;
            return( dist );
        }
    }
    else if( outside_ea2 && v_da2 > 0.0 && v_da2 < len_da2 )
    {
        d = v_da2 / len_da2;
        px = a2x + d * da2x;
        py = a2y + d * da2y;
        pz = a2z + d * da2z;
        SUB( d, b0, p );

        if( DOT( d, db0 ) >= -tolerance && DOT( d, db2 ) <= tolerance )
        {
            weights[2] = 1.0 - d;
            weights[0] = d;
            weights[3] = 1.0;
            dist = DOT( d, d );
            *tri_case = A0 | A2 | B0;
            return( dist );
        }
    }
    else if( v_da0 <= 0.0 && v_da2 >= len_da2 )
    {
        if( b0_dot_db0 - a0_dot_db0 >= -tolerance &&
            DOT( b0, db2 ) - a0_dot_db2 <= tolerance )
        {
            weights[0] = 1.0;
            weights[3] = 1.0;
            SUB( d, b0, a0 );
            dist = DOT( d, d );
            *tri_case = A0 | B0;
            return( dist );
        }
    }
    else if( v_da1 <= 0.0 && v_da0 >= len_da0 )
    {
        if( b0_dot_db0 - a1_dot_db0 >= -tolerance &&
            DOT( b0, db2 ) - a1_dot_db2 <= tolerance )
        {
            weights[1] = 1.0;
            weights[3] = 1.0;
            SUB( d, b0, a1 );
            dist = DOT( d, d );
            *tri_case = A1 | B0;
            return( dist );
        }
    }
    else
    {
        if( b0_dot_db0 - a2_dot_db0 >= -tolerance &&
            DOT( b0, db2 ) - a2_dot_db2 <= tolerance )
        {
            weights[2] = 1.0;
            weights[3] = 1.0;
            SUB( d, b0, a2 );
            dist = DOT( d, d );
            *tri_case = A2 | B0;
            return( dist );
        }
    }

    /*--- test b1 against the triangle a0, a1, a2 */

    v_da0 = b1_dot_da0 - a0_dot_da0;
    v_da1 = b1_dot_da1 - a1_dot_da1;
    v_da2 = b1_dot_da2 - a2_dot_da2;

    v_ea0 = DOT( b1, ea0 ) - a0_dot_ea0;
    v_ea1 = DOT( b1, ea1 ) - a1_dot_ea1;
    v_ea2 = DOT( b1, ea2 ) - a2_dot_ea2;

    outside_ea0 = (v_ea0 <= 0.0);
    outside_ea1 = (v_ea1 <= 0.0);
    outside_ea2 = (v_ea2 <= 0.0);

    if( !outside_ea0 && !outside_ea1 && !outside_ea2 )
    {
        d = DOT( b1, an ) - a0_dot_an;
        if( d * db1_dot_an >= -tolerance && d * db0_dot_an <= tolerance )
        {
            x_dot_y = DOT( a2, da0 ) - a0_dot_da0;

            bottom = len_da0 * len_da2 - x_dot_y * x_dot_y;

            if( bottom != 0.0 )
            {
                x_dot_v = v_da0;
                y_dot_v = -(DOT( b1, da2 ) - DOT( a0, da2 ));

                x_pos = (x_dot_v * len_da2 - x_dot_y * y_dot_v) / bottom;
                y_pos = (y_dot_v * len_da0 - x_dot_y * x_dot_v) / bottom;

                weights[0] = 1.0 - x_pos - y_pos;
                weights[1] = x_pos;
                weights[2] = y_pos;
                weights[4] = 1.0;

                dist = d * d / len_an;

                *tri_case = A0 | A1 | A2 | B1;
                return( dist );
            }
        }
    }
    else if( outside_ea0 && v_da0 > 0.0 && v_da0 < len_da0 )
    {
        d = v_da0 / len_da0;
        px = a0x + d * da0x;
        py = a0y + d * da0y;
        pz = a0z + d * da0z;
        SUB( d, b1, p );

        if( DOT( d, db1 ) >= -tolerance && DOT( d, db0 ) <= tolerance )
        {
            weights[0] = 1.0 - d;
            weights[1] = d;
            weights[4] = 1.0;
            dist = DOT( d, d );
            *tri_case = A0 | A1 | B1;
            return( dist );
        }
    }
    else if( outside_ea1 && v_da1 > 0.0 && v_da1 < len_da1 )
    {
        d = v_da1 / len_da1;
        px = a1x + d * da1x;
        py = a1y + d * da1y;
        pz = a1z + d * da1z;
        SUB( d, b1, p );

        if( DOT( d, db1 ) >= -tolerance && DOT( d, db0 ) <= tolerance )
        {
            weights[1] = 1.0 - d;
            weights[2] = d;
            weights[4] = 1.0;
            dist = DOT( d, d );
            *tri_case = A1 | A2 | B1;
            return( dist );
        }
    }
    else if( outside_ea2 && v_da2 > 0.0 && v_da2 < len_da2 )
    {
        d = v_da2 / len_da2;
        px = a2x + d * da2x;
        py = a2y + d * da2y;
        pz = a2z + d * da2z;
        SUB( d, b1, p );

        if( DOT( d, db1 ) >= -tolerance && DOT( d, db0 ) <= tolerance )
        {
            weights[2] = 1.0 - d;
            weights[0] = d;
            weights[4] = 1.0;
            dist = DOT( d, d );
            *tri_case = A0 | A2 | B1;
            return( dist );
        }
    }
    else if( v_da0 <= 0.0 && v_da2 >= len_da2 )
    {
        if( b1_dot_db1 - a0_dot_db1 >= -tolerance &&
            DOT( b1, db0 ) - a0_dot_db0 <= tolerance )
        {
            weights[0] = 1.0;
            weights[4] = 1.0;
            SUB( d, b1, a0 );
            dist = DOT( d, d );
            *tri_case = A0 | B1;
            return( dist );
        }
    }
    else if( v_da1 <= 0.0 && v_da0 >= len_da0 )
    {
        if( b1_dot_db1 - a1_dot_db1 >= -tolerance &&
            DOT( b1, db0 ) - a1_dot_db0 <= tolerance )
        {
            weights[1] = 1.0;
            weights[4] = 1.0;
            SUB( d, b1, a1 );
            dist = DOT( d, d );
            *tri_case = A1 | B1;
            return( dist );
        }
    }
    else
    {
        if( b1_dot_db1 - a2_dot_db1 >= -tolerance &&
            DOT( b1, db0 ) - a2_dot_db0 <= tolerance )
        {
            weights[2] = 1.0;
            weights[4] = 1.0;
            SUB( d, b1, a2 );
            dist = DOT( d, d );
            *tri_case = A2 | B1;
            return( dist );
        }
    }

    /*--- test b2 against the triangle a0, a1, a2 */

    v_da0 = b2_dot_da0 - a0_dot_da0;
    v_da1 = b2_dot_da1 - a1_dot_da1;
    v_da2 = b2_dot_da2 - a2_dot_da2;

    v_ea0 = DOT( b2, ea0 ) - a0_dot_ea0;
    v_ea1 = DOT( b2, ea1 ) - a1_dot_ea1;
    v_ea2 = DOT( b2, ea2 ) - a2_dot_ea2;

    outside_ea0 = (v_ea0 <= 0.0);
    outside_ea1 = (v_ea1 <= 0.0);
    outside_ea2 = (v_ea2 <= 0.0);

    if( !outside_ea0 && !outside_ea1 && !outside_ea2 )
    {
        d = DOT( b2, an ) - a0_dot_an;
        if( d * db2_dot_an >= -tolerance && d * db1_dot_an <= tolerance )
        {
            x_dot_y = DOT( a2, da0 ) - a0_dot_da0;

            bottom = len_da0 * len_da2 - x_dot_y * x_dot_y;

            if( bottom != 0.0 )
            {
                x_dot_v = v_da0;
                y_dot_v = -(DOT( b2, da2 ) - DOT( a0, da2 ));

                x_pos = (x_dot_v * len_da2 - x_dot_y * y_dot_v) / bottom;
                y_pos = (y_dot_v * len_da0 - x_dot_y * x_dot_v) / bottom;

                weights[0] = 1.0 - x_pos - y_pos;
                weights[1] = x_pos;
                weights[2] = y_pos;
                weights[5] = 1.0;

                dist = d * d / len_an;

                *tri_case = A0 | A1 | A2 | B2;
                return( dist );
            }
        }
    }
    else if( outside_ea0 && v_da0 > 0.0 && v_da0 < len_da0 )
    {
        d = v_da0 / len_da0;
        px = a0x + d * da0x;
        py = a0y + d * da0y;
        pz = a0z + d * da0z;
        SUB( d, b2, p );

        if( DOT( d, db2 ) >= -tolerance && DOT( d, db1 ) <= tolerance )
        {
            weights[0] = 1.0 - d;
            weights[1] = d;
            weights[5] = 1.0;
            dist = DOT( d, d );
            *tri_case = A0 | A1 | B2;
            return( dist );
        }
    }
    else if( outside_ea1 && v_da1 > 0.0 && v_da1 < len_da1 )
    {
        d = v_da1 / len_da1;
        px = a1x + d * da1x;
        py = a1y + d * da1y;
        pz = a1z + d * da1z;
        SUB( d, b2, p );

        if( DOT( d, db2 ) >= -tolerance && DOT( d, db1 ) <= tolerance )
        {
            weights[1] = 1.0 - d;
            weights[2] = d;
            weights[5] = 1.0;
            dist = DOT( d, d );
            *tri_case = A1 | A2 | B2;
            return( dist );
        }
    }
    else if( outside_ea2 && v_da2 > 0.0 && v_da2 < len_da2 )
    {
        d = v_da2 / len_da2;
        px = a2x + d * da2x;
        py = a2y + d * da2y;
        pz = a2z + d * da2z;
        SUB( d, b2, p );

        if( DOT( d, db2 ) >= -tolerance && DOT( d, db1 ) <= tolerance )
        {
            weights[2] = 1.0 - d;
            weights[0] = d;
            weights[5] = 1.0;
            dist = DOT( d, d );
            *tri_case = A0 | A2 | B2;
            return( dist );
        }
    }
    else if( v_da0 <= 0.0 && v_da2 >= len_da2 )
    {
        if( b2_dot_db2 - a0_dot_db2 >= -tolerance &&
            DOT( b2, db1 ) - a0_dot_db1 <= tolerance )
        {
            weights[0] = 1.0;
            weights[5] = 1.0;
            SUB( d, b2, a0 );
            dist = DOT( d, d );
            *tri_case = A0 | B2;
            return( dist );
        }
    }
    else if( v_da1 <= 0.0 && v_da0 >= len_da0 )
    {
        if( b2_dot_db2 - a1_dot_db2 >= -tolerance &&
            DOT( b2, db1 ) - a1_dot_db1 <= tolerance )
        {
            weights[1] = 1.0;
            weights[5] = 1.0;
            SUB( d, b2, a1 );
            dist = DOT( d, d );
            *tri_case = A1 | B2;
            return( dist );
        }
    }
    else
    {
        if( b2_dot_db2 - a2_dot_db2 >= -tolerance &&
            DOT( b2, db1 ) - a2_dot_db1 <= tolerance )
        {
            weights[2] = 1.0;
            weights[5] = 1.0;
            SUB( d, b2, a2 );
            dist = DOT( d, d );
            *tri_case = A2 | B2;
            return( dist );
        }
    }

    /*--- test a0 against the triangle b0, b1, b2 */

    v_db0 = a0_dot_db0 - b0_dot_db0;
    v_db1 = a0_dot_db1 - b1_dot_db1;
    v_db2 = a0_dot_db2 - b2_dot_db2;

    v_eb0 = DOT( a0, eb0 ) - b0_dot_eb0;
    v_eb1 = DOT( a0, eb1 ) - b1_dot_eb1;
    v_eb2 = DOT( a0, eb2 ) - b2_dot_eb2;

    outside_eb0 = (v_eb0 <= 0.0);
    outside_eb1 = (v_eb1 <= 0.0);
    outside_eb2 = (v_eb2 <= 0.0);

    if( !outside_eb0 && !outside_eb1 && !outside_eb2 )
    {
        d = DOT( a0, bn ) - b0_dot_bn;
        if( d * da0_dot_bn >= -tolerance && d * da2_dot_bn <= tolerance )
        {
            x_dot_y = DOT( b2, db0 ) - b0_dot_db0;

            bottom = len_db0 * len_db2 - x_dot_y * x_dot_y;

            if( bottom != 0.0 )
            {
                x_dot_v = v_db0;
                y_dot_v = -(a0_dot_db2 - DOT( b0, db2 ));

                x_pos = (x_dot_v * len_db2 - x_dot_y * y_dot_v) / bottom;
                y_pos = (y_dot_v * len_db0 - x_dot_y * x_dot_v) / bottom;

                weights[0] = 1.0;
                weights[3] = 1.0 - x_pos - y_pos;
                weights[4] = x_pos;
                weights[5] = y_pos;

                dist = d * d / len_bn;

                *tri_case = A0 | B0 | B1 | B2;

                return( dist );
            }
        }
    }
    else if( outside_eb0 && v_db0 > 0.0 && v_db0 < len_db0 )
    {
        d = v_db0 / len_db0;
        px = b0x + d * db0x;
        py = b0y + d * db0y;
        pz = b0z + d * db0z;
        SUB( d, a0, p );

        if( DOT( d, da0 ) >= -tolerance && DOT( d, da2 ) <= tolerance )
        {
            weights[0] = 1.0;
            weights[3] = 1.0 - d;
            weights[4] = d;
            dist = DOT( d, d );
            *tri_case = A0 | B0 | B1;
            return( dist );
        }
    }
    else if( outside_eb1 && v_db1 > 0.0 && v_db1 < len_db1 )
    {
        d = v_db1 / len_db1;
        px = b1x + d * db1x;
        py = b1y + d * db1y;
        pz = b1z + d * db1z;
        SUB( d, a0, p );

        if( DOT( d, da0 ) >= -tolerance && DOT( d, da2 ) <= tolerance )
        {
            weights[0] = 1.0;
            weights[4] = 1.0 - d;
            weights[5] = d;
            dist = DOT( d, d );
            *tri_case = A0 | B1 | B2;
            return( dist );
        }
    }
    else if( outside_eb2 && v_db2 > 0.0 && v_db2 < len_db2 )
    {
        d = v_db2 / len_db2;
        px = b2x + d * db2x;
        py = b2y + d * db2y;
        pz = b2z + d * db2z;
        SUB( d, a0, p );

        if( DOT( d, da0 ) >= -tolerance && DOT( d, da2 ) <= tolerance )
        {
            weights[0] = 1.0;
            weights[5] = 1.0 - d;
            weights[3] = d;
            dist = DOT( d, d );
            *tri_case = A0 | B0 | B2;
            return( dist );
        }
    }

    /*--- test a1 against the triangle b0, b1, b2 */

    v_db0 = a1_dot_db0 - b0_dot_db0;
    v_db1 = a1_dot_db1 - b1_dot_db1;
    v_db2 = a1_dot_db2 - b2_dot_db2;

    v_eb0 = DOT( a1, eb0 ) - b0_dot_eb0;
    v_eb1 = DOT( a1, eb1 ) - b1_dot_eb1;
    v_eb2 = DOT( a1, eb2 ) - b2_dot_eb2;

    outside_eb0 = (v_eb0 <= 0.0);
    outside_eb1 = (v_eb1 <= 0.0);
    outside_eb2 = (v_eb2 <= 0.0);

    if( !outside_eb0 && !outside_eb1 && !outside_eb2 )
    {
        d = DOT( a1, bn ) - b0_dot_bn;
        if( d * da1_dot_bn >= -tolerance && d * da0_dot_bn <= tolerance )
        {
            x_dot_y = DOT( b2, db0 ) - b0_dot_db0;

            bottom = len_db0 * len_db2 - x_dot_y * x_dot_y;

            if( bottom != 0.0 )
            {
                x_dot_v = v_db0;
                y_dot_v = -(a1_dot_db2 - DOT(b0,db2));

                x_pos = (x_dot_v * len_db2 - x_dot_y * y_dot_v) / bottom;
                y_pos = (y_dot_v * len_db0 - x_dot_y * x_dot_v) / bottom;

                weights[1] = 1.0;
                weights[3] = 1.0 - x_pos - y_pos;
                weights[4] = x_pos;
                weights[5] = y_pos;

                dist = d * d / len_bn;

                *tri_case = A1 | B0 | B1 | B2;
                return( dist );
            }
        }
    }
    else if( outside_eb0 && v_db0 > 0.0 && v_db0 < len_db0 )
    {
        d = v_db0 / len_db0;
        px = b0x + d * db0x;
        py = b0y + d * db0y;
        pz = b0z + d * db0z;
        SUB( d, a1, p );

        if( DOT( d, da1 ) >= -tolerance && DOT( d, da0 ) <= tolerance )
        {
            weights[1] = 1.0;
            weights[3] = 1.0 - d;
            weights[4] = d;
            dist = DOT( d, d );
            *tri_case = A1 | B0 | B1;
            return( dist );
        }
    }
    else if( outside_eb1 && v_db1 > 0.0 && v_db1 < len_db1 )
    {
        d = v_db1 / len_db1;
        px = b1x + d * db1x;
        py = b1y + d * db1y;
        pz = b1z + d * db1z;
        SUB( d, a1, p );

        if( DOT( d, da1 ) >= -tolerance && DOT( d, da0 ) <= tolerance )
        {
            weights[1] = 1.0;
            weights[4] = 1.0 - d;
            weights[5] = d;
            dist = DOT( d, d );
            *tri_case = A1 | B1 | B2;
            return( dist );
        }
    }
    else if( outside_eb2 && v_db2 > 0.0 && v_db2 < len_db2 )
    {
        d = v_db2 / len_db2;
        px = b2x + d * db2x;
        py = b2y + d * db2y;
        pz = b2z + d * db2z;
        SUB( d, a1, p );

        if( DOT( d, da1 ) >= -tolerance && DOT( d, da0 ) <= tolerance )
        {
            weights[1] = 1.0;
            weights[5] = 1.0 - d;
            weights[3] = d;
            dist = DOT( d, d );
            *tri_case = A1 | B0 | B2;
            return( dist );
        }
    }

    /*--- test a2 against the triangle b0, b1, b2 */

    v_db0 = a2_dot_db0 - b0_dot_db0;
    v_db1 = a2_dot_db1 - b1_dot_db1;
    v_db2 = a2_dot_db2 - b2_dot_db2;

    v_eb0 = DOT( a2, eb0 ) - b0_dot_eb0;
    v_eb1 = DOT( a2, eb1 ) - b1_dot_eb1;
    v_eb2 = DOT( a2, eb2 ) - b2_dot_eb2;

    outside_eb0 = (v_eb0 <= 0.0);
    outside_eb1 = (v_eb1 <= 0.0);
    outside_eb2 = (v_eb2 <= 0.0);

    if( !outside_eb0 && !outside_eb1 && !outside_eb2 )
    {
        d = DOT( a2, bn ) - b0_dot_bn;
        if( d * da2_dot_bn >= -tolerance && d * da1_dot_bn <= tolerance )
        {
            x_dot_y = DOT( b2, db0 ) - b0_dot_db0;

            bottom = len_db0 * len_db2 - x_dot_y * x_dot_y;

            if( bottom != 0.0 )
            {
                x_dot_v = v_db0;
                y_dot_v = -(a2_dot_db2 - DOT( b0, db2 ));

                x_pos = (x_dot_v * len_db2 - x_dot_y * y_dot_v) / bottom;
                y_pos = (y_dot_v * len_db0 - x_dot_y * x_dot_v) / bottom;

                weights[2] = 1.0;
                weights[3] = 1.0 - x_pos - y_pos;
                weights[4] = x_pos;
                weights[5] = y_pos;

                dist = d * d / len_bn;

                *tri_case = A2 | B0 | B1 | B2;
                return( dist );
            }
        }
    }
    else if( outside_eb0 && v_db0 > 0.0 && v_db0 < len_db0 )
    {
        d = v_db0 / len_db0;
        px = b0x + d * db0x;
        py = b0y + d * db0y;
        pz = b0z + d * db0z;
        SUB( d, a2, p );

        if( DOT( d, da2 ) >= -tolerance && DOT( d, da1 ) <= tolerance )
        {
            weights[2] = 1.0;
            weights[3] = 1.0 - d;
            weights[4] = d;
            dist = DOT( d, d );
            *tri_case = A2 | B0 | B1;
            return( dist );
        }
    }
    else if( outside_eb1 && v_db1 > 0.0 && v_db1 < len_db1 )
    {
        d = v_db1 / len_db1;
        px = b1x + d * db1x;
        py = b1y + d * db1y;
        pz = b1z + d * db1z;
        SUB( d, a2, p );

        if( DOT( d, da2 ) >= -tolerance && DOT( d, da1 ) <= tolerance )
        {
            weights[2] = 1.0;
            weights[4] = 1.0 - d;
            weights[5] = d;
            dist = DOT( d, d );
            *tri_case = A2 | B1 | B2;
            return( dist );
        }
    }
    else if( outside_eb2 && v_db2 > 0.0 && v_db2 < len_db2 )
    {
        d = v_db2 / len_db2;
        px = b2x + d * db2x;
        py = b2y + d * db2y;
        pz = b2z + d * db2z;
        SUB( d, a2, p );

        if( DOT( d, da2 ) >= -tolerance && DOT( d, da1 ) <= tolerance )
        {
            weights[2] = 1.0;
            weights[5] = 1.0 - d;
            weights[3] = d;
            dist = DOT( d, d );
            *tri_case = A2 | B0 | B2;
            return( dist );
        }
    }

    /*--- a0-a1 vs b0-b1 */

    c = a0_dot_da0 - b0_dot_da0;
    d = a1_dot_db0 - a0_dot_db0;
    f = a0_dot_db0 - b0_dot_db0;

    denom = d * d - len_da0 * len_db0;
    s = len_db0 * c - f * d;
    t = d * c - len_da0 * f;

    if( denom > 0.0 && s > 0.0 && s < denom && t > 0.0 && t < denom ||
        denom < 0.0 && s < 0.0 && s > denom && t < 0.0 && t > denom )
    {
        s /= denom;
        t /= denom;

        dx = b0x + t * db0x - (a0x + s * da0x);
        dy = b0y + t * db0y - (a0y + s * da0y);
        dz = b0z + t * db0z - (a0z + s * da0z);

        if( DOT( d, ea0 ) <= tolerance && DOT( d, eb0 ) >= -tolerance )
        {
            weights[0] = 1.0 - s;
            weights[1] = s;
            weights[3] = 1.0 - t;
            weights[4] = t;
            dist = dx * dx + dy * dy + dz * dz;
            *tri_case = A0 | A1 | B0 | B1;
            return( dist );
        }
    }

    /*--- a0-a1 vs b1-b2 */

    c = a0_dot_da0 - b1_dot_da0;
    d = a1_dot_db1 - a0_dot_db1;
    f = a0_dot_db1 - b1_dot_db1;

    denom = d * d - len_da0 * len_db1;
    s = len_db1 * c - f * d;
    t = d * c - len_da0 * f;

    if( denom > 0.0 && s > 0.0 && s < denom && t > 0.0 && t < denom ||
        denom < 0.0 && s < 0.0 && s > denom && t < 0.0 && t > denom )
    {
        s /= denom;
        t /= denom;

        dx = b1x + t * db1x - (a0x + s * da0x);
        dy = b1y + t * db1y - (a0y + s * da0y);
        dz = b1z + t * db1z - (a0z + s * da0z);

        if( DOT( d, ea0 ) <= tolerance && DOT( d, eb1 ) >= -tolerance )
        {
            weights[0] = 1.0 - s;
            weights[1] = s;
            weights[4] = 1.0 - t;
            weights[5] = t;
            dist = dx * dx + dy * dy + dz * dz;
            *tri_case = A0 | A1 | B1 | B2;
            return( dist );
        }
    }

    /*--- a0-a1 vs b2-b0 */

    c = a0_dot_da0 - b2_dot_da0;
    d = a1_dot_db2 - a0_dot_db2;
    f = a0_dot_db2 - b2_dot_db2;

    denom = d * d - len_da0 * len_db2;
    s = len_db2 * c - f * d;
    t = d * c - len_da0 * f;

    if( denom > 0.0 && s > 0.0 && s < denom && t > 0.0 && t < denom ||
        denom < 0.0 && s < 0.0 && s > denom && t < 0.0 && t > denom )
    {
        s /= denom;
        t /= denom;

        dx = b2x + t * db2x - (a0x + s * da0x);
        dy = b2y + t * db2y - (a0y + s * da0y);
        dz = b2z + t * db2z - (a0z + s * da0z);

        if( DOT( d, ea0 ) <= tolerance && DOT( d, eb2 ) >= -tolerance )
        {
            weights[0] = 1.0 - s;
            weights[1] = s;
            weights[3] = t;
            weights[5] = 1.0 - t;
            dist = dx * dx + dy * dy + dz * dz;
            *tri_case = A0 | A1 | B0 | B2;
            return( dist );
        }
    }

    /*--- a1-a2 vs b0-b1 */

    c = a1_dot_da1 - b0_dot_da1;
    d = a2_dot_db0 - a1_dot_db0;
    f = a1_dot_db0 - b0_dot_db0;

    denom = d * d - len_da1 * len_db0;
    s = len_db0 * c - f * d;
    t = d * c - len_da1 * f;

    if( denom > 0.0 && s > 0.0 && s < denom && t > 0.0 && t < denom ||
        denom < 0.0 && s < 0.0 && s > denom && t < 0.0 && t > denom )
    {
        s /= denom;
        t /= denom;

        dx = b0x + t * db0x - (a1x + s * da1x);
        dy = b0y + t * db0y - (a1y + s * da1y);
        dz = b0z + t * db0z - (a1z + s * da1z);

        if( DOT( d, ea1 ) <= tolerance && DOT( d, eb0 ) >= -tolerance )
        {
            weights[1] = 1.0 - s;
            weights[2] = s;
            weights[3] = 1.0 - t;
            weights[4] = t;
            dist = dx * dx + dy * dy + dz * dz;
            *tri_case = A1 | A2 | B0 | B1;
            return( dist );
        }
    }

    /*--- a1-a2 vs b1-b2 */

    c = a1_dot_da1 - b1_dot_da1;
    d = a2_dot_db1 - a1_dot_db1;
    f = a1_dot_db1 - b1_dot_db1;

    denom = d * d - len_da1 * len_db1;
    s = len_db1 * c - f * d;
    t = d * c - len_da1 * f;

    if( denom > 0.0 && s > 0.0 && s < denom && t > 0.0 && t < denom ||
        denom < 0.0 && s < 0.0 && s > denom && t < 0.0 && t > denom )
    {
        s /= denom;
        t /= denom;

        dx = b1x + t * db1x - (a1x + s * da1x);
        dy = b1y + t * db1y - (a1y + s * da1y);
        dz = b1z + t * db1z - (a1z + s * da1z);

        if( DOT( d, ea1 ) <= tolerance && DOT( d, eb1 ) >= -tolerance )
        {
            weights[1] = 1.0 - s;
            weights[2] = s;
            weights[4] = 1.0 - t;
            weights[5] = t;
            dist = dx * dx + dy * dy + dz * dz;
            *tri_case = A1 | A2 | B1 | B2;
            return( dist );
        }
    }

    /*--- a1-a2 vs b2-b0 */

    c = a1_dot_da1 - b2_dot_da1;
    d = a2_dot_db2 - a1_dot_db2;
    f = a1_dot_db2 - b2_dot_db2;

    denom = d * d - len_da1 * len_db2;
    s = len_db2 * c - f * d;
    t = d * c - len_da1 * f;

    if( denom > 0.0 && s > 0.0 && s < denom && t > 0.0 && t < denom ||
        denom < 0.0 && s < 0.0 && s > denom && t < 0.0 && t > denom )
    {
        s /= denom;
        t /= denom;

        dx = b2x + t * db2x - (a1x + s * da1x);
        dy = b2y + t * db2y - (a1y + s * da1y);
        dz = b2z + t * db2z - (a1z + s * da1z);

        if( DOT( d, ea1 ) <= tolerance && DOT( d, eb2 ) >= -tolerance )
        {
            weights[1] = 1.0 - s;
            weights[2] = s;
            weights[3] = t;
            weights[5] = 1.0 - t;
            dist = dx * dx + dy * dy + dz * dz;
            *tri_case = A1 | A2 | B0 | B2;
            return( dist );
        }
    }

    /*--- a2-a0 vs b0-b1 */

    c = a2_dot_da2 - b0_dot_da2;
    d = a0_dot_db0 - a2_dot_db0;
    f = a2_dot_db0 - b0_dot_db0;

    denom = d * d - len_da2 * len_db0;
    s = len_db0 * c - f * d;
    t = d * c - len_da2 * f;

    if( denom > 0.0 && s > 0.0 && s < denom && t > 0.0 && t < denom ||
        denom < 0.0 && s < 0.0 && s > denom && t < 0.0 && t > denom )
    {
        s /= denom;
        t /= denom;

        dx = b0x + t * db0x - (a2x + s * da2x);
        dy = b0y + t * db0y - (a2y + s * da2y);
        dz = b0z + t * db0z - (a2z + s * da2z);

        if( DOT( d, ea2 ) <= tolerance && DOT( d, eb0 ) >= -tolerance )
        {
            weights[0] = s;
            weights[2] = 1.0 - s;
            weights[3] = 1.0 - t;
            weights[4] = t;
            dist = dx * dx + dy * dy + dz * dz;
            *tri_case = A0 | A2 | B0 | B1;
            return( dist );
        }
    }

    /*--- a2-a0 vs b1-b2 */

    c = a2_dot_da2 - b1_dot_da2;
    d = a0_dot_db1 - a2_dot_db1;
    f = a2_dot_db1 - b1_dot_db1;

    denom = d * d - len_da2 * len_db1;
    s = len_db1 * c - f * d;
    t = d * c - len_da2 * f;

    if( denom > 0.0 && s > 0.0 && s < denom && t > 0.0 && t < denom ||
        denom < 0.0 && s < 0.0 && s > denom && t < 0.0 && t > denom )
    {
        s /= denom;
        t /= denom;

        dx = b1x + t * db1x - (a2x + s * da2x);
        dy = b1y + t * db1y - (a2y + s * da2y);
        dz = b1z + t * db1z - (a2z + s * da2z);

        if( DOT( d, ea2 ) <= tolerance && DOT( d, eb1 ) >= -tolerance )
        {
            weights[0] = s;
            weights[2] = 1.0 - s;
            weights[4] = 1.0 - t;
            weights[5] = t;
            dist = dx * dx + dy * dy + dz * dz;
            *tri_case = A0 | A2 | B1 | B2;
            return( dist );
        }
    }

    /*--- a2-a0 vs b2-b0 */

    c = a2_dot_da2 - b2_dot_da2;
    d = a0_dot_db2 - a2_dot_db2;
    f = a2_dot_db2 - b2_dot_db2;

    denom = d * d - len_da2 * len_db2;
    s = len_db2 * c - f * d;
    t = d * c - len_da2 * f;

    if( denom > 0.0 && s > 0.0 && s < denom && t > 0.0 && t < denom ||
        denom < 0.0 && s < 0.0 && s > denom && t < 0.0 && t > denom )
    {
        s /= denom;
        t /= denom;

        dx = b2x + t * db2x - (a2x + s * da2x);
        dy = b2y + t * db2y - (a2y + s * da2y);
        dz = b2z + t * db2z - (a2z + s * da2z);

        if( DOT( d, ea2 ) <= tolerance && DOT( d, eb2 ) >= -tolerance )
        {
            weights[0] = s;
            weights[2] = 1.0 - s;
            weights[3] = t;
            weights[5] = 1.0 - t;
            dist = dx * dx + dy * dy + dz * dz;
            *tri_case = A0 | A2 | B0 | B2;
            return( dist );
        }
    }

    /*--- otherwise, they must intersect */

    weights[0] = 0.0;
    weights[1] = 0.0;
    weights[2] = 0.0;
    weights[3] = 0.0;
    weights[4] = 0.0;
    weights[5] = 0.0;
    *tri_case = 0;

    return( 0.0 );
}

Real   Weights[6];

#ifdef STATS
static  int  vv = 0;
static  int  ve = 0;
static  int  ee = 0;
static  int  vp = 0;
static  int  other = 0;
static  int  count = 0;
#endif  /* STATS */


public  Real  sq_triangle_triangle_dist(
    Real   a0[],
    Real   a1[],
    Real   a2[],
    Real   b0[],
    Real   b1[],
    Real   b2[],
    int    *which_case )
{
    int   i, tri_case;
    Real  dist;
#ifdef STATS
    int   n1, n2;
#endif  /* STATS */

    if( which_case != NULL && *which_case > 0 &&
        test_distance_with_known_case( a0, a1, a2, b0, b1, b2, *which_case,
                                       &dist ) )
    {
        return( dist );
    }

    dist = sq_triangle_triangle_dist_and_weights( a0, a1, a2, b0, b1, b2,
                                                  Weights, &tri_case );

    if( which_case != NULL )
        *which_case = tri_case;

#ifdef STATS
    n1 = 0;
    for_less( i, 0, 3 )
        if( Weights[i] != 0.0 )  ++n1;

    n2 = 0;
    for_less( i, 3, 6 )
        if( Weights[i] != 0.0 )  ++n2;

    if( n1 == 1 && n2 == 1 )
        ++vv;
    else if( n1 == 1 && n2 == 2 || n1 == 2 && n2 == 1 )
        ++ve;
    else if( n1 == 2 && n2 == 2 )
        ++ee;
    else if( n1 == 1 && n2 == 3 || n1 == 3 && n2 == 1 )
        ++vp;
    else
        ++other;

    ++count;
    if( count % 1000000 == 0 )
    {
        print( "vv %5.1f    ve %5.1f     ee %5.1f   vp %5.1f   other  %5.1f\n",
                100.0 * (Real) vv / (Real) count,
                100.0 * (Real) ve / (Real) count,
                100.0 * (Real) ee / (Real) count,
                100.0 * (Real) vp / (Real) count,
                100.0 * (Real) other / (Real) count );
    }
#endif  /* STATS */

    return( dist );
}

private  void  compute_point_point_derivative(
    Real   p0[],
    Real   deriv_0[],
    Real   p1[],
    Real   deriv_1[] )
{
    Real  dx, dy, dz;

    dx = p1[0] - p0[0];
    dy = p1[1] - p0[1];
    dz = p1[2] - p0[2];

    deriv_0[0] = 2.0 * dx * -1.0;
    deriv_0[1] = 2.0 * dy * -1.0;
    deriv_0[2] = 2.0 * dz * -1.0;

    deriv_1[0] = 2.0 * dx * 1.0;
    deriv_1[1] = 2.0 * dy * 1.0;
    deriv_1[2] = 2.0 * dz * 1.0;
}

private  void  compute_point_edge_derivative(
    Real   point[],
    Real   deriv_point[],
    Real   edge0[],
    Real   edge0_deriv[],
    Real   edge1[],
    Real   edge1_deriv[] )
{
    Real     x0, y0, z0, x1, y1, z1, px, py, pz;
    Real     v01x, v01y, v01z;
    Real     pv0x, pv0y, pv0z;
    Real     rratio;
    Real     rv_dot_e, re_dot_e;

    x0 = edge0[0];
    y0 = edge0[1];
    z0 = edge0[2];
    x1 = edge1[0];
    y1 = edge1[1];
    z1 = edge1[2];

    px = point[0];
    py = point[1];
    pz = point[2];

    v01x = x1 - x0;
    v01y = y1 - y0;
    v01z = z1 - z0;

    pv0x = px - x0;
    pv0y = py - y0;
    pv0z = pz - z0;

    rv_dot_e = pv0x * v01x + pv0y * v01y + pv0z * v01z;
    re_dot_e = v01x * v01x + v01y * v01y + v01z * v01z;

    rratio = rv_dot_e / re_dot_e;

    edge0_deriv[0] = 2.0 * pv0x * (-1.0) -
                     2.0 * rv_dot_e / re_dot_e * (-pv0x - v01x) -
                     rv_dot_e * rv_dot_e / re_dot_e / re_dot_e * 2.0 *v01x;
    edge0_deriv[1] = 2.0 * pv0y * (-1.0) -
                     2.0 * rv_dot_e / re_dot_e * (-pv0y - v01y) -
                     rv_dot_e * rv_dot_e / re_dot_e / re_dot_e * 2.0 *v01y;
    edge0_deriv[2] = 2.0 * pv0z * (-1.0) -
                     2.0 * rv_dot_e / re_dot_e * (-pv0z - v01z) -
                     rv_dot_e * rv_dot_e / re_dot_e / re_dot_e * 2.0 *v01z;

    edge1_deriv[0] = -2.0 * rratio * pv0x +
                     rv_dot_e * rratio / re_dot_e * 2.0 * v01x;
    edge1_deriv[1] = -2.0 * rratio * pv0y +
                     rv_dot_e * rratio / re_dot_e * 2.0 * v01y;
    edge1_deriv[2] = -2.0 * rratio * pv0z +
                     rv_dot_e * rratio / re_dot_e * 2.0 * v01z;

    deriv_point[0] = 2.0*pv0x -2.0*rv_dot_e * v01x / re_dot_e;
    deriv_point[1] = 2.0*pv0y -2.0*rv_dot_e * v01y / re_dot_e;
    deriv_point[2] = 2.0*pv0z -2.0*rv_dot_e * v01z / re_dot_e;
}

private  void  compute_point_plane_derivative(
    Real   point[],
    Real   deriv_point[],
    Real   *tri_points[],
    Real   *tri_derivs[] )
{
    Real     x0, y0, z0, x1, y1, z1, x2, y2, z2, px, py, pz;
    Real     v01x, v01y, v01z;
    Real     v12x, v12y, v12z;
    Real     v20x, v20y, v20z;
    Real     pv0x, pv0y, pv0z;
    Real     pv1x, pv1y, pv1z;
    Real     pv2x, pv2y, pv2z;
    Real     nx, ny, nz;
    Real     v_dot_n, n_dot_n;
    Real     vn_x0, vn_y0, vn_z0, vn_x1, vn_y1, vn_z1;
    Real     vn_x2, vn_y2, vn_z2, vn_x, vn_y, vn_z;
    Real     nn_x0, nn_y0, nn_z0, nn_x1, nn_y1, nn_z1;
    Real     nn_x2, nn_y2, nn_z2;

    x0 = tri_points[0][0];
    y0 = tri_points[0][1];
    z0 = tri_points[0][2];
    x1 = tri_points[1][0];
    y1 = tri_points[1][1];
    z1 = tri_points[1][2];
    x2 = tri_points[2][0];
    y2 = tri_points[2][1];
    z2 = tri_points[2][2];

    px = point[0];
    py = point[1];
    pz = point[2];

    v01x = x1 - x0;
    v01y = y1 - y0;
    v01z = z1 - z0;
    v12x = x2 - x1;
    v12y = y2 - y1;
    v12z = z2 - z1;
    v20x = x0 - x2;
    v20y = y0 - y2;
    v20z = z0 - z2;

    nx = v20y * v01z - v20z * v01y;
    ny = v20z * v01x - v20x * v01z;
    nz = v20x * v01y - v20y * v01x;

    pv0x = px - x0;
    pv0y = py - y0;
    pv0z = pz - z0;
    pv1x = px - x1;
    pv1y = py - y1;
    pv1z = pz - z1;
    pv2x = px - x2;
    pv2y = py - y2;
    pv2z = pz - z2;

    v_dot_n = pv0x * nx + pv0y * ny + pv0z * nz;
    n_dot_n = nx * nx + ny * ny + nz * nz;

    vn_x0 = -nx + pv0y * v12z + pv0z * -v12y;
    vn_y0 = -ny + pv0z * v12x + pv0x * -v12z;
    vn_z0 = -nz + pv0x * v12y + pv0y * -v12x;

    vn_x1 = -nx + pv1y * v20z + pv1z * -v20y;
    vn_y1 = -ny + pv1z * v20x + pv1x * -v20z;
    vn_z1 = -nz + pv1x * v20y + pv1y * -v20x;

    vn_x2 = -nx + pv2y * v01z + pv2z * -v01y;
    vn_y2 = -ny + pv2z * v01x + pv2x * -v01z;
    vn_z2 = -nz + pv2x * v01y + pv2y * -v01x;

    vn_x = v20y * v01z - v20z * v01y;
    vn_y = v20z * v01x - v20x * v01z;
    vn_z = v20x * v01y - v20y * v01x;

    nn_x0 = 2.0 * ny * v12z - 2.0 * nz * v12y;
    nn_y0 = 2.0 * nz * v12x - 2.0 * nx * v12z;
    nn_z0 = 2.0 * nx * v12y - 2.0 * ny * v12x;

    nn_x1 = 2.0 * ny * v20z - 2.0 * nz * v20y;
    nn_y1 = 2.0 * nz * v20x - 2.0 * nx * v20z;
    nn_z1 = 2.0 * nx * v20y - 2.0 * ny * v20x;

    nn_x2 = 2.0 * ny * v01z - 2.0 * nz * v01y;
    nn_y2 = 2.0 * nz * v01x - 2.0 * nx * v01z;
    nn_z2 = 2.0 * nx * v01y - 2.0 * ny * v01x;

    tri_derivs[0][0] = 2.0 * v_dot_n / n_dot_n * vn_x0 -
                       v_dot_n * v_dot_n / n_dot_n / n_dot_n * nn_x0;
    tri_derivs[0][1] = 2.0 * v_dot_n / n_dot_n * vn_y0 -
                       v_dot_n * v_dot_n / n_dot_n / n_dot_n * nn_y0;
    tri_derivs[0][2] = 2.0 * v_dot_n / n_dot_n * vn_z0 -
                       v_dot_n * v_dot_n / n_dot_n / n_dot_n * nn_z0;

    tri_derivs[1][0] = 2.0 * v_dot_n / n_dot_n * vn_x1 -
                       v_dot_n * v_dot_n / n_dot_n / n_dot_n * nn_x1;
    tri_derivs[1][1] = 2.0 * v_dot_n / n_dot_n * vn_y1 -
                       v_dot_n * v_dot_n / n_dot_n / n_dot_n * nn_y1;
    tri_derivs[1][2] = 2.0 * v_dot_n / n_dot_n * vn_z1 -
                       v_dot_n * v_dot_n / n_dot_n / n_dot_n * nn_z1;

    tri_derivs[2][0] = 2.0 * v_dot_n / n_dot_n * vn_x2 -
                       v_dot_n * v_dot_n / n_dot_n / n_dot_n * nn_x2;
    tri_derivs[2][1] = 2.0 * v_dot_n / n_dot_n * vn_y2 -
                       v_dot_n * v_dot_n / n_dot_n / n_dot_n * nn_y2;
    tri_derivs[2][2] = 2.0 * v_dot_n / n_dot_n * vn_z2 -
                       v_dot_n * v_dot_n / n_dot_n / n_dot_n * nn_z2;

    deriv_point[0] = 2.0 * v_dot_n / n_dot_n * vn_x;
    deriv_point[1] = 2.0 * v_dot_n / n_dot_n * vn_y;
    deriv_point[2] = 2.0 * v_dot_n / n_dot_n * vn_z;
}

private  void  compute_edge_edge_derivative(
    Real   a0[],
    Real   deriv_a0[],
    Real   a1[],
    Real   deriv_a1[],
    Real   b0[],
    Real   deriv_b0[],
    Real   b1[],
    Real   deriv_b1[] )
{
    int   dim, d;
    Real  *deriv[4];
    Real  d_v_dot_d[4][N_DIMENSIONS];
    Real  d_d_dot_d[4][N_DIMENSIONS];
    Real  v_dot_d, d_dot_d;
    Real  vx, vy, vz, v0x, v0y, v0z, v1x, v1y, v1z, nx, ny, nz;

    deriv[0] = deriv_a0;
    deriv[1] = deriv_a1;
    deriv[2] = deriv_b0;
    deriv[3] = deriv_b1;

    vx = a0[0] - b0[0];
    vy = a0[1] - b0[1];
    vz = a0[2] - b0[2];

    v0x = a1[0] - a0[0];
    v0y = a1[1] - a0[1];
    v0z = a1[2] - a0[2];

    v1x = b1[0] - b0[0];
    v1y = b1[1] - b0[1];
    v1z = b1[2] - b0[2];

    nx = v0y * v1z - v0z * v1y;
    ny = v0z * v1x - v0x * v1z;
    nz = v0x * v1y - v0y * v1x;

    v_dot_d = vx * nx + vy * ny + vz * nz;
    d_dot_d = nx * nx + ny * ny + nz * nz;

    d_v_dot_d[0][0] = nx + vy * v1z - vz * v1y;
    d_v_dot_d[0][1] = ny + vz * v1x - vx * v1z;
    d_v_dot_d[0][2] = nz + vx * v1y - vy * v1x;

    d_v_dot_d[1][0] = vz * v1y - vy * v1z;
    d_v_dot_d[1][1] = vx * v1z - vz * v1x;
    d_v_dot_d[1][2] = vy * v1x - vx * v1y;

    d_v_dot_d[2][0] = -nx + vz * v0y - vy * v0z;
    d_v_dot_d[2][1] = -ny + vx * v0z - vz * v0x;
    d_v_dot_d[2][2] = -nz + vy * v0x - vx * v0y;

    d_v_dot_d[3][0] = vy * v0z - vz * v0y;
    d_v_dot_d[3][1] = vz * v0x - vx * v0z;
    d_v_dot_d[3][2] = vx * v0y - vy * v0x;

    d_d_dot_d[0][0] = 2.0 * ny * v1z - 2.0 * nz * v1y;
    d_d_dot_d[0][1] = 2.0 * nz * v1x - 2.0 * nx * v1z;
    d_d_dot_d[0][2] = 2.0 * nx * v1y - 2.0 * ny * v1x;

    d_d_dot_d[1][0] = - 2.0 * ny * v1z + 2.0 * nz * v1y;
    d_d_dot_d[1][1] = - 2.0 * nz * v1x + 2.0 * nx * v1z;
    d_d_dot_d[1][2] = - 2.0 * nx * v1y + 2.0 * ny * v1x;

    d_d_dot_d[2][0] = - 2.0 * ny * v0z + 2.0 * nz * v0y;
    d_d_dot_d[2][1] = - 2.0 * nz * v0x + 2.0 * nx * v0z;
    d_d_dot_d[2][2] = - 2.0 * nx * v0y + 2.0 * ny * v0x;

    d_d_dot_d[3][0] = 2.0 * ny * v0z - 2.0 * nz * v0y;
    d_d_dot_d[3][1] = 2.0 * nz * v0x - 2.0 * nx * v0z;
    d_d_dot_d[3][2] = 2.0 * nx * v0y - 2.0 * ny * v0x;

    for_less( dim, 0, N_DIMENSIONS )
    {
        for_less( d, 0, 4 )
        {
            deriv[d][dim] = 2.0 * v_dot_d / d_dot_d * d_v_dot_d[d][dim] -
                            v_dot_d * v_dot_d / d_dot_d / d_dot_d *
                            d_d_dot_d[d][dim];
        }
    }

#ifdef DEBUG
 {
    Real test_dist, s, t, v_dot_v, true_dist;

    v_dot_v = vx * vx + vy * vy + vz * vz;
    true_dist =  line_line_dist( a0, a1, b0, b1, &s, &t );
    test_dist = v_dot_d * v_dot_d / d_dot_d;

    if( !numerically_close( true_dist, test_dist, 1.0e-3 ) )
        handle_internal_error( "true dist test dist" );
 }
#endif
}

private  void  compute_edge_plane_derivative(
    Real   a0[],
    Real   deriv_a0[],
    Real   a1[],
    Real   deriv_a1[],
    Real   *tri_points[],
    Real   *tri_derivs[] )
{
    int   dim;

    /*  nothing to do, edge must intersect plane */

    for_less( dim, 0, 3 )
    {
        deriv_a0[dim] = 0.0;
        deriv_a1[dim] = 0.0;
        tri_derivs[0][dim] = 0.0;
        tri_derivs[1][dim] = 0.0;
        tri_derivs[2][dim] = 0.0;
    }
}

private  void  compute_edge_derivative(
    Real   weight0,
    Real   edge0[],
    Real   edge0_deriv[],
    Real   weight1,
    Real   edge1[],
    Real   edge1_deriv[],
    Real   weights[],
    Real   *tri_points[],
    Real   *tri_derivs[] )
{
    int    w, n_non_null_b, which_b[3];

    if( !numerically_close( weight0 + weight1, 1.0, 1.0e-5 ) )
        handle_internal_error( "compute_edge_derivative: weights:" );

    n_non_null_b = 0;
    for_less( w, 0, 3 )
    {
        if( weights[w] != 0.0 )
        {
            which_b[n_non_null_b] = w;
            ++n_non_null_b;
        }
    }

    if( n_non_null_b == 0 )
        handle_internal_error( "compute_edge_derivative" );
    else if( n_non_null_b == 1 )
    {
        if( weights[which_b[0]] != 1.0 )
            handle_internal_error( "compute_edge_derivative" );

        compute_point_edge_derivative( tri_points[which_b[0]],
                                       tri_derivs[which_b[0]],
                                       edge0, edge0_deriv,
                                       edge1, edge1_deriv );
    }
    else if( n_non_null_b == 2 )
    {
        compute_edge_edge_derivative( edge0, edge0_deriv,
                                      edge1, edge1_deriv,
                                      tri_points[which_b[0]],
                                      tri_derivs[which_b[0]],
                                      tri_points[which_b[1]],
                                      tri_derivs[which_b[1]] );
    }
    else
    {
        compute_edge_plane_derivative( edge0, edge0_deriv,
                                       edge1, edge1_deriv,
                                       tri_points, tri_derivs );
    }
}

private  void  compute_point_derivative(
    Real   point[],
    Real   point_deriv[],
    Real   weights[],
    Real   *tri_points[],
    Real   *tri_derivs[] )
{
    int    w, n_non_null_b, which_b[3];

    n_non_null_b = 0;
    for_less( w, 0, 3 )
    {
        if( weights[w] != 0.0 )
        {
            which_b[n_non_null_b] = w;
            ++n_non_null_b;
        }
    }

    if( n_non_null_b == 0 )
        handle_internal_error( "compute_point_derivative" );
    else if( n_non_null_b == 1 )
    {
        if( weights[which_b[0]] != 1.0 )
            handle_internal_error( "compute_point_derivative" );

        compute_point_point_derivative( point, point_deriv,
                                        tri_points[which_b[0]],
                                        tri_derivs[which_b[0]] );
    }
    else if( n_non_null_b == 2 )
    {
        compute_point_edge_derivative( point, point_deriv,
                                       tri_points[which_b[0]],
                                       tri_derivs[which_b[0]],
                                       tri_points[which_b[1]],
                                       tri_derivs[which_b[1]] );
    }
    else
    {
        compute_point_plane_derivative( point, point_deriv,
                                        tri_points, tri_derivs );
    }
}

public  void  sq_triangle_triangle_dist_deriv(
    Real   a0[],
    Real   a1[],
    Real   a2[],
    Real   b0[],
    Real   b1[],
    Real   b2[],
    Real   deriv_a0[],
    Real   deriv_a1[],
    Real   deriv_a2[],
    Real   deriv_b0[],
    Real   deriv_b1[],
    Real   deriv_b2[] )
{
    int    dim, w, n_non_null_a, n_non_null_b, which_a[3], which_b[3];
    int    tri_case;
    Real   weights[6];
    Real   *points[6], *derivs[6];

    (void) sq_triangle_triangle_dist_and_weights( a0, a1, a2, b0, b1, b2,
                                                  weights, &tri_case );

    points[0] = a0;
    points[1] = a1;
    points[2] = a2;
    points[3] = b0;
    points[4] = b1;
    points[5] = b2;
    derivs[0] = deriv_a0;
    derivs[1] = deriv_a1;
    derivs[2] = deriv_a2;
    derivs[3] = deriv_b0;
    derivs[4] = deriv_b1;
    derivs[5] = deriv_b2;

    for_less( w, 0, 6 )
    for_less( dim, 0, N_DIMENSIONS )
        derivs[w][dim] = 0.0;

    n_non_null_a = 0;
    for_less( w, 0, 3 )
    {
        if( weights[w] != 0.0 )
        {
            which_a[n_non_null_a] = w;
            ++n_non_null_a;
        }
    }

    if( n_non_null_a == 1 )
    {
        if( weights[which_a[0]] != 1.0 )
            handle_internal_error( "compute_point_derivative: a" );
        compute_point_derivative( points[which_a[0]], derivs[which_a[0]],
                                  &weights[3], &points[3], &derivs[3] );
    }
    else if( n_non_null_a == 2 )
    {
        compute_edge_derivative( weights[which_a[0]],
                                 points[which_a[0]], derivs[which_a[0]],
                                 weights[which_a[1]],
                                 points[which_a[1]], derivs[which_a[1]],
                                 &weights[3], &points[3], &derivs[3] );
    }
    else if( n_non_null_a == 3 )
    {
        n_non_null_b = 0;
        for_less( w, 3, 6 )
        {
            if( weights[w] != 0.0 )
            {
                which_b[n_non_null_b] = w;
                ++n_non_null_b;
            }
        }

        if( n_non_null_b == 1 )
        {
            compute_point_plane_derivative(
                                     points[which_b[0]], derivs[which_b[0]],
                                     &points[0], &derivs[0] );
        }
        else if( n_non_null_b == 2 )
        {
            compute_edge_plane_derivative( 
                                points[which_b[0]], derivs[which_b[0]],
                                points[which_b[1]], derivs[which_b[1]],
                                &points[0], &derivs[0] );
        }
        else
            handle_internal_error( "sq_triangle_triangle_dist_deriv: b" );
    }
    else if( n_non_null_a == 0 )
    {
    }
}


private  Real  sq_triangle_point_dist_and_weights(
    Real            a0[],
    Real            a1[],
    Real            a2[],
    Real            p[],
    Real            weights[] )
{
    Real      a0x, a0y, a0z, a1x, a1y, a1z, a2x, a2y, a2z;
    Real      da0x, da0y, da0z, da1x, da1y, da1z, da2x, da2y, da2z;
    Real      anx, any, anz;
    Real      x_dot_y, x_dot_v, y_dot_v;
    Real      bottom, x_pos, y_pos;
    Real      dx, dy, dz, ix, iy, iz;
    Real      d, dist;
    Real      len_an;
    Real      v_da0, v_da1, v_da2, v_ea0, v_ea1, v_ea2;
    BOOLEAN   outside_ea0, outside_ea1, outside_ea2;
    Real      ea0x, ea0y, ea0z, ea1x, ea1y, ea1z, ea2x, ea2y, ea2z;
    Real      a0_dot_da0, a1_dot_da1, a2_dot_da2;
    Real      len_da0, len_da1, len_da2;
    Real      px, py, pz;

    weights[0] = 0.0;
    weights[1] = 0.0;
    weights[2] = 0.0;

    a0x = a0[X];     a0y = a0[Y];     a0z = a0[Z];
    a1x = a1[X];     a1y = a1[Y];     a1z = a1[Z];
    a2x = a2[X];     a2y = a2[Y];     a2z = a2[Z];

    px = p[X];     py = p[Y];     pz = p[Z];

    SUB( da0, a1, a0 );
    SUB( da1, a2, a1 );
    SUB( da2, a0, a2 );

    CROSS( an, da2, da0 );

    CROSS( ea0, an, da0 );
    CROSS( ea1, an, da1 );
    CROSS( ea2, an, da2 );

    a0_dot_da0 = DOT( a0, da0 );
    a1_dot_da1 = DOT( a1, da1 );
    a2_dot_da2 = DOT( a2, da2 );

    len_da0 = DOT( da0, da0 );
    len_da1 = DOT( da1, da1 );
    len_da2 = DOT( da2, da2 );

    len_an = DOT( an, an );

    /*--- test b0 against the triangle a0, a1, a2 */

    v_da0 = DOT( p, da0 ) - a0_dot_da0;
    v_da1 = DOT( p, da1 ) - a1_dot_da1;
    v_da2 = DOT( p, da2 ) - a2_dot_da2;

    v_ea0 = DOT( p, ea0 ) - DOT( a0, ea0 );
    v_ea1 = DOT( p, ea1 ) - DOT( a1, ea1 );
    v_ea2 = DOT( p, ea2 ) - DOT( a2, ea2 );

    outside_ea0 = (v_ea0 <= 0.0);
    outside_ea1 = (v_ea1 <= 0.0);
    outside_ea2 = (v_ea2 <= 0.0);

    if( !outside_ea0 && !outside_ea1 && !outside_ea2 )
    {
        d = DOT( p, an ) - DOT( a0, an );
        x_dot_y = DOT( a2, da0 ) - a0_dot_da0;
        bottom = len_da0 * len_da2 - x_dot_y * x_dot_y;

        if( bottom != 0.0 )
        {
            x_dot_v = v_da0;
            y_dot_v = -(DOT( p, da2 ) - DOT( a0, da2 ));

            x_pos = (x_dot_v * len_da2 - x_dot_y * y_dot_v) / bottom;
            y_pos = (y_dot_v * len_da0 - x_dot_y * x_dot_v) / bottom;

            weights[0] = 1.0 - x_pos - y_pos;
            weights[1] = x_pos;
            weights[2] = y_pos;

            dist = d * d / len_an;

            return( dist );
        }
    }
    else if( outside_ea0 && v_da0 > 0.0 && v_da0 < len_da0 )
    {
        d = v_da0 / len_da0;
        ix = a0x + d * da0x;
        iy = a0y + d * da0y;
        iz = a0z + d * da0z;
        SUB( d, p, i );

        weights[0] = 1.0 - d;
        weights[1] = d;
        dist = DOT( d, d );
        return( dist );
    }
    else if( outside_ea1 && v_da1 > 0.0 && v_da1 < len_da1 )
    {
        d = v_da1 / len_da1;
        ix = a1x + d * da1x;
        iy = a1y + d * da1y;
        iz = a1z + d * da1z;
        SUB( d, p, i );

        weights[1] = 1.0 - d;
        weights[2] = d;
        dist = DOT( d, d );
        return( dist );
    }
    else if( outside_ea2 && v_da2 > 0.0 && v_da2 < len_da2 )
    {
        d = v_da2 / len_da2;
        ix = a2x + d * da2x;
        iy = a2y + d * da2y;
        iz = a2z + d * da2z;
        SUB( d, p, i );

        weights[2] = 1.0 - d;
        weights[0] = d;
        dist = DOT( d, d );
        return( dist );
    }
    else if( v_da0 <= 0.0 && v_da2 >= len_da2 )
    {
        weights[0] = 1.0;
        SUB( d, p, a0 );
        dist = DOT( d, d );
        return( dist );
    }
    else if( v_da1 <= 0.0 && v_da0 >= len_da0 )
    {
        weights[1] = 1.0;
        SUB( d, p, a1 );
        dist = DOT( d, d );
        return( dist );
    }
    else
    {
        weights[2] = 1.0;
        SUB( d, p, a2 );
        dist = DOT( d, d );
        return( dist );
    }

    return( 0.0 );
}

Real   Weights[6];


public  Real  sq_triangle_point_dist(
    Real   a0[],
    Real   a1[],
    Real   a2[],
    Real   point[] )
{
    return( sq_triangle_point_dist_and_weights( a0, a1, a2, point, Weights ) );
}

public  void  sq_triangle_point_dist_deriv(
    Real   a0[],
    Real   a1[],
    Real   a2[],
    Real   point[],
    Real   deriv_a0[],
    Real   deriv_a1[],
    Real   deriv_a2[],
    Real   deriv_point[] )
{
    int    dim, w;
    Real   weights[3];
    Real   *points[3], *derivs[4];

    (void) sq_triangle_point_dist_and_weights( a0, a1, a2, point, weights );

    points[0] = a0;
    points[1] = a1;
    points[2] = a2;

    derivs[0] = deriv_a0;
    derivs[1] = deriv_a1;
    derivs[2] = deriv_a2;
    derivs[3] = deriv_point;

    for_less( w, 0, 4 )
    for_less( dim, 0, N_DIMENSIONS )
        derivs[w][dim] = 0.0;

    compute_point_derivative( point, deriv_point,
                              weights, points, derivs );
}
