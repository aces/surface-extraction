#include  <fit_3d.h>

#define MIN_AND_MAX3( min, max, t1, t2, t3 ) \
    if( (t1) < (t2) )                                               \
    {                                                               \
        if( (t3) > (t2) )                                           \
        {                                                           \
            (min) = (t1);                                           \
            (max) = (t3);                                           \
        }                                                           \
        else if( (t3) > (t1) )                                      \
        {                                                           \
            (min) = (t1);                                           \
            (max) = (t2);                                           \
        }                                                           \
        else                                                        \
        {                                                           \
            (min) = (t3);                                           \
            (max) = (t2);                                           \
        }                                                           \
    }                                                               \
    else if( (t3) > (t1) )                                          \
    {                                                               \
        (min) = (t2);                                               \
        (max) = (t3);                                               \
    }                                                               \
    else if( (t3) > (t2) )                                          \
    {                                                               \
        (min) = (t2);                                               \
        (max) = (t1);                                               \
    }                                                               \
    else                                                            \
    {                                                               \
        (min) = (t3);                                               \
        (max) = (t1);                                               \
    }

public  Real sq_triangle_triangle_dist_estimate(
    Real    a0[],
    Real    a1[],
    Real    a2[],
    Real    b0[],
    Real    b1[],
    Real    b2[],
    Real    search_distance_sq )
{
    Real   a0x, a0y, a0z, a1x, a1y, a1z, a2x, a2y, a2z;
    Real   b0x, b0y, b0z, b1x, b1y, b1z, b2x, b2y, b2z;
    Real   a01x, a01y, a01z, a20x, a20y, a20z, a12x, a12y, a12z;
    Real   b01x, b01y, b01z, b20x, b20y, b20z, b12x, b12y, b12z;
    Real   nx, ny, nz, vx, vy, vz, dot0, dot1, dot2;
    Real   min_dir1, max_dir1, min_dir2, max_dir2, min_dir3;
    Real   dist, delta, mag_dir, min_b, max_b, min_a, max_a;
    Real   dist0, dist1, dist2, dist3, dista, distb;

    a0x = a0[0];
    a0y = a0[1];
    a0z = a0[2];
    a1x = a1[0];
    a1y = a1[1];
    a1z = a1[2];
    a2x = a2[0];
    a2y = a2[1];
    a2z = a2[2];

    b0x = b0[0];
    b0y = b0[1];
    b0z = b0[2];
    b1x = b1[0];
    b1y = b1[1];
    b1z = b1[2];
    b2x = b2[0];
    b2y = b2[1];
    b2z = b2[2];

#ifdef SPHERE_TEST
{
    Real   cx1, cy1, cz1, cx2, cy2, cz2, dist1, dist2, dist3;
    Real   dx, dy, dz, radius1, radius2;
    Real   sumx1, sumy1, sumz1, sumxx1, sumyy1, sumzz1;

    cx1 = (a0x + a1x + a2x) / 3.0;
    cy1 = (a0y + a1y + a2y) / 3.0;
    cz1 = (a0z + a1z + a2z) / 3.0;

/*
    dist1 = a0x * a0x + a0y * a0y + a0z * a0z -
            2.0 * (cx1 * a0x + cy1 * a0y + cz1 * a0z);
    dist2 = a1x * a1x + a1y * a1y + a1z * a1z -
            2.0 * (cx1 * a1x + cy1 * a1y + cz1 * a1z);
    dist3 = a2x * a2x + a2y * a2y + a2z * a2z -
            2.0 * (cx1 * a2x + cy1 * a2y + cz1 * a2z);
*/

    dx = a0x - cx1; dy = a0y - cy1; dz = a0z - cz1;
    dist1 = dx * dx + dy * dy + dz * dz;

    dx = a1x - cx1; dy = a1y - cy1; dz = a1z - cz1;
    dist2 = dx * dx + dy * dy + dz * dz;

    dx = a2x - cx1; dy = a2y - cy1; dz = a2z - cz1;
    dist3 = dx * dx + dy * dy + dz * dz;

    radius1 = MAX3( dist1, dist2, dist3 );
    radius1 = sqrt( radius1 + cx1 * cx1 + cy1 * cy1 + cz1 * cz1 );

    cx2 = (b0x + b1x + b2x) / 3.0;
    cy2 = (b0y + b1y + b2y) / 3.0;
    cz2 = (b0z + b1z + b2z) / 3.0;

    dx = b0x - cx2; dy = b0y - cy2; dz = b0z - cz2;
    dist1 = dx * dx + dy * dy + dz * dz;

    dx = b1x - cx2; dy = b1y - cy2; dz = b1z - cz2;
    dist2 = dx * dx + dy * dy + dz * dz;

    dx = b2x - cx2; dy = b2y - cy2; dz = b2z - cz2;
    dist3 = dx * dx + dy * dy + dz * dz;

    radius2 = MAX3( dist1, dist2, dist3 );
    radius2 = sqrt( radius2 );

    dx = cx1 - cx2;
    dy = cy1 - cy2;
    dz = cz1 - cz2;

    dist0 = dx * dx + dy * dy + dz * dz;
    dist0 = sqrt( dist0 ) - radius1 - radius2;
    if( dist0 > 0.0 && dist0 * dist0 >= search_distance_sq )
        return( dist0 * dist0 );
}
#endif

    /*--- test a to b */

    a01x = a1x - a0x;
    a01y = a1y - a0y;
    a01z = a1z - a0z;

    a12x = a2x - a1x;
    a12y = a2y - a1y;
    a12z = a2z - a1z;

    a20x = a0x - a2x;
    a20y = a0y - a2y;
    a20z = a0z - a2z;

    nx = a01y * (-a20z) - a01z * (-a20y);
    ny = a01z * (-a20x) - a01x * (-a20z);
    nz = a01x * (-a20y) - a01y * (-a20x);

    /*--- do dir3 distance */

    dot0 = b0x * nx + b0y * ny + b0z * nz;
    dot1 = b1x * nx + b1y * ny + b1z * nz;
    dot2 = b2x * nx + b2y * ny + b2z * nz;
    min_dir3 = a0x * nx + a0y * ny + a0z * nz;

    MIN_AND_MAX3( min_b, max_b, dot0, dot1, dot2 );

    if( min_b > min_dir3 )
    {
        delta = min_b - min_dir3;
        mag_dir = nx * nx + ny * ny + nz * nz;
        dist3 = delta * delta / mag_dir;
    }
    else if( max_b < min_dir3 )
    {
        delta = min_dir3 - max_b;
        mag_dir = nx * nx + ny * ny + nz * nz;
        dist3 = delta * delta / mag_dir;
    }
    else
        dist3 = 0.0;

    if( dist3 >= search_distance_sq )
        return( dist3 );

    /*--- get dist0 */

    dist0 = dist3;

    vx = ny * a01z - nz * a01y;
    vy = nz * a01x - nx * a01z;
    vz = nx * a01y - ny * a01x;

    dot0 = a0x * a01x + a0y * a01y + a0z * a01z;
    dot1 = a1x * a01x + a1y * a01y + a1z * a01z;
    dot2 = a2x * a01x + a2y * a01y + a2z * a01z;

    min_dir1 = MIN( dot0, dot2 );
    max_dir1 = MAX( dot1, dot2 );

    min_dir2 = a0x * vx + a0y * vy + a0z * vz;
    max_dir2 = a2x * vx + a2y * vy + a2z * vz;

    /*--- do dir1 distance */

    dot0 = b0x * a01x + b0y * a01y + b0z * a01z;
    dot1 = b1x * a01x + b1y * a01y + b1z * a01z;
    dot2 = b2x * a01x + b2y * a01y + b2z * a01z;

    MIN_AND_MAX3( min_b, max_b, dot0, dot1, dot2 );

    if( min_b > max_dir1 )
    {
        delta = min_b - max_dir1;
        mag_dir = a01x * a01x + a01y * a01y + a01z * a01z;
        dist0 += delta * delta / mag_dir;
    }
    else if( max_b < min_dir1 )
    {
        delta = min_dir1 - max_b;
        mag_dir = a01x * a01x + a01y * a01y + a01z * a01z;
        dist0 += delta * delta / mag_dir;
    }

    /*--- do dir2 distance */

    dot0 = b0x * vx + b0y * vy + b0z * vz;
    dot1 = b1x * vx + b1y * vy + b1z * vz;
    dot2 = b2x * vx + b2y * vy + b2z * vz;

    MIN_AND_MAX3( min_b, max_b, dot0, dot1, dot2 );

    if( min_b > max_dir2 )
    {
        delta = min_b - max_dir2;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist0 += delta * delta / mag_dir;
    }
    else if( max_b < min_dir2 )
    {
        delta = min_dir2 - max_b;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist0 += delta * delta / mag_dir;
    }

    if( dist0 >= search_distance_sq )
        return( dist0 );

    /*--- get dist1 */

    dist1 = dist3;

    vx = ny * a12z - nz * a12y;
    vy = nz * a12x - nx * a12z;
    vz = nx * a12y - ny * a12x;

    dot0 = a0x * a12x + a0y * a12y + a0z * a12z;
    dot1 = a1x * a12x + a1y * a12y + a1z * a12z;
    dot2 = a2x * a12x + a2y * a12y + a2z * a12z;

    min_dir1 = MIN( dot1, dot0 );
    max_dir1 = MAX( dot2, dot0 );

    min_dir2 = a1x * vx + a1y * vy + a1z * vz;
    max_dir2 = a0x * vx + a0y * vy + a0z * vz;

    /*--- do dir1 distance */

    dot0 = b0x * a12x + b0y * a12y + b0z * a12z;
    dot1 = b1x * a12x + b1y * a12y + b1z * a12z;
    dot2 = b2x * a12x + b2y * a12y + b2z * a12z;

    MIN_AND_MAX3( min_b, max_b, dot0, dot1, dot2 );

    if( min_b > max_dir1 )
    {
        delta = min_b - max_dir1;
        mag_dir = a12x * a12x + a12y * a12y + a12z * a12z;
        dist1 += delta * delta / mag_dir;
    }
    else if( max_b < min_dir1 )
    {
        delta = min_dir1 - max_b;
        mag_dir = a12x * a12x + a12y * a12y + a12z * a12z;
        dist1 += delta * delta / mag_dir;
    }

    /*--- do dir2 distance */

    dot0 = b0x * vx + b0y * vy + b0z * vz;
    dot1 = b1x * vx + b1y * vy + b1z * vz;
    dot2 = b2x * vx + b2y * vy + b2z * vz;

    MIN_AND_MAX3( min_b, max_b, dot0, dot1, dot2 );

    if( min_b > max_dir2 )
    {
        delta = min_b - max_dir2;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist1 += delta * delta / mag_dir;
    }
    else if( max_b < min_dir2 )
    {
        delta = min_dir2 - max_b;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist1 += delta * delta / mag_dir;
    }

    if( dist1 >= search_distance_sq )
        return( dist1 );

    /*--- get dist2 */

    dist2 = dist3;

    vx = ny * a20z - nz * a20y;
    vy = nz * a20x - nx * a20z;
    vz = nx * a20y - ny * a20x;

    dot0 = a0x * a20x + a0y * a20y + a0z * a20z;
    dot1 = a1x * a20x + a1y * a20y + a1z * a20z;
    dot2 = a2x * a20x + a2y * a20y + a2z * a20z;

    min_dir1 = MIN( dot2, dot1 );
    max_dir1 = MAX( dot0, dot1 );

    min_dir2 = a2x * vx + a2y * vy + a2z * vz;
    max_dir2 = a1x * vx + a1y * vy + a1z * vz;

    /*--- do dir1 distance */

    dot0 = b0x * a20x + b0y * a20y + b0z * a20z;
    dot1 = b1x * a20x + b1y * a20y + b1z * a20z;
    dot2 = b2x * a20x + b2y * a20y + b2z * a20z;

    MIN_AND_MAX3( min_b, max_b, dot0, dot1, dot2 );

    if( min_b > max_dir1 )
    {
        delta = min_b - max_dir1;
        mag_dir = a20x * a20x + a20y * a20y + a20z * a20z;
        dist2 += delta * delta / mag_dir;
    }
    else if( max_b < min_dir1 )
    {
        delta = min_dir1 - max_b;
        mag_dir = a20x * a20x + a20y * a20y + a20z * a20z;
        dist2 += delta * delta / mag_dir;
    }

    /*--- do dir2 distance */

    dot0 = b0x * vx + b0y * vy + b0z * vz;
    dot1 = b1x * vx + b1y * vy + b1z * vz;
    dot2 = b2x * vx + b2y * vy + b2z * vz;

    MIN_AND_MAX3( min_b, max_b, dot0, dot1, dot2 );

    if( min_b > max_dir2 )
    {
        delta = min_b - max_dir2;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist2 += delta * delta / mag_dir;
    }
    else if( max_b < min_dir2 )
    {
        delta = min_dir2 - max_b;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist2 += delta * delta / mag_dir;
    }

    dista = MAX3( dist0, dist1, dist2 );

    /*--- test a to b */

    b01x = b1x - b0x;
    b01y = b1y - b0y;
    b01z = b1z - b0z;

    b12x = b2x - b1x;
    b12y = b2y - b1y;
    b12z = b2z - b1z;

    b20x = b0x - b2x;
    b20y = b0y - b2y;
    b20z = b0z - b2z;

    nx = b01y * (-b20z) - b01z * (-b20y);
    ny = b01z * (-b20x) - b01x * (-b20z);
    nz = b01x * (-b20y) - b01y * (-b20x);

    /*--- do dir3 distance */

    dot0 = a0x * nx + a0y * ny + a0z * nz;
    dot1 = a1x * nx + a1y * ny + a1z * nz;
    dot2 = a2x * nx + a2y * ny + a2z * nz;
    min_dir3 = b0x * nx + b0y * ny + b0z * nz;

    MIN_AND_MAX3( min_a, max_a, dot0, dot1, dot2 );

    if( min_a > min_dir3 )
    {
        delta = min_a - min_dir3;
        mag_dir = nx * nx + ny * ny + nz * nz;
        dist3 = delta * delta / mag_dir;
    }
    else if( max_a < min_dir3 )
    {
        delta = min_dir3 - max_a;
        mag_dir = nx * nx + ny * ny + nz * nz;
        dist3 = delta * delta / mag_dir;
    }
    else
        dist3 = 0.0;

    if( dist3 >= search_distance_sq )
        return( dist3 );

    return( 0.0 );

/*---  profiling indicates that it is faster just to stop here, rather than
       to repeat the first half of this function with the roles of a and b
       inverted */

    /*--- get dist0 */

    dist0 = dist3;

    vx = ny * b01z - nz * b01y;
    vy = nz * b01x - nx * b01z;
    vz = nx * b01y - ny * b01x;

    dot0 = b0x * b01x + b0y * b01y + b0z * b01z;
    dot1 = b1x * b01x + b1y * b01y + b1z * b01z;
    dot2 = b2x * b01x + b2y * b01y + b2z * b01z;

    min_dir1 = MIN( dot0, dot2 );
    max_dir1 = MAX( dot1, dot2 );

    min_dir2 = b0x * vx + b0y * vy + b0z * vz;
    max_dir2 = b2x * vx + b2y * vy + b2z * vz;

    /*--- do dir1 distance */

    dot0 = a0x * b01x + a0y * b01y + a0z * b01z;
    dot1 = a1x * b01x + a1y * b01y + a1z * b01z;
    dot2 = a2x * b01x + a2y * b01y + a2z * b01z;

    MIN_AND_MAX3( min_a, max_a, dot0, dot1, dot2 );

    if( min_a > max_dir1 )
    {
        delta = min_a - max_dir1;
        mag_dir = b01x * b01x + b01y * b01y + b01z * b01z;
        dist0 += delta * delta / mag_dir;
    }
    else if( max_a < min_dir1 )
    {
        delta = min_dir1 - max_a;
        mag_dir = b01x * b01x + b01y * b01y + b01z * b01z;
        dist0 += delta * delta / mag_dir;
    }

    /*--- do dir2 distance */

    dot0 = a0x * vx + a0y * vy + a0z * vz;
    dot1 = a1x * vx + a1y * vy + a1z * vz;
    dot2 = a2x * vx + a2y * vy + a2z * vz;

    MIN_AND_MAX3( min_a, max_a, dot0, dot1, dot2 );

    if( min_a > max_dir2 )
    {
        delta = min_a - max_dir2;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist0 += delta * delta / mag_dir;
    }
    else if( max_a < min_dir2 )
    {
        delta = min_dir2 - max_a;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist0 += delta * delta / mag_dir;
    }

    if( dist0 >= search_distance_sq )
        return( dist0 );

    /*--- get dist1 */

    dist1 = dist3;

    vx = ny * b12z - nz * b12y;
    vy = nz * b12x - nx * b12z;
    vz = nx * b12y - ny * b12x;

    dot0 = b0x * b12x + b0y * b12y + b0z * b12z;
    dot1 = b1x * b12x + b1y * b12y + b1z * b12z;
    dot2 = b2x * b12x + b2y * b12y + b2z * b12z;

    min_dir1 = MIN( dot1, dot0 );
    max_dir1 = MAX( dot2, dot0 );

    min_dir2 = b1x * vx + b1y * vy + b1z * vz;
    max_dir2 = b0x * vx + b0y * vy + b0z * vz;

    /*--- do dir1 distance */

    dot0 = a0x * b12x + a0y * b12y + a0z * b12z;
    dot1 = a1x * b12x + a1y * b12y + a1z * b12z;
    dot2 = a2x * b12x + a2y * b12y + a2z * b12z;

    MIN_AND_MAX3( min_a, max_a, dot0, dot1, dot2 );

    if( min_a > max_dir1 )
    {
        delta = min_a - max_dir1;
        mag_dir = b12x * b12x + b12y * b12y + b12z * b12z;
        dist1 += delta * delta / mag_dir;
    }
    else if( max_a < min_dir1 )
    {
        delta = min_dir1 - max_a;
        mag_dir = b12x * b12x + b12y * b12y + b12z * b12z;
        dist1 += delta * delta / mag_dir;
    }

    /*--- do dir2 distance */

    dot0 = a0x * vx + a0y * vy + a0z * vz;
    dot1 = a1x * vx + a1y * vy + a1z * vz;
    dot2 = a2x * vx + a2y * vy + a2z * vz;

    MIN_AND_MAX3( min_a, max_a, dot0, dot1, dot2 );

    if( min_a > max_dir2 )
    {
        delta = min_a - max_dir2;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist1 += delta * delta / mag_dir;
    }
    else if( max_a < min_dir2 )
    {
        delta = min_dir2 - max_a;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist1 += delta * delta / mag_dir;
    }

    if( dist1 >= search_distance_sq )
        return( dist1 );

    /*--- get dist2 */

    dist2 = dist3;

    vx = ny * b20z - nz * b20y;
    vy = nz * b20x - nx * b20z;
    vz = nx * b20y - ny * b20x;

    dot0 = b0x * b20x + b0y * b20y + b0z * b20z;
    dot1 = b1x * b20x + b1y * b20y + b1z * b20z;
    dot2 = b2x * b20x + b2y * b20y + b2z * b20z;

    min_dir1 = MIN( dot2, dot1 );
    max_dir1 = MAX( dot0, dot1 );

    min_dir2 = b2x * vx + b2y * vy + b2z * vz;
    max_dir2 = b1x * vx + b1y * vy + b1z * vz;

    /*--- do dir1 distance */

    dot0 = a0x * b20x + a0y * b20y + a0z * b20z;
    dot1 = a1x * b20x + a1y * b20y + a1z * b20z;
    dot2 = a2x * b20x + a2y * b20y + a2z * b20z;

    MIN_AND_MAX3( min_a, max_a, dot0, dot1, dot2 );

    if( min_a > max_dir1 )
    {
        delta = min_a - max_dir1;
        mag_dir = b20x * b20x + b20y * b20y + b20z * b20z;
        dist2 += delta * delta / mag_dir;
    }
    else if( max_a < min_dir1 )
    {
        delta = min_dir1 - max_a;
        mag_dir = b20x * b20x + b20y * b20y + b20z * b20z;
        dist2 += delta * delta / mag_dir;
    }

    /*--- do dir2 distance */

    dot0 = a0x * vx + a0y * vy + a0z * vz;
    dot1 = a1x * vx + a1y * vy + a1z * vz;
    dot2 = a2x * vx + a2y * vy + a2z * vz;

    MIN_AND_MAX3( min_a, max_a, dot0, dot1, dot2 );

    if( min_a > max_dir2 )
    {
        delta = min_a - max_dir2;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist2 += delta * delta / mag_dir;
    }
    else if( max_a < min_dir2 )
    {
        delta = min_dir2 - max_a;
        mag_dir = vx * vx + vy * vy + vz * vz;
        dist2 += delta * delta / mag_dir;
    }

    distb = MAX3( dist0, dist1, dist2 );

    dist = MAX( dista, distb );

    return( dist );
}

#define  MAX_THRESHOLD_N_POLYGONS   200
#define  THRESHOLD_N_POLYGONS   60

#ifdef PRINT_DIST
#define PRINT_DIST
static Real sum_dist;
#endif

#ifdef DEBUG
#define DEBUG
static int n_tests1 = 0;
static int n_less1 = 0;
static int n_tests2 = 0;
static int n_less2 = 0;


private  void  increment_test_more1( void )
{
    ++n_tests1;
    if( n_tests1 % 10000 == 0 )
    {
        print( "tri-tri tests 1: %d / %d\n", n_less1, n_tests1 );
    }
}

private  void  increment_test_less1( void )
{
    ++n_less1;
    increment_test_more1();
}

private  void  increment_test_more2( void )
{
    ++n_tests2;
    if( n_tests2 % 10000 == 0 )
    {
        print( "tri-tri tests 2: %d / %d\n", n_less2, n_tests2 );
    }
}

private  void  increment_test_less2( void )
{
    ++n_less2;
    increment_test_more2();
}

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
    Real               parameters[],
    poly_info_struct   poly_info[],
    int                **list_of_polys,
    unsigned int       **poly_classes,
    int                *n_alloced,
    int                offset_index,
    int                n_polygons,
    int                next_offset_index,
    Real               line_dir[],
    int                *n_pairs_ptr,
    int                *p1s_ptr[],
    int                *p2s_ptr[],
    unsigned char      *cases_ptr[],
    float              *min_line_dists_ptr[],
    int                n_points,
    Smallest_int       used_flags[],
    int                n_neighbours[],
    int                *neighbours[],
    Real               *closest_distance )
{
    int               parms1[3], parms2[3];
    Real              mid_plane, min_distance_sq, search_distance_sq;
    Real              movement, movement_sq, dx, dy, dz, d, t_min, dist;
    poly_info_struct  *info_ptr;
    poly_info_struct  *static_compare[MAX_THRESHOLD_N_POLYGONS];
    poly_info_struct  **compare;
    int               which1, which2, n1, n2;
    float             mid_plane_plus, mid_plane_minus;
    Real              new_min[N_DIMENSIONS], new_max[N_DIMENSIONS];
    float             x1_min, x1_max, y1_min, y1_max, z1_min, z1_max;
    float             fl_search_distance;
    float             x1_low, x1_high, y1_low, y1_high, z1_low, z1_high;
    int               split_axis, dim, *list, poly1, p1, p2, poly;
    unsigned int      *classes, *src_classes, class_bit;
    unsigned int      *new_classes1, *new_classes2;
    int               i, i1, i2, diff;
    int               *src_list, *new_list1, *new_list2;
    int               new_n_polygons1, new_n_polygons2;
    int               n11, n12, n21, n22;
    int               parm_p1, parm_p2, parm_n11, parm_n12, parm_n21, parm_n22;
    BOOLEAN           on_left, on_right;
    unsigned int      cl;
    int               which_case;
    int               count[17], c1, c2, start;
    int               static_list[MAX_THRESHOLD_N_POLYGONS];
    int               *sorted_list;
    int               n_to_compare;
    static int        threshold_n_polygons = THRESHOLD_N_POLYGONS;
    static BOOLEAN    first = TRUE;

    if( first )
    {
        first = FALSE;
        if( getenv( "THRESHOLD_N_POLYGONS1" ) == NULL ||
            sscanf( getenv("THRESHOLD_N_POLYGONS1" ), "%d",
                    &threshold_n_polygons ) != 1  ||
            threshold_n_polygons > MAX_THRESHOLD_N_POLYGONS )
            threshold_n_polygons = THRESHOLD_N_POLYGONS;
    }

    if( n_polygons == 0 )
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
        n_polygons >= threshold_n_polygons )
    {
        if( next_offset_index + 2 * n_polygons > *n_alloced )
        {
            *n_alloced = next_offset_index + 2 * n_polygons;
            REALLOC( *list_of_polys, *n_alloced );
            REALLOC( *poly_classes, *n_alloced );
        }

        mid_plane = (max_range[split_axis] + min_range[split_axis]) / 2.0;
        mid_plane_plus = (float) (mid_plane + search_distance / 2.0);
        mid_plane_minus = (float) (mid_plane - search_distance / 2.0);

        /*--- do everything to left */

        src_list = &(*list_of_polys)[offset_index];
        src_classes = &(*poly_classes)[offset_index];

        new_list1 = &(*list_of_polys)[next_offset_index];
        new_classes1 = &(*poly_classes)[next_offset_index];
        new_n_polygons1 = 0;

        new_list2 = &(*list_of_polys)[next_offset_index+n_polygons];
        new_classes2 = &(*poly_classes)[next_offset_index+n_polygons];
        new_n_polygons2 = 0;

        class_bit = (1u << split_axis);

        for_less( i, 0, n_polygons )
        {
            poly = src_list[i];
            on_left = poly_info[poly].low_limits[split_axis] <= mid_plane_plus;
            on_right = poly_info[poly].high_limits[split_axis] >=
                                        mid_plane_minus;
            if( on_left )
            {
                new_list1[new_n_polygons1] = poly;
                new_classes1[new_n_polygons1] = src_classes[i];
                ++new_n_polygons1;
            }
            if( on_right )
            {
                new_list2[new_n_polygons2] = poly;
                if( on_left )
                    new_classes2[new_n_polygons2] = (src_classes[i] |class_bit);
                else
                    new_classes2[new_n_polygons2] = src_classes[i];
                ++new_n_polygons2;
            }
        }

        diff = n_polygons - new_n_polygons1;
        for_less( i, 0, new_n_polygons2 )
        {
            new_list2[i-diff] = new_list2[i];
            new_classes2[i-diff] = new_classes2[i];
        }

        if( new_n_polygons1 < n_polygons && new_n_polygons2 < n_polygons )
        {
            new_max[0] = max_range[0];
            new_max[1] = max_range[1];
            new_max[2] = max_range[2];
            new_max[split_axis] = (Real) mid_plane_plus;

            recursive_find_close_pairs( min_range, new_max,
                        min_distance, search_distance, parameters,
                        poly_info, list_of_polys, poly_classes,
                        n_alloced, next_offset_index, new_n_polygons1,
                        next_offset_index + new_n_polygons1 + new_n_polygons2,
                        line_dir, n_pairs_ptr, p1s_ptr, p2s_ptr, cases_ptr,
                        min_line_dists_ptr, n_points, used_flags,
                        n_neighbours, neighbours, closest_distance );

            new_min[0] = min_range[0];
            new_min[1] = min_range[1];
            new_min[2] = min_range[2];
            new_min[split_axis] = (Real) mid_plane_minus;

            recursive_find_close_pairs( new_min, max_range,
                        min_distance, search_distance, parameters,
                        poly_info, list_of_polys, poly_classes,
                        n_alloced, next_offset_index + new_n_polygons1,
                        new_n_polygons2,
                        next_offset_index + new_n_polygons1 + new_n_polygons2,
                        line_dir, n_pairs_ptr, p1s_ptr, p2s_ptr, cases_ptr,
                        min_line_dists_ptr, n_points, used_flags,
                        n_neighbours, neighbours, closest_distance );

            return;
        }
    }

    min_distance_sq = min_distance * min_distance;
    search_distance_sq = search_distance * search_distance;
    fl_search_distance = (float) search_distance;
    list = &(*list_of_polys)[offset_index];
    classes = &(*poly_classes)[offset_index];

    if( n_polygons > threshold_n_polygons )
    {
        ALLOC( sorted_list, n_polygons );
        ALLOC( compare, n_polygons );
    }
    else
    {
        sorted_list = static_list;
        compare = static_compare;
    }
        
    for_less( i1, 0, 16 )
        count[i1] = 0;

    for_less( i1, 0, n_polygons )
        ++count[classes[i1]];

    for_less( i1, 1, 16 )
        count[i1] += count[i1-1];

    for_less( i1, 0, n_polygons )
    {
        cl = classes[i1];
        --count[cl];
        sorted_list[count[cl]] = list[i1];
    }

    count[16] = n_polygons;

    for_less( c1, 0, 16 )
    {
        if( count[c1] == count[c1+1] )
            continue;

        n_to_compare = 0;

        for_less( c2, c1, 16 )
        {
            if( (c1 & c2) != 0 || count[c2] == count[c2+1] )
                continue;

            for_less( i2, count[c2], count[c2+1] )
            {
                compare[n_to_compare] = &poly_info[sorted_list[i2]];
                ++n_to_compare;
            }
        }

        for_less( i1, count[c1], count[c1+1] )
        {
            poly1 = sorted_list[i1];
            p1 = poly_info[poly1].p1;
            n11 = poly_info[poly1].n11;
            n12 = poly_info[poly1].n12;

            used_flags[p1] = TRUE;
            used_flags[n11] = TRUE;
            used_flags[n12] = TRUE;

            parm_p1 = IJ( p1, 0, 3 );
            parm_n11 = IJ( n11, 0, 3 );
            parm_n12 = IJ( n12, 0, 3 );

            x1_min = poly_info[poly1].low_limits[X];
            x1_max = poly_info[poly1].high_limits[X];
            y1_min = poly_info[poly1].low_limits[Y];
            y1_max = poly_info[poly1].high_limits[Y];
            z1_min = poly_info[poly1].low_limits[Z];
            z1_max = poly_info[poly1].high_limits[Z];

            x1_low = x1_min - fl_search_distance;
            x1_high = x1_max + fl_search_distance;
            y1_low = y1_min - fl_search_distance;
            y1_high = y1_max + fl_search_distance;
            z1_low = z1_min - fl_search_distance;
            z1_high = z1_max + fl_search_distance;

            if( c1 == 0 )
                start = i1 + 1; 
            else
                start = 0; 

            for_less( i2, start, n_to_compare )
            {
                info_ptr = compare[i2];

                if( info_ptr->low_limits[X] > x1_high ||
                    info_ptr->high_limits[X] < x1_low ||
                    info_ptr->low_limits[Y] > y1_high ||
                    info_ptr->high_limits[Y] < y1_low ||
                    info_ptr->low_limits[Z] > z1_high ||
                    info_ptr->high_limits[Z] < z1_low ||
                    used_flags[ p2 = info_ptr->p1 ] ||
                    used_flags[ n21 = info_ptr->n11 ] ||
                    used_flags[ n22 = info_ptr->n12 ] )
                {
                    continue;
                }

                parm_p2 = 3 * p2;
                parm_n21 = 3 * n21;
                parm_n22 = 3 * n22;

                if( sq_triangle_triangle_dist_estimate(
                                                &parameters[parm_p1],
                                                &parameters[parm_n11],
                                                &parameters[parm_n12],
                                                &parameters[parm_p2],
                                                &parameters[parm_n21],
                                                &parameters[parm_n22],
                                                search_distance_sq ) >=
                                                search_distance_sq )
                {
                    continue;
                }

                which_case = 0;
                dist = sq_triangle_triangle_dist( &parameters[parm_p1],
                                                  &parameters[parm_n11],
                                                  &parameters[parm_n12],
                                                  &parameters[parm_p2],
                                                  &parameters[parm_n21],
                                                  &parameters[parm_n22],
                                                  &which_case );

                if( dist >= search_distance_sq )
                {
#ifdef DEBUG
                    increment_test_more1();
#endif
                    continue;
                }
#ifdef DEBUG
                increment_test_less1();
#endif

#ifdef PRINT_DIST
sum_dist += dist;
#endif

                if( *closest_distance < 0.0 || dist < *closest_distance )
                    *closest_distance = dist;

                if( dist <= min_distance_sq || line_dir == NULL )
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
                        dx = line_dir[parms1[which1]+0] -
                             line_dir[parms2[which2]+0];
                        dy = line_dir[parms1[which1]+1] -
                             line_dir[parms2[which2]+1];
                        dz = line_dir[parms1[which1]+2] -
                             line_dir[parms2[which2]+2];
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

                for_less( n1, 0, n_neighbours[p1] )
                {
                    if( neighbours[p1][n1] == n11 &&
                        neighbours[p1][(n1+1)%n_neighbours[p1]] == n12 )
                        break;
                }
                if( n1 >= n_neighbours[p1] )
                    handle_internal_error( "n1 >= n_neighbours[p1]" );

                for_less( n2, 0, n_neighbours[p2] )
                {
                    if( neighbours[p2][n2] == n21 &&
                        neighbours[p2][(n2+1)%n_neighbours[p2]] == n22 )
                        break;
                }
                if( n2 >= n_neighbours[p2] )
                    handle_internal_error( "n2 >= n_neighbours[p2]" );

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
                (*p1s_ptr)[*n_pairs_ptr] = IJ( n1, p1, n_points );
                (*p2s_ptr)[*n_pairs_ptr] = IJ( n2, p2, n_points );
                (*cases_ptr)[*n_pairs_ptr] = (unsigned char) which_case;
                ++(*n_pairs_ptr);
            }

            used_flags[p1] = FALSE;
            used_flags[n11] = FALSE;
            used_flags[n12] = FALSE;
        }
    }

    if( n_polygons > threshold_n_polygons )
    {
        FREE( sorted_list );
        FREE( compare );
    }
}

private  void   find_self_intersect_candidates(
    Real               min_distance,
    int                *n_pairs_ptr,
    int                *p1s_ptr[],
    int                *p2s_ptr[],
    unsigned char      *cases_ptr[],
    float              *min_line_dists_ptr[],
    int                n_neighbours[],
    int                *neighbours[],
    int                n_points,
    Real               parameters[],
    Real               max_movement,
    Real               line_dir[],
    Real               *closest_distance )
{
    poly_info_struct  *poly_info;
    int               poly1, n_polygons;
    int               p1, dim, p, n;
    int               n11, n12;
    int               *list_of_polys, n_alloced;
    unsigned int      *poly_classes;
    Real              min_range[N_DIMENSIONS], max_range[N_DIMENSIONS];
    Real              search_distance;
    Smallest_int      *used_flags;

#ifdef PRINT_DIST
sum_dist = 0.0;
#endif

    search_distance = min_distance + max_movement;

    n_polygons = 0;
    for_less( p, 0, n_points )
    {
        for_less( n, 0, n_neighbours[p] )
        {
            n11 = neighbours[p][n];
            n12 = neighbours[p][(n+1) % n_neighbours[p]];
            if( n11 > p && n12 > p )
                ++n_polygons;
        }
    }

    ALLOC( poly_info, n_polygons );
    poly1 = 0;
    for_less( p, 0, n_points )
    {
        for_less( n, 0, n_neighbours[p] )
        {
            n11 = neighbours[p][n];
            n12 = neighbours[p][(n+1) % n_neighbours[p]];
            if( n11 > p && n12 > p )
            {
                poly_info[poly1].p1 = p;
                poly_info[poly1].n11 = n11;
                poly_info[poly1].n12 = n12;

                for_less( dim, 0, N_DIMENSIONS )
                {
                    poly_info[poly1].low_limits[dim] = (float)
                                            MIN3( parameters[IJ(p,dim,3)],
                                                  parameters[IJ(n11,dim,3)],
                                                  parameters[IJ(n12,dim,3)] );
                    poly_info[poly1].high_limits[dim] = (float)
                                            MAX3( parameters[IJ(p,dim,3)],
                                                  parameters[IJ(n11,dim,3)],
                                                  parameters[IJ(n12,dim,3)] );
                }

                ++poly1;
            }
        }
    }

    for_less( dim, 0, N_DIMENSIONS )
    {
        min_range[dim] = 0.0;
        max_range[dim] = 0.0;
    }

    for_less( p, 0, n_points )
    {
        for_less( dim, 0, N_DIMENSIONS )
        {
            if( p == 0 || parameters[IJ(p,dim,3)] < min_range[dim] )
                min_range[dim] = parameters[IJ(p,dim,3)];
            if( p == 0 || parameters[IJ(p,dim,3)] > max_range[dim] )
                max_range[dim] = parameters[IJ(p,dim,3)];
        }
    }

    n_alloced = 3 * n_polygons;
    ALLOC( list_of_polys, n_alloced );
    ALLOC( poly_classes, n_alloced );

    for_less( poly1, 0, n_polygons )
    {
        list_of_polys[poly1] = poly1;
        p1 = poly_info[poly1].p1;
        n11 = poly_info[poly1].n11;
        n12 = poly_info[poly1].n12;

        if( line_dir[IJ(p1,0,3)] == 0.0 &&
            line_dir[IJ(p1,1,3)] == 0.0 &&
            line_dir[IJ(p1,2,3)] == 0.0 &&
            line_dir[IJ(n11,0,3)] == 0.0 &&
            line_dir[IJ(n11,1,3)] == 0.0 &&
            line_dir[IJ(n11,2,3)] == 0.0 &&
            line_dir[IJ(n12,0,3)] == 0.0 &&
            line_dir[IJ(n12,1,3)] == 0.0 &&
            line_dir[IJ(n12,2,3)] == 0.0 )
        {
            poly_classes[poly1] = 8;
        }
        else
        {
            poly_classes[poly1] = 0;
        }
    }

    ALLOC( used_flags, n_points );
    for_less( p1, 0, n_points )
        used_flags[p1] = FALSE;

    *closest_distance = -1.0;

    *n_pairs_ptr = 0;
    *p1s_ptr = NULL;
    *p2s_ptr = NULL;
    *cases_ptr = NULL;
    *min_line_dists_ptr = NULL;
    recursive_find_close_pairs( min_range, max_range,
                        min_distance, search_distance, parameters,
                        poly_info, &list_of_polys, &poly_classes,
                        &n_alloced, 0, n_polygons, n_polygons,
                        line_dir, n_pairs_ptr, p1s_ptr, p2s_ptr, cases_ptr,
                        min_line_dists_ptr, n_points, used_flags,
                        n_neighbours, neighbours, closest_distance );
    FREE( used_flags );

    if( n_alloced > 0 )
    {
        FREE( list_of_polys );
        FREE( poly_classes );
    }

    FREE( poly_info );

    if( *closest_distance > 0.0 )
        *closest_distance = sqrt( *closest_distance );

#ifdef PRINT_DIST
print( "Sum dist: %.15g\n", sum_dist );
#endif
}

public  Real  recompute_self_intersects(
    Deform_struct                 *deform,
    int                           grid_size,
    int                           start_parameter[],
    Real                          parameters[],
    Real                          max_movement,
    Real                          line_dir[],
    self_intersect_lookup_struct  **si_lookup )
{
    int                    i, w, surface;
    Real                   *this_parms, *this_line, closest_distance;
    self_intersect_struct  *self;
    Real                   min_distance, close;

    closest_distance = -1.0;

    for_less( surface, 0, deform->n_surfaces )
    {
        this_parms = &parameters[start_parameter[surface]];
        if( line_dir == NULL )
            this_line = NULL;
        else
            this_line = &line_dir[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_self_intersects )
        {
            self = &deform->surfaces[surface].self_intersects[i];

            min_distance = 0.0;
            for_less( w, 0, self->n_weights )
                min_distance = MAX( min_distance, self->min_distances[w] );

            find_self_intersect_candidates(
                       min_distance,
                       &si_lookup[surface][i].n_pairs,
                       &si_lookup[surface][i].p1s,
                       &si_lookup[surface][i].p2s,
                       &si_lookup[surface][i].cases,
                       &si_lookup[surface][i].min_line_dists,
                       deform->surfaces[surface].surface.n_neighbours,
                       deform->surfaces[surface].surface.neighbours,
                       deform->surfaces[surface].surface.n_points,
                       this_parms, max_movement, this_line, &close );

            if( close >= 0.0 &&
                (closest_distance < 0.0 || close < closest_distance) )
                closest_distance = close;
        }
    }

    return( closest_distance );
}

public  void  initialize_self_intersect_lookup(
    self_intersect_lookup_struct  *si_lookup )
{
    si_lookup->n_pairs = 0;
    si_lookup->p1s = NULL;
    si_lookup->p2s = NULL;
    si_lookup->cases = NULL;
    si_lookup->min_line_dists = NULL;
}

public  void  delete_self_intersect_lookup(
    self_intersect_lookup_struct  *si_lookup )
{
    if( si_lookup->n_pairs > 0 )
    {
        FREE( si_lookup->p1s );
        FREE( si_lookup->p2s );
        FREE( si_lookup->cases );
        FREE( si_lookup->min_line_dists );
        si_lookup->n_pairs = 0;
    }
}

public  int  get_n_self_intersect_candidate(
    self_intersect_lookup_struct  *si_lookup )
{
    return( si_lookup->n_pairs );
}

public  BOOLEAN   test_self_intersect_candidate(
    self_intersect_lookup_struct  *si_lookup,
    BOOLEAN                       use_tri_tri_dist,
    int                           which,
    Real                          dist_from_computed_self_intersect,
    Real                          max_distance_sq,
    int                           n_points,
    Real                          parameters[],
    Smallest_int                  active_flags[],
    int                           n_neighbours[],
    int                           *neighbours[],
    Real                          *dist_sq )
{
    int   p1, n1, p2, n2, n11, n12, n21, n22, prev_case, which_case;

    if( dist_from_computed_self_intersect> 0.0 &&
        (Real) si_lookup->min_line_dists[which] >
        dist_from_computed_self_intersect)
    {
        return( FALSE );
    }

    p1 = si_lookup->p1s[which] % n_points;
    n1 = si_lookup->p1s[which] / n_points;
    p2 = si_lookup->p2s[which] % n_points;
    n2 = si_lookup->p2s[which] / n_points;

    prev_case = si_lookup->cases[which];

    n11 = neighbours[p1][n1];
    n12 = neighbours[p1][(n1+1)%n_neighbours[p1]];
    n21 = neighbours[p2][n2];
    n22 = neighbours[p2][(n2+1)%n_neighbours[p2]];

    if( active_flags != NULL &&
        !active_flags[p1] && !active_flags[n11] && !active_flags[n12] &&
        !active_flags[p2] && !active_flags[n21] && !active_flags[n22] )
    {
        return( FALSE );
    }

    which_case = prev_case;
    *dist_sq = sq_triangle_triangle_dist( &parameters[p1*3],
                                          &parameters[n11*3],
                                          &parameters[n12*3],
                                          &parameters[p2*3],
                                          &parameters[n21*3],
                                          &parameters[n22*3], &which_case );

    if( which_case != prev_case )
        si_lookup->cases[which] = (unsigned char) which_case;

#ifdef DEBUG
    if( *dist_sq > max_distance_sq )
        increment_test_more2();
    else
        increment_test_less2();
#endif

    return( TRUE );
}

private  void   create_self_intersect_deriv_info_tritri(
    self_intersect_lookup_struct  *si_lookup,
    int                           which,
    Real                          dist_sq,
    int                           n_points,
    Real                          parameters[],
    int                           n_neighbours[],
    int                           *neighbours[],
    self_intersect_deriv_struct   *deriv )
{
    int   p1, n1, p2, n2, n11, n12, n21, n22;

    p1 = si_lookup->p1s[which] % n_points;
    n1 = si_lookup->p1s[which] / n_points;
    p2 = si_lookup->p2s[which] % n_points;
    n2 = si_lookup->p2s[which] / n_points;

    n11 = neighbours[p1][n1];
    n12 = neighbours[p1][(n1+1)%n_neighbours[p1]];
    n21 = neighbours[p2][n2];
    n22 = neighbours[p2][(n2+1)%n_neighbours[p2]];

    if( dist_sq > 0.0 )
        deriv->dist = sqrt( dist_sq );
    else
        deriv->dist = sqrt( dist_sq );

    sq_triangle_triangle_dist_deriv( &parameters[IJ(p1,0,3)],
                                     &parameters[IJ(n11,0,3)],
                                     &parameters[IJ(n12,0,3)],
                                     &parameters[IJ(p2,0,3)],
                                     &parameters[IJ(n21,0,3)],
                                     &parameters[IJ(n22,0,3)],
                                     deriv->tri_tri_deriv[0],
                                     deriv->tri_tri_deriv[1],
                                     deriv->tri_tri_deriv[2],
                                     deriv->tri_tri_deriv[3],
                                     deriv->tri_tri_deriv[4],
                                     deriv->tri_tri_deriv[5] );
}

private  void   get_self_intersect_deriv_tritri(
    BOOLEAN                       use_square_flag,
    self_intersect_lookup_struct  *si_lookup,
    int                           which,
    Real                          min_distance,
    Real                          weight,
    int                           n_points,
    Real                          deriv[],
    int                           n_neighbours[],
    int                           *neighbours[],
    self_intersect_deriv_struct   *deriv_info )
{
    int   p1, n1, p2, n2, n11, n12, n21, n22;
    Real  diff, factor;

    p1 = si_lookup->p1s[which] % n_points;
    n1 = si_lookup->p1s[which] / n_points;
    p2 = si_lookup->p2s[which] % n_points;
    n2 = si_lookup->p2s[which] / n_points;

    n11 = neighbours[p1][n1];
    n12 = neighbours[p1][(n1+1)%n_neighbours[p1]];
    n21 = neighbours[p2][n2];
    n22 = neighbours[p2][(n2+1)%n_neighbours[p2]];

    diff = min_distance - deriv_info->dist;

    if( deriv_info->dist == 0.0 )
        factor = 0.0;
#ifdef USE_CUBE
#define USE_CUBE
    else
    factor = -3.0 * diff * diff * weight / 2.0 / deriv_info->dist / min_distance;
#else
    else if( use_square_flag )
    {
        factor = -2.0 * diff * weight / 2.0 / deriv_info->dist;
    }
    else
    {
        if( diff < 0.0 )
            factor = 1.0 / 2.0 / deriv_info->dist;
        else
            factor = -1.0 / 2.0 / deriv_info->dist;
    }
#endif

    deriv[IJ(p1,0,3)] += factor * deriv_info->tri_tri_deriv[0][0];
    deriv[IJ(p1,1,3)] += factor * deriv_info->tri_tri_deriv[0][1];
    deriv[IJ(p1,2,3)] += factor * deriv_info->tri_tri_deriv[0][2];

    deriv[IJ(n11,0,3)] += factor * deriv_info->tri_tri_deriv[1][0];
    deriv[IJ(n11,1,3)] += factor * deriv_info->tri_tri_deriv[1][1];
    deriv[IJ(n11,2,3)] += factor * deriv_info->tri_tri_deriv[1][2];

    deriv[IJ(n12,0,3)] += factor * deriv_info->tri_tri_deriv[2][0];
    deriv[IJ(n12,1,3)] += factor * deriv_info->tri_tri_deriv[2][1];
    deriv[IJ(n12,2,3)] += factor * deriv_info->tri_tri_deriv[2][2];

    deriv[IJ(p2,0,3)] += factor * deriv_info->tri_tri_deriv[3][0];
    deriv[IJ(p2,1,3)] += factor * deriv_info->tri_tri_deriv[3][1];
    deriv[IJ(p2,2,3)] += factor * deriv_info->tri_tri_deriv[3][2];

    deriv[IJ(n21,0,3)] += factor * deriv_info->tri_tri_deriv[4][0];
    deriv[IJ(n21,1,3)] += factor * deriv_info->tri_tri_deriv[4][1];
    deriv[IJ(n21,2,3)] += factor * deriv_info->tri_tri_deriv[4][2];

    deriv[IJ(n22,0,3)] += factor * deriv_info->tri_tri_deriv[5][0];
    deriv[IJ(n22,1,3)] += factor * deriv_info->tri_tri_deriv[5][1];
    deriv[IJ(n22,2,3)] += factor * deriv_info->tri_tri_deriv[5][2];
}

/*--------------- */

public  void   create_self_intersect_deriv_info(
    self_intersect_lookup_struct  *si_lookup,
    BOOLEAN                       use_tri_tri_dist,
    int                           which,
    Real                          dist_sq,
    int                           n_points,
    Real                          parameters[],
    int                           n_neighbours[],
    int                           *neighbours[],
    self_intersect_deriv_struct   *deriv )
{
    create_self_intersect_deriv_info_tritri(
                       si_lookup, which, dist_sq,
                       n_points, parameters, 
                       n_neighbours, neighbours, deriv );
}

public  void   get_self_intersect_deriv(
    self_intersect_lookup_struct  *si_lookup,
    BOOLEAN                       use_tri_tri_dist,
    BOOLEAN                       use_square_flag,
    int                           which,
    Real                          dist_sq,
    Real                          min_distance,
    Real                          weight,
    int                           n_points,
    Real                          parameters[],
    Real                          deriv[],
    int                           n_neighbours[],
    int                           *neighbours[],
    self_intersect_deriv_struct   *deriv_info )
{
    get_self_intersect_deriv_tritri( use_square_flag,
                       si_lookup, which, min_distance, weight,
                       n_points, deriv, n_neighbours, neighbours, deriv_info );
}
