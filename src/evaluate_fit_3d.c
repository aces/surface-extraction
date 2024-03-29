/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <time.h>
#include <fit_3d.h>

#define ADAPTIVE_RATIO 1
#define NO_ANCHOR 0
#define BOUNDARY_DECREASE 1

#define EVAL_NEAREST -1
#define EVAL_LINEAR 0
#define EVAL_CUBIC 2

private  Real   evaluate_intersect_wm_fit(
    object_struct **objects,
    Real          weight,
    Real          parameter[],
    int           *intersected,
    int           n_intersected,
    int           start_point,
    int           end_point
)
{
  Real         fit = 0.0f;
  int          n_intersections;
  Point        origin;
  Vector       direction;
  Real         dist, *distances;
  int          obj_index, point, p_index;
  int          i;

  if( weight > 0 )
  {
    fill_Vector(direction, 0, 1, 1);

    for( point=0; point<n_intersected && n_intersected<1000; point++ )

    {
      p_index = IJ( intersected[point], 0, 3 );

      fill_Point( origin, parameter[p_index+0], 
                  parameter[p_index+1], parameter[p_index+2] );

      n_intersections = intersect_ray_with_object( &origin, &direction,
                                               objects[0], &obj_index, &dist,
                                               &distances );

      // inside of the wm surface if number of intersections is odd
      if( n_intersections % 2 != 0 )
      {
        fit += 1 * weight;
      }

      FREE( distances );
    }
  }

  return fit;
}

private  Real   evaluate_intersect_wm_fit_deriv(
    object_struct **objects,
    Real          weight,
    Real          parameter[],
    Volume        wm_masked,
    int           *intersected,
    int           *n_intersected,
    int           n_neighbours[],
    int           *neighbours[],
    int           start_point,
    int           end_point,
    Real          deriv[]
)
{
  int             point, p_index;
  int             neigh, n, obj_index, n_intersections, max_neighbours;
  Point           *neigh_points, origin;
  Vector          normal, direction;
  Real            dist, *distances;
  Real            mask_value;

  if( weight > 0 )
  {
    *n_intersected = 0;
    max_neighbours = 0;
    for_less( n, start_point, end_point ){
      max_neighbours = MAX( max_neighbours, n_neighbours[n] );
    }
    ALLOC( neigh_points, max_neighbours);

    fill_Vector(direction, 1, 0, 0);

    for_less( point, start_point, end_point )
    {
      p_index = IJ( point, 0, 3 );

      evaluate_volume_in_world( wm_masked,
                                parameter[p_index+0],
                                parameter[p_index+1],
                                parameter[p_index+2],
                                EVAL_LINEAR, FALSE,  // isoval=0.5 on 0:1 mask, so need linear
                                0.0, &mask_value,
                                NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL );

      // When the vertex is suspected to be in the wm surface
      if( mask_value > 0 )
      {
        fill_Point( origin, parameter[p_index+0], 
                    parameter[p_index+1], parameter[p_index+2] );

        n_intersections = intersect_ray_with_object( &origin, &direction,
                                               objects[0], &obj_index, &dist,
                                               &distances );

        // if n_intersections is odd, that point is inside of the polygon
        if( n_intersections%2 != 0 )
        {
          for_less( n, 0, n_neighbours[point] ){
            neigh = neighbours[point][n];
            fill_Point( neigh_points[n],
                      parameter[IJ(neigh,0,3)],
                      parameter[IJ(neigh,1,3)],
                      parameter[IJ(neigh,2,3)] );
          }
          find_polygon_normal( n_neighbours[point], neigh_points, &normal);

          deriv[p_index+0] += -normal.coords[0] * weight;
          deriv[p_index+1] += -normal.coords[1] * weight;
          deriv[p_index+2] += -normal.coords[2] * weight;

          intersected[(*n_intersected)++] = point;
          //printf("%d: (%g,%g,%g) L_value=%g\n",point,parameter[p_index+0],parameter[p_index+1], parameter[p_index+2],laplace_value);
        }

        if( n_intersections > 0 )
        {
          FREE( distances );
        }
      }
    }

    FREE( neigh_points );
  }
}

// Now allow oversampling.
// This function is valid only for "to > from" (forward direction).

private  Real   evaluate_laplacian_fit(
    Real         weight,
    Volume       laplacian_map,
    Real         to,
    int          type,
    int          n_neighbours[],
    int        * neighbours[],
    Real         sampling,
    Real         parameter[],
    int          start_point,
    int          end_point)
{
  int         point, p_index, n;
  Real        fit=0.0f;
  Real        value;
  Real        dxyz[3];
  Point      *neigh_points;
  Vector      normal;

  if( weight > 0 ) {
    int max_neighbours = 0;
    for_less( n, start_point, end_point ){
      max_neighbours = MAX( max_neighbours, n_neighbours[n] );
    }
    ALLOC( neigh_points, max_neighbours);

    for_less( point, start_point, end_point ) {

      p_index = IJ( point, 0, 3 );

      // NOTE: It's probably ok to evaluate dxyz with EVAL_LINEAR
      //       since we are using the derivative only to establish
      //       the in/out direction with the surface normal. CL.
      evaluate_volume_in_world( laplacian_map,
                                parameter[p_index+0],
                                parameter[p_index+1],
                                parameter[p_index+2],
                                EVAL_LINEAR, FALSE,
                                0.0, &value,
                                dxyz, dxyz+1, dxyz+2,
                                NULL, NULL, NULL, NULL, NULL, NULL );

      // Make sure the normal direction is consistent with
      // the image value and gradient direction.

      if( type == 0 && value < to ) {
        for_less( n, 0, n_neighbours[point] ){
          int neigh = neighbours[point][n];
          fill_Point( neigh_points[n],
                      parameter[IJ(neigh,0,3)],
                      parameter[IJ(neigh,1,3)],
                      parameter[IJ(neigh,2,3)] );
        }
        find_polygon_normal( n_neighbours[point], neigh_points, &normal);
        Real dot = dxyz[0] * normal.coords[0] +
                   dxyz[1] * normal.coords[1] +
                   dxyz[2] * normal.coords[2];

        // this will negate the "reverse" contribution
        if( dot < 0.0 ) value = to;   
      }
      fit += (value-to)*(value-to);

    }

    // Do oversampling. This is useful, but slows down the code.

    if( sampling > 0 ) {

      int  n, neigh1, neigh2, p0, p1, p2, w1, w2;
      Real alpha1, alpha2, N0, N1, N2, xp, yp, zp;

      weight /= ((sampling+1.0)*(sampling+1.0));

      for( point = start_point; point < end_point; point++ ) {

        // WARNING: This is wrong in deep narrow sulci with "thin"
        //          surface since the normal vector changes rapidly.
        //          We must use the exact normal vector at each 
        //          sampling point, otherwise the direction of
        //          attraction at the sampling point will be wrong.
        //          CL.

        // Find the normal vector at this point. It's approximate only.

        for_less( n, 0, n_neighbours[point] ){
          int neigh = neighbours[point][n];
          fill_Point( neigh_points[n],
                      parameter[IJ(neigh,0,3)],
                      parameter[IJ(neigh,1,3)],
                      parameter[IJ(neigh,2,3)] );
        }
        find_polygon_normal( n_neighbours[point], neigh_points, &normal);

        for_less( n, 0, n_neighbours[point] ) {
          neigh1 = neighbours[point][n];
          neigh2 = neighbours[point][(n+1)%n_neighbours[point]];

          if( point > neigh1 || point > neigh2 ) continue;

          p0 = 3*point;
          p1 = 3*neigh1;
          p2 = 3*neigh2;

          for( w1 = 0; w1 <= sampling; w1++ ) {

            alpha1 = (Real) w1 / (Real) (sampling+1);

            for( w2 = 1; w2 <= sampling - w1 + 1; w2++ ) {
              if( ( w1 == 0 && w2 == sampling + 1 ) ||
                  ( w2 == sampling - w1 + 1 && neigh1 < neigh2 ) ) {
                continue;
              }

              alpha2 = (Real) w2 / (Real) (sampling - w1 +1);
              N0 = (1.0 - alpha1) * (1.0 - alpha2);
              N1 = alpha1;
              N2 = (1.0 - alpha1) * alpha2;
              xp = N0 * parameter[p0] + N1 * parameter[p1] + N2 * parameter[p2];
              yp = N0 * parameter[p0+1] + N1 * parameter[p1+1] + N2 * parameter[p2+1];
              zp = N0 * parameter[p0+2] + N1 * parameter[p1+2] + N2 * parameter[p2+2];

              // NOTE: It's probably ok to evaluate dxyz with EVAL_LINEAR
              //       since we are using the derivative only to establish
              //       the in/out direction with the surface normal. CL.

              evaluate_volume_in_world( laplacian_map, xp, yp, zp,
                                        EVAL_LINEAR, FALSE, 0.0, &value,
                                        dxyz, dxyz+1, dxyz+2,
                                        NULL, NULL, NULL, NULL, NULL, NULL );

              // this will negate the "reverse" contribution if we
              // are moving in the wrong direction
              if( type == 0 && value < to ) {
                Real dot = dxyz[0] * normal.coords[0] +
                           dxyz[1] * normal.coords[1] +
                           dxyz[2] * normal.coords[2];
                if( dot < 0.0 ) value = to;   
              }

              fit += (value-to)*(value-to);
            }
          }
        }
      }
    }

    fit *= 0.5 * weight;
    FREE( neigh_points );
  }

  return fit;
}

//
// Note: In the old code, the evaluation of the derivative would check
//       the tissue classification value against a threshold=1.5. However
//       the cls_correct file does not contain pve_csf values, thus it
//       reports value=2.0 (gray) in a sulcus, with phi >= 10 (uses
//       pve_csf to build Laplacian field). So to get the good gray
//       surface boundary, use the information from phi only to decide
//       to move inwards or outwards and ignore the cls value. There is
//       another confusing case where the white surface cuts through the
//       ventricules, with cls values < 1.5 but with possibly phi<10.
//       No matter what we do in this case, go inwards or outwards,
//       we'll never get it really right. So go outwards by default.
//       The old code would not assign a derivative to this node. CL.
//

private  Real   evaluate_laplacian_fit_deriv(
    Real         weight,
    Volume       laplacian_map,
    Volume       volume,
    Real         from,
    Real         to,
    int          type,
    Real         parameter[],
    int          start_point,
    int          end_point,
    int          n_neighbours[],
    int        * neighbours[],
    Real         sampling,
    Real         deriv[]) {

    Real         value1, volume_value;
    Real         dxyz[3];
    int          point, p_index, n;
    Real         factor = 1e0;
    Point        *neigh_points;
    Vector       normal, normal_p1, normal_p2, normal_p3;

    int          count_res = 0;
    Real         phi_res = 0.0;
    Real         phi_min = MAX( to, from );
    Real         phi_max = MIN( to, from );

    if( weight > 0 ) {

      int max_neighbours = 0;
      for_less( n, start_point, end_point ){
        max_neighbours = MAX( max_neighbours, n_neighbours[n] );
      }
      ALLOC( neigh_points, max_neighbours);

      if( sampling > 0 ) {
        weight /= (Real)((sampling+1)*(sampling+1));
      }
      count_res = end_point - start_point;
      for_less( point, start_point, end_point ) {
        p_index = IJ( point, 0, 3 );

        /* IMPORTANT: use CUBIC for derivatives since a symmetric
           stencil is used, whereas an upwind stencil is used for
           LINEAR, thus causing a bias in the derivatives. CL. */

        evaluate_volume_in_world( laplacian_map,
                                  parameter[p_index+0],
                                  parameter[p_index+1],
                                  parameter[p_index+2],
                                  EVAL_CUBIC, FALSE,
                                  0.0, &value1,
                                  dxyz, dxyz+1, dxyz+2,
                                  NULL, NULL, NULL, NULL, NULL, NULL );

        /* IMPORTANT: use LINEAR for the function values since CUBIC
           can cause over/undershoots and lead to missing the correct
           boundary of "to". CL. */

        evaluate_volume_in_world( laplacian_map,
                                  parameter[p_index+0],
                                  parameter[p_index+1],
                                  parameter[p_index+2],
                                  EVAL_LINEAR, FALSE,
                                  0.0, &value1,
                                  NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL );

        phi_res += ABS( to - value1 );
        if( value1 < phi_min ) phi_min = value1;
        if( value1 > phi_max ) phi_max = value1;

        factor = weight * (value1 - to);

        // Make sure the normal direction is consistent with
        // the image value and gradient direction.

        if( type == 0 && value1 < to ) {
          for_less( n, 0, n_neighbours[point] ){
            int neigh = neighbours[point][n];
            fill_Point( neigh_points[n],
                        parameter[IJ(neigh,0,3)],
                        parameter[IJ(neigh,1,3)],
                        parameter[IJ(neigh,2,3)] );
          }
          find_polygon_normal( n_neighbours[point], neigh_points, &normal);
          Real dot = dxyz[0] * normal.coords[0] +
                     dxyz[1] * normal.coords[1] +
                     dxyz[2] * normal.coords[2];
          // if( dot < 0.0 ) factor = -factor;
          if( dot < 0.0 ) factor = 0.0;
        }
        deriv[p_index+0] += dxyz[0]*factor;
        deriv[p_index+1] += dxyz[1]*factor;
        deriv[p_index+2] += dxyz[2]*factor;
      }

      // Do oversampling. This is useful, but slows down the code.

      if( sampling > 0 ) {

        int  n, nn, neigh1, neigh2, p0, p1, p2, w1, w2;
        Real dx0, dy0, dz0, dx1, dy1, dz1, nx, ny, nz, mag;
        Real alpha1, alpha2, N0, N1, N2, xp, yp, zp;

        for( point = start_point; point < end_point; point++ ) {

          for_less( n, 0, n_neighbours[point] ) {
            neigh1 = neighbours[point][n];
            neigh2 = neighbours[point][(n+1)%n_neighbours[point]];

            if( point > neigh1 || point > neigh2 ) continue;

            p0 = 3*point;
            p1 = 3*neigh1;
            p2 = 3*neigh2;

            for( w1 = 0; w1 <= sampling; w1++ ) {

              alpha1 = (Real) w1 / (Real) (sampling+1);

              for( w2 = 1; w2 <= sampling - w1 + 1; w2++ ) {
                if( ( w1 == 0 && w2 == sampling + 1 ) ||
                    ( w2 == sampling - w1 + 1 && neigh1 < neigh2 ) ) {
                  continue;
                }

                alpha2 = (Real) w2 / (Real) (sampling - w1 +1);
                N0 = (1.0 - alpha1) * (1.0 - alpha2);
                N1 = alpha1;
                N2 = (1.0 - alpha1) * alpha2;
                xp = N0 * parameter[p0] + N1 * parameter[p1] + N2 * parameter[p2];
                yp = N0 * parameter[p0+1] + N1 * parameter[p1+1] + N2 * parameter[p2+1];
                zp = N0 * parameter[p0+2] + N1 * parameter[p1+2] + N2 * parameter[p2+2];

                /* IMPORTANT: use CUBIC for derivatives since a symmetric
                   stencil is used, whereas an upwind stencil is used for
                   LINEAR, thus causing a bias in the derivatives. CL. */

                evaluate_volume_in_world( laplacian_map, xp, yp, zp,
                                          EVAL_CUBIC, FALSE, 0.0, &value1,
                                          dxyz, dxyz+1, dxyz+2,
                                          NULL, NULL, NULL, NULL, NULL, NULL );

                /* IMPORTANT: use LINEAR for the function values since CUBIC
                   can cause over/undershoots and lead to missing the correct
                   boundary of "to". CL. */

                evaluate_volume_in_world( laplacian_map, xp, yp, zp,
                                          EVAL_LINEAR, FALSE, 0.0, &value1,
                                          NULL, NULL, NULL,
                                          NULL, NULL, NULL, NULL, NULL, NULL );

                count_res++;
                phi_res += ABS( to - value1 );
                if( value1 < phi_min ) phi_min = value1;
                if( value1 > phi_max ) phi_max = value1;

                factor = weight * (value1 - to);

                // Make sure the normal direction is consistent with
                // the image value and gradient direction.

                if( type == 0 && value1 < to ) {

                  for_less( nn, 0, n_neighbours[point] ){
                    int neigh = neighbours[point][nn];
                    fill_Point( neigh_points[nn],
                                parameter[IJ(neigh,0,3)],
                                parameter[IJ(neigh,1,3)],
                                parameter[IJ(neigh,2,3)] );
                  }
                  find_polygon_normal( n_neighbours[point], neigh_points, 
                                       &normal_p1);

                  for_less( nn, 0, n_neighbours[neigh1] ){
                    int neigh = neighbours[neigh1][nn];
                    fill_Point( neigh_points[nn],
                                parameter[IJ(neigh,0,3)],
                                parameter[IJ(neigh,1,3)],
                                parameter[IJ(neigh,2,3)] );
                  }
                  find_polygon_normal( n_neighbours[neigh1], neigh_points, 
                                       &normal_p2);

                  for_less( nn, 0, n_neighbours[neigh2] ){
                    int neigh = neighbours[neigh2][nn];
                    fill_Point( neigh_points[nn],
                                parameter[IJ(neigh,0,3)],
                                parameter[IJ(neigh,1,3)],
                                parameter[IJ(neigh,2,3)] );
                  }
                  find_polygon_normal( n_neighbours[neigh2], neigh_points, 
                                       &normal_p3);

                  for( nn = 0; nn < 3; nn++ ) {
                    normal.coords[nn] = N0 * normal_p1.coords[nn] + 
                                        N1 * normal_p2.coords[nn] + 
                                        N2 * normal_p3.coords[nn];
                  }
                  Real dot = dxyz[0] * normal.coords[0] +
                             dxyz[1] * normal.coords[1] +
                             dxyz[2] * normal.coords[2];
                  // if( dot < 0.0 ) factor = -factor;
                  if( dot < 0.0 ) factor = 0.0;
                }

                deriv[p0+0] += dxyz[0]*factor*N0;
                deriv[p0+1] += dxyz[1]*factor*N0;
                deriv[p0+2] += dxyz[2]*factor*N0;
                deriv[p1+0] += dxyz[0]*factor*N1;
                deriv[p1+1] += dxyz[1]*factor*N1;
                deriv[p1+2] += dxyz[2]*factor*N1;
                deriv[p2+0] += dxyz[0]*factor*N2;
                deriv[p2+1] += dxyz[1]*factor*N2;
                deriv[p2+2] += dxyz[2]*factor*N2;
              }
            }
          }
        }
      }

      FREE( neigh_points );
      printf( "phi_res(%d) = %g  min = %g  max = %g\n", (end_point-start_point),
              phi_res/(float)(count_res), phi_min, phi_max );
    }
}

private  Real   evaluate_volume_fit(
    Real         weight,
    Real         max_weight,
    int          direction,
    Volume       volume,
    Volume       laplacian,
    Real         threshold,
    Real         parameter[],
    int          start_point,
    int          end_point,
    int          n_points)
{
    int          point;
    int          p_index;
    Real         value;
    Real         fit1=0.0, fit2=0.0;
    //clock_t      start, end, s_start, s_end;

    //start = clock();
    if( weight > 0 || max_weight > 0 )
    {
      for_less( point, start_point, end_point )
      {
        p_index = IJ( point, 0, 3 );
        // for the WM surface
        if( direction < 0 )
        {
          evaluate_volume_in_world( volume,
                                parameter[p_index+0],
                                parameter[p_index+1],
                                parameter[p_index+2],
                                EVAL_LINEAR, FALSE,
                                0.0, &value,
                                NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL );
          fit1 += value/6 * max_weight;
        }
        // for the GM surface
        else
          // If the vertex is in the CSF region
          //if( (value <= threshold && wm_value>0 ) )
          //if( (value>=threshold && direction<0) || (value<=threshold && direction>0) )
        {
          //dx = RPoint_x(a->anchor_point) - parameter[p_index+0];
          //dy = RPoint_y(a->anchor_point) - parameter[p_index+1];
          //dz = RPoint_z(a->anchor_point) - parameter[p_index+2];
          /*dx = RPoint_x(anchor_points[point]) - parameter[p_index+0];
          dy = RPoint_y(anchor_points[point]) - parameter[p_index+1];
          dz = RPoint_z(anchor_points[point]) - parameter[p_index+2];
          dist = sqrt(dx*dx + dy*dy + dz*dz);
          fit1 += dist * max_weight;*/
          evaluate_volume_in_world( laplacian,
                                parameter[p_index+0],
                                parameter[p_index+1],
                                parameter[p_index+2],
                                EVAL_LINEAR, FALSE,
                                0.0, &value,
                                NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL );
          if(value<0){
            fit1 += -value/2 * max_weight;
          }
        }
      }
    }
    //end = clock();
    //printf("vol=%f,%f\n",((double)(end - start))/CLOCKS_PER_SEC, fit1);

    return fit1 + fit2;
}

private  Real   evaluate_volume_fit_deriv(
    Real         weight,
    Real         max_weight,
    int          direction,
    Volume       volume,
    Volume       laplacian,
    Real         threshold,
    Real         parameter[],
    int          n_neighbours[],
    int          *neighbours[],
    int          start_point,
    int          end_point,
    int          n_points,
    Real         deriv[] )
{
    int          point;
    int          p_index;
    Real         value, wm_value;
    anchor_point_struct *a;
    Real         dx, dy, dz, dist;
    Real         diff;
    Real         factor=1;
    int          neigh, n, max_neighbours;
    Point        *neigh_points;
    Vector       normal;
    Real         dxyz[3];

    if( weight > 0 || max_weight > 0 ){
      max_neighbours = 0;
      for_less( n, start_point, end_point ){
        max_neighbours = MAX( max_neighbours, n_neighbours[n] );
      }
      ALLOC( neigh_points, max_neighbours);

      for_less( point, start_point, end_point )
      {
        //a = &anchor_points[point];
        p_index = IJ( point, 0, 3 );
        if( direction < 0 )
        {
          evaluate_volume_in_world( volume,
                                parameter[p_index+0],
                                parameter[p_index+1],
                                parameter[p_index+2],
                                EVAL_LINEAR, FALSE,
                                0.0, &value,
                                NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL );
          if( value <= threshold )
          {
            for_less( n, 0, n_neighbours[point] ){
              neigh = neighbours[point][n];
              fill_Point( neigh_points[n],
                          parameter[IJ(neigh,0,3)],
                          parameter[IJ(neigh,1,3)],
                          parameter[IJ(neigh,2,3)] );
            }
            find_polygon_normal( n_neighbours[point], neigh_points, &normal);
            deriv[p_index+0] +=  -normal.coords[0] * max_weight * factor;
            deriv[p_index+1] +=  -normal.coords[1] * max_weight * factor;
            deriv[p_index+2] +=  -normal.coords[2] * max_weight * factor;
          }
        }
        else
        //if( ((value>=threshold && direction<0) || (value<=threshold && direction>0)) && wm_value==0 )
        {
          evaluate_volume_in_world( laplacian,
                                    parameter[p_index+0],
                                    parameter[p_index+1],
                                    parameter[p_index+2],
                                    EVAL_LINEAR, FALSE,
                                    0.0, &value,
                                    NULL, NULL, NULL,
                                    NULL, NULL, NULL, NULL, NULL, NULL );
          if( value < 0 )
          {
            for_less( n, 0, n_neighbours[point] ){
              neigh = neighbours[point][n];
              fill_Point( neigh_points[n],
                          parameter[IJ(neigh,0,3)],
                          parameter[IJ(neigh,1,3)],
                          parameter[IJ(neigh,2,3)] );
            }
            find_polygon_normal( n_neighbours[point], neigh_points, &normal);
            deriv[p_index+0] +=  normal.coords[0] * max_weight * factor;
            deriv[p_index+1] +=  normal.coords[1] * max_weight * factor;
            deriv[p_index+2] +=  normal.coords[2] * max_weight * factor;
          }
        }
      }
      FREE( neigh_points );
    }
}

private  Real   evaluate_boundary_search_fit(
    Real         image_weight_in,
    Real         image_weight_out,
    Real         max_inward,
    Real         max_outward,
    Real         max_dist_threshold,
    Real         max_dist_weight,
    int          oversample,
    Smallest_int boundary_flags[],
    Point        boundary_points[],
    int          start_point,
    int          end_point,
    int          n_points,
    Real         parameters[],
    Smallest_int active_flags[],
    int          n_neighbours[],
    int          *neighbours[],
    Real         *weights,
    Real         max_weight_value,
    Real         adaptive_ratio )
{
    int      point, p_index;
    Real     x, y, z, dx, dy, dz, dist_sq, fit2, weight;
    Real     dist, max_dist, max_diff2;
    Real     max_dist_sq;
    Real     adaptive_weight;

    if( max_dist_weight <= 0.0 || max_dist_threshold <= 0.0 )
        return( 0.0 );

    if( max_dist_weight > 0.0 && max_dist_threshold > 0.0 )
    {
        max_dist_sq = max_dist_threshold * max_dist_threshold;
        max_dist = MAX( max_outward, max_inward ) - max_dist_threshold;
        if( max_dist > 0.0 )
            max_diff2 = MAX( image_weight_out, image_weight_in ) *
                        max_dist * max_dist;
        else
            max_diff2 = 0.0;
    }
    else
    {
        max_dist_sq = 0.0;
        max_diff2 = 0.0;
    }

    fit2 = 0.0;

    if( active_flags ) {
      for_less( point, start_point, end_point ) {
        if( !active_flags[point] ) continue;

        if( boundary_flags[point] == BOUNDARY_NOT_FOUND ) {
          // Added by June Sic Kim 7/12/2002
	  if(adaptive_ratio==0){
	    fit2 += max_diff2;
	  }
	  else{
/*	    adaptive_weight = (max_weight_value-weights[point])/(max_weight_value/2);
            fit2 += max_diff2 * adaptive_weight * 
               (adaptive_weight>=1?ADAPTIVE_RATIO_BOUNDARY:1/ADAPTIVE_RATIO_BOUNDARY);
*/          adaptive_weight = (max_weight_value-weights[point])/max_weight_value;
	    fit2 += max_diff2 * adaptive_weight * adaptive_ratio;
	  }
            continue;
        } else if( boundary_flags[point] == BOUNDARY_IS_OUTSIDE )
            weight = image_weight_out;
        else
            weight = image_weight_in;

        p_index = IJ( point, 0, 3 );

        x = parameters[p_index+0];
        y = parameters[p_index+1];
        z = parameters[p_index+2];

        dx = x - RPoint_x( boundary_points[point] );
        dy = y - RPoint_y( boundary_points[point] );
        dz = z - RPoint_z( boundary_points[point] );

        dist_sq = dx * dx + dy * dy + dz * dz;

        if( max_dist_sq > 0.0 && dist_sq > max_dist_sq ) {
            dist = sqrt( dist_sq ) - max_dist_threshold;
	    if(adaptive_ratio==0){
	      fit2 += weight * dist * dist;
	    }
	    else{
/*	      adaptive_weight = (max_weight_value-weights[point])/(max_weight_value/2);
	      fit2 += weight * dist * dist * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_BOUNDARY:1/ADAPTIVE_RATIO_BOUNDARY);
*/            adaptive_weight = (max_weight_value-weights[point])/max_weight_value;
	      fit2 += weight * dist * dist * adaptive_weight * adaptive_ratio;
	    }
        }
      }
    }

    if( oversample > 0 ) {
        handle_internal_error( "gotta implement this" );
#ifdef NOT_YET
        ind = n_points;

        if( start_point > 0 )
            handle_internal_error( "gotta implement this" );

        for_less( point, start_point, end_point ) {
            p_index = IJ( point, 0, 3 );
            x = parameters[p_index+0];
            y = parameters[p_index+1];
            z = parameters[p_index+2];

            neigh_ptr = neighbours[point];
            n_neighs = n_neighbours[point];

            for_less( n, 0, n_neighs )
            {
                neigh1 = neigh_ptr[n];
                neigh2 = neigh_ptr[(n+1) % n_neighs];

                if( point > neigh1 || point > neigh2 )
                    continue;

                n1_index = IJ( neigh1, 0, 3 );
                n2_index = IJ( neigh2, 0, 3 );

                x1 = parameters[n1_index+0];
                y1 = parameters[n1_index+1];
                z1 = parameters[n1_index+2];

                x2 = parameters[n2_index+0];
                y2 = parameters[n2_index+1];
                z2 = parameters[n2_index+2];

                delta_x = (x2 - x) / (Real) (oversample + 1);
                delta_y = (y2 - y) / (Real) (oversample + 1);
                delta_z = (z2 - z) / (Real) (oversample + 1);

                x_int = x;
                y_int = y;
                z_int = z;

                for_less( w, 0, oversample )
                {
                    x_int += delta_x;
                    y_int += delta_y;
                    z_int += delta_z;

                    if( boundary_flags[ind] == BOUNDARY_NOT_FOUND )
                    {
		      if(adaptive_ratio==0){
			fit2+= max_diff2;
		      }
		      else{
/*			adaptive_weight = (max_weight_value-weights[point])/(max_weight_value/2);
                        fit2 += max_diff2 * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_BOUNDARY:1/ADAPTIVE_RATIO_BOUNDARY);
*/                      adaptive_weight = (max_weight_value-weights[point])/max_weight_value;
			fit2 += max_diff2 * adaptive_weight * adaptive_ratio;
                        ++ind;
		      }
                        continue;
                    }
                    else if( boundary_flags[ind] == BOUNDARY_IS_OUTSIDE )
                        weight = image_weight_out;
                    else
                        weight = image_weight_in;

                    bound_point = &boundary_points[ind];

                    dx = x_int - RPoint_x( *bound_point );
                    dy = y_int - RPoint_y( *bound_point );
                    dz = z_int - RPoint_z( *bound_point );

                    dist_sq = dx * dx + dy * dy + dz * dz;

                    ++ind;

                    if( max_dist_sq > 0.0 && dist_sq > max_dist_sq )
                    {
                        dist = sqrt( dist_sq ) - max_dist_threshold;
			if(adaptive_ratio==0){
			  fit2 += weight * dist * dist;
			}
			else{
/*			  adaptive_weight = (max_weight_value-weights[point])/(max_weight_value/2);
			  fit2 += weight * dist * dist * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_BOUNDARY:1/ADAPTIVE_RATIO_BOUNDARY);
*/                        adaptive_weight = (max_weight_value-weights[point])/max_weight_value;
			  fit2 += weight * dist * dist * adaptive_weight * adaptive_ratio;
			}
		    }
                }
            }
        }
#endif
    }

    return( max_dist_weight * fit2 );
}

private  void   evaluate_boundary_search_fit_deriv(
    Real         image_weight_in,
    Real         image_weight_out,
    Real         max_dist_threshold,
    Real         max_dist_weight,
    int          oversample,
    Smallest_int boundary_flags[],
    Point        boundary_points[],
    int          start_point,
    int          end_point,
    int          n_points,
    Real         parameters[],
    int          n_neighbours[],
    int          *neighbours[],
    Real         deriv[],
    Real         *weights,
    Real         max_weight_value,
    Real         adaptive_ratio )
{
    int    point, p_index;
    Real   x, y, z, dx, dy, dz, weight;
    Real   factor, diff, dist, dist_sq, max_dist_sq;
    Real   adaptive_weight;

    if( image_weight_in < 0.0 && image_weight_out < 0.0 ||
        max_dist_weight <= 0.0 || max_dist_threshold < 0.0 )
        return;

    max_dist_sq = max_dist_threshold * max_dist_threshold;

    for_less( point, start_point, end_point )
    {
        if( boundary_flags[point] == BOUNDARY_NOT_FOUND )
            continue;
        else if( boundary_flags[point] == BOUNDARY_IS_OUTSIDE )
            weight = image_weight_out;
        else
            weight = image_weight_in;

        p_index = IJ( point, 0, 3 );

        x = parameters[p_index+0];
        y = parameters[p_index+1];
        z = parameters[p_index+2];

        dx = x - RPoint_x( boundary_points[point] );
        dy = y - RPoint_y( boundary_points[point] );
        dz = z - RPoint_z( boundary_points[point] );

        dist_sq = dx * dx + dy * dy + dz * dz;
        if( dist_sq > max_dist_sq ) {
            dist = sqrt( dist_sq );
            diff = dist - max_dist_threshold;
            // Modified by June Sic Kim at 9/12/2002
	    if(adaptive_ratio==0){
	      factor = weight * max_dist_weight * 2.0 * diff / dist;
	    } else{
/*	      adaptive_weight = (max_weight_value-weights[point])/(max_weight_value/2);
	      factor = weight * max_dist_weight * 2.0 * diff / dist * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_BOUNDARY:1/ADAPTIVE_RATIO_BOUNDARY);
*/            adaptive_weight = (max_weight_value-weights[point])/(max_weight_value);
	      factor = weight * max_dist_weight * 2.0 * diff / dist * adaptive_weight * adaptive_ratio;
	    }
            deriv[p_index+0] += factor * dx;
            deriv[p_index+1] += factor * dy;
            deriv[p_index+2] += factor * dz;
        }
    }

    if( oversample > 0 ) {
        handle_internal_error( "not yet" );
#ifdef NOT_YET
        ind = n_points;
        for_less( point, 0, start_point ) {
            for_less( n, 0, n_neighbours[point] ) {
                neigh = neighbours[point][n];
                if( THIS_IS_UNIQUE_EDGE( point, neigh ) )
                    ind += oversample;
            }
        }

        for_less( point, start_point, end_point ) {
            p_index = IJ( point, 0, 3 );
            x = parameters[p_index+0];
            y = parameters[p_index+1];
            z = parameters[p_index+2];

            for_less( n, 0, n_neighbours[point] ) {
                neigh = neighbours[point][n];
                if( !THIS_IS_UNIQUE_EDGE( point, neigh ) )
                    continue;

                n_index = IJ( neigh, 0, 3 );
                x2 = parameters[n_index+0];
                y2 = parameters[n_index+1];
                z2 = parameters[n_index+2];

                for_less( w, 0, oversample ) {
                    if( boundary_flags[ind] == BOUNDARY_NOT_FOUND ) {
                        ++ind;
                        continue;
                    } else if( boundary_flags[ind] == BOUNDARY_IS_OUTSIDE )
                        weight = image_weight_out;
                    else
                        weight = image_weight_in;

                    alpha = (Real) (w+1) / (Real) (oversample+1);

                    dx = (1.0 - alpha) * x + alpha * x2 -
                         RPoint_x( boundary_points[ind] );
                    dy = (1.0 - alpha) * y + alpha * y2 -
                         RPoint_y( boundary_points[ind] );
                    dz = (1.0 - alpha) * z + alpha * z2 -
                         RPoint_z( boundary_points[ind] );

                    dist_sq = dx * dx + dy * dy + dz * dz;
                    if( dist_sq > max_dist_sq ) {
                        dist = sqrt( dist_sq );
                        diff = dist - max_dist_threshold;
                        // Modified by June Sic Kim at 9/12/2002
			if(adaptive_ratio==0){
			  factor = weight * max_dist_weight * 2.0 * diff / dist;
			}
			else{
/*			  adaptive_weight = (max_weight_value-weights[point])/(max_weight_value/2);
			  factor = weight * max_dist_weight * 2.0 * diff / dist * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_BOUNDARY:1/ADAPTIVE_RATIO_BOUNDARY);
*/                        adaptive_weight = (max_weight_value-weights[point])/max_weight_value;
			  factor = weight * max_dist_weight * 2.0 * diff / dist * adaptive_weight * adaptive_ratio;
			}
                        deriv[p_index+0] += factor * (1.0 - alpha) * dx;
                        deriv[p_index+1] += factor * (1.0 - alpha) * dy;
                        deriv[p_index+2] += factor * (1.0 - alpha) * dz;
                        deriv[n_index+0] += factor * alpha * dx;
                        deriv[n_index+1] += factor * alpha * dy;
                        deriv[n_index+2] += factor * alpha * dz;
                    }

                    ++ind;
                }
            }
        }
#endif
    }
}

private  Real    trilinear_interpolate(
    Volume               volume,
    voxel_coef_struct    *lookup,
    Real                 voxel[] )
{
    int    c, i, j, k, *sizes, dx, dy, dz;
    Real   x, y, z, u, v, w, u1, v1, w1, value;
    Real   coefs[8];

    sizes = &volume->array.sizes[0];

    x = voxel[0];
    y = voxel[1];
    z = voxel[2];

    if( x >= 0.0 && x < (Real) sizes[0]-1.0 &&
        y >= 0.0 && y < (Real) sizes[1]-1.0 &&
        z >= 0.0 && z < (Real) sizes[2]-1.0 )
    {
        i = (int) x;
        j = (int) y;
        k = (int) z;

        lookup_volume_coeficients( lookup, i, j, k, coefs );
    }
    else
    {
        i = FLOOR( x );
        j = FLOOR( y );
        k = FLOOR( z );

        c = 0;
        for_less( dx, 0, 2 )
        for_less( dy, 0, 2 )
        for_less( dz, 0, 2 )
        {
            if( i + dx >= 0 && i + dx < sizes[0] &&
                j + dy >= 0 && j + dy < sizes[1] &&
                k + dz >= 0 && k + dz < sizes[2] )
            {
                coefs[c] = lookup->voxel_to_real_values[lookup->voxel_ptr[i+dx]
                                                          [j+dy][k+dz]];
            }
            else
                coefs[c] = 0.0;
            ++c;
        }
    }

    u = x - (Real) i;
    u1 = 1.0 - u;
    v = y - (Real) j;
    v1 = 1.0 - v;
    w = z - (Real) k;
    w1 = 1.0 - w;

    /*--- if the value is desired, interpolate in 1D to get the value */

    value = u1 * (v1 * (w1 * coefs[0] + w * coefs[1]) +
                   v * (w1 * coefs[2] + w * coefs[3])) +
             u * (v1 * (w1 * coefs[4] + w * coefs[5]) +
                   v * (w1 * coefs[6] + w * coefs[7]));

    return( value );
}

private  Real    trilinear_interpolate_with_deriv(
    Volume               volume,
    BOOLEAN              inside_flag,
    Real                 trilin_sizes[],
    voxel_coef_struct    *lookup,
    Real                 voxel[],
    Real                 deriv[] )
{
    int    c, i, j, k, *sizes, dx, dy, dz;
    Real   x, y, z, u, v, w, value;
    Real   coefs[8];
    Real   du00, du01, du10, du11, c00, c01, c10, c11, c0, c1, du0, du1;
    Real   dv0, dv1, dw;

    x = voxel[0];
    y = voxel[1];
    z = voxel[2];

    if( inside_flag || x >= 0.0 && x < trilin_sizes[0] &&
                       y >= 0.0 && y < trilin_sizes[1] &&
                       z >= 0.0 && z < trilin_sizes[2] )
    {
        i = (int) x;
        j = (int) y;
        k = (int) z;

        lookup_volume_coeficients( lookup, i, j, k, coefs );
    }
    else
    {
        i = FLOOR( x );
        j = FLOOR( y );
        k = FLOOR( z );

        sizes = &volume->array.sizes[0];

        c = 0;
        for_less( dx, 0, 2 )
        for_less( dy, 0, 2 )
        for_less( dz, 0, 2 )
        {
            if( i + dx >= 0 && i + dx < sizes[0] &&
                j + dy >= 0 && j + dy < sizes[1] &&
                k + dz >= 0 && k + dz < sizes[2] )
            {
                coefs[c] = lookup->voxel_to_real_values[lookup->voxel_ptr[i+dx]
                                                          [j+dy][k+dz]];
            }
            else
                coefs[c] = 0.0;
            ++c;
        }
    }

    u = x - (Real) i;
    v = y - (Real) j;
    w = z - (Real) k;

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

    /*--- if the value is desired, interpolate in 1D to get the value */

    value = c0 + w * dw;

    du0 = INTERPOLATE( v, du00, du10 );
    du1 = INTERPOLATE( v, du01, du11 );

    deriv[X] = INTERPOLATE( w, du0, du1 );
    deriv[Y] = INTERPOLATE( w, dv0, dv1 );
    deriv[Z] = dw;

    return( value );
}

private  Real   differential_weights(
    Real   diff,
    Real   min_diff,
    Real   max_diff,
    Real   diff_offset,
    Real   ln_start_weight )
{
    Real   fit;

    if( diff <= min_diff )
        diff = min_diff - diff;
    else
        diff = diff - max_diff;

    if( diff > diff_offset )
        fit = diff * diff;
    else
        fit = exp( ln_start_weight + (-ln_start_weight) *
                   diff / diff_offset ) * diff * diff;

    return( fit );
}

private  Real   differential_weights_deriv(
    Real   diff,
    Real   min_diff,
    Real   max_diff,
    Real   diff_offset,
    Real   ln_start_weight )
{
    Real   deriv, sign, ln_weight, weight;

    if( diff <= min_diff )
    {
        diff = min_diff - diff;
        sign = -1.0;
    }
    else
    {
        diff = diff - max_diff;
        sign = 1.0;
    }

    if( diff > diff_offset )
    {
        deriv = sign * 2.0 * diff;
    }
    else
    {
        ln_weight = ln_start_weight + (-ln_start_weight) * diff/
                    diff_offset;
        weight = exp( ln_weight );
        deriv = weight * sign * 2.0 * diff +
                weight * sign * (-ln_start_weight) / diff_offset *
                             diff * diff;
    }

    return( deriv );
}

typedef  struct
{
    Real      alpha1;
    Real      alpha2;
}
weighting_struct;

private  Real   evaluate_gradient_fit(
    Volume               volume,
    voxel_coef_struct    *lookup,
    int                  continuity,
    Real                 image_weight,
    Real                 threshold,
    Real                 min_diff,
    Real                 max_diff,
    Real                 max_diff_weight,
    Real                 differential_offset,
    Real                 differential_ratio,
    int                  n_points,
    Real                 parameters[],
    Smallest_int         active_flags[] ) {

    int        p1, p1_index;
    Real       x, y, z, value, diff, fit1, fit2;
    Real       ln_start_weight;
    Real       wtv00, wtv01, wtv02, wtv03;
    Real       wtv10, wtv11, wtv12, wtv13;
    Real       wtv20, wtv21, wtv22, wtv23;
    Real       voxel[3], voxel_000[3], voxel_100[3], voxel_010[3], voxel_001[3];

    if( image_weight < 0.0 )
        return( 0.0 );

    fit1 = 0.0;
    fit2 = 0.0;

    if( active_flags != NULL )
        handle_internal_error( "Need to implement this" );

    convert_world_to_voxel( volume, 0.0, 0.0, 0.0, voxel_000 );
    convert_world_to_voxel( volume, 1.0, 0.0, 0.0, voxel_100 );
    convert_world_to_voxel( volume, 0.0, 1.0, 0.0, voxel_010 );
    convert_world_to_voxel( volume, 0.0, 0.0, 1.0, voxel_001 );

    wtv00 = voxel_100[0] - voxel_000[0];
    wtv01 = voxel_010[0] - voxel_000[0];
    wtv02 = voxel_001[0] - voxel_000[0];
    wtv03 = voxel_000[0];

    wtv10 = voxel_100[1] - voxel_000[1];
    wtv11 = voxel_010[1] - voxel_000[1];
    wtv12 = voxel_001[1] - voxel_000[1];
    wtv13 = voxel_000[1];

    wtv20 = voxel_100[2] - voxel_000[2];
    wtv21 = voxel_010[2] - voxel_000[2];
    wtv22 = voxel_001[2] - voxel_000[2];
    wtv23 = voxel_000[2];

    ln_start_weight = log( differential_ratio );

    for_less( p1, 0, n_points ) {
        p1_index = IJ( p1, 0, 3 );

        x = parameters[p1_index+0];
        y = parameters[p1_index+1];
        z = parameters[p1_index+2];

        voxel[0] = wtv00 * x + wtv01 * y + wtv02 * z + wtv03;
        voxel[1] = wtv10 * x + wtv11 * y + wtv12 * z + wtv13;
        voxel[2] = wtv20 * x + wtv21 * y + wtv22 * z + wtv23;

        if( continuity == 0 ) {
            value = trilinear_interpolate( volume, lookup, voxel );
        } else {
            (void) evaluate_volume( volume, voxel, NULL, continuity,
                                    FALSE, 0.0, &value, NULL, NULL );
        }

        //////////////////////////////////////////////////////////////////
        // modified by June S. Kim at 7/MAY/2004
        if( value < threshold ) {
          diff = value - threshold;
          fit1 += diff * diff;
        }
    }

    return fit1*image_weight + fit2;
}

private Real evaluate_gw_gradient_fit( Volume t1,
                                       Real weight,
                                       int oversample,
                                       int * mask,
                                       int start_point,
                                       int end_point,
                                       int   n_ngh[],
                                       int * ngh[],
                                       Real * coords,
                                       Real * t1grad ) {

  int          n, point;
  Real         fit = 0.0f;

  if( weight > 0 && t1grad ) {
    Real value;
    for( point = start_point; point < end_point; point++ ) {
      if( !mask[point-start_point] ) {
        evaluate_volume_in_world( t1, coords[3*point+0],
                                  coords[3*point+1], coords[3*point+2],
                                  EVAL_LINEAR, FALSE,
                                  0.0, &value,
                                  NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL );
        fit += ( value / t1grad[point-start_point] - 1.0 ) *
               ( value / t1grad[point-start_point] - 1.0 );
      }
    }

    // Do oversampling. There is only one level of oversampling on edges.

    if( oversample ) {

      for( point = start_point; point < end_point; point++ ) {
        if( mask[point-start_point] ) continue;

        for_less( n, 0, n_ngh[point] ) {
          int neigh1 = ngh[point][n];
          if( point > neigh1 ) continue;

          Real xmid = 0.5 * ( coords[3*point] + coords[3*neigh1] );
          Real ymid = 0.5 * ( coords[3*point+1] + coords[3*neigh1+1] );
          Real zmid = 0.5 * ( coords[3*point+2] + coords[3*neigh1+2] );
          Real t1mid = 0.5 * ( t1grad[point-start_point] + t1grad[neigh1-start_point] );

          evaluate_volume_in_world( t1, xmid, ymid, zmid,
                                    EVAL_LINEAR, FALSE,
                                    0.0, &value,
                                    NULL, NULL, NULL,
                                    NULL, NULL, NULL, NULL, NULL, NULL );
          fit += ( value / t1mid - 1.0 ) * ( value / t1mid - 1.0 );
        }
      }
      weight *= 0.25;   // this is very very close to the correct factor N/(4*N-6).
    }

    fit *= 0.5 * weight;
  }
  return( fit );
}


private  Real   evaluate_image_value_fit(
    Volume               volume,
    voxel_coef_struct    *lookup,
    int                  continuity,
    Real                 image_weight,
    Real                 threshold,
    Real                 min_diff,
    Real                 max_diff,
    Real                 max_diff_weight,
    Real                 differential_offset,
    Real                 differential_ratio,
    int                  oversample,
    int                  n_points,
    int                  n_neighbours[],
    int                  *neighbours[],
    Real                 parameters[],
    Smallest_int         active_flags[] )
{
    int        p1, p1_index, n1_index, n2_index, neigh1, neigh2, w1, w2, n;
    int        i, j, k, sizes[N_DIMENSIONS];
    int        n_neighs, *neighs;
    Real       x, y, z, value, diff, fit1, fit2;
    Real       ln_start_weight;
    Real       x1, y1, z1, x2, y2, z2, x3, y3, z3;
    Real       xv1, yv1, zv1, xv2, yv2, zv2, xv3, yv3, zv3;
    Real       dx12, dy12, dz12, dx13, dy13, dz13;
    Real       alpha1, alpha2;
    Real       wtv00, wtv01, wtv02, wtv03;
    Real       wtv10, wtv11, wtv12, wtv13;
    Real       wtv20, wtv21, wtv22, wtv23;
    Real       voxel[3], voxel_000[3], voxel_100[3], voxel_010[3];
    Real       voxel_001[3];
    Real       u, v, w, u1, v1, ww1;
    Real       real_scale, real_translation;
    unsigned char  *ptr;
    BOOLEAN    inside;
    int               n_weights, n_weights1, n_weights2;
    weighting_struct  tt, *weights1, *weights2, *weights;
    unsigned  char  ***voxel_ptr;
    int             offset1, offset2, offset3, offset4;
    int             offset5, offset6, offset7;
    BOOLEAN    checking_differential, fast_interp;

    if( image_weight < 0.0 )
        return( 0.0 );

    get_volume_sizes( volume, sizes );

    fit1 = 0.0;
    fit2 = 0.0;

    if( active_flags != NULL )
        handle_internal_error( "Need to implement this" );

    convert_world_to_voxel( volume, 0.0, 0.0, 0.0, voxel_000 );
    convert_world_to_voxel( volume, 1.0, 0.0, 0.0, voxel_100 );
    convert_world_to_voxel( volume, 0.0, 1.0, 0.0, voxel_010 );
    convert_world_to_voxel( volume, 0.0, 0.0, 1.0, voxel_001 );

    wtv00 = voxel_100[0] - voxel_000[0];
    wtv01 = voxel_010[0] - voxel_000[0];
    wtv02 = voxel_001[0] - voxel_000[0];
    wtv03 = voxel_000[0];

    wtv10 = voxel_100[1] - voxel_000[1];
    wtv11 = voxel_010[1] - voxel_000[1];
    wtv12 = voxel_001[1] - voxel_000[1];
    wtv13 = voxel_000[1];

    wtv20 = voxel_100[2] - voxel_000[2];
    wtv21 = voxel_010[2] - voxel_000[2];
    wtv22 = voxel_001[2] - voxel_000[2];
    wtv23 = voxel_000[2];

    ln_start_weight = log( differential_ratio );

    checking_differential = (min_diff < max_diff && max_diff_weight > 0.0);

    real_translation = convert_voxel_to_value( volume, 0.0 );
    real_scale = convert_voxel_to_value( volume, 1.0 ) - real_translation;
    real_translation -= threshold;

    for_less( p1, 0, n_points ) {
        p1_index = IJ( p1, 0, 3 );

        x = parameters[p1_index+0];
        y = parameters[p1_index+1];
        z = parameters[p1_index+2];

        voxel[0] = wtv00 * x + wtv01 * y + wtv02 * z + wtv03;
        voxel[1] = wtv10 * x + wtv11 * y + wtv12 * z + wtv13;
        voxel[2] = wtv20 * x + wtv21 * y + wtv22 * z + wtv23;

        if( continuity == 0 )
            value = trilinear_interpolate( volume, lookup, voxel );
        else
        {
            (void) evaluate_volume( volume, voxel, NULL, continuity,
                                    FALSE, 0.0, &value, NULL, NULL );
        }

        diff = value - threshold;

        fit1 += diff * diff;
        if( !checking_differential || diff >= min_diff && diff <= max_diff )
            continue;

        fit2 += differential_weights( diff, min_diff, max_diff,
                                      differential_offset,
                                      ln_start_weight );
    }

    if( oversample > 0 )
    {
        voxel_ptr = lookup->voxel_ptr;
        offset1 = lookup->offset1;
        offset2 = lookup->offset2;
        offset3 = lookup->offset3;
        offset4 = lookup->offset4;
        offset5 = lookup->offset5;
        offset6 = lookup->offset6;
        offset7 = lookup->offset7;

        inside = TRUE;
        for_less( p1, 0, n_points )
        {
            p1_index = IJ(p1,0,3);
            x1 = parameters[p1_index+0];
            y1 = parameters[p1_index+1];
            z1 = parameters[p1_index+2];

            xv1 = wtv00 * x1 + wtv01 * y1 + wtv02 * z1 + wtv03;
            yv1 = wtv10 * x1 + wtv11 * y1 + wtv12 * z1 + wtv13;
            zv1 = wtv20 * x1 + wtv21 * y1 + wtv22 * z1 + wtv23;

            if( xv1 < 0.0 || xv1 >= (Real) sizes[0]-1.0 ||
                yv1 < 0.0 || yv1 >= (Real) sizes[1]-1.0 ||
                zv1 < 0.0 || zv1 >= (Real) sizes[2]-1.0 )
            {
                inside = FALSE;

                break;
            }
        }

        fast_interp = (continuity == 0 && inside);

        n_weights1 = 0;
        n_weights2 = 0;
        weights1 = NULL;
        weights2 = NULL;

        for_inclusive( w1, 0, oversample )
        {
            tt.alpha1 = (Real) w1 / (Real) (oversample+1);
            for_inclusive( w2, 1, oversample - w1 + 1 )
            {
                tt.alpha2 = (1.0 - tt.alpha1) *
                            (Real) w2 / (Real) (oversample-w1+1);


                if( w2 < oversample + 1 )
                {
                    ADD_ELEMENT_TO_ARRAY( weights1, n_weights1, tt,
                                          DEFAULT_CHUNK_SIZE );
                }

                if( w2 < oversample - w1 + 1 )
                {
                    ADD_ELEMENT_TO_ARRAY( weights2, n_weights2, tt,
                                          DEFAULT_CHUNK_SIZE );
                }
            }
        }

        for_less( p1, 0, n_points )
        {
            p1_index = IJ(p1,0,3);
            x1 = parameters[p1_index+0];
            y1 = parameters[p1_index+1];
            z1 = parameters[p1_index+2];

            xv1 = wtv00 * x1 + wtv01 * y1 + wtv02 * z1 + wtv03;
            yv1 = wtv10 * x1 + wtv11 * y1 + wtv12 * z1 + wtv13;
            zv1 = wtv20 * x1 + wtv21 * y1 + wtv22 * z1 + wtv23;

            n_neighs = n_neighbours[p1];
            neighs = neighbours[p1];
            for_less( n, 0, n_neighs )
            {
                neigh1 = neighs[n];
                neigh2 = neighs[(n+1)%n_neighs];

                if( p1 > neigh1 || p1 > neigh2 )
                    continue;

                n1_index = IJ(neigh1,0,3);
                n2_index = IJ(neigh2,0,3);

                x2 = parameters[n1_index+0];
                y2 = parameters[n1_index+1];
                z2 = parameters[n1_index+2];
                x3 = parameters[n2_index+0];
                y3 = parameters[n2_index+1];
                z3 = parameters[n2_index+2];

                xv2 = wtv00 * x2 + wtv01 * y2 + wtv02 * z2 + wtv03;
                yv2 = wtv10 * x2 + wtv11 * y2 + wtv12 * z2 + wtv13;
                zv2 = wtv20 * x2 + wtv21 * y2 + wtv22 * z2 + wtv23;

                xv3 = wtv00 * x3 + wtv01 * y3 + wtv02 * z3 + wtv03;
                yv3 = wtv10 * x3 + wtv11 * y3 + wtv12 * z3 + wtv13;
                zv3 = wtv20 * x3 + wtv21 * y3 + wtv22 * z3 + wtv23;

                dx12 = xv2 - xv1;
                dy12 = yv2 - yv1;
                dz12 = zv2 - zv1;

                dx13 = xv3 - xv1;
                dy13 = yv3 - yv1;
                dz13 = zv3 - zv1;

                if( neigh1 < neigh2 )
                {
                    n_weights = n_weights1;
                    weights = weights1;
                }
                else
                {
                    n_weights = n_weights2;
                    weights = weights2;
                }

                if( fast_interp )
                {
                for_less( w1, 0, n_weights )
                {
                    alpha1 = weights[w1].alpha1;
                    alpha2 = weights[w1].alpha2;

                    x = xv1 + alpha1 * dx12 + alpha2 * dx13;
                    y = yv1 + alpha1 * dy12 + alpha2 * dy13;
                    z = zv1 + alpha1 * dz12 + alpha2 * dz13;

                    i = (int) x;
                    j = (int) y;
                    k = (int) z;

                    u = x - (Real) i;
                    u1 = 1.0 - u;
                    v = y - (Real) j;
                    v1 = 1.0 - v;
                    w = z - (Real) k;
                    ww1 = 1.0 - w;

                    ptr = &voxel_ptr[i][j][k];

                    diff = real_scale *
                            (u1 * (v1 * (ww1 * (Real) *ptr +
                               w * (Real) ptr[offset1]) +
                         v * (ww1 * (Real) ptr[offset2] +
                              w * (Real) ptr[offset3])) +

                         u*(
                          v1 * (ww1 * (Real) ptr[offset4]+
                               w * (Real) ptr[offset5]) +
                         v * (ww1 * (Real) ptr[offset6] +
                              w * (Real) ptr[offset7]))) +
                             real_translation;

                    fit1 += diff * diff;

                    if( checking_differential &&
                        (diff < min_diff || diff > max_diff) )
                    {
                        fit2 += differential_weights( diff, min_diff, max_diff,
                                                      differential_offset,
                                                      ln_start_weight );
                    }
                }
                }
                else
                {
                for_less( w1, 0, n_weights )
                {
                    alpha1 = weights[w1].alpha1;
                    alpha2 = weights[w1].alpha2;

                    voxel[0] = xv1 + alpha1 * dx12 + alpha2 * dx13;
                    voxel[1] = yv1 + alpha1 * dy12 + alpha2 * dy13;
                    voxel[2] = zv1 + alpha1 * dz12 + alpha2 * dz13;

                    if( continuity == 0 )
                    {
                        value = trilinear_interpolate( volume, lookup,
                                                       voxel );
                    }
                    else
                    {
                        (void) evaluate_volume( volume, voxel, NULL,
                                                continuity, FALSE, 0.0,
                                                &value, NULL, NULL);
                    }

                    diff = value - threshold;

                    fit1 += diff * diff;

                    if( checking_differential &&
                        (diff < min_diff || diff > max_diff) )
                    {
                        fit2 += differential_weights( diff, min_diff, max_diff,
                                                      differential_offset,
                                                      ln_start_weight );
                    }
                }
                }
            }
        }

        FREE( weights1 );
        FREE( weights2 );
    }

    fit1 *= image_weight;
    fit2 *= max_diff_weight;

    return( fit1 + fit2 );
}

private  void   evaluate_gradient_fit_deriv(
    Volume               volume,
    voxel_coef_struct    *lookup,
    int                  continuity,
    Real                 image_weight,
    Real                 threshold,
    Real                 differential_ratio,
    int                  n_points,
    Real                 parameters[],
    Real                 deriv[] )
{
    int     p1, p1_index;
    Real    x, y, z, value, diff, dx, dy, dz;
    Real    value_deriv[N_DIMENSIONS], *value_deriv_ptr[1];
    Real    wtv00, wtv01, wtv02, wtv03;
    Real    wtv10, wtv11, wtv12, wtv13;
    Real    wtv20, wtv21, wtv22, wtv23;
    Real    voxel[3], voxel_000[3], voxel_100[3], voxel_010[3], voxel_001[3];
    Real    trilin_sizes[N_DIMENSIONS];
    int     sizes[N_DIMENSIONS];
    Real    xv1, yv1, zv1, xv2, yv2, zv2, xv3, yv3, zv3;
    BOOLEAN inside;

    if( image_weight < 0.0 )
        return;

    Real    two_image_weight = image_weight + image_weight;

    get_volume_sizes( volume, sizes );
    trilin_sizes[0] = (Real) sizes[0] - 1.0;
    trilin_sizes[1] = (Real) sizes[1] - 1.0;
    trilin_sizes[2] = (Real) sizes[2] - 1.0;

    convert_world_to_voxel( volume, 0.0, 0.0, 0.0, voxel_000 );
    convert_world_to_voxel( volume, 1.0, 0.0, 0.0, voxel_100 );
    convert_world_to_voxel( volume, 0.0, 1.0, 0.0, voxel_010 );
    convert_world_to_voxel( volume, 0.0, 0.0, 1.0, voxel_001 );

    wtv00 = voxel_100[0] - voxel_000[0];
    wtv01 = voxel_010[0] - voxel_000[0];
    wtv02 = voxel_001[0] - voxel_000[0];
    wtv03 = voxel_000[0];

    wtv10 = voxel_100[1] - voxel_000[1];
    wtv11 = voxel_010[1] - voxel_000[1];
    wtv12 = voxel_001[1] - voxel_000[1];
    wtv13 = voxel_000[1];

    wtv20 = voxel_100[2] - voxel_000[2];
    wtv21 = voxel_010[2] - voxel_000[2];
    wtv22 = voxel_001[2] - voxel_000[2];
    wtv23 = voxel_000[2];

    inside = TRUE;
    for_less( p1, 0, n_points ) {
        p1_index = IJ(p1,0,3);
        x = parameters[p1_index+0];
        y = parameters[p1_index+1];
        z = parameters[p1_index+2];

        xv1 = wtv00 * x + wtv01 * y + wtv02 * z + wtv03;
        yv1 = wtv10 * x + wtv11 * y + wtv12 * z + wtv13;
        zv1 = wtv20 * x + wtv21 * y + wtv22 * z + wtv23;

        if( xv1 < 0.0 || xv1 >= (Real) sizes[0]-1.0 ||
            yv1 < 0.0 || yv1 >= (Real) sizes[1]-1.0 ||
            zv1 < 0.0 || zv1 >= (Real) sizes[2]-1.0 ) {
            inside = FALSE;

            break;
        }
    }

    for_less( p1, 0, n_points ) {
        p1_index = IJ( p1, 0, 3 );

        x = parameters[p1_index+0];
        y = parameters[p1_index+1];
        z = parameters[p1_index+2];

        voxel[0] = wtv00 * x + wtv01 * y + wtv02 * z + wtv03;
        voxel[1] = wtv10 * x + wtv11 * y + wtv12 * z + wtv13;
        voxel[2] = wtv20 * x + wtv21 * y + wtv22 * z + wtv23;

        if( continuity == 0 ) {
            value = trilinear_interpolate_with_deriv( volume, inside,
                                                      trilin_sizes,
                                                      lookup, voxel,
                                                      value_deriv );
        } else {
            value_deriv_ptr[0] = value_deriv;
            (void) evaluate_volume( volume, voxel, NULL, continuity,
                                    FALSE, 0.0, &value, value_deriv_ptr, NULL );
        }

        convert_voxel_normal_vector_to_world( volume, value_deriv,
                                              &dx, &dy, &dz );

        if( value < threshold ) {
          diff = two_image_weight * ( value - threshold );
          deriv[p1_index+0] += dx * diff;
          deriv[p1_index+1] += dy * diff;
          deriv[p1_index+2] += dz * diff;
        }
    }
}


private void evaluate_gw_gradient_fit_deriv( Volume t1,
                                             Real weight,
                                             int oversample,
                                             int image_type,
                                             int * mask,
                                             int start_point,
                                             int end_point,
                                             int   n_ngh[],
                                             int * ngh[],
                                             Real * coords,
                                             Real * t1grad,
                                             Real deriv[] ) {

  enum { UNKNOWN, T1, HISTO, T2 };

  int          n, point;

  if( weight > 0 && t1grad ) {
    if( oversample ) {
      weight *= 0.25;   // this is very very close to the correct factor N/(4*N-6).
    }

    int max_neighbours = 0;
    for( point = start_point; point < end_point; point++ ) {
      if( n_ngh[point] > max_neighbours ) max_neighbours = n_ngh[point];
    }
    Point * neigh_points = (Point *)malloc( max_neighbours * sizeof( Point ) );

    for( point = start_point; point < end_point; point++ ) {
      if( !mask[point-start_point] ) {
        Real value, dx, dy, dz;
        evaluate_volume_in_world( t1, coords[3*point+0],
                                  coords[3*point+1], coords[3*point+2],
                                  EVAL_LINEAR, FALSE,
                                  0.0, &value, &dx, &dy, &dz,
                                  NULL, NULL, NULL, NULL, NULL, NULL );

        // Make sure that the expected direction of the t1-gradient
        // is in a coherent direction with the outward normal vector.
        // For example, on in-vivo MRI, with CSF < GM < WM, the
        // t1 gradient points inwards whereas the normal points
        // outwards. The converse is true for histology. In case
        // there appears to be some incompatibility, do nothing.

        for( n = 0; n < n_ngh[point]; n++ ) {
          int neigh = ngh[point][n];
          fill_Point( neigh_points[n], coords[3*neigh+0],
                      coords[3*neigh+1], coords[3*neigh+2] );
        }
        Vector  normal;
        find_polygon_normal( n_ngh[point], neigh_points, &normal );

        Real nx = normal.coords[0];
        Real ny = normal.coords[1];
        Real nz = normal.coords[2];

        if( image_type == T1 ) {
          if( nx*dx + ny*dy + nz*dz > 0.0 ) {
            dx = 0.0;
            dy = 0.0;
            dz = 0.0;
          }
        } else if( image_type == HISTO ) {
          if( nx*dx + ny*dy + nz*dz < 0.0 ) {
            dx = 0.0;
            dy = 0.0;
            dz = 0.0;
          }
        }

        Real factor = weight * ( value / t1grad[point-start_point] - 1.0 ) / 
                      t1grad[point-start_point];
        deriv[3*point+0] += factor * dx;
        deriv[3*point+1] += factor * dy;
        deriv[3*point+2] += factor * dz;
      }

    }

    if( oversample ) {

      for( point = start_point; point < end_point; point++ ) {
        if( mask[point-start_point] ) continue;

        for_less( n, 0, n_ngh[point] ) {
          int neigh1 = ngh[point][n];
          if( point > neigh1 ) continue;

          Real xmid = 0.5 * ( coords[3*point] + coords[3*neigh1] );
          Real ymid = 0.5 * ( coords[3*point+1] + coords[3*neigh1+1] );
          Real zmid = 0.5 * ( coords[3*point+2] + coords[3*neigh1+2] );
          Real t1mid = 0.5 * ( t1grad[point-start_point] + t1grad[neigh1-start_point] );

          Real value, dx, dy, dz;
          evaluate_volume_in_world( t1, xmid, ymid, zmid,
                                    EVAL_LINEAR, FALSE,
                                    0.0, &value, &dx, &dy, &dz,
                                    NULL, NULL, NULL, NULL, NULL, NULL );

          // This is an approximate vector based on the triangular face.
          // Should consider triangles on each side of edge, not only one side.
          int neigh2 = ngh[point][(n+1)%n_ngh[point]];
          fill_Point( neigh_points[0], coords[3*point+0],
                      coords[3*point+1], coords[3*point+2] );
          fill_Point( neigh_points[1], coords[3*neigh1+0],
                      coords[3*neigh1+1], coords[3*neigh1+2] );
          fill_Point( neigh_points[2], coords[3*neigh2+0],
                      coords[3*neigh2+1], coords[3*neigh2+2] );

          Vector  normal;
          find_polygon_normal( 3, neigh_points, &normal );

          Real nx = normal.coords[0];
          Real ny = normal.coords[1];
          Real nz = normal.coords[2];

          if( image_type == T1 ) {
            if( nx*dx + ny*dy + nz*dz > 0.0 ) {
              dx = 0.0;
              dy = 0.0;
              dz = 0.0;
            }
          } else if( image_type == HISTO ) {
            if( nx*dx + ny*dy + nz*dz < 0.0 ) {
              dx = 0.0;
              dy = 0.0;
              dz = 0.0;
            }
          }

          // half to each vertex of the edge
          Real factor = 0.5 * weight * ( value / t1mid - 1.0 ) / t1mid;
          deriv[3*point+0] += factor * dx;
          deriv[3*point+1] += factor * dy;
          deriv[3*point+2] += factor * dz;
          deriv[3*neigh1+0] += factor * dx;
          deriv[3*neigh1+1] += factor * dy;
          deriv[3*neigh1+2] += factor * dz;

        }
      }
    }
    free( neigh_points );
  }
}


private  void   evaluate_image_value_fit_deriv(
    Volume               volume,
    voxel_coef_struct    *lookup,
    int                  continuity,
    Real                 image_weight,
    Real                 threshold,
    Real                 min_diff,
    Real                 max_diff,
    Real                 max_diff_weight,
    Real                 differential_offset,
    Real                 differential_ratio,
    int                  oversample,
    int                  n_points,
    int                  n_neighbours[],
    int                  *neighbours[],
    Real                 parameters[],
    Real                 deriv[] )
{
    int     p1, p1_index, n1_index, n2_index, neigh1, neigh2, w1, w2, n;
    Real    x, y, z, value, diff, dx, dy, dz;
    Real    x1, y1, z1, x2, y2, z2, x3, y3, z3, derivative;
    Real    weight1, weight2, weight3, alpha1, alpha2, f1, f2, f3;
    Real    ln_start_weight, value_deriv[N_DIMENSIONS], *value_deriv_ptr[1];
    Real    wtv00, wtv01, wtv02, wtv03;
    Real    wtv10, wtv11, wtv12, wtv13;
    Real    wtv20, wtv21, wtv22, wtv23;
    Real    voxel[3], voxel_000[3], voxel_100[3], voxel_010[3];
    Real    voxel_001[3];
    Real              trilin_sizes[N_DIMENSIONS];
    int               n_weights, n_weights1, n_weights2, sizes[N_DIMENSIONS];
    weighting_struct  tt, *t_ptr, *weights1, *weights2, *weights;
    Real              dx12, dy12, dz12, dx13, dy13, dz13;
    Real              xv1, yv1, zv1, xv2, yv2, zv2, xv3, yv3, zv3;
    BOOLEAN           inside;

    if( image_weight < 0.0 )
        return;

    get_volume_sizes( volume, sizes );
    trilin_sizes[0] = (Real) sizes[0] - 1.0;
    trilin_sizes[1] = (Real) sizes[1] - 1.0;
    trilin_sizes[2] = (Real) sizes[2] - 1.0;

    convert_world_to_voxel( volume, 0.0, 0.0, 0.0, voxel_000 );
    convert_world_to_voxel( volume, 1.0, 0.0, 0.0, voxel_100 );
    convert_world_to_voxel( volume, 0.0, 1.0, 0.0, voxel_010 );
    convert_world_to_voxel( volume, 0.0, 0.0, 1.0, voxel_001 );

    wtv00 = voxel_100[0] - voxel_000[0];
    wtv01 = voxel_010[0] - voxel_000[0];
    wtv02 = voxel_001[0] - voxel_000[0];
    wtv03 = voxel_000[0];

    wtv10 = voxel_100[1] - voxel_000[1];
    wtv11 = voxel_010[1] - voxel_000[1];
    wtv12 = voxel_001[1] - voxel_000[1];
    wtv13 = voxel_000[1];

    wtv20 = voxel_100[2] - voxel_000[2];
    wtv21 = voxel_010[2] - voxel_000[2];
    wtv22 = voxel_001[2] - voxel_000[2];
    wtv23 = voxel_000[2];

    ln_start_weight = log( differential_ratio );

    inside = TRUE;
    for_less( p1, 0, n_points )
    {
        p1_index = IJ(p1,0,3);
        x1 = parameters[p1_index+0];
        y1 = parameters[p1_index+1];
        z1 = parameters[p1_index+2];

        xv1 = wtv00 * x1 + wtv01 * y1 + wtv02 * z1 + wtv03;
        yv1 = wtv10 * x1 + wtv11 * y1 + wtv12 * z1 + wtv13;
        zv1 = wtv20 * x1 + wtv21 * y1 + wtv22 * z1 + wtv23;

        if( xv1 < 0.0 || xv1 >= (Real) sizes[0]-1.0 ||
            yv1 < 0.0 || yv1 >= (Real) sizes[1]-1.0 ||
            zv1 < 0.0 || zv1 >= (Real) sizes[2]-1.0 )
        {
            inside = FALSE;

            break;
        }
    }

    for_less( p1, 0, n_points )
    {
        p1_index = IJ( p1, 0, 3 );

        x = parameters[p1_index+0];
        y = parameters[p1_index+1];
        z = parameters[p1_index+2];

        voxel[0] = wtv00 * x + wtv01 * y + wtv02 * z + wtv03;
        voxel[1] = wtv10 * x + wtv11 * y + wtv12 * z + wtv13;
        voxel[2] = wtv20 * x + wtv21 * y + wtv22 * z + wtv23;

        if( continuity == 0 )
        {
            value = trilinear_interpolate_with_deriv( volume, inside,
                                                      trilin_sizes,
                                                      lookup, voxel,
                                                      value_deriv );
        }
        else
        {
            value_deriv_ptr[0] = value_deriv;
            (void) evaluate_volume( volume, voxel, NULL, continuity,
                                    FALSE, 0.0, &value, value_deriv_ptr, NULL );
        }

        convert_voxel_normal_vector_to_world( volume, value_deriv,
                                              &dx, &dy, &dz );

        diff = value - threshold;

        deriv[p1_index+0] += image_weight * 2.0 * dx * diff;
        deriv[p1_index+1] += image_weight * 2.0 * dy * diff;
        deriv[p1_index+2] += image_weight * 2.0 * dz * diff;

        if( min_diff < max_diff && max_diff_weight > 0.0 &&
            (diff < min_diff || diff > max_diff) )
        {
#ifdef OLD
            if( diff < min_diff )
                diff = diff - min_diff;
            else
                diff = diff - max_diff;
#endif

            derivative = differential_weights_deriv( diff, min_diff, max_diff,
                                                     differential_offset,
                                                     ln_start_weight );

            deriv[p1_index+0] += max_diff_weight * derivative * dx;
            deriv[p1_index+1] += max_diff_weight * derivative * dy;
            deriv[p1_index+2] += max_diff_weight * derivative * dz;
        }
    }

    if( oversample <= 0 )
        return;

    n_weights1 = 0;
    n_weights2 = 0;
    weights1 = NULL;
    weights2 = NULL;

    for_inclusive( w1, 0, oversample )
    {
        tt.alpha1 = (Real) w1 / (Real) (oversample+1);
        for_inclusive( w2, 1, oversample - w1 + 1 )
        {
            tt.alpha2 = (1.0 - tt.alpha1) *
                        (Real) w2 / (Real) (oversample-w1+1);


            if( w2 < oversample + 1 )
            {
                ADD_ELEMENT_TO_ARRAY( weights1, n_weights1, tt,
                                      DEFAULT_CHUNK_SIZE );
            }

            if( w2 < oversample - w1 + 1 )
            {
                ADD_ELEMENT_TO_ARRAY( weights2, n_weights2, tt,
                                      DEFAULT_CHUNK_SIZE );
            }
        }
    }

    for_less( p1, 0, n_points )
    {
        for_less( n, 0, n_neighbours[p1] )
        {
            neigh1 = neighbours[p1][n];
            neigh2 = neighbours[p1][(n+1)%n_neighbours[p1]];

            if( p1 > neigh1 || p1 > neigh2 )
                continue;

            p1_index = IJ(p1,0,3);
            n1_index = IJ(neigh1,0,3);
            n2_index = IJ(neigh2,0,3);

            x1 = parameters[p1_index+0];
            y1 = parameters[p1_index+1];
            z1 = parameters[p1_index+2];
            x2 = parameters[n1_index+0];
            y2 = parameters[n1_index+1];
            z2 = parameters[n1_index+2];
            x3 = parameters[n2_index+0];
            y3 = parameters[n2_index+1];
            z3 = parameters[n2_index+2];

            xv1 = wtv00 * x1 + wtv01 * y1 + wtv02 * z1 + wtv03;
            yv1 = wtv10 * x1 + wtv11 * y1 + wtv12 * z1 + wtv13;
            zv1 = wtv20 * x1 + wtv21 * y1 + wtv22 * z1 + wtv23;

            xv2 = wtv00 * x2 + wtv01 * y2 + wtv02 * z2 + wtv03;
            yv2 = wtv10 * x2 + wtv11 * y2 + wtv12 * z2 + wtv13;
            zv2 = wtv20 * x2 + wtv21 * y2 + wtv22 * z2 + wtv23;

            xv3 = wtv00 * x3 + wtv01 * y3 + wtv02 * z3 + wtv03;
            yv3 = wtv10 * x3 + wtv11 * y3 + wtv12 * z3 + wtv13;
            zv3 = wtv20 * x3 + wtv21 * y3 + wtv22 * z3 + wtv23;

            dx12 = xv2 - xv1;
            dy12 = yv2 - yv1;
            dz12 = zv2 - zv1;

            dx13 = xv3 - xv1;
            dy13 = yv3 - yv1;
            dz13 = zv3 - zv1;

            if( neigh1 < neigh2 )
            {
                n_weights = n_weights1;
                weights = weights1;
            }
            else
            {
                n_weights = n_weights2;
                weights = weights2;
            }

            for_less( w1, 0, n_weights )
            {
                t_ptr = &weights[w1];
                alpha1 = t_ptr->alpha1;
                alpha2 = t_ptr->alpha2;

                voxel[0] = xv1 + alpha1 * dx12 + alpha2 * dx13;
                voxel[1] = yv1 + alpha1 * dy12 + alpha2 * dy13;
                voxel[2] = zv1 + alpha1 * dz12 + alpha2 * dz13;

                if( continuity == 0 )
                {
                    value = trilinear_interpolate_with_deriv( volume, inside,
                                   trilin_sizes, lookup, voxel, value_deriv );
                }
                else
                {
                    value_deriv_ptr[0] = value_deriv;
                    (void) evaluate_volume( volume, voxel, NULL, continuity,
                                    FALSE, 0.0, &value, value_deriv_ptr, NULL );
                }

                convert_voxel_normal_vector_to_world( volume, value_deriv,
                                                      &dx, &dy, &dz );

                diff = value - threshold;

                weight1 = (1.0 - alpha1 - alpha2);
                weight2 = alpha1;
                weight3 = alpha2;

                f1 = image_weight * weight1 * 2.0 * diff;
                f2 = image_weight * weight2 * 2.0 * diff;
                f3 = image_weight * weight3 * 2.0 * diff;

                deriv[p1_index+0] += f1 * dx;
                deriv[p1_index+1] += f1 * dy;
                deriv[p1_index+2] += f1 * dz;

                deriv[n1_index+0] += f2 * dx;
                deriv[n1_index+1] += f2 * dy;
                deriv[n1_index+2] += f2 * dz;

                deriv[n2_index+0] += f3 * dx;
                deriv[n2_index+1] += f3 * dy;
                deriv[n2_index+2] += f3 * dz;

                if( min_diff >= max_diff || max_diff_weight <= 0.0 ||
                    diff >= min_diff && diff <= max_diff )
                    continue;

                derivative = differential_weights_deriv( diff, min_diff,
                                      max_diff, differential_offset,
                                      ln_start_weight );

                f1 = max_diff_weight * weight1 * derivative;
                f2 = max_diff_weight * weight2 * derivative;
                f3 = max_diff_weight * weight3 * derivative;

                deriv[p1_index+0] += f1 * dx;
                deriv[p1_index+1] += f1 * dy;
                deriv[p1_index+2] += f1 * dz;

                deriv[n1_index+0] += f2 * dx;
                deriv[n1_index+1] += f2 * dy;
                deriv[n1_index+2] += f2 * dz;

                deriv[n2_index+0] += f3 * dx;
                deriv[n2_index+1] += f3 * dy;
                deriv[n2_index+2] += f3 * dz;
            }
        }
    }

    FREE( weights1 );
    FREE( weights2 );
}

private  Real   evaluate_stretch_fit(
    Real          stretch_weight,
    int           n_points,
    Real          parameters[],
    int           n_neighbours[],
    int           *neighbours[] )
{
    int     point, neigh, ind, n, ind0, ind1, ind2, ind3, *neigh_ptr, n_neighs;
    Real    fit1, fit2, dx, dy, dz, len, x1, y1, z1, elen[3];

    if( stretch_weight < 0.0 )
        return( 0.0 );

    fit1 = 0.0;
    fit2 = 0.0;

    Real * areas;
    ALLOC( areas, n_points );
    for_less( point, 0, n_points ) areas[point] = 0.0;

    for_less( point, 0, n_points ) {
      ind0 = IJ(point,0,3);
      for_less( n, 0, n_neighbours[point] ) {
        int n1 = neighbours[point][n];
        int n2 = neighbours[point][(n+1)%n_neighbours[point]];
        if( !( point < n1 && point < n2 ) ) continue;
        ind1 = IJ(n1,0,3);
        ind2 = IJ(n2,0,3);
        elen[0] = sqrt( ( parameters[ind0+0] - parameters[ind1+0] ) *
                        ( parameters[ind0+0] - parameters[ind1+0] ) +
                        ( parameters[ind0+1] - parameters[ind1+1] ) *
                        ( parameters[ind0+1] - parameters[ind1+1] ) +
                        ( parameters[ind0+2] - parameters[ind1+2] ) *
                        ( parameters[ind0+2] - parameters[ind1+2] ) );
        elen[1] = sqrt( ( parameters[ind2+0] - parameters[ind1+0] ) *
                        ( parameters[ind2+0] - parameters[ind1+0] ) +
                        ( parameters[ind2+1] - parameters[ind1+1] ) *
                        ( parameters[ind2+1] - parameters[ind1+1] ) +
                        ( parameters[ind2+2] - parameters[ind1+2] ) *
                        ( parameters[ind2+2] - parameters[ind1+2] ) );
        elen[2] = sqrt( ( parameters[ind0+0] - parameters[ind2+0] ) *
                        ( parameters[ind0+0] - parameters[ind2+0] ) +
                        ( parameters[ind0+1] - parameters[ind2+1] ) *
                        ( parameters[ind0+1] - parameters[ind2+1] ) +
                        ( parameters[ind0+2] - parameters[ind2+2] ) *
                        ( parameters[ind0+2] - parameters[ind2+2] ) );
        Real s = 0.5 * ( elen[0] + elen[1] + elen[2] );
        Real local_area = sqrt( fabs( s * ( s - elen[0] ) * ( s - elen[1] ) * ( s - elen[2] ) ) +
                                1.0e-20 );
        areas[point] += local_area;
        areas[n1] += local_area;
        areas[n2] += local_area;
      }
    }

    // Starting edge
    ind = 0;

    for_less( point, 0, n_points ) {
        ind0 = IJ(point,0,3);
        x1 = parameters[ind0+0];
        y1 = parameters[ind0+1];
        z1 = parameters[ind0+2];

        neigh_ptr = neighbours[point];
        n_neighs = n_neighbours[point];

        Real sum_areas = 0.0;
        dx = 0.0; dy = 0.0; dz = 0.0;
        for_less( n, 0, n_neighs ) {
            neigh = neigh_ptr[n];
            ind1 = IJ(neigh,0,3);
            dx += areas[neigh] * parameters[ind1+0];
            dy += areas[neigh] * parameters[ind1+1];
            dz += areas[neigh] * parameters[ind1+2];
            sum_areas += areas[neigh];
        }

        dx = dx / sum_areas - x1;
        dy = dy / sum_areas - y1;
        dz = dz / sum_areas - z1;
        fit1 += dx * dx + dy * dy + dz * dz;
    }
    FREE( areas );

    return( 0.5 * fit1 * stretch_weight );
}

private  void   evaluate_stretch_fit_deriv(
    Real          stretch_weight,
    int           n_points,
    Real          parameters[],
    int           n_neighbours[],
    int           *neighbours[],
    Real          deriv[] )
{
    int    point, neigh, ind, n, ind0, ind1, ind2, ind3, j;
    Real   dx, dy, dz;

    if( stretch_weight < 0.0 )
        return;

    // Starting edge

    ind = 0;

    Real   elen[3];
    Real * areas;
    ALLOC( areas, n_points );
    for_less( point, 0, n_points ) areas[point] = 0.0;

    for_less( point, 0, n_points ) {
      ind0 = IJ(point,0,3);
      for_less( n, 0, n_neighbours[point] ) {
        int n1 = neighbours[point][n];
        int n2 = neighbours[point][(n+1)%n_neighbours[point]];
        if( !( point < n1 && point < n2 ) ) continue;
        ind1 = IJ(n1,0,3);
        ind2 = IJ(n2,0,3);
        elen[0] = sqrt( ( parameters[ind0+0] - parameters[ind1+0] ) *
                        ( parameters[ind0+0] - parameters[ind1+0] ) +
                        ( parameters[ind0+1] - parameters[ind1+1] ) *
                        ( parameters[ind0+1] - parameters[ind1+1] ) +
                        ( parameters[ind0+2] - parameters[ind1+2] ) *
                        ( parameters[ind0+2] - parameters[ind1+2] ) );
        elen[1] = sqrt( ( parameters[ind2+0] - parameters[ind1+0] ) *
                        ( parameters[ind2+0] - parameters[ind1+0] ) +
                        ( parameters[ind2+1] - parameters[ind1+1] ) *
                        ( parameters[ind2+1] - parameters[ind1+1] ) +
                        ( parameters[ind2+2] - parameters[ind1+2] ) *
                        ( parameters[ind2+2] - parameters[ind1+2] ) );
        elen[2] = sqrt( ( parameters[ind0+0] - parameters[ind2+0] ) *
                        ( parameters[ind0+0] - parameters[ind2+0] ) +
                        ( parameters[ind0+1] - parameters[ind2+1] ) *
                        ( parameters[ind0+1] - parameters[ind2+1] ) +
                        ( parameters[ind0+2] - parameters[ind2+2] ) *
                        ( parameters[ind0+2] - parameters[ind2+2] ) );
        Real s = 0.5 * ( elen[0] + elen[1] + elen[2] );
        Real local_area = sqrt( fabs( s * ( s - elen[0] ) * ( s - elen[1] ) * ( s - elen[2] ) ) +
                                1.0e-20 );
        areas[point] += local_area;
        areas[n1] += local_area;
        areas[n2] += local_area;
      }
    }

    Real * delta, * norm_areas;
    ALLOC( delta, 3 * n_points ); 
    ALLOC( norm_areas, n_points ); 

    for_less( point, 0, n_points ) {
      ind0 = IJ(point,0,3);
      for( j = 0; j < 3; j++ ) delta[ind0+j] = 0.0;
      Real sum_areas = 0.0;
      for_less( n, 0, n_neighbours[point] ) {
        neigh = neighbours[point][n];
        ind1 = IJ(neigh,0,3);
        for( j = 0; j < 3; j++ ) {
          delta[ind0+j] += areas[neigh] * parameters[ind1+j];
        }
        sum_areas += areas[neigh];
      }
      norm_areas[point] = areas[point] / sum_areas;
      for( j = 0; j < 3; j++ ) {
        delta[ind0+j] /= sum_areas;
        delta[ind0+j] -= parameters[ind0+j];
      }
    }
    FREE( areas );
    for_less( point, 0, n_points ) {
        ind0 = IJ(point,0,3);

        deriv[ind0+0] -= stretch_weight * delta[ind0+0];
        deriv[ind0+1] -= stretch_weight * delta[ind0+1];
        deriv[ind0+2] -= stretch_weight * delta[ind0+2];
        for_less( n, 0, n_neighbours[point] ) {
            neigh = neighbours[point][n];
            ind1 = IJ(neigh,0,3);
            deriv[ind0+0] += stretch_weight * delta[ind1+0] * norm_areas[neigh];
            deriv[ind0+1] += stretch_weight * delta[ind1+1] * norm_areas[neigh];
            deriv[ind0+2] += stretch_weight * delta[ind1+2] * norm_areas[neigh];
        }
    }

    FREE( delta );
    FREE( norm_areas );
}

public  Real  get_base_length(
    int    n_neighbours,
    Point  neighbours[],
    Point  *centroid )
{
    int    p;
    Real   len, x, y, z, dx, dy, dz;

    x = RPoint_x(*centroid);
    y = RPoint_y(*centroid);
    z = RPoint_z(*centroid);
    len = 0.0;
    for_less( p, 0, n_neighbours )
    {
        dx = RPoint_x(neighbours[p]) - x;
        dy = RPoint_y(neighbours[p]) - y;
        dz = RPoint_z(neighbours[p]) - z;
        len += dx * dx + dy * dy + dz * dz;
    }

    return( sqrt( len / (Real) n_neighbours ) );
}

public  void  get_neighbours_normal(
    int      n_points,
    Point    points[],
    Vector   *normal )
{
    int     i, next_i;
    Vector  v1, v2;
    Real    vx, vy, vz, x, y, z, tx, ty, tz;

    vx = 0.0;
    vy = 0.0;
    vz = 0.0;

    tx = RPoint_x(points[0]);
    ty = RPoint_y(points[0]);
    tz = RPoint_z(points[0]);

    for_less( i, 0, n_points )
    {
        next_i = (i + 1) % n_points;

        x = tx;
        y = ty;
        z = tz;
        tx = RPoint_x(points[next_i]);
        ty = RPoint_y(points[next_i]);
        tz = RPoint_z(points[next_i]);

        vx += y * tz - ty * z;
        vy += z * tx - tz * x;
        vz += x * ty - tx * y;
    }

    /*--- if result is null, try to find one vertex for which a normal can
          be computed */

    if( vx == 0.0 && vy == 0.0 && vz == 0.0 )
    {
        for_less( i, 0, n_points )
        {
            SUB_POINTS( v1, points[(i+1)%n_points], points[i] );
            SUB_POINTS( v2, points[(i-1)%n_points], points[i] );
            CROSS_VECTORS( *normal, v1, v2 );
            if( !null_Vector( normal ) )
                break;
        }
    }
    else
        fill_Point( *normal, vx, vy, vz );
}

public  void  compute_centroid(
    int     n_points,
    Point   points[],
    Point   *centroid )
{
    int    i;
    Real   x, y, z, factor;

    x = 0.0;
    y = 0.0;
    z = 0.0;

    for_less( i, 0, n_points )
    {
        x += RPoint_x(points[i]);
        y += RPoint_y(points[i]);
        z += RPoint_z(points[i]);
    }

    factor = 1.0 / (Real) n_points;
    fill_Point( *centroid, x * factor, y * factor, z * factor );
}

private  Real   evaluate_curvature_fit(
    Real          curvature_weight,
    Real          max_curvature_weight,
    Real          min_curvature,
    Real          max_curvature,
    int           start_point,
    int           end_point,
    Real          parameters[],
    Smallest_int  active_flags[],
    int           n_neighbours[],
    int           *neighbours[],
    float         curvatures[] )
{
    int     point, neigh, n, ind0, ind1, max_neighbours;
    int     *neigh_ptr, n_neighs;
    Real    fit1, fit2, normal_len;
    Real    curv;
    Point   *neigh_points;
    Real    x, y, z, nx, ny, nz, factor, cx, cy, cz, tx, ty, tz;
    Real    x_this, y_this, z_this, ss;

    if( curvature_weight < 0.0 )
        return( 0.0 );

    fit1 = 0.0;
    fit2 = 0.0;
    max_neighbours = 0;
    for_less( point, start_point, end_point )
        max_neighbours = MAX( max_neighbours, n_neighbours[point] );

    ALLOC( neigh_points, max_neighbours );

    for_less( point, start_point, end_point )
    {
        neigh_ptr = neighbours[point];
        n_neighs = n_neighbours[point];

        if( active_flags != NULL && !active_flags[point] )
            continue;

        ind0 = IJ(point,0,3);
        x_this = parameters[ind0+0];
        y_this = parameters[ind0+1];
        z_this = parameters[ind0+2];

        if( n_neighs == 6 )
        {
            ss = 0.0;
            cx = 0.0;
            cy = 0.0;
            cz = 0.0;
            nx = 0.0;
            ny = 0.0;
            nz = 0.0;
            neigh = neigh_ptr[5];
            ind1 = IJ(neigh,0,3);
            tx = parameters[ind1+0];
            ty = parameters[ind1+1];
            tz = parameters[ind1+2];

            x = tx;
            y = ty;
            z = tz;
            neigh = neigh_ptr[0];
            ind1 = IJ(neigh,0,3);
            tx = parameters[ind1+0];
            ty = parameters[ind1+1];
            tz = parameters[ind1+2];
            cx += tx;
            cy += ty;
            cz += tz;
            ss += tx * tx + ty * ty + tz * tz;
            nx += y * tz - ty * z;
            ny += z * tx - tz * x;
            nz += x * ty - tx * y;

            x = tx;
            y = ty;
            z = tz;
            neigh = neigh_ptr[1];
            ind1 = IJ(neigh,0,3);
            tx = parameters[ind1+0];
            ty = parameters[ind1+1];
            tz = parameters[ind1+2];
            cx += tx;
            cy += ty;
            cz += tz;
            ss += tx * tx + ty * ty + tz * tz;
            nx += y * tz - ty * z;
            ny += z * tx - tz * x;
            nz += x * ty - tx * y;

            x = tx;
            y = ty;
            z = tz;
            neigh = neigh_ptr[2];
            ind1 = IJ(neigh,0,3);
            tx = parameters[ind1+0];
            ty = parameters[ind1+1];
            tz = parameters[ind1+2];
            cx += tx;
            cy += ty;
            cz += tz;
            ss += tx * tx + ty * ty + tz * tz;
            nx += y * tz - ty * z;
            ny += z * tx - tz * x;
            nz += x * ty - tx * y;

            x = tx;
            y = ty;
            z = tz;
            neigh = neigh_ptr[3];
            ind1 = IJ(neigh,0,3);
            tx = parameters[ind1+0];
            ty = parameters[ind1+1];
            tz = parameters[ind1+2];
            cx += tx;
            cy += ty;
            cz += tz;
            ss += tx * tx + ty * ty + tz * tz;
            nx += y * tz - ty * z;
            ny += z * tx - tz * x;
            nz += x * ty - tx * y;

            x = tx;
            y = ty;
            z = tz;
            neigh = neigh_ptr[4];
            ind1 = IJ(neigh,0,3);
            tx = parameters[ind1+0];
            ty = parameters[ind1+1];
            tz = parameters[ind1+2];
            cx += tx;
            cy += ty;
            cz += tz;
            ss += tx * tx + ty * ty + tz * tz;
            nx += y * tz - ty * z;
            ny += z * tx - tz * x;
            nz += x * ty - tx * y;

            x = tx;
            y = ty;
            z = tz;
            neigh = neigh_ptr[5];
            ind1 = IJ(neigh,0,3);
            tx = parameters[ind1+0];
            ty = parameters[ind1+1];
            tz = parameters[ind1+2];
            cx += tx;
            cy += ty;
            cz += tz;
            ss += tx * tx + ty * ty + tz * tz;
            nx += y * tz - ty * z;
            ny += z * tx - tz * x;
            nz += x * ty - tx * y;

            factor = 1.0 / 6.0;

            cx *= factor;
            cy *= factor;
            cz *= factor;
        }
        else
        {
            ss = 0.0;
            cx = 0.0;
            cy = 0.0;
            cz = 0.0;
            nx = 0.0;
            ny = 0.0;
            nz = 0.0;
            neigh = neigh_ptr[n_neighs-1];
            ind1 = IJ(neigh,0,3);
            tx = parameters[ind1+0];
            ty = parameters[ind1+1];
            tz = parameters[ind1+2];

            for_less( n, 0, n_neighs )
            {
                x = tx;
                y = ty;
                z = tz;

                neigh = neigh_ptr[n];
                ind1 = IJ(neigh,0,3);

                tx = parameters[ind1+0];
                ty = parameters[ind1+1];
                tz = parameters[ind1+2];
                cx += tx;
                cy += ty;
                cz += tz;

                ss += tx * tx + ty * ty + tz * tz;

                nx += y * tz - ty * z;
                ny += z * tx - tz * x;
                nz += x * ty - tx * y;
            }

            factor = 1.0 / (Real) n_neighs;

            cx *= factor;
            cy *= factor;
            cz *= factor;
        }

#ifdef USE_CURV_SQRT
        normal_len = sqrt( (nx*nx+ny*ny+nz*nz) *
                           (ss * factor - cx*cx - cy*cy - cz*cz) );
#else
        normal_len = nx*nx+ny*ny+nz*nz;
#endif

        curv = ((x_this-cx)*nx + (y_this-cy)*ny + (z_this-cz)*nz) /
               normal_len - (Real) curvatures[point];

        fit1 += curv * curv;

        if( min_curvature >= max_curvature || max_curvature_weight <= 0.0 )
            continue;

        if( curv < min_curvature )
        {
            curv = curv - min_curvature;
            fit2 += curv * curv;
        }
        else if( curv > max_curvature )
        {
            curv = curv - max_curvature;
            fit2 += curv * curv;
        }
    }

    FREE( neigh_points );

    return( fit1 * curvature_weight + fit2 * max_curvature_weight );
}

private  void   evaluate_curvature_fit_deriv(
    Real          curvature_weight,
    Real          max_curvature_weight,
    Real          min_curvature,
    Real          max_curvature,
    int           n_points,
    Real          parameters[],
    int           n_neighbours[],
    int           *neighbours[],
    float         curvatures[],
    Real          deriv[] )
{
    int     point, neigh, n, ind0, ind1, nn, dim, a1, a2;
    Real    factor, curv, normal_len;
    Real    x, y, z, nx, ny, nz, cx, cy, cz, tx, ty, tz;
    Real    x_this, y_this, z_this, ss, d_ss;
    Real    d_left, d_right, d_curv, d_nxnynz, d_sscxcycz;
    Real    d_centroid[N_DIMENSIONS], d_normal[N_DIMENSIONS];

    if( curvature_weight < 0.0 )
        return;

    for_less( point, 0, n_points )
    {
        ind0 = IJ(point,0,3);

        x_this = parameters[ind0+0];
        y_this = parameters[ind0+1];
        z_this = parameters[ind0+2];

        ss = 0.0;
        cx = 0.0;
        cy = 0.0;
        cz = 0.0;
        nx = 0.0;
        ny = 0.0;
        nz = 0.0;
        nn = n_neighbours[point];
        neigh = neighbours[point][nn-1];
        ind1 = IJ(neigh,0,3);
        tx = parameters[ind1+0];
        ty = parameters[ind1+1];
        tz = parameters[ind1+2];

        for_less( n, 0, nn )
        {
            x = tx;
            y = ty;
            z = tz;

            neigh = neighbours[point][n];
            ind1 = IJ(neigh,0,3);

            tx = parameters[ind1+0];
            ty = parameters[ind1+1];
            tz = parameters[ind1+2];
            cx += tx;
            cy += ty;
            cz += tz;

            ss += tx * tx + ty * ty + tz * tz;

            nx += y * tz - ty * z;
            ny += z * tx - tz * x;
            nz += x * ty - tx * y;
        }

        factor = 1.0 / (Real) nn;

        cx *= factor;
        cy *= factor;
        cz *= factor;

#ifdef USE_CURV_SQRT
        normal_len = sqrt( (nx*nx+ny*ny+nz*nz) *
                           (ss * factor - cx*cx - cy*cy - cz*cz) );
#else
        normal_len = nx*nx+ny*ny+nz*nz;
#endif

        curv = ((x_this-cx)*nx + (y_this-cy)*ny + (z_this-cz)*nz) /
               normal_len - (Real) curvatures[point];

        factor = curvature_weight * 2.0 * curv;

        deriv[ind0+0] += factor * nx / normal_len;
        deriv[ind0+1] += factor * ny / normal_len;
        deriv[ind0+2] += factor * nz / normal_len;

        /*--- now for other derivatives */

        for_less( n, 0, nn )
        {
            neigh = neighbours[point][n];
            ind1 = IJ(neigh,0,3);
            for_less( dim, 0, 3 )
            {
                a1 = (dim+1) % N_DIMENSIONS;
                a2 = (dim+2) % N_DIMENSIONS;

                d_centroid[dim] = 1.0 / (Real) nn;
                d_centroid[a1] = 0.0;
                d_centroid[a2] = 0.0;

                d_normal[dim] = 0.0;
                d_normal[a1] =
                          parameters[IJ(neighbours[point][(n-1+nn)%nn],a2,3)] -
                          parameters[IJ(neighbours[point][(n+1)%nn],a2,3)];
                d_normal[a2] =
                          parameters[IJ(neighbours[point][(n+1)%nn],a1,3)] -
                          parameters[IJ(neighbours[point][(n-1+nn)%nn],a1,3)];

                d_left = -d_centroid[X] * nx + (x_this - cx) * d_normal[X]
                         -d_centroid[Y] * ny + (y_this - cy) * d_normal[Y]
                         -d_centroid[Z] * nz + (z_this - cz) * d_normal[Z];

                d_nxnynz = 2.0 * nx * d_normal[X] +
                           2.0 * ny * d_normal[Y] +
                           2.0 * nz * d_normal[Z];

                d_ss = 2.0 * parameters[ind1+dim] / (Real) nn;

                d_sscxcycz = d_ss - 2.0 * cx * d_centroid[X] -
                                    2.0 * cy * d_centroid[Y] -
                                    2.0 * cz * d_centroid[Z];

                d_right = -0.5 / normal_len / normal_len / normal_len *
                          (d_nxnynz * (ss / (Real) nn - cx*cx - cy*cy - cz*cz) +
                           (nx*nx+ny*ny+nz*nz) * d_sscxcycz);

                d_curv = d_left / normal_len +
                         ((x_this-cx)*nx + (y_this-cy)*ny + (z_this-cz)*nz) *
                         d_right;

                deriv[ind1+dim] += factor * d_curv;
            }
        }

        if( min_curvature >= max_curvature || max_curvature_weight <= 0.0 ||
            (min_curvature <= curv && curv <= max_curvature) )
            continue;

        if( curv < min_curvature )
            curv = curv - min_curvature;
        else
            curv = curv + min_curvature;

        factor = max_curvature_weight * 2.0 * curv;

        deriv[ind0+0] += factor * nx / normal_len;
        deriv[ind0+1] += factor * ny / normal_len;
        deriv[ind0+2] += factor * nz / normal_len;

        /*--- now for other derivatives */

        for_less( n, 0, nn )
        {
            neigh = neighbours[point][n];
            ind1 = IJ(neigh,0,3);
            for_less( dim, 0, 3 )
            {
                a1 = (dim+1) % N_DIMENSIONS;
                a2 = (dim+2) % N_DIMENSIONS;

                d_centroid[dim] = 1.0 / (Real) nn;
                d_centroid[a1] = 0.0;
                d_centroid[a2] = 0.0;

                d_normal[dim] = 0.0;
                d_normal[a1] =
                          parameters[IJ(neighbours[point][(n-1+nn)%nn],a2,3)] -
                          parameters[IJ(neighbours[point][(n+1)%nn],a2,3)];
                d_normal[a2] =
                          parameters[IJ(neighbours[point][(n+1)%nn],a1,3)] -
                          parameters[IJ(neighbours[point][(n-1+nn)%nn],a1,3)];

                d_left = -d_centroid[X] * nx + (x_this - cx) * d_normal[X]
                         -d_centroid[Y] * ny + (y_this - cy) * d_normal[Y]
                         -d_centroid[Z] * nz + (z_this - cz) * d_normal[Z];

                d_nxnynz = 2.0 * nx * d_normal[X] +
                           2.0 * ny * d_normal[Y] +
                           2.0 * nz * d_normal[Z];

                d_ss = 2.0 * parameters[ind1+dim] / (Real) nn;

                d_sscxcycz = d_ss - 2.0 * cx * d_centroid[X] -
                                    2.0 * cy * d_centroid[Y] -
                                    2.0 * cz * d_centroid[Z];

                d_right = -0.5 / normal_len / normal_len / normal_len *
                          (d_nxnynz * (ss / (Real) nn - cx*cx - cy*cy - cz*cz) +
                           (nx*nx+ny*ny+nz*nz) * d_sscxcycz);

                d_curv = d_left / normal_len +
                         ((x_this-cx)*nx + (y_this-cy)*ny + (z_this-cz)*nz) *
                         d_right;

                deriv[ind1+dim] += factor * d_curv;
            }
        }
    }
}

#ifdef DEBUG
private  void  test_get(
    Real   p0[],
    Real   p1[],
    Real   p2[],
    Real   p3[],
    Real   *x,
    Real   *y )
{
    Real   ax, ay, az, bx, by, bz, cx, cy, cz;
    Real   abx, aby, abz, len_ab2, len, len2, xn, yn, hx, hy, hz, len_b2;
    Real   len1;

    ax = p1[X] - p0[X];
    ay = p1[Y] - p0[Y];
    az = p1[Z] - p0[Z];

    bx = p2[X] - p0[X];
    by = p2[Y] - p0[Y];
    bz = p2[Z] - p0[Z];

    cx = p3[X] - p0[X];
    cy = p3[Y] - p0[Y];
    cz = p3[Z] - p0[Z];

    abx = ay * bz - az * by;
    aby = az * bx - ax * bz;
    abz = ax * by - ay * bx;

    hx = by * abz - bz * aby;
    hy = bz * abx - bx * abz;
    hz = bx * aby - by * abx;

    len1 = abx * abx + aby * aby + abz * abz;
    len2 = hx * hx + hy * hy + hz * hz;

    xn = (cx * hx + cy * hy + cz * hz) / sqrt( len2 );
    yn = (cx * abx + cy * aby + cz * abz) / sqrt( len1 );

    len = xn * xn + yn * yn;

    if( len == 0.0 )
    {
        *x = -1.0;
        *y = 0.0;
        return;
    }

    len = sqrt( len );

    *x = xn / len;
    *y = yn / len;
}
#endif

public  void  get_bending_xy(
    Real   p0[],
    Real   p1[],
    Real   p2[],
    Real   p3[],
    Real   *x,
    Real   *y )
{
    Real   ax, ay, az, bx, by, bz, cx, cy, cz;
    Real   abx, aby, abz, len, len2, xn, yn, hx, hy, hz, len_b2;

    ax = p1[X] - p0[X];
    ay = p1[Y] - p0[Y];
    az = p1[Z] - p0[Z];

    bx = p2[X] - p0[X];
    by = p2[Y] - p0[Y];
    bz = p2[Z] - p0[Z];

    cx = p3[X] - p0[X];
    cy = p3[Y] - p0[Y];
    cz = p3[Z] - p0[Z];

    abx = ay * bz - az * by;
    aby = az * bx - ax * bz;
    abz = ax * by - ay * bx;

    hx = by * abz - bz * aby;
    hy = bz * abx - bx * abz;
    hz = bx * aby - by * abx;

    len_b2 = bx * bx + by * by + bz * bz;

    xn = cx * hx + cy * hy + cz * hz; 
    yn = (cx * abx + cy * aby + cz * abz) * sqrt( len_b2 );

    len2 = xn * xn + yn * yn;
    if( len2 != 0.0 )
    {
        len = sqrt( len2 );
        *x = xn / len;
        *y = yn / len;
    }
    else
    {
        *x = -1.0;
        *y = 0.0;
    }

#ifdef TEST_DEBUG
#define TEST_DEBUG
{
    Real  test_x, test_y;

    test_get( p0, p1, p2, p3, &test_x, &test_y );

    if( !numerically_close( test_x, *x, 1.0e-3 ) ||
        !numerically_close( test_y, *y, 1.0e-3 ) )
    {
        print( " %g %g %g %g\n", test_x, test_y, *x, *y );
        test_get( p0, p1, p2, p3, &test_x, &test_y );
    }
}
#endif
}

private  void  get_bending_xy_and_deriv(
    Real   x0,
    Real   y0,
    Real   z0,
    Real   x1,
    Real   y1,
    Real   z1,
    Real   x2,
    Real   y2,
    Real   z2,
    Real   x3,
    Real   y3,
    Real   z3,
    Real   *x,
    Real   *y,
    Real   x_deriv[],
    Real   y_deriv[] )
{
    int    d;
    Real   ax, ay, az, bx, by, bz, cx, cy, cz;
    Real   abx, aby, abz, len_b2, len;
    Real   d_ax[12], d_ay[12], d_az[12];
    Real   d_bx[12], d_by[12], d_bz[12];
    Real   d_cx[12], d_cy[12], d_cz[12];
    Real   d_abx, d_aby, d_abz, d_len_b2, d_len, len_b;
    Real   first, d_first;
    Real   xn, yn, d_xn, d_yn, len2, d_len2;
    Real   hx, hy, hz, d_hx, d_hy, d_hz;

    for_less( d, 0, 12 )
    {
        d_ax[d] = 0.0;
        d_ay[d] = 0.0;
        d_az[d] = 0.0;
        d_bx[d] = 0.0;
        d_by[d] = 0.0;
        d_bz[d] = 0.0;
        d_cx[d] = 0.0;
        d_cy[d] = 0.0;
        d_cz[d] = 0.0;
    }

    ax = x1 - x0;
    ay = y1 - y0;
    az = z1 - z0;

    bx = x2 - x0;
    by = y2 - y0;
    bz = z2 - z0;

    cx = x3 - x0;
    cy = y3 - y0;
    cz = z3 - z0;

    d_ax[0] = -1.0;
    d_ay[1] = -1.0;
    d_az[2] = -1.0;

    d_bx[0] = -1.0;
    d_by[1] = -1.0;
    d_bz[2] = -1.0;

    d_cx[0] = -1.0;
    d_cy[1] = -1.0;
    d_cz[2] = -1.0;

    d_ax[3] = 1.0;
    d_ay[4] = 1.0;
    d_az[5] = 1.0;

    d_bx[6] = 1.0;
    d_by[7] = 1.0;
    d_bz[8] = 1.0;

    d_cx[9] = 1.0;
    d_cy[10] = 1.0;
    d_cz[11] = 1.0;

    abx = ay * bz - az * by;
    aby = az * bx - ax * bz;
    abz = ax * by - ay * bx;

    hx = by * abz - bz * aby;
    hy = bz * abx - bx * abz;
    hz = bx * aby - by * abx;

    len_b2 = bx * bx + by * by + bz * bz;
    if( len_b2 == 0.0 )
    {
        handle_internal_error( "len_b2" );
    }
    len_b = sqrt( len_b2 );

    xn = cx * hx + cy * hy + cz * hz; 
    first = cx * abx + cy * aby + cz * abz;
    yn = first * len_b;

    len2 = xn * xn + yn * yn;
    if( len2 != 0.0 )
    {
        len = sqrt( len2 );
        *x = xn / len;
        *y = yn / len;
    }
    else
    {
        print_error( "get_bending_xy_and_deriv(): flattened out edge" );
        *x = -1.0;
        *y = 0.0;
        for_less( d, 0, 12 )
        {
            x_deriv[d] = 0.0;
            y_deriv[d] = 0.0;
        }
        return;
    }

    for_less( d, 0, 12 )
    {
        d_abx = ay * d_bz[d] + d_ay[d] * bz - az * d_by[d] - d_az[d] * by;
        d_aby = az * d_bx[d] + d_az[d] * bx - ax * d_bz[d] - d_ax[d] * bz;
        d_abz = ax * d_by[d] + d_ax[d] * by - ay * d_bx[d] - d_ay[d] * bx;

        d_hx = by * d_abz + d_by[d] * abz - bz * d_aby - d_bz[d] * aby;
        d_hy = bz * d_abx + d_bz[d] * abx - bx * d_abz - d_bx[d] * abz;
        d_hz = bx * d_aby + d_bx[d] * aby - by * d_abx - d_by[d] * abx;

        d_len_b2 = 2.0 * bx * d_bx[d] + 2.0 * by * d_by[d] + 2.0 * bz * d_bz[d];

        d_first = cx * d_abx + d_cx[d] * abx +
                  cy * d_aby + d_cy[d] * aby +
                  cz * d_abz + d_cz[d] * abz;

        d_xn = cx * d_hx + d_cx[d] * hx +
               cy * d_hy + d_cy[d] * hy +
               cz * d_hz + d_cz[d] * hz;

        d_yn = first * 0.5 / len_b * d_len_b2 + d_first * len_b;

        d_len2 = 2.0 * xn * d_xn + 2.0 * yn * d_yn;
        d_len = 0.5 / len * d_len2;
        x_deriv[d] = d_xn / len - xn / len / len * d_len;
        y_deriv[d] = d_yn / len - yn / len / len * d_len;
    }
}

#ifdef old
{
    int    which;
    Real   ax, ay, az, bx, by, bz, cx, cy, cz;
    Real   a_dot_a, a_dot_b, b_dot_b, a_dot_c, b_dot_c;
    Real   abx, aby, abz, len_ab2, len_aab2, len, d_len;
    Real   d_len_ab2, d_len_aab2, len_sq;
    Real   dx, dy, dy_top, d_sqrt, factor;
    Real   d_a_dot_a[12], d_a_dot_b[12], d_a_dot_c[12];
    Real   d_b_dot_b[12], d_b_dot_c[12];
    Real   d_abx[12], d_aby[12], d_abz[12];
    Real   d_cx[12], d_cy[12], d_cz[12];

    ax = x1 - x0;
    ay = y1 - y0;
    az = z1 - z0;
    bx = x2 - x0;
    by = y2 - y0;
    bz = z2 - z0;
    cx = x3 - x0;
    cy = y3 - y0;
    cz = z3 - z0;

    a_dot_a = ax * ax + ay * ay + az * az;
    a_dot_b = ax * bx + ay * by + az * bz;
    a_dot_c = ax * cx + ay * cy + az * cz;
    b_dot_b = bx * bx + by * by + bz * bz;
    b_dot_c = bx * cx + by * cy + bz * cz;

    abx = ay * bz - az * by;
    aby = az * bx - ax * bz;
    abz = ax * by - ay * bx;

    for_less( which, 0, 12 )
    {
        d_cx[which] = 0.0;
        d_cy[which] = 0.0;
        d_cz[which] = 0.0;
    }

    d_cx[0] = -1.0;
    d_cx[9] = 1.0;
    d_cy[1] = -1.0;
    d_cy[10] = 1.0;
    d_cz[2] = -1.0;
    d_cz[11] = 1.0;

    d_a_dot_a[0] = -2.0*x1+2.0*x0;
    d_a_dot_a[1] = -2.0*y1+2.0*y0;
    d_a_dot_a[2] = -2.0*z1+2.0*z0;
    d_a_dot_a[3] = 2.0*x1-2.0*x0;
    d_a_dot_a[4] = 2.0*y1-2.0*y0;
    d_a_dot_a[5] = 2.0*z1-2.0*z0;
    d_a_dot_a[6] = 0.0;
    d_a_dot_a[7] = 0.0;
    d_a_dot_a[8] = 0.0;
    d_a_dot_a[9] = 0.0;
    d_a_dot_a[10] = 0.0;
    d_a_dot_a[11] = 0.0;

    d_a_dot_b[0] = -x2+2.0*x0-x1;
    d_a_dot_b[1] = -y2+2.0*y0-y1;
    d_a_dot_b[2] = -z2+2.0*z0-z1;
    d_a_dot_b[3] = x2-x0;
    d_a_dot_b[4] = y2-y0;
    d_a_dot_b[5] = z2-z0;
    d_a_dot_b[6] = x1-x0;
    d_a_dot_b[7] = y1-y0;
    d_a_dot_b[8] = z1-z0;
    d_a_dot_b[9] = 0.0;
    d_a_dot_b[10] = 0.0;
    d_a_dot_b[11] = 0.0;

    d_a_dot_c[0] = -x3+2.0*x0-x1;
    d_a_dot_c[1] = -y3+2.0*y0-y1;
    d_a_dot_c[2] = -z3+2.0*z0-z1;
    d_a_dot_c[3] = x3-x0;
    d_a_dot_c[4] = y3-y0;
    d_a_dot_c[5] = z3-z0;
    d_a_dot_c[6] = 0.0;
    d_a_dot_c[7] = 0.0;
    d_a_dot_c[8] = 0.0;
    d_a_dot_c[9] = x1-x0;
    d_a_dot_c[10] = y1-y0;
    d_a_dot_c[11] = z1-z0;

    d_b_dot_b[0] = -2.0*x2+2.0*x0;
    d_b_dot_b[1] = -2.0*y2+2.0*y0;
    d_b_dot_b[2] = -2.0*z2+2.0*z0;
    d_b_dot_b[3] = 0.0;
    d_b_dot_b[4] = 0.0;
    d_b_dot_b[5] = 0.0;
    d_b_dot_b[6] = 2.0*x2-2.0*x0;
    d_b_dot_b[7] = 2.0*y2-2.0*y0;
    d_b_dot_b[8] = 2.0*z2-2.0*z0;
    d_b_dot_b[9] = 0.0;
    d_b_dot_b[10] = 0.0;
    d_b_dot_b[11] = 0.0;

    d_b_dot_c[0] = -x3+2.0*x0-x2;
    d_b_dot_c[1] = -y3+2.0*y0-y2;
    d_b_dot_c[2] = -z3+2.0*z0-z2;
    d_b_dot_c[3] = 0.0;
    d_b_dot_c[4] = 0.0;
    d_b_dot_c[5] = 0.0;
    d_b_dot_c[6] = x3-x0;
    d_b_dot_c[7] = y3-y0;
    d_b_dot_c[8] = z3-z0;
    d_b_dot_c[9] = x2-x0;
    d_b_dot_c[10] = y2-y0;
    d_b_dot_c[11] = z2-z0;

    d_abx[0] = 0.0;
    d_abx[1] = -z2+z1;
    d_abx[2] = -y1+y2;
    d_abx[3] = 0.0;
    d_abx[4] = z2-z0;
    d_abx[5] = -y2+y0;
    d_abx[6] = 0.0;
    d_abx[7] = -z1+z0;
    d_abx[8] = y1-y0;
    d_abx[9] = 0.0;
    d_abx[10] = 0.0;
    d_abx[11] = 0.0;

    d_aby[0] = -z1+z2;
    d_aby[1] = 0.0;
    d_aby[2] = -x2+x1;
    d_aby[3] = -z2+z0;
    d_aby[4] = 0.0;
    d_aby[5] = x2-x0;
    d_aby[6] = z1-z0;
    d_aby[7] = 0.0;
    d_aby[8] = -x1+x0;
    d_aby[9] = 0.0;
    d_aby[10] = 0.0;
    d_aby[11] = 0.0;

    d_abz[0] = -y2+y1;
    d_abz[1] = -x1+x2;
    d_abz[2] = 0.0;
    d_abz[3] = y2-y0;
    d_abz[4] = -x2+x0;
    d_abz[5] = 0.0;
    d_abz[6] = -y1+y0;
    d_abz[7] = x1-x0;
    d_abz[8] = 0.0;
    d_abz[9] = 0.0;
    d_abz[10] = 0.0;
    d_abz[11] = 0.0;

    len_ab2 = abx * abx + aby * aby + abz * abz;

    if( len_ab2 == 0.0 )
    {
        print( "get_bending_xy() : len_ab2 == 0.0\n" );
        len_ab2 = 1.0;
    }

    len_aab2 = a_dot_a * (a_dot_a * b_dot_b - a_dot_b * a_dot_b);

    if( len_aab2 == 0.0 )
    {
        print( "get_bending_xy() : len_aab2 == 0.0\n" );
        len_aab2 = 1.0;
    }

    factor = sqrt( len_aab2 / len_ab2 );

    *x = (a_dot_b * a_dot_c - a_dot_a * b_dot_c);
    *y = (cx * abx + cy * aby + cz * abz) * factor;

    len_sq = *x * *x + *y * *y;
    if( len_sq == 0.0 )  len_sq = 1.0;

    len = sqrt( len_sq );

    for_less( which, 0, 12 )
    {
        dx = d_a_dot_b[which] * a_dot_c + a_dot_b * d_a_dot_c[which] -
             (d_a_dot_a[which] * b_dot_c + a_dot_a * d_b_dot_c[which]);

        dy_top = d_cx[which] * abx + cx * d_abx[which] +
                 d_cy[which] * aby + cy * d_aby[which] +
                 d_cz[which] * abz + cz * d_abz[which];

        d_len_ab2 = 2.0 * abx * d_abx[which] +
                    2.0 * aby * d_aby[which] +
                    2.0 * abz * d_abz[which];

        d_len_aab2 = d_a_dot_a[which] * (a_dot_a * b_dot_b - a_dot_b * a_dot_b)
                     + a_dot_a *
                     (d_a_dot_a[which] * b_dot_b + a_dot_a * d_b_dot_b[which] -
                      2.0 * a_dot_b * d_a_dot_b[which]);

        d_sqrt = 0.5 / factor *
                 (d_len_aab2 / len_ab2 -
                  len_aab2 / len_ab2 / len_ab2 * d_len_ab2);

        dy = dy_top * factor + (cx * abx + cy * aby + cz * abz) * d_sqrt;

        d_len = 2.0 * *x * dx + 2.0 * *y * dy;;

        x_deriv[which] = dx / len - 0.5 * *x / len / len_sq * d_len;
        y_deriv[which] = dy / len - 0.5 * *y / len / len_sq * d_len;
    }

    *x /= len;
    *y /= len;
}
#endif

#ifdef OLD
private  void  get_bending_xy_and_deriv(
    Real   x0,
    Real   y0,
    Real   z0,
    Real   x1,
    Real   y1,
    Real   z1,
    Real   x2,
    Real   y2,
    Real   z2,
    Real   x3,
    Real   y3,
    Real   z3,
    Real   *x,
    Real   *y,
    Real   x_deriv[],
    Real   y_deriv[] )
{
    int             which;
    Real            xs[12], save, x_prev, y_prev, x_next, y_next;
    Real            used_delta, max_dx, max_dy, dx, dy;
    static BOOLEAN  first = TRUE;
    static Real     delta;

    if( first )
    {
        first = FALSE;
        if( getenv("DELTA") == NULL ||
            sscanf( getenv("DELTA"), "%lf", &delta ) != 1 || delta <= 0.0 )
            delta = 1.0e-5;
    }

    get_bending_xy( x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y );

    max_dx = delta * FABS( *x );
    max_dy = delta * FABS( *y );

    xs[0] =  x0;
    xs[1] =  y0;
    xs[2] =  z0;
    xs[3] =  x1;
    xs[4] =  y1;
    xs[5] =  z1;
    xs[6] =  x2;
    xs[7] =  y2;
    xs[8] =  z2;
    xs[9] =  x3;
    xs[10] = y3;
    xs[11] = z3;

    for_less( which, 0, 12 )
    {
        save = xs[which];

        used_delta = 1e-4;

        do
        {
            used_delta /= 2.0;

            xs[which] = save - used_delta;
            get_bending_xy( xs[0], xs[1], xs[2], xs[3], xs[4], xs[5], xs[6], xs[7],
                            xs[8], xs[9], xs[10], xs[11], &x_prev, &y_prev );

            xs[which] = save + used_delta;
            get_bending_xy( xs[0], xs[1], xs[2], xs[3], xs[4], xs[5], xs[6], xs[7],
                            xs[8], xs[9], xs[10], xs[11], &x_next, &y_next );

            dx = FABS( x_next - x_prev );
            dy = FABS( y_next - y_prev );
        }
        while( dx > max_dx || dy > max_dy );

        xs[which] = save;

        x_deriv[which] = (x_next - x_prev) / 2.0 / used_delta;
        y_deriv[which] = (y_next - y_prev) / 2.0 / used_delta;
    }
}
#endif

private  Real  convert_xy_to_curvature(
    Real     x,
    Real     y )
{
    Real   dx, curvature;

    dx = x - -1.0;
    curvature = (dx * dx + y * y) / 4.0;

    if( y >= 0.0 )
        curvature = -curvature;

    return( curvature );
}

private  Real  convert_xy_to_curvature_deriv(
    Real     x,
    Real     y,
    Real     *x_deriv,
    Real     *y_deriv )
{
    Real   dx, curvature;

    dx = x - -1.0;
    curvature = (dx * dx + y * y) / 4.0;

    *x_deriv = 2.0 * dx / 4.0;
    *y_deriv = 2.0 * y / 4.0;

    if( y >= 0.0 )
    {
        curvature = -curvature;
        *x_deriv = - *x_deriv;
        *y_deriv = - *y_deriv;
    }

    return( curvature );
}

private  Real   evaluate_bend_fit(
    Real          bend_weight,
    Real          max_bend_weight,
    Real          min_bend,
    Real          max_bend,
    int           n_points,
    Real          parameters[],
    Smallest_int  active_flags[],
    int           n_neighbours[],
    int           *neighbours[],
    float         model_x[],
    float         model_y[] )
{
    int     point, neigh, ind, n, ind0, ind1, ind2, ind3, *neigh_ptr, n_neighs;
    int     next_neigh, prev_neigh;
    Real    x, y, curvature, model_curvature, diff, delta;
    Real    fit1, fit2, dx, dy, len;
    Real    *this_ptr, *current_ptr, *prev_ptr, *next_ptr;
    BOOLEAN active;

    if( bend_weight < 0.0 )
        return( 0.0 );

    fit1 = 0.0;
    fit2 = 0.0;
    ind = 0;

    if( active_flags == NULL )
    {
        for_less( point, 0, n_points )
        {
            this_ptr = &parameters[3*point];

            neigh_ptr = neighbours[point];
            n_neighs = n_neighbours[point];

            current_ptr = &parameters[3*neigh_ptr[n_neighs-1]];
            next_ptr = &parameters[3*neigh_ptr[0]];

            for_less( n, 0, n_neighs )
            {
                prev_ptr = current_ptr;
                current_ptr = next_ptr;
                next_ptr = &parameters[3*neigh_ptr[(n+1)%n_neighs]];

                if( current_ptr < this_ptr )
                    continue;

                get_bending_xy( this_ptr, prev_ptr, current_ptr, next_ptr,
                                &x, &y );

                curvature = convert_xy_to_curvature( x, y );
                model_curvature = convert_xy_to_curvature( (Real) model_x[ind],
                                                        (Real) model_y[ind] );
                ++ind;

                delta = curvature - model_curvature;

                fit1 += delta * delta;

                if( min_bend >= max_bend || max_bend_weight <= 0.0 ||
                    delta >= min_bend && delta <= max_bend )
                    continue;

                if( delta < min_bend )
                    diff = delta - min_bend;
                else
                    diff = delta - max_bend;

                fit2 += diff * diff;
            }
        }
    }
    else
    {
        handle_internal_error( "Not implemented yet\n" );

        for_less( point, 0, n_points )
        {
            ind0 = IJ(point,0,3);

            active = active_flags == NULL || active_flags[point];

            neigh_ptr = neighbours[point];
            n_neighs = n_neighbours[point];

            neigh = neigh_ptr[n_neighs-1];
            next_neigh = neigh_ptr[0];

            for_less( n, 0, n_neighs )
            {
                prev_neigh = neigh;
                neigh = next_neigh;
                next_neigh = neigh_ptr[(n+1)%n_neighs];

                if( !THIS_IS_UNIQUE_EDGE( point, neigh ) )
                    continue;

                if( !active && !active_flags[neigh] &&
                    !active_flags[next_neigh] && !active_flags[prev_neigh] )
                {
                    ++ind;
                    continue;
                }

                ind1 = IJ(prev_neigh,0,3);
                ind2 = IJ(neigh,0,3);
                ind3 = IJ(next_neigh,0,3);

                get_bending_xy( &parameters[ind0],
                                &parameters[ind1],
                                &parameters[ind2],
                                &parameters[ind3],
                                &x, &y );

                dx = x - (Real) model_x[ind];
                dy = y - (Real) model_y[ind];
                ++ind;

                len = dx * dx + dy * dy;

                fit1 += len;

                if( min_bend >= max_bend || max_bend_weight <= 0.0 )
                    continue;

                if( len > max_bend * max_bend )
                {
                    if( len != 0.0 )
                        len = sqrt( len );

                    len = len - max_bend;

                    fit2 += len * len;
                }
            }
        }
    }

    return( fit1 * bend_weight + fit2 * max_bend_weight );
}

private  void   evaluate_bend_fit_deriv(
    Real          bend_weight,
    Real          max_bend_weight,
    Real          min_bend,
    Real          max_bend,
    int           n_points,
    Real          parameters[],
    int           n_neighbours[],
    int           *neighbours[],
    float         model_x[],
    float         model_y[],
    Real          deriv[] )
{
    int     point, neigh, ind, n, ind0, ind1, ind2, ind3, *neigh_ptr, n_neighs;
    int     next_neigh, prev_neigh;
    Real    x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y;
    Real    x_factor, y_factor;
    Real    x_deriv[12], y_deriv[12];
    Real    delta, curvature, model_curvature, d_curv_x, d_curv_y;

    if( bend_weight < 0.0 )
        return;

    ind = 0;

    for_less( point, 0, n_points )
    {
        ind0 = IJ(point,0,3);
        x0 = parameters[ind0+0];
        y0 = parameters[ind0+1];
        z0 = parameters[ind0+2];

        neigh_ptr = neighbours[point];
        n_neighs = n_neighbours[point];

        for_less( n, 0, n_neighs )
        {
            neigh = neigh_ptr[n];
            next_neigh = neigh_ptr[(n+1)%n_neighs];
            prev_neigh = neigh_ptr[(n-1+n_neighs)%n_neighs];

            if( !THIS_IS_UNIQUE_EDGE( point, neigh ) )
                continue;

            ind1 = IJ(prev_neigh,0,3);
            x1 = parameters[ind1+0];
            y1 = parameters[ind1+1];
            z1 = parameters[ind1+2];

            ind2 = IJ(neigh,0,3);
            x2 = parameters[ind2+0];
            y2 = parameters[ind2+1];
            z2 = parameters[ind2+2];

            ind3 = IJ(next_neigh,0,3);
            x3 = parameters[ind3+0];
            y3 = parameters[ind3+1];
            z3 = parameters[ind3+2];

            get_bending_xy_and_deriv( x0, y0, z0, x1, y1, z1, x2, y2, z2,
                                      x3, y3, z3, &x, &y,
                                      x_deriv, y_deriv );

            curvature = convert_xy_to_curvature_deriv( x, y,
                                                       &d_curv_x, &d_curv_y );
            model_curvature = convert_xy_to_curvature( (Real) model_x[ind],
                                                       (Real) model_y[ind] );
            ++ind;

            delta = curvature - model_curvature;

            x_factor = bend_weight * 2.0 * delta * d_curv_x;
            y_factor = bend_weight * 2.0 * delta * d_curv_y;

            deriv[ind0+0] += x_factor * x_deriv[0] + y_factor * y_deriv[0];
            deriv[ind0+1] += x_factor * x_deriv[1] + y_factor * y_deriv[1];
            deriv[ind0+2] += x_factor * x_deriv[2] + y_factor * y_deriv[2];

            deriv[ind1+0] += x_factor * x_deriv[3] + y_factor * y_deriv[3];
            deriv[ind1+1] += x_factor * x_deriv[4] + y_factor * y_deriv[4];
            deriv[ind1+2] += x_factor * x_deriv[5] + y_factor * y_deriv[5];

            deriv[ind2+0] += x_factor * x_deriv[6] + y_factor * y_deriv[6];
            deriv[ind2+1] += x_factor * x_deriv[7] + y_factor * y_deriv[7];
            deriv[ind2+2] += x_factor * x_deriv[8] + y_factor * y_deriv[8];

            deriv[ind3+0] += x_factor * x_deriv[9] + y_factor * y_deriv[9];
            deriv[ind3+1] += x_factor * x_deriv[10]+ y_factor * y_deriv[10];
            deriv[ind3+2] += x_factor * x_deriv[11]+ y_factor * y_deriv[11];

            if( min_bend >= max_bend || max_bend_weight <= 0.0 ||
                delta >= min_bend && delta <= max_bend )
                continue;

            if( delta < min_bend )
                delta = delta - min_bend;
            else
                delta = delta - max_bend;
                 
            x_factor = bend_weight * 2.0 * delta * d_curv_x;
            y_factor = bend_weight * 2.0 * delta * d_curv_y;

            deriv[ind0+0] += x_factor * x_deriv[0] + y_factor * y_deriv[0];
            deriv[ind0+1] += x_factor * x_deriv[1] + y_factor * y_deriv[1];
            deriv[ind0+2] += x_factor * x_deriv[2] + y_factor * y_deriv[2];

            deriv[ind1+0] += x_factor * x_deriv[3] + y_factor * y_deriv[3];
            deriv[ind1+1] += x_factor * x_deriv[4] + y_factor * y_deriv[4];
            deriv[ind1+2] += x_factor * x_deriv[5] + y_factor * y_deriv[5];

            deriv[ind2+0] += x_factor * x_deriv[6] + y_factor * y_deriv[6];
            deriv[ind2+1] += x_factor * x_deriv[7] + y_factor * y_deriv[7];
            deriv[ind2+2] += x_factor * x_deriv[8] + y_factor * y_deriv[8];

            deriv[ind3+0] += x_factor * x_deriv[9] + y_factor * y_deriv[9];
            deriv[ind3+1] += x_factor * x_deriv[10]+ y_factor * y_deriv[10];
            deriv[ind3+2] += x_factor * x_deriv[11]+ y_factor * y_deriv[11];
        }
    }
}

private  Real  get_nth_root(
    Real   x,
    int    n )
{
    return( exp( log(x) / (Real) n ) );
}

private  Real   evaluate_inter_surface_fit(
    Real               weight, 
    Real               max_dist_weight1,
    Real               max_dist_weight2,
    int                n_weight_steps,
    int                n_connections,
    connection_struct  connections[],
    Real               parameters1[],
    Smallest_int       active_flags1[],
    Real               parameters2[],
    Smallest_int       active_flags2[] )
{
    int                con, p1_index, p2_index, i;
    Real               fit1, fit2, dist_threshold;
    Real               dist, scale_factor, weight_factor;
    Real               x1, y1, z1, x2, y2, z2, dx, dy, dz, off;
    connection_struct  *c;
    BOOLEAN            using_max_dist;

    if( weight < 0.0 )
        return( 0.0 );

    fit1 = 0.0;
    fit2 = 0.0;

    if( max_dist_weight1 > 0.0 )
    {
        scale_factor = get_nth_root( max_dist_weight2 / max_dist_weight1,
                                     n_weight_steps-1 );
        using_max_dist = TRUE;
    }
    else
        using_max_dist = FALSE;

    for_less( con, 0, n_connections )
    {
        c = &connections[con];
        if( active_flags1 != NULL && !active_flags1[c->surface1_point] &&
            active_flags2 != NULL && !active_flags2[c->surface2_point] )
            continue;

        p1_index = IJ( c->surface1_point, 0, 3 );
        x1 = parameters1[p1_index+0];
        y1 = parameters1[p1_index+1];
        z1 = parameters1[p1_index+2];

        p2_index = IJ( c->surface2_point, 0, 3 );
        x2 = parameters2[p2_index+0];
        y2 = parameters2[p2_index+1];
        z2 = parameters2[p2_index+2];

        dx = x2 - x1;
        dy = y2 - y1;
        dz = z2 - z1;

        dist = dx * dx + dy * dy + dz * dz;
        if( dist > 0.0 )
            dist = sqrt( dist );

        dist -= c->desired_distance;
        fit1 += dist * dist;

        if( using_max_dist && c->min_distance1 < c->max_distance1 )
        {
            if( dist < c->min_distance1 )
            {
                weight_factor = max_dist_weight1;
                for_less( i, 0, n_weight_steps )
                {
                    dist_threshold = INTERPOLATE( (Real) (i) /
                                                  (Real) (n_weight_steps-1),
                                                  c->min_distance1,
                                                  c->min_distance2 );
                                     
                    if( dist >= dist_threshold )
                        break;

                    off = dist - dist_threshold;
                    fit2 += weight_factor * off * off;
                    weight_factor *= scale_factor;
                }
            }
            else if( dist > c->max_distance1 )
            {
                weight_factor = max_dist_weight1;
                for_less( i, 0, n_weight_steps )
                {
                    dist_threshold = INTERPOLATE( (Real) (i) /
                                                  (Real) (n_weight_steps-1),
                                                  c->max_distance1,
                                                  c->max_distance2 );
                                     
                    if( dist <= dist_threshold )
                        break;

                    off = dist - dist_threshold;
                    fit2 += weight_factor * off * off;
                    weight_factor *= scale_factor;
                }
            }
        }
    }

    return( fit1 * weight + fit2 );
}

private  void   evaluate_inter_surface_fit_deriv(
    Real               weight, 
    Real               max_dist_weight1,
    Real               max_dist_weight2,
    int                n_weight_steps,
    int                n_connections,
    connection_struct  connections[],
    int                n_neighbours[],
    int                *neighbours[],
    Real               parameters1[],
    Real               parameters2[],
    Real               deriv1[],
    Real               deriv2[] )
{
    int                con, p1_index, p2_index, i, n;
    int                *neigh_ptr, n_neighs, point, neigh, ind1;
    Real               dist, diff, factor, dist_threshold;
    Real               x1, y1, z1, x2, y2, z2, dx, dy, dz, off;
    Real               scale_factor, weight_factor;
    Real               x, y, z, tx, ty, tz, dot;
    Real               len, cx, cy, cz, nx, ny, nz;
    BOOLEAN            using_max_dist;
    connection_struct  *c;

    if( weight < 0.0 )
        return;

    if( max_dist_weight1 > 0.0 )
    {
        scale_factor = get_nth_root( max_dist_weight2 / max_dist_weight1,
                                     n_weight_steps-1 );
        using_max_dist = TRUE;
    }
    else
        using_max_dist = FALSE;

    for_less( con, 0, n_connections )
    {
        c = &connections[con];
        p1_index = IJ( c->surface1_point, 0, 3 );
        x1 = parameters1[p1_index+0];
        y1 = parameters1[p1_index+1];
        z1 = parameters1[p1_index+2];

        p2_index = IJ( c->surface2_point, 0, 3 );
        x2 = parameters2[p2_index+0];
        y2 = parameters2[p2_index+1];
        z2 = parameters2[p2_index+2];

        dx = x2 - x1;
        dy = y2 - y1;
        dz = z2 - z1;

        dist = dx * dx + dy * dy + dz * dz;

        if( dist <= 0.0 )
        {
            if( c->desired_distance == 0.0 )
                continue;

            factor = 0.0;
            cx = 0.0;
            cy = 0.0;
            cz = 0.0;
            nx = 0.0;
            ny = 0.0;
            nz = 0.0;

            point = c->surface1_point;
            neigh_ptr = neighbours[point];
            n_neighs = n_neighbours[point];

            neigh = neigh_ptr[n_neighs-1];
            ind1 = IJ(neigh,0,3);
            tx = parameters1[ind1+0];
            ty = parameters1[ind1+1];
            tz = parameters1[ind1+2];

            for_less( n, 0, n_neighs )
            {
                x = tx;
                y = ty;
                z = tz;

                neigh = neigh_ptr[n];
                ind1 = IJ(neigh,0,3);

                tx = parameters1[ind1+0];
                ty = parameters1[ind1+1];
                tz = parameters1[ind1+2];
                cx += tx;
                cy += ty;
                cz += tz;

                nx += y * tz - ty * z;
                ny += z * tx - tz * x;
                nz += x * ty - tx * y;
            }

            factor = 1.0 / (Real) n_neighs;

            cx *= factor;
            cy *= factor;
            cz *= factor;

            dx = x1 - cx;
            dy = y1 - cy;
            dz = z1 - cz;

            dot = dx * nx + dy * ny + dz * nz;
            if( dot < 0.0 )
            {
                dx = -dx;
                dy = -dy;
                dz = -dz;
            }
            else if( dot == 0.0 )
            {
                dx = nx;
                dy = ny;
                dz = nz;
            }

            len = dx * dx + dy * dy + dz * dz;
            if( len > 0.0 )
            {
                len = sqrt( len );
                factor = weight / len;
                deriv1[p1_index+0] += factor * dx;
                deriv1[p1_index+1] += factor * dy;
                deriv1[p1_index+2] += factor * dz;
                deriv2[p2_index+0] += factor * -dx;
                deriv2[p2_index+1] += factor * -dy;
                deriv2[p2_index+2] += factor * -dz;
            }
        }
        else
        {
            dist = sqrt( dist );
            diff = dist - c->desired_distance;
            factor = weight * diff * 2.0 / dist;
            deriv1[p1_index+0] += factor * -dx;
            deriv1[p1_index+1] += factor * -dy;
            deriv1[p1_index+2] += factor * -dz;
            deriv2[p2_index+0] += factor * dx;
            deriv2[p2_index+1] += factor * dy;
            deriv2[p2_index+2] += factor * dz;
        }

        if( dist > 0.0 &&
            using_max_dist && c->min_distance1 < c->max_distance1 )
        {
            if( diff < c->min_distance1 )
            {
                factor = 0.0;
                weight_factor = max_dist_weight1;
                for_less( i, 0, n_weight_steps )
                {
                    dist_threshold = INTERPOLATE( (Real) (i) /
                                                  (Real) (n_weight_steps-1),
                                                  c->min_distance1,
                                                  c->min_distance2 );
                                     
                    if( diff >= dist_threshold )
                        break;

                    off = diff - dist_threshold;
                    factor += weight_factor * off * 2.0 / dist;
                    weight_factor *= scale_factor;
                }
                deriv1[p1_index+0] += factor * -dx;
                deriv1[p1_index+1] += factor * -dy;
                deriv1[p1_index+2] += factor * -dz;
                deriv2[p2_index+0] += factor * dx;
                deriv2[p2_index+1] += factor * dy;
                deriv2[p2_index+2] += factor * dz;
            }
            else if( diff > c->max_distance1 )
            {
                factor = 0.0;
                weight_factor = max_dist_weight1;
                for_less( i, 0, n_weight_steps )
                {
                    dist_threshold = INTERPOLATE( (Real) (i) /
                                                  (Real) (n_weight_steps-1),
                                                  c->max_distance1,
                                                  c->max_distance2 );
                                     
                    if( diff <= dist_threshold )
                        break;

                    off = diff - dist_threshold;
                    factor += weight_factor * off * 2.0 / dist;
                    weight_factor *= scale_factor;
                }

                deriv1[p1_index+0] += factor * -dx;
                deriv1[p1_index+1] += factor * -dy;
                deriv1[p1_index+2] += factor * -dz;
                deriv2[p2_index+0] += factor * dx;
                deriv2[p2_index+1] += factor * dy;
                deriv2[p2_index+2] += factor * dz;
            }
        }
    }
}

private  Real   evaluate_oversampled_inter_surface_fit(
    Real               weight, 
    Real               max_dist_weight1,
    Real               max_dist_weight2,
    int                n_weight_steps,
    int                n_connections,
    connection_struct  connections[],
    int                oversample,
    int                n_points,
    int                n_neighbours[],
    int                *neighbours[],
    Real               parameters1[],
    Smallest_int       active_flags1[],
    Real               parameters2[],
    Smallest_int       active_flags2[] )
{
    int                p1, p1_index, n1_index, n2_index, neigh1, neigh2, i;
    int                w1, w2, n, n_on_line, end, n_neighs, *neighs;
    Real               fit1, fit2, alpha2, one_minus_alpha1, off;
    Real               weight1, weight2, weight3, dist_threshold;
    Real               dist, min_dist1, min_dist2;
    Real               max_dist1, max_dist2;
    Real               dist1, dist2, dist3, desired_dist;
    Real               dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3;
    Real               dx, dy, dz, real_n_on_line;
    Real               scale_factor, weight_factor;
    Real               min_dist1a, min_dist1b, min_dist2a, min_dist2b;
    Real               min_dist3a, min_dist3b;
    Real               max_dist1a, max_dist1b, max_dist2a, max_dist2b;
    Real               max_dist3a, max_dist3b;
    BOOLEAN            using_max_dist;

    if( weight < 0.0 )
        return( 0.0 );

    fit1 = 0.0;
    fit2 = 0.0;

    if( max_dist_weight1 > 0.0 )
    {
        scale_factor = get_nth_root( max_dist_weight2 / max_dist_weight1,
                                     n_weight_steps-1 );
        using_max_dist = TRUE;
    }
    else
        using_max_dist = FALSE;

    for_less( p1, 0, n_points )
    {
        p1_index = IJ(p1,0,3);
        dist1 = connections[p1].desired_distance;
        min_dist1a = connections[p1].min_distance1;
        min_dist1b = connections[p1].min_distance2;
        max_dist1a = connections[p1].max_distance1;
        max_dist1b = connections[p1].max_distance2;

        n_neighs = n_neighbours[p1];
        neighs = neighbours[p1];

        neigh2 = neighs[n_neighs-1];
        for_less( n, 0, n_neighs )
        {
            neigh1 = neigh2;
            neigh2 = neighs[n];

        if( p1 > neigh1 || p1 > neigh2 )
            continue;

        if( active_flags1 != NULL  && !active_flags1[p1] &&
            !active_flags1[neigh1] && !active_flags1[neigh2] &&
            !active_flags2[p1]     && !active_flags2[neigh1] &&
            !active_flags2[neigh2] )
            continue;
            
        n1_index = IJ(neigh1,0,3);
        n2_index = IJ(neigh2,0,3);

        dx1 = parameters1[p1_index+0] - parameters2[p1_index+0];
        dx2 = parameters1[n1_index+0] - parameters2[n1_index+0];
        dx3 = parameters1[n2_index+0] - parameters2[n2_index+0];

        dy1 = parameters1[p1_index+1] - parameters2[p1_index+1];
        dy2 = parameters1[n1_index+1] - parameters2[n1_index+1];
        dy3 = parameters1[n2_index+1] - parameters2[n2_index+1];

        dz1 = parameters1[p1_index+2] - parameters2[p1_index+2];
        dz2 = parameters1[n1_index+2] - parameters2[n1_index+2];
        dz3 = parameters1[n2_index+2] - parameters2[n2_index+2];

        dist2 = connections[neigh1].desired_distance;
        dist3 = connections[neigh2].desired_distance;

        min_dist2a = connections[neigh1].min_distance1;
        min_dist2b = connections[neigh1].min_distance2;
        min_dist3a = connections[neigh2].min_distance1;
        min_dist3b = connections[neigh2].min_distance2;

        max_dist2a = connections[neigh1].max_distance1;
        max_dist2b = connections[neigh1].max_distance2;
        max_dist3a = connections[neigh2].max_distance1;
        max_dist3b = connections[neigh2].max_distance2;

        for_inclusive( w1, 0, oversample )
        {
            weight2 = (Real) w1 / (Real) (oversample+1);
            one_minus_alpha1 = 1.0 - weight2;

            n_on_line = oversample - w1 + 1;
            real_n_on_line = (Real) n_on_line;

            if( neigh1 > neigh2 || w1 == 0 )
                end = n_on_line-1;
            else
                end = n_on_line;

            for_inclusive( w2, 1, end )
            {
                alpha2 = (Real) w2 / real_n_on_line;

                weight1 = one_minus_alpha1 * (1.0 - alpha2);
                weight3 = one_minus_alpha1 * alpha2;

                dx = weight1 * dx1 + weight2 * dx2 + weight3 * dx3;
                dy = weight1 * dy1 + weight2 * dy2 + weight3 * dy3;
                dz = weight1 * dz1 + weight2 * dz2 + weight3 * dz3;

                dist = dx * dx + dy * dy + dz * dz;
                if( dist > 0.0 )
                    dist = sqrt( dist );

                desired_dist = weight1 * dist1 + weight2 * dist2 +
                               weight3 * dist3;
                min_dist1 = weight1 * min_dist1a + weight2 * min_dist2a +
                            weight3 * min_dist3a;
                min_dist2 = weight1 * min_dist1b + weight2 * min_dist2b +
                            weight3 * min_dist3b;
                max_dist1 = weight1 * max_dist1a + weight2 * max_dist2a +
                            weight3 * max_dist3a;
                max_dist2 = weight1 * max_dist1b + weight2 * max_dist2b +
                            weight3 * max_dist3b;

                dist -= desired_dist;
                fit1 += dist * dist;

                if( using_max_dist && min_dist1 < max_dist1 )
                {
                    if( dist < min_dist1 )
                    {
                        weight_factor = max_dist_weight1;
                        for_less( i, 0, n_weight_steps )
                        {
                            dist_threshold = INTERPOLATE( (Real) (i) /
                                                 (Real) (n_weight_steps-1),
                                                 min_dist1, min_dist2 );
                                     
                            if( dist >= dist_threshold )
                                break;

                            off = dist - dist_threshold;
                            fit2 += weight_factor * off * off;
                            weight_factor *= scale_factor;
                        }
                    }
                    else if( dist > max_dist1 )
                    {
                        weight_factor = max_dist_weight1;
                        for_less( i, 0, n_weight_steps )
                        {
                            dist_threshold = INTERPOLATE( (Real) (i) /
                                                  (Real) (n_weight_steps-1),
                                                  max_dist1, max_dist2 );
                                     
                            if( dist <= dist_threshold )
                                break;

                            off = dist - dist_threshold;
                            fit2 += weight_factor * off * off;
                            weight_factor *= scale_factor;
                        }
                    }
                }
            }
        }
        }
    }

    return( fit1 * weight + fit2 );
}

private  void   evaluate_oversampled_inter_surface_fit_deriv(
    Real               weight, 
    Real               max_dist_weight1,
    Real               max_dist_weight2,
    int                n_weight_steps,
    int                n_connections,
    connection_struct  connections[],
    int                oversample,
    int                n_points,
    int                n_neighbours[],
    int                *neighbours[],
    Real               parameters1[],
    Real               parameters2[],
    Real               deriv1[],
    Real               deriv2[] )
{
    int                p1, p1_index, n1_index, n2_index, neigh1, neigh2, i;
    int                w1, w2, n;
    Real               alpha1, alpha2, factor, diff, off, dist_threshold;
    Real               weight1, weight2, weight3;
    Real               dist, min_dist1, min_dist2;
    Real               max_dist1, max_dist2;
    Real               dist1, dist2, dist3, desired_dist;
    Real               x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z, dx, dy, dz;
    Real               a1, b1, c1, a2, b2, c2, a3, b3, c3, a, b, c;
    Real               scale_factor, weight_factor;
    Real               min_dist1a, min_dist1b, min_dist2a, min_dist2b;
    Real               min_dist3a, min_dist3b;
    Real               max_dist1a, max_dist1b, max_dist2a, max_dist2b;
    Real               max_dist3a, max_dist3b;
    BOOLEAN            using_max_dist;

    if( weight < 0.0 )
        return;

    if( max_dist_weight1 > 0.0 )
    {
        scale_factor = get_nth_root( max_dist_weight2 / max_dist_weight1,
                                     n_weight_steps-1 );
        using_max_dist = TRUE;
    }
    else
        using_max_dist = FALSE;

    for_less( p1, 0, n_points )
    {
        for_less( n, 0, n_neighbours[p1] )
        {
        neigh1 = neighbours[p1][n];
        neigh2 = neighbours[p1][(n+1)%n_neighbours[p1]];

        if( p1 > neigh1 || p1 > neigh2 )
            continue;

        p1_index = IJ(p1,0,3);
        n1_index = IJ(neigh1,0,3);
        n2_index = IJ(neigh2,0,3);

        x1 = parameters1[p1_index+0];
        y1 = parameters1[p1_index+1];
        z1 = parameters1[p1_index+2];
        x2 = parameters1[n1_index+0];
        y2 = parameters1[n1_index+1];
        z2 = parameters1[n1_index+2];
        x3 = parameters1[n2_index+0];
        y3 = parameters1[n2_index+1];
        z3 = parameters1[n2_index+2];

        a1 = parameters2[p1_index+0];
        b1 = parameters2[p1_index+1];
        c1 = parameters2[p1_index+2];
        a2 = parameters2[n1_index+0];
        b2 = parameters2[n1_index+1];
        c2 = parameters2[n1_index+2];
        a3 = parameters2[n2_index+0];
        b3 = parameters2[n2_index+1];
        c3 = parameters2[n2_index+2];

        dist1 = connections[p1].desired_distance;
        dist2 = connections[neigh1].desired_distance;
        dist3 = connections[neigh2].desired_distance;

        min_dist1a = connections[p1].min_distance1;
        min_dist1b = connections[p1].min_distance2;
        min_dist2a = connections[neigh1].min_distance1;
        min_dist2b = connections[neigh1].min_distance2;
        min_dist3a = connections[neigh2].min_distance1;
        min_dist3b = connections[neigh2].min_distance2;

        max_dist1a = connections[p1].max_distance1;
        max_dist1b = connections[p1].max_distance2;
        max_dist2a = connections[neigh1].max_distance1;
        max_dist2b = connections[neigh1].max_distance2;
        max_dist3a = connections[neigh2].max_distance1;
        max_dist3b = connections[neigh2].max_distance2;

        for_inclusive( w1, 0, oversample )
        {
            alpha1 = (Real) w1 / (Real) (oversample+1);
            weight2 = alpha1;

            for_inclusive( w2, 1, oversample - w1 + 1 )
            {
                if( w1 == 0 && w2 == oversample + 1 ||
                    w2 == oversample - w1 + 1 && neigh1 > neigh2 )
                {
                    continue;
                }

                alpha2 = (Real) w2 / (Real) (oversample - w1 + 1);

                weight1 = (1.0 - alpha1) * (1.0 - alpha2);
                weight3 = (1.0 - alpha1) * alpha2;

                x = weight1 * x1 + weight2 * x2 + weight3 * x3;
                y = weight1 * y1 + weight2 * y2 + weight3 * y3;
                z = weight1 * z1 + weight2 * z2 + weight3 * z3;

                a = weight1 * a1 + weight2 * a2 + weight3 * a3;
                b = weight1 * b1 + weight2 * b2 + weight3 * b3;
                c = weight1 * c1 + weight2 * c2 + weight3 * c3;

                dx = x - a;
                dy = y - b;
                dz = z - c;

                dist = dx * dx + dy * dy + dz * dz;
                if( dist <= 0.0 )
                    continue;

                desired_dist = weight1 * dist1 + weight2 * dist2 +
                               weight3 * dist3;
                min_dist1 = weight1 * min_dist1a + weight2 * min_dist2a +
                            weight3 * min_dist3a;
                min_dist2 = weight1 * min_dist1b + weight2 * min_dist2b +
                            weight3 * min_dist3b;
                max_dist1 = weight1 * max_dist1a + weight2 * max_dist2a +
                            weight3 * max_dist3a;
                max_dist2 = weight1 * max_dist1b + weight2 * max_dist2b +
                            weight3 * max_dist3b;

                diff = dist - desired_dist;

                factor = weight * diff * 2.0 / dist;

                deriv1[p1_index+0] += factor * dx * weight1;
                deriv1[p1_index+1] += factor * dy * weight1;
                deriv1[p1_index+2] += factor * dz * weight1;
                deriv1[n1_index+0] += factor * dx * weight2;
                deriv1[n1_index+1] += factor * dy * weight2;
                deriv1[n1_index+2] += factor * dz * weight2;
                deriv1[n2_index+0] += factor * dx * weight3;
                deriv1[n2_index+1] += factor * dy * weight3;
                deriv1[n2_index+2] += factor * dz * weight3;

                deriv2[p1_index+0] += factor * -dx * weight1;
                deriv2[p1_index+1] += factor * -dy * weight1;
                deriv2[p1_index+2] += factor * -dz * weight1;
                deriv2[n1_index+0] += factor * -dx * weight2;
                deriv2[n1_index+1] += factor * -dy * weight2;
                deriv2[n1_index+2] += factor * -dz * weight2;
                deriv2[n2_index+0] += factor * -dx * weight3;
                deriv2[n2_index+1] += factor * -dy * weight3;
                deriv2[n2_index+2] += factor * -dz * weight3;

                if( using_max_dist && min_dist1 < max_dist1 &&
                    (diff < min_dist1 || diff > max_dist1) )
                {
                    factor = 0.0;
                    if( diff < min_dist1 )
                    {
                        weight_factor = max_dist_weight1;
                        for_less( i, 0, n_weight_steps )
                        {
                            dist_threshold = INTERPOLATE( (Real) (i) /
                                                 (Real) (n_weight_steps-1),
                                                 min_dist1, min_dist2 );

                            if( diff >= dist_threshold )
                                break;

                            off = diff - dist_threshold;
                            factor += weight_factor * off * 2.0 / dist;
                            weight_factor *= scale_factor;
                        }
                    }
                    else
                    {
                        weight_factor = max_dist_weight1;
                        for_less( i, 0, n_weight_steps )
                        {
                            dist_threshold = INTERPOLATE( (Real) (i) /
                                                  (Real) (n_weight_steps-1),
                                                  max_dist1, max_dist2 );

                            if( diff <= dist_threshold )
                                break;

                            off = diff - dist_threshold;
                            factor += weight_factor * off * 2.0 / dist;
                            weight_factor *= scale_factor;
                        }
                    }

                    deriv1[p1_index+0] += factor * dx * weight1;
                    deriv1[p1_index+1] += factor * dy * weight1;
                    deriv1[p1_index+2] += factor * dz * weight1;
                    deriv1[n1_index+0] += factor * dx * weight2;
                    deriv1[n1_index+1] += factor * dy * weight2;
                    deriv1[n1_index+2] += factor * dz * weight2;
                    deriv1[n2_index+0] += factor * dx * weight3;
                    deriv1[n2_index+1] += factor * dy * weight3;
                    deriv1[n2_index+2] += factor * dz * weight3;

                    deriv2[p1_index+0] += factor * -dx * weight1;
                    deriv2[p1_index+1] += factor * -dy * weight1;
                    deriv2[p1_index+2] += factor * -dz * weight1;
                    deriv2[n1_index+0] += factor * -dx * weight2;
                    deriv2[n1_index+1] += factor * -dy * weight2;
                    deriv2[n1_index+2] += factor * -dz * weight2;
                    deriv2[n2_index+0] += factor * -dx * weight3;
                    deriv2[n2_index+1] += factor * -dy * weight3;
                    deriv2[n2_index+2] += factor * -dz * weight3;
                }
            }
        }
        }
    }
}

private  Real   evaluate_anchor_fit(
    Real                  adaptive_ratio,
    Real                  weight, 
    Real                  max_dist_weight,
    int                   n_anchor_points,
    anchor_point_struct   anchor_points[],
    Real                  parameters[],
    Smallest_int          active_flags[],
    Real                  max_weight_value )
{
    int                  anchor, p_index;
    Real                 fit1, fit2;
    Real                 dist;
    Real                 x1, y1, z1, x2, y2, z2, dx, dy, dz;
    anchor_point_struct  *a;
    Real                 adaptive_weight;

    if( weight < 0.0 )
        return( 0.0 );

    fit1 = 0.0;
    fit2 = 0.0;

    for_less( anchor, 0, n_anchor_points )
    {
        a = &anchor_points[anchor];

        if( active_flags != NULL && !active_flags[a->surface_point] )
            continue;

        p_index = IJ( a->surface_point, 0, 3 );
        x1 = parameters[p_index+0];
        y1 = parameters[p_index+1];
        z1 = parameters[p_index+2];

        x2 = RPoint_x( a->anchor_point );
        y2 = RPoint_y( a->anchor_point );
        z2 = RPoint_z( a->anchor_point );

        dx = x2 - x1;
        dy = y2 - y1;
        dz = z2 - z1;

        dist = dx * dx + dy * dy + dz * dz;
        if( dist >= 0.0 )
            dist = sqrt( dist );
        else
            dist = 1.0;

        dist -= a->desired_distance;

        // Added by June Sic Kim 7/12/2002
	if(adaptive_ratio==0){
	  fit1 += dist * dist;
	}
	else{
/*	  adaptive_weight = (a->weight / (max_weight_value/2));
	  fit1 += dist * dist * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_ANCHOR:1/ADAPTIVE_RATIO_ANCHOR);
*/        adaptive_weight = (a->weight / max_weight_value/2);
	  fit1 += dist * dist * adaptive_weight * adaptive_ratio;
	}

        if( a->min_distance < a->max_distance && max_dist_weight > 0.0 )
        {
            if( dist < a->min_distance )
            {
                dist = dist - a->min_distance;
		if(adaptive_ratio==0){
		  fit2 += dist * dist;
		}
		else{
/*		  adaptive_weight = (a->weight / (max_weight_value/2));
		  fit2 += dist * dist * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_ANCHOR:1/ADAPTIVE_RATIO_ANCHOR);
*/                adaptive_weight = (a->weight / max_weight_value/2);
		  fit2 += dist * dist * adaptive_weight * adaptive_ratio;
		}
            }
            else if( dist > a->max_distance )
            {
                dist = dist - a->max_distance;
		if(adaptive_ratio==0){
		  fit2 += dist * dist;
		}
		else{
/*		  adaptive_weight = (a->weight / (max_weight_value/2));
		  fit2 += dist * dist * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_ANCHOR:1/ADAPTIVE_RATIO_ANCHOR);
*/                adaptive_weight = (a->weight / max_weight_value/2);
		  fit2 += dist * dist * adaptive_weight * adaptive_ratio;
		}
            }
        }
    }

    return( fit1 * weight + fit2 * max_dist_weight );
}

private  void   evaluate_anchor_fit_deriv(
    Real                   adaptive_ratio,
    Real                   weight, 
    Real                   max_dist_weight,
    int                    n_anchor_points,
    anchor_point_struct    anchor_points[],
    Real                   parameters[],
    Real                   max_weight_value,
    Real                   deriv[] )
{
    int                  anchor, p_index;
    Real                 dist, diff, factor;
    Real                 x1, y1, z1, x2, y2, z2, dx, dy, dz;
    anchor_point_struct  *a;
    Real                 adaptive_weight;

    if( weight < 0.0 )
        return;

    for_less( anchor, 0, n_anchor_points )
    {
        a = &anchor_points[anchor];
        p_index = IJ( a->surface_point, 0, 3 );
        x1 = parameters[p_index+0];
        y1 = parameters[p_index+1];
        z1 = parameters[p_index+2];

        x2 = RPoint_x( a->anchor_point );
        y2 = RPoint_y( a->anchor_point );
        z2 = RPoint_z( a->anchor_point );

        dx = x2 - x1;
        dy = y2 - y1;
        dz = z2 - z1;

        dist = dx * dx + dy * dy + dz * dz;
        if( dist > 0.0 )
            dist = sqrt( dist );
        else
            dist = 1.0;

        diff = dist - a->desired_distance;

        // Added by June Sic Kim at 9/12/2002
	if(adaptive_ratio==0){
	  factor = weight * diff * 2.0 / dist;
	}
	else{
/*	  adaptive_weight = (a->weight / (max_weight_value/2));
	  factor = weight * diff * 2.0 / dist * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_ANCHOR:1/ADAPTIVE_RATIO_ANCHOR);
*/        adaptive_weight = (a->weight / max_weight_value);
	  factor = weight * diff * 2.0 / dist * adaptive_weight * adaptive_ratio;
	}

        deriv[p_index+0] += factor * -dx;
        deriv[p_index+1] += factor * -dy;
        deriv[p_index+2] += factor * -dz;

        if( a->min_distance < a->max_distance && max_dist_weight > 0.0 &&
            (diff < a->min_distance || diff > a->max_distance) )
        {
            if( diff < a->min_distance )
                diff = diff - a->min_distance;
            else
                diff = diff - a->max_distance;

            // Modified by June Sic Kim at 9/12/2002
	    if(adaptive_ratio==0){
	      factor = max_dist_weight * diff * 2.0 / dist;
	    }
	    else{
/*	      adaptive_weight = (a->weight / (max_weight_value/2));
	      factor = max_dist_weight * diff * 2.0 / dist * adaptive_weight * (adaptive_weight>=1?ADAPTIVE_RATIO_ANCHOR:1/ADAPTIVE_RATIO_ANCHOR);
*/            adaptive_weight = (a->weight / max_weight_value);
	      factor = max_dist_weight * diff * 2.0 / dist * adaptive_weight * adaptive_ratio;
	    }

            deriv[p_index+0] += factor * -dx;
            deriv[p_index+1] += factor * -dy;
            deriv[p_index+2] += factor * -dz;
        }
    }
}

private  Real   evaluate_weight_point_fit(
    Real                  weight, 
    Real                  max_dist_weight,
    int                   n_weight_points,
    weight_struct         weight_points[],
    Real                  parameters[],
    Smallest_int          active_flags[] )
{
    int                  w_index, p_index, p;
    Real                 fit1, fit2, this_weight;
    Real                 dist;
    Real                 x1, y1, z1, x2, y2, z2, dx, dy, dz;
    weight_struct        *w;

    if( weight < 0.0 )
        return( 0.0 );

    fit1 = 0.0;
    fit2 = 0.0;

    for_less( w_index, 0, n_weight_points )
    {
        w = &weight_points[w_index];

        if( active_flags != NULL )
        {
            for_less( p, 0, w->n_surface_points )
                if( !active_flags[w->surface_points[p]] )
                    continue;
        }

        x1 = 0.0;
        y1 = 0.0;
        z1 = 0.0;
        for_less( p, 0, w->n_surface_points )
        {
            p_index = IJ( w->surface_points[p], 0, 3 );
            this_weight = w->surface_weights[p];

            x1 += this_weight * parameters[p_index+0];
            y1 += this_weight * parameters[p_index+1];
            z1 += this_weight * parameters[p_index+2];
        }
            
        x2 = RPoint_x( w->anchor_point );
        y2 = RPoint_y( w->anchor_point );
        z2 = RPoint_z( w->anchor_point );

        dx = x1 - x2;
        dy = y1 - y2;
        dz = z1 - z2;

        dist = dx * dx + dy * dy + dz * dz;
        if( dist > 0.0 )
            dist = sqrt( dist );

        dist -= w->desired_distance;
        fit1 += dist * dist;

        if( w->min_distance <= w->max_distance && max_dist_weight > 0.0 )
        {
            if( dist < w->min_distance )
            {
                dist = dist - w->min_distance;
                fit2 += dist * dist;
            }
            else if( dist > w->max_distance )
            {
                dist = dist - w->max_distance;
                fit2 += dist * dist;
            }
        }
    }

    return( fit1 * weight + fit2 * max_dist_weight );
}

private  void   evaluate_weight_point_fit_deriv(
    Real                   weight, 
    Real                   max_dist_weight,
    int                    n_weight_points,
    weight_struct          weight_points[],
    Real                   parameters[],
    Real                   deriv[] )
{
    int                  w_index, p_index, p;
    Real                 this_weight;
    Real                 dist, diff, factor;
    Real                 x1, y1, z1, x2, y2, z2, dx, dy, dz;
    weight_struct        *w;

    if( weight < 0.0 )
        return;

    for_less( w_index, 0, n_weight_points )
    {
        w = &weight_points[w_index];

        x1 = 0.0;
        y1 = 0.0;
        z1 = 0.0;
        for_less( p, 0, w->n_surface_points )
        {
            p_index = IJ( w->surface_points[p], 0, 3 );
            this_weight = w->surface_weights[p];

            x1 += this_weight * parameters[p_index+0];
            y1 += this_weight * parameters[p_index+1];
            z1 += this_weight * parameters[p_index+2];
        }
            
        x2 = RPoint_x( w->anchor_point );
        y2 = RPoint_y( w->anchor_point );
        z2 = RPoint_z( w->anchor_point );

        dx = x1 - x2;
        dy = y1 - y2;
        dz = z1 - z2;

        dist = dx * dx + dy * dy + dz * dz;
        if( dist > 0.0 )
            dist = sqrt( dist );

        diff = dist - w->desired_distance;

        factor = weight * diff * 2.0 / dist;

        for_less( p, 0, w->n_surface_points )
        {
            p_index = IJ( w->surface_points[p], 0, 3 );
            this_weight = w->surface_weights[p];

            deriv[p_index+0] += this_weight * factor * dx;
            deriv[p_index+1] += this_weight * factor * dy;
            deriv[p_index+2] += this_weight * factor * dz;
        }

        if( w->min_distance <= w->max_distance && max_dist_weight > 0.0 &&
            (diff < w->min_distance || diff > w->max_distance) )
        {
            if( diff < w->min_distance )
                diff = diff - w->min_distance;
            else
                diff = diff - w->max_distance;

            factor = max_dist_weight * diff * 2.0 / dist;

            for_less( p, 0, w->n_surface_points )
            {
                p_index = IJ( w->surface_points[p], 0, 3 );
                this_weight = w->surface_weights[p];

                deriv[p_index+0] += this_weight * factor * dx;
                deriv[p_index+1] += this_weight * factor * dy;
                deriv[p_index+2] += this_weight * factor * dz;
            }
        }
    }
}

private  Real   evaluate_self_intersect_fit(
    int                           n_weights,
    Real                          weights[],
    Real                          min_distances[],
    self_intersect_lookup_struct  *si_lookup,
    BOOLEAN                       use_tri_tri_dist,
    BOOLEAN                       use_square_flag,
    Real                          dist_from_computed_self_intersect,
    int                           n_triangles,
    int                           triangles[],
    Real                          parameters[],
    Smallest_int                  active_flags[],
    Real                          *closest_dist ) {

    int     i, k, w, n_candidates;
    Real    fit, *fits=0, dist, dist_sq, diff, *sq_min_distances=0;
    Real    max_distance_sq;

    Real max_min_dist_sq = min_distances[0];
    for( w = 1; w < n_weights; w++ ) {
      if( min_distances[w] > max_min_dist_sq ) max_min_dist_sq = min_distances[w];
    }
    max_min_dist_sq *= max_min_dist_sq;

    Real * tri_distance;
    ALLOC( tri_distance, n_triangles );
    for( i = 0; i < n_triangles; i++ ) {
      tri_distance[i] = max_min_dist_sq;
    }

    n_candidates = get_n_self_intersect_candidate( si_lookup );

    for_less( i, 0, n_candidates ) {

      int poly1 = si_lookup->p1s[i];
      int p0 = triangles[3*poly1];
      int n01 = triangles[3*poly1+1];
      int n02 = triangles[3*poly1+2];
      int poly2 = si_lookup->p2s[i];
      int p1 = triangles[3*poly2];
      int n11 = triangles[3*poly2+1];
      int n12 = triangles[3*poly2+2];

      if( !test_self_intersect_candidate( p0, n01, n02, p1, n11, n12,
                                          &si_lookup->cases[i],
                                          si_lookup->min_line_dists[i],
                                          dist_from_computed_self_intersect,
                                          parameters,
                                          active_flags, &dist_sq ) ) {
        continue;
      }

      // Check for intersection. This should be a small tolerance to
      // avoid rounding errors.

      if( dist_sq == 0.0 ) {
        fit = 1.0e+30;
        break;
      }
      if( dist_sq < tri_distance[poly1] ) tri_distance[poly1] = dist_sq;
      if( dist_sq < tri_distance[poly2] ) tri_distance[poly2] = dist_sq;
    }

    if( i == n_candidates ) {
      ALLOC( fits, n_weights );
      for_less( w, 0, n_weights ) {
        fits[w] = 0.0;
      }

      *closest_dist = 1.0e+20;

      for( i = 0; i < n_triangles; i++ ) {

        if( tri_distance[i] >= max_min_dist_sq ) {
          tri_distance[i] = 0.0;                     // do not account in fit evaluation
          continue;
        } else {
          tri_distance[i] = sqrt( tri_distance[i] );
          if( tri_distance[i] < *closest_dist ) *closest_dist = tri_distance[i];
        }

        for( w = 0; w < n_weights; w++ ) {
          if( tri_distance[i] >= min_distances[w] ){
            continue;
          }

          diff = min_distances[w] - tri_distance[i];

#ifdef USE_CUBE
#define USE_CUBE
          fits[w] += diff * diff * diff / min_distances[w];
#else
          if( use_square_flag ) {
            fits[w] += diff * diff;
          } else {
            fits[w] += FABS( diff );
          }
#endif
        }
      }

      fit = 0.0;
      for_less( w, 0, n_weights ) {
        fit += weights[w] * fits[w];
      }

      FREE( fits );
    }

    FREE( tri_distance );

    return( fit );
}

private  void   evaluate_self_intersect_fit_deriv(
    int                n_weights,
    Real               weights[],
    Real               min_distances[],
    self_intersect_lookup_struct  *si_lookup,
    BOOLEAN            use_tri_tri_dist,
    BOOLEAN            use_square_flag,
    int                triangles[],
    Real               parameters[],
    Real               deriv[] ) {

    int        i, w, n_candidates;
    Real       dist, dist_sq, *sq_min_distances, max_distance_sq;
    BOOLEAN    deriv_done;
    Real       tri_tri_deriv[6][3];

    ALLOC( sq_min_distances, n_weights );

    max_distance_sq = 0.0;
    for_less( w, 0, n_weights ) {
        sq_min_distances[w] = min_distances[w] * min_distances[w];
        max_distance_sq = MAX( max_distance_sq, sq_min_distances[w] );
    }

    n_candidates = get_n_self_intersect_candidate( si_lookup );

    for_less( i, 0, n_candidates ) {

        int poly1 = si_lookup->p1s[i];
        int p0 = triangles[3*poly1];
        int n01 = triangles[3*poly1+1];
        int n02 = triangles[3*poly1+2];
        int poly2 = si_lookup->p2s[i];
        int p1 = triangles[3*poly2];
        int n11 = triangles[3*poly2+1];
        int n12 = triangles[3*poly2+2];

        dist_sq = si_lookup->min_line_dists[i];

        deriv_done = FALSE;
        Real factor = 0.0;

        for_less( w, 0, n_weights ) {
            if( dist_sq >= sq_min_distances[w] )
                continue;

            if( !deriv_done ) {
                deriv_done = TRUE;
                dist = sqrt( dist_sq );
                sq_triangle_triangle_dist_deriv( &parameters[3*p0],
                                                 &parameters[3*n01],
                                                 &parameters[3*n02],
                                                 &parameters[3*p1],
                                                 &parameters[3*n11],
                                                 &parameters[3*n12],
                                                 tri_tri_deriv[0],
                                                 tri_tri_deriv[1],
                                                 tri_tri_deriv[2],
                                                 tri_tri_deriv[3],
                                                 tri_tri_deriv[4],
                                                 tri_tri_deriv[5],
                                                 si_lookup->cases[i] );
            }

            factor += get_self_intersect_deriv_factor( 
                                      use_square_flag, min_distances[w],
                                      weights[w], dist );
        }

        if( deriv_done ) {
          get_self_intersect_deriv( p0, n01, n02, p1, n11, n12, factor,
                                    deriv, tri_tri_deriv );
        }

    }

    FREE( sq_min_distances );
}

typedef struct
{
    Real  min_distance;
    Real  weight;
    Real  a;
    Real  b;
    Real  c;
} piecewise_weight_struct;

private  Real   evaluate_surf_surf_fit(
    int                           n_weights,
    Real                          weights[],
    Real                          min_distances[],
    surf_surf_lookup_struct       *ss_lookup,
    Real                          dist_from_computed_self_intersect,
    int                           triangle1[],
    Real                          parameters1[],
    Smallest_int                  active_flags1[],
    int                           triangle2[],
    Real                          parameters2[],
    Smallest_int                  active_flags2[],
    Real                          *closest_dist )
{
    int     i, w, w1, n_candidates, best;
    int     p1, n1, p2, n2, n11, n12, n21, n22;
    Real    fit, dist, dist_sq, sq_min_dist;
    Real    max_distance_sq, closest;
    Real    *parameters_p1, *parameters_n11, *parameters_n12;
    int     prev_p1s, p1s, p2s;
    piecewise_weight_struct  *sorted, tmp;
    Real    intersect_dist, tmp_dist;
    int     intersect_num = 0;

    closest = -1.0;

    max_distance_sq = 0.0;
    for_less( w, 0, n_weights ) {
        sq_min_dist = min_distances[w] * min_distances[w];
        max_distance_sq = MAX( max_distance_sq, sq_min_dist );
    }

    ALLOC( sorted, n_weights );

    for_less( w, 0, n_weights ) {
        sorted[w].min_distance = min_distances[w];
        sorted[w].weight = weights[w];
    }

    for_less( w, 0, n_weights - 1 ) {
        best = w;
        for_less( w1, w+1, n_weights ) {
            if( sorted[w1].min_distance < sorted[best].min_distance )
                best = w1;
        }

        tmp = sorted[best];
        sorted[best] = sorted[w];
        sorted[w] = tmp;
    }

    for_down( w, n_weights-1, 0 ) {
        if( w == n_weights-1 ) {
            sorted[w].a = 0.0;
            sorted[w].b = 0.0;
            sorted[w].c = 0.0;
        } else {
            sorted[w].a = sorted[w+1].a;
            sorted[w].b = sorted[w+1].b;
            sorted[w].c = sorted[w+1].c;
        }

        sorted[w].a += sorted[w].weight;
        sorted[w].b += -2.0 * sorted[w].weight * sorted[w].min_distance;
        sorted[w].c += sorted[w].weight * sorted[w].min_distance *
                       sorted[w].min_distance;
    }

    n_candidates = get_n_surf_surf_candidate( ss_lookup );

    prev_p1s = -1;
    fit = 0.0;

    for_less( i, 0, n_candidates ) {
        if( dist_from_computed_self_intersect > 0.0 &&
            (Real) ss_lookup->min_line_dists[i] >
            dist_from_computed_self_intersect ) {
            continue;
        }

        p1s = ss_lookup->p1s[i];

        if( p1s != prev_p1s ) {
            prev_p1s = p1s;
            p1 = triangle1[3*p1s];
            n11 = triangle1[3*p1s+1];
            n12 = triangle1[3*p1s+2];
            parameters_p1 = &parameters1[p1*3];
            parameters_n11 = &parameters1[n11*3];
            parameters_n12 = &parameters1[n12*3];
        }

        p2s = ss_lookup->p2s[i];
        p2 = triangle2[3*p2s];
        n21 = triangle2[3*p2s+1];
        n22 = triangle2[3*p2s+2];

        if( active_flags1 != NULL &&
            !active_flags1[p1] && !active_flags1[n11] && !active_flags1[n12] &&
            !active_flags2[p2] && !active_flags2[n21] && !active_flags2[n22] ) {
            continue;
        }

        dist_sq = sq_triangle_triangle_dist( parameters_p1,
                                             parameters_n11,
                                             parameters_n12,
                                             &parameters2[p2*3],
                                             &parameters2[n21*3],
                                             &parameters2[n22*3],
                                             &ss_lookup->cases[i] );

        if( dist_sq > max_distance_sq )
            continue;

        if( closest < 0.0 || dist_sq < closest )
            closest = dist_sq;

        if( dist_sq > 0.0 ) {
            dist = sqrt( dist_sq );
            w = 0;
            while( dist > sorted[w].min_distance && w < n_weights-1 ) {
                ++w;
            }
            fit += sorted[w].c + dist * (sorted[w].b + dist * sorted[w].a);
        } else {
          //fit += 1e-4;//2e-3;
          intersect_dist = 0;
          tmp_dist = sq_triangle_point_dist(&parameters2[p2*3],
                                            &parameters2[n21*3],
                                            &parameters2[n22*3], 
                                            parameters_p1 );
          intersect_dist = (tmp_dist>intersect_dist)?tmp_dist:intersect_dist;
          /*tmp_dist = sq_triangle_point_dist(&parameters2[p2*3],
                                            &parameters2[n21*3],
                                            &parameters2[n22*3], 
                                            parameters_n11 );
          intersect_dist = (tmp_dist>intersect_dist)?tmp_dist:intersect_dist;
          tmp_dist = sq_triangle_point_dist(&parameters2[p2*3],
                                            &parameters2[n21*3],
                                            &parameters2[n22*3], 
                                            parameters_n12 );
          intersect_dist = (tmp_dist>intersect_dist)?tmp_dist:intersect_dist;
          if(intersect_dist!=0)printf("%g\t", intersect_dist);*/
          fit += intersect_dist*1e0;
          intersect_num++;
        }
    }

    FREE( sorted );
    //printf("intersect num = %d\n",intersect_num);
    if( closest > 0.0 )
        closest = sqrt( closest );

    *closest_dist = closest;

    return( fit );
}

private  void   evaluate_surf_surf_fit_deriv(
    int                      n_weights,
    Real                     weights[],
    Real                     min_distances[],
    surf_surf_lookup_struct  *ss_lookup,
    int                      n_points1,
    int                      triangle1[],
    Real                     parameters1[],
    Real                     deriv1[],
    int                      n_points2,
    int                      triangle2[],
    Real                     parameters2[],
    Real                     deriv2[] ) {

    int                     i, w, n_candidates;
    Real                    dist, dist_sq, *sq_min_distances, max_distance_sq;
    BOOLEAN                 deriv_done;
    Real                    tri_tri_deriv[6][3];

    ALLOC( sq_min_distances, n_weights );

    max_distance_sq = 0.0;
    for_less( w, 0, n_weights ) {
        sq_min_distances[w] = min_distances[w] * min_distances[w];
        max_distance_sq = MAX( max_distance_sq, sq_min_distances[w] );
    }

    n_candidates = get_n_surf_surf_candidate( ss_lookup );

    for_less( i, 0, n_candidates ) {

        int p1s = ss_lookup->p1s[i];
        int p1 = triangle1[3*p1s];
        int n11 = triangle1[3*p1s+1];
        int n12 = triangle1[3*p1s+2];

        int p2s = ss_lookup->p2s[i];
        int p2 = triangle2[3*p2s];
        int n21 = triangle2[3*p2s+1];
        int n22 = triangle2[3*p2s+2];

        dist_sq = sq_triangle_triangle_dist( &parameters1[p1*3],
                                             &parameters1[n11*3],
                                             &parameters1[n12*3],
                                             &parameters2[p2*3],
                                             &parameters2[n21*3],
                                             &parameters2[n22*3],
                                             &ss_lookup->cases[i] );

        deriv_done = FALSE;
        Real factor = 0.0;

        for_less( w, 0, n_weights ) {
            if( dist_sq >= sq_min_distances[w] )
                continue;

            if( !deriv_done ) {
                deriv_done = TRUE;

                if( dist_sq > 0.0 ) {
                  dist = sqrt( dist_sq );
                } else {
                  dist = 0.0;
                }

                dist = sqrt( dist_sq );
                sq_triangle_triangle_dist_deriv( &parameters1[3*p1],
                                                 &parameters1[3*n11],
                                                 &parameters1[3*n12],
                                                 &parameters2[3*p2],
                                                 &parameters2[3*n21],
                                                 &parameters2[3*n22],
                                                 tri_tri_deriv[0],
                                                 tri_tri_deriv[1],
                                                 tri_tri_deriv[2],
                                                 tri_tri_deriv[3],
                                                 tri_tri_deriv[4],
                                                 tri_tri_deriv[5],
                                                 ss_lookup->cases[i] );
            }

            factor += get_surf_surf_deriv_factor( min_distances[w], weights[w],
                                                  dist );
        }

        if( deriv_done ) {
          get_surf_surf_deriv( p1, n11, n12, p2, n21, n22, factor,
                               deriv1, deriv2, tri_tri_deriv );
        }


    }
    memset(deriv2, 0, sizeof(Real)*n_points2*3);

    FREE( sq_min_distances );
}

public  void   compute_boundary_line_coefficients(
    int                           n_parameters,
    Real                          parameters[],
    Real                          constant,
    Real                          linear[],
    Real                          square[],
    int                           n_cross_terms[],
    int                           *cross_parms[],
    Real                          *cross_terms[],
    Real                          line_dir[],
    Real                          coefs[] )
{
    int                    p, n, c;
    Real                   weight;

    coefs[0] = 0.0;
    coefs[1] = 0.0;
    coefs[2] = 0.0;

    if( linear == NULL )
        return;

    coefs[0] = constant;

    for_less( p, 0, n_parameters )
    {
        /*--- linear term */

        coefs[0] += linear[p] * parameters[p];
        coefs[1] += linear[p] * line_dir[p];

        /*--- square term */

        coefs[0] += square[p] * parameters[p] * parameters[p];
        coefs[1] += square[p] * 2.0 * parameters[p] * line_dir[p];
        coefs[2] += square[p] * line_dir[p] * line_dir[p];

        for_less( c, 0, n_cross_terms[p] )
        {
            n = cross_parms[p][c];
            weight = cross_terms[p][c];
            coefs[0] += weight * parameters[p] * parameters[n];
            coefs[1] += weight * (parameters[p] * line_dir[n] +
                                  parameters[n] * line_dir[p]);
            coefs[2] += weight * line_dir[p] * line_dir[n];
        }
    }
}

private  int   private_evaluate_fit(
    int                           which,
    Deform_struct                 *deform,
    int                           start_parameter[],
    Real                          parameters[],
    Smallest_int                  active_flags[],
    Smallest_int                  evaluate_flags[],
    Real                          boundary_coefs[],
    Real                          t_dist,
    Smallest_int                  boundary_flags[],
    Point                         boundary_points[],
    Real                          max_value,
    Real                          dist_from_computed_self_intersect,
    self_intersect_lookup_struct  **si_lookup,
    surf_surf_lookup_struct       *ss_lookup,
    Real                          *fit,
    fit_eval_struct               *eval )
{
    int                    i, surface, ind, count;
    Real                   f, fv, closest;
    surface_bound_struct   *bound;
    stretch_struct         *stretch;
    curvature_struct       *curvature;
    bend_struct            *bend;
    gradient_struct        *gradient;
    gw_gradient_struct     *gw_gradient;
    surface_value_struct   *value;
    inter_surface_struct   *inter;
    self_intersect_struct  *self;
    surf_surf_struct       *surf;
    anchor_struct          *anchor;
    weight_point_struct    *weight_point;
    Real                   *this_parms;
    Smallest_int           *this_active, *this_evaluate;

    count = 0;

    for_less( i, 0, deform->n_inter_surfaces ) {
        inter = &deform->inter_surfaces[i];

        if( which == -2 || which == count ) {
            f = evaluate_inter_surface_fit(
               inter->weight,
               inter->max_weight1,
               inter->max_weight2,
               inter->n_weight_steps,
               inter->n_connections, inter->connections,
               &parameters[start_parameter[inter->surface_index1]],
               (active_flags == NULL) ? NULL :
                  &active_flags[start_parameter[inter->surface_index1]/3],
               &parameters[start_parameter[inter->surface_index2]],
               (active_flags == NULL) ? NULL :
                  &active_flags[start_parameter[inter->surface_index2]/3] );

            eval->inter_surface_fit += f;
            *fit += f;

            if( which == -2 && max_value > 0.0 && *fit > max_value )
                return( 0 );
        }
        ++count;

        if( inter->oversample > 0 ) {
            if( which == -2 || which == count ) {
                f = evaluate_oversampled_inter_surface_fit(
                   inter->weight,
                   inter->max_weight1,
                   inter->max_weight2,
                   inter->n_weight_steps,
                   inter->n_connections, inter->connections,
                   inter->oversample,
                   deform->surfaces[inter->surface_index1].surface.n_points,
                   deform->surfaces[inter->surface_index1].surface.n_neighbours,
                   deform->surfaces[inter->surface_index1].surface.neighbours,
                   &parameters[start_parameter[inter->surface_index1]],
                   (active_flags == NULL) ? NULL :
                      &active_flags[start_parameter[inter->surface_index1]/3],
                   &parameters[start_parameter[inter->surface_index2]],
                   (active_flags == NULL) ? NULL :
                      &active_flags[start_parameter[inter->surface_index2]/3] );

                eval->inter_surface_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                    return( 0 );
            }
            ++count;
        }
    }

    for_less( i, 0, deform->n_surf_surfs ) {
        surf = &deform->surf_surfs[i];

        if( which == -2 || which == count ) {
            f = evaluate_surf_surf_fit(
                surf->n_weights, surf->weights, surf->min_distances,
                &ss_lookup[i],
                dist_from_computed_self_intersect,
                deform->surfaces[surf->surface_index1].surface.triangles,
                &parameters[start_parameter[surf->surface_index1]],
                (active_flags == NULL) ? NULL :
                      &active_flags[start_parameter[surf->surface_index1]/3],
                deform->surfaces[surf->surface_index2].surface.triangles,
                &parameters[start_parameter[surf->surface_index2]],
                (active_flags == NULL) ? NULL :
                      &active_flags[start_parameter[surf->surface_index2]/3],
                &closest );

            if( closest >= 0.0 && (closest < eval->closest_surf_surf ||
                                   eval->closest_surf_surf < 0.0 ) )
                eval->closest_surf_surf = closest;

            eval->surf_surf_fit += f;
            *fit += f;

            if( which == -2 && max_value > 0.0 && *fit > max_value )
                return( 0 );
        }
        ++count;
    }

    for_less( surface, 0, deform->n_surfaces>0?1:0 ) {

        this_parms = &parameters[start_parameter[surface]];
        this_active = (active_flags == NULL) ? NULL :
                         &active_flags[start_parameter[surface]/N_DIMENSIONS];

        for_less( i, 0, deform->surfaces[surface].n_self_intersects ) {
            self = &deform->surfaces[surface].self_intersects[i];
            if( which == -2 || which == count ) {

                f = evaluate_self_intersect_fit(
                    self->n_weights, self->weights,
                    self->min_distances,
                    &si_lookup[surface][i],
                    self->use_tri_tri_dist,
                    self->square_flag,
                    dist_from_computed_self_intersect,
                    deform->surfaces[surface].surface.n_polygons,
                    deform->surfaces[surface].surface.triangles,
                    this_parms, this_active, &closest );

                if( closest >= 0.0 && (closest < eval->closest_self_intersect ||
                                   eval->closest_self_intersect < 0.0 ) )
                    eval->closest_self_intersect = closest;

                eval->self_intersect_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                    return( 0 );
            }
            ++count;
        }

        //////////////////////////////////////////////////////////////////
        // prevent from intersecting with the WM surface
        for_less( i, 0, deform->n_intersect_wm ) {
          if( which == -2 || which == count ) {
            f = evaluate_intersect_wm_fit( deform->intersect_wm->objects,
                           deform->intersect_wm->weight, this_parms,
                           deform->intersect_wm->intersected,
                           deform->intersect_wm->n_intersected,
                           0, deform->surfaces[surface].surface.n_points);
            eval->surf_surf_fit += f;
            *fit += f;
	  }
	  ++count;
        }

        //////////////////////////////////////////////////////////////////
        // LAPLACIAN constraint
        for_less( i, 0, deform->surfaces[surface].n_laplacian ) {
          //anchor = &deform->surfaces[surface].anchors[0];
          if( which == -2 || which == count ) {
            f = evaluate_laplacian_fit(
                           deform->surfaces[surface].laplacian->weight,
                           deform->surfaces[surface].laplacian->volume,
                           deform->surfaces[surface].laplacian->to_value,
                           deform->surfaces[surface].laplacian->type,
                           deform->surfaces[surface].surface.n_neighbours,
                           deform->surfaces[surface].surface.neighbours,
                           deform->surfaces[surface].laplacian->oversample,
                           this_parms,
                           0,
                           deform->surfaces[surface].surface.n_points);
            eval->laplacian_fit += f;
            *fit += f;
          }
          ++count;
        }

        for_less( i, 0, deform->surfaces[surface].n_volume )
        {
            if( which == -2 || which == count ) {
              //anchor = &deform->surfaces[surface].anchors[0];
              //bound = &deform->surfaces[surface].bound[0];
              bound = &deform->surfaces[surface].bound[i];   // modified by Claude
              fv = evaluate_volume_fit( deform->surfaces[surface].volume->weight,
                             deform->surfaces[surface].volume->max_weight,
                             deform->surfaces[surface].volume->direction,
                             bound->volume,
                             deform->surfaces[surface].laplacian->volume,
                             bound->threshold,
                             this_parms,
                             0,
                             deform->surfaces[surface].surface.n_points,
                             deform->surfaces[surface].surface.n_points);
              eval->volume_fit += fv;
              *fit += fv;
              if( which == -2 && max_value > 0.0 && *fit > max_value )
                  return( 0 );
            }
            count++;
        }          

        for_less( i, 0, deform->surfaces[surface].n_anchors )
        {
            if( which == -2 || which == count ) {
              /////////////////////////////////////////////////////////
              /////// Modified by June ////////////////////////////////
              if( NO_ANCHOR != 1 ) {
                anchor = &deform->surfaces[surface].anchors[i];
                f = evaluate_anchor_fit( 0,
                       anchor->weight*1, anchor->max_dist_weight*1,
                       anchor->n_anchor_points, anchor->anchor_points,
                       this_parms, this_active, anchor->max_weight_value );
        
                eval->anchor_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                  return( 0 );
              }
            }
            ++count;
        }

        for_less( i, 0, deform->surfaces[surface].n_weight_points ) {

            if( which == -2 || which == count ) {
                weight_point = &deform->surfaces[surface].weight_points[i];
                f = evaluate_weight_point_fit(
                       weight_point->weight, weight_point->max_dist_weight,
                       weight_point->n_weight_points,
                       weight_point->weight_points,
                       this_parms, this_active );

                eval->weight_point_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                    return( 0 );
            }

            ++count;
        }
    }

    ind = 0;
    for_less( surface, 0, deform->n_surfaces>0?1:0 ) {
        this_parms = &parameters[start_parameter[surface]];
        this_evaluate = (evaluate_flags == NULL) ? NULL :
                         &evaluate_flags[start_parameter[surface]/N_DIMENSIONS];
        this_active = (active_flags == NULL) ? NULL :
                         &active_flags[start_parameter[surface]/N_DIMENSIONS];

        if( which == -2 || which == count ) {
            f = boundary_coefs[0] + t_dist *
                      (boundary_coefs[1] + t_dist * boundary_coefs[2]);
            //////////////////////////////////////////////////////////
            // Modified by June at 22/11
            // reduce boundary_search_fit term to 0.5 times
            if(BOUNDARY_DECREASE == 1){
	      eval->boundary_fit += (f*1.0);
	      *fit += (f*1.0);
            }
            //////////////////////////////////////////////////////////
        }

        ++count;

        for_less( i, 0, deform->surfaces[surface].n_bound ) {
            bound = &deform->surfaces[surface].bound[i];

            if( bound->max_dist_weight > 0.0 ) {
                if( which == -2 || which == count ) {
                    f = evaluate_boundary_search_fit( bound->image_weight_in,
                                          bound->image_weight_out,
                                          bound->max_inward,
                                          bound->max_outward,
                                          bound->max_dist_threshold,
                                          bound->max_dist_weight,
                                          bound->oversample,
                                          &boundary_flags[ind],
                                          &boundary_points[ind],
                                          0,
                       deform->surfaces[surface].surface.n_points,
                       deform->surfaces[surface].surface.n_points,
                       this_parms, this_evaluate,
                       deform->surfaces[surface].surface.n_neighbours,
                       deform->surfaces[surface].surface.neighbours, 
                       bound->weights,
                       bound->max_weight_value,
                       deform->surfaces[surface].volume->adaptive_boundary_ratio);

                    ind += deform->surfaces[surface].surface.n_points +
                           bound->oversample * deform->surfaces[surface].surface.n_edges +
                           ( bound->oversample * ( bound->oversample - 1 ) ) / 2 *
                           deform->surfaces[surface].surface.n_polygons;

                    eval->boundary_fit += f;
                    *fit += f;

                    if( which == -2 && max_value > 0.0 && *fit > max_value )
                        return( 0 );
                }
                ++count;
            }
        }

        for_less( i, 0, deform->surfaces[surface].n_gradient ) {
            gradient = &deform->surfaces[surface].gradient[i];

            if( which == -2 || which == count ) {
                f = evaluate_gradient_fit( gradient->volume,
                                          gradient->voxel_lookup,
                                          gradient->continuity,
                                          gradient->image_weight,
                                          gradient->threshold,
                                          gradient->min_diff,
                                          gradient->max_diff,
                                          gradient->max_diff_weight,
                                          gradient->differential_offset,
                                          gradient->differential_ratio,
                   deform->surfaces[surface].surface.n_points,
                   this_parms, this_active);

                eval->gradient_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                    return( 0 );
            }
            ++count;
        }

        for_less( i, 0, deform->surfaces[surface].n_gw_gradient ) {
            gw_gradient = &deform->surfaces[surface].gw_gradient[i];

            if( which == -2 || which == count ) {
                f = evaluate_gw_gradient_fit( gw_gradient->t1,
                                 gw_gradient->image_weight, 
                                 gw_gradient->oversample,
                                 gw_gradient->mask, 0,
                                 deform->surfaces[surface].surface.n_points,
                                 deform->surfaces[surface].surface.n_neighbours,
                                 deform->surfaces[surface].surface.neighbours,
                                 this_parms, gw_gradient->t1grad );

                eval->gw_gradient_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                    return( 0 );
            }
            ++count;
        }

        for_less( i, 0, deform->surfaces[surface].n_value ) {

            if( which == -2 || which == count ) {
                value = &deform->surfaces[surface].value[i];
                f = evaluate_image_value_fit( value->volume,
                                          value->voxel_lookup,
                                          value->continuity,
                                          value->image_weight,
                                          value->threshold,
                                          value->min_diff,
                                          value->max_diff,
                                          value->max_diff_weight,
                                          value->differential_offset,
                                          value->differential_ratio,
                                          value->oversample,
                   deform->surfaces[surface].surface.n_points,
                   deform->surfaces[surface].surface.n_neighbours,
                   deform->surfaces[surface].surface.neighbours,
                   this_parms, this_active);

                eval->value_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                    return( 0 );
            }
            ++count;
        }

        for_less( i, 0, deform->surfaces[surface].n_stretch ) {

            if( which == -2 || which == count ) {
                stretch = &deform->surfaces[surface].stretch[i];
                f = evaluate_stretch_fit( stretch->stretch_weight,
                                          stretch->n_points,
                                          this_parms,
                                          stretch->n_neighbours,
                                          stretch->neighbours );

                eval->stretch_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                    return( 0 );
            }
            ++count;
        }

        for_less( i, 0, deform->surfaces[surface].n_curvature ) {

            if( which == -2 || which == count ) {
                curvature = &deform->surfaces[surface].curvature[i];
                f = evaluate_curvature_fit( curvature->curvature_weight,
                                        curvature->max_curvature_weight,
                                        curvature->min_curvature,
                                        curvature->max_curvature,
                                        0,
                                        curvature->n_points,
                                        this_parms, this_evaluate,
                                        curvature->n_neighbours,
                                        curvature->neighbours,
                                        curvature->curvatures );

                eval->curvature_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                    return( 0 );
            }
            ++count;
        }

        for_less( i, 0, deform->surfaces[surface].n_bend ) {

            if( which == -2 || which == count ) {
                bend = &deform->surfaces[surface].bend[i];
                f = evaluate_bend_fit( bend->bend_weight,
                                       bend->max_bend_weight,
                                       bend->min_bend,
                                       bend->max_bend,
                                       bend->n_points,
                                       this_parms, this_active,
                                       bend->n_neighbours,
                                       bend->neighbours,
                                       bend->model_x, bend->model_y );

                eval->bend_fit += f;
                *fit += f;

                if( which == -2 && max_value > 0.0 && *fit > max_value )
                    return( 0 );
            }
            ++count;
        }
    }

    return( count );
}

public  Real   evaluate_fit(
    Deform_struct                 *deform,
    int                           start_parameter[],
    Real                          parameters[],
    Smallest_int                  active_flags[],
    Smallest_int                  evaluate_flags[],
    Real                          boundary_coefs[],
    Real                          t_dist,
    Smallest_int                  boundary_flags[],
    Point                         boundary_points[],
    Real                          max_value,
    Real                          dist_from_computed_self_intersect,
    self_intersect_lookup_struct  **si_lookup,
    surf_surf_lookup_struct       *ss_lookup,
    fit_eval_struct               *fit_info )
{
    Real                   fit;
    fit_eval_struct        eval;

    eval.boundary_fit = 0.0;
    eval.gradient_fit = 0.0;
    eval.gw_gradient_fit = 0.0;
    eval.value_fit = 0.0;
    eval.stretch_fit = 0.0;
    eval.curvature_fit = 0.0;
    eval.bend_fit = 0.0;
    eval.self_intersect_fit = 0.0;
    eval.closest_self_intersect = -1.0;
    eval.surf_surf_fit = 0.0;
    eval.closest_surf_surf = -1.0;
    eval.inter_surface_fit = 0.0;
    eval.anchor_fit = 0.0;
    eval.weight_point_fit = 0.0;
    eval.volume_fit = 0.0;
    eval.laplacian_fit = 0.0;

    fit = 0.0;

    (void) private_evaluate_fit( -2, deform, start_parameter, parameters,
                                 active_flags, evaluate_flags,
                                 boundary_coefs, t_dist,
                                 boundary_flags, boundary_points,
                                 max_value,
                                 dist_from_computed_self_intersect,
                                 si_lookup, ss_lookup,
                                 &fit, &eval );

    if( fit_info != NULL )
        *fit_info = eval;

    return( fit );
}

private  Real  compute_deriv_mag(
    int           n_parameters,
    Real          deriv[] )
{
    int     i;
    Real    sum;

    sum = 0.0;
    for_less( i, 0, n_parameters )
        sum += deriv[i] * deriv[i];

    if( sum > 0.0 )
        sum = sqrt( sum );

    return( sum );
}

private  void  remove_self_intersect_components(
    Real          remove_ratio,
    int           n_parameters,
    Real          deriv[],
    Real          si_deriv[] )
{
    int    p, p_index;
    Real   si_dot_si, d_dot_si, amount;

    for_less( p, 0, n_parameters / 3 )
    {
        p_index = IJ( p, 0, 3 );

        deriv[p_index+0] += si_deriv[p_index+0];
        deriv[p_index+1] += si_deriv[p_index+1];
        deriv[p_index+2] += si_deriv[p_index+2];

        d_dot_si = deriv[p_index+0] * si_deriv[p_index+0] +
                   deriv[p_index+1] * si_deriv[p_index+1] +
                   deriv[p_index+2] * si_deriv[p_index+2];
        si_dot_si = si_deriv[p_index+0] * si_deriv[p_index+0] +
                    si_deriv[p_index+1] * si_deriv[p_index+1] +
                    si_deriv[p_index+2] * si_deriv[p_index+2];

        if( d_dot_si < 0.0 && si_dot_si > 0.0 )
        {
            amount = remove_ratio * (d_dot_si / si_dot_si);

            deriv[p_index+0] -= amount * si_deriv[p_index+0];
            deriv[p_index+1] -= amount * si_deriv[p_index+1];
            deriv[p_index+2] -= amount * si_deriv[p_index+2];
        }
    }
}

public  void   evaluate_fit_deriv(
    Deform_struct                 *deform,
    int                           n_parameters,
    int                           start_parameter[],
    Real                          parameters[],
    Smallest_int                  boundary_flags[],
    Point                         boundary_points[],
    Real                          linear[],
    Real                          square[],
    int                           n_cross_terms[],
    int                           *cross_parms[],
    Real                          *cross_terms[],
    self_intersect_lookup_struct  **si_lookup,
    surf_surf_lookup_struct       *ss_lookup,
    Real                          full_deriv[],
    fit_eval_struct               *fit_info )
{
    int                    p, i, surface, ind;
    stretch_struct         *stretch;
    curvature_struct       *curvature;
    bend_struct            *bend;
    surface_bound_struct   *bound;
    gradient_struct        *gradient;
    gw_gradient_struct     *gw_gradient;
    surface_value_struct   *value;
    inter_surface_struct   *inter;
    self_intersect_struct  *self;
    surf_surf_struct       *surf;
    anchor_struct          *anchor;
    weight_point_struct    *weight_point;
    Real                   *this_parms, *this_deriv;
    Real                   *derivative;
    Real                   remove_ratio;
    BOOLEAN                deriv_modified;
    fit_eval_struct        eval;

    eval.boundary_fit = 0.0;
    eval.gradient_fit = 0.0;
    eval.gw_gradient_fit = 0.0;
    eval.value_fit = 0.0;
    eval.stretch_fit = 0.0;
    eval.curvature_fit = 0.0;
    eval.bend_fit = 0.0;
    eval.self_intersect_fit = 0.0;
    eval.closest_self_intersect = -1.0;
    eval.surf_surf_fit = 0.0;
    eval.closest_surf_surf = -1.0;
    eval.inter_surface_fit = 0.0;
    eval.anchor_fit = 0.0;
    eval.weight_point_fit = 0.0;
    eval.volume_fit = 0.0;
    eval.laplacian_fit = 0.0;

    ALLOC( derivative, n_parameters );

    deriv_modified = FALSE;
    for_less( p, 0, n_parameters ) {
      full_deriv[p] = 0.0;
      derivative[p] = 0.0;
    }

    for_less( surface, 0, deform->n_surfaces ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_gradient ) {
            deriv_modified = TRUE;
            gradient = &deform->surfaces[surface].gradient[i];

            evaluate_gradient_fit_deriv( gradient->volume,
                                         gradient->voxel_lookup,
                                         gradient->continuity,
                                         gradient->image_weight,
                                         gradient->threshold,
                                         gradient->differential_ratio,
                   deform->surfaces[surface].surface.n_points,
                   this_parms, this_deriv );
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.gradient_fit = compute_deriv_mag( n_parameters, derivative );
      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    for_less( surface, 0, deform->n_surfaces ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_gw_gradient ) {
            deriv_modified = TRUE;
            gw_gradient = &deform->surfaces[surface].gw_gradient[i];
            if( gw_gradient->t1grad == NULL ) {
              gw_gradient->t1grad = 
                (Real*)malloc( deform->surfaces[surface].surface.n_points *
                               sizeof( Real ) );
            }
            if( gw_gradient->mask == NULL ) {
              gw_gradient->mask = 
                (int*)malloc( deform->surfaces[surface].surface.n_points *
                               sizeof( int ) );
              FILE * fp = fopen( gw_gradient->mask_filename, "r+t" );
              int point;
              for( point = 0; point < deform->surfaces[surface].surface.n_points; point++ ) {
                fscanf( fp, "%d", &gw_gradient->mask[point] );
              }
              fclose( fp );
	    }

            if( gw_gradient->update == 0 ) {
              adjust_gm_wm_gradient( gw_gradient->t1, gw_gradient->cls, 
                                     gw_gradient->search_distance, 
                                     gw_gradient->search_increment,
                                     0,
                                     deform->surfaces[surface].surface.n_points,
                                     deform->surfaces[surface].surface.n_neighbours,
                                     deform->surfaces[surface].surface.neighbours,
                                     this_parms,
                                     gw_gradient->mask,
                                     gw_gradient->t1grad,
                                     &gw_gradient->image_type );
            }
            gw_gradient->update = ( gw_gradient->update + 1 ) % 2;
            evaluate_gw_gradient_fit_deriv( gw_gradient->t1,
                                   gw_gradient->image_weight,
                                   gw_gradient->oversample,
                                   gw_gradient->image_type,
                                   gw_gradient->mask, 0,
                                   deform->surfaces[surface].surface.n_points,
                                   deform->surfaces[surface].surface.n_neighbours,
                                   deform->surfaces[surface].surface.neighbours,
                                   this_parms, gw_gradient->t1grad,
                                   this_deriv );
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.gw_gradient_fit = compute_deriv_mag( n_parameters, derivative );
      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    for_less( surface, 0, deform->n_surfaces ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_value ) {
            deriv_modified = TRUE;
            value = &deform->surfaces[surface].value[i];

            evaluate_image_value_fit_deriv( value->volume,
                                            value->voxel_lookup,
                                            value->continuity,
                                            value->image_weight,
                                            value->threshold,
                                            value->min_diff,
                                            value->max_diff,
                                            value->max_diff_weight,
                                            value->differential_offset,
                                            value->differential_ratio,
                                            value->oversample,
                   deform->surfaces[surface].surface.n_points,
                   deform->surfaces[surface].surface.n_neighbours,
                   deform->surfaces[surface].surface.neighbours,
                   this_parms, this_deriv );
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.value_fit = compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    // This is the stretch constraint for the mesh smoothness.

    for_less( surface, 0, deform->n_surfaces>0?1:0 ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_stretch ) {
            deriv_modified = TRUE;
            stretch = &deform->surfaces[surface].stretch[i];

            evaluate_stretch_fit_deriv( stretch->stretch_weight,
                                        stretch->n_points, 
                                        this_parms,
                                        stretch->n_neighbours,
                                        stretch->neighbours,
                                        this_deriv );
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.stretch_fit = compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    // This evaluates the gradient due to boundary terms (for white matter).
    // Because it uses ray tracing to find the nearest white voxel in the
    // normal direction (which is slow), the terms of the fit constraint
    // (a second order polynomial) are stored in linear, square, cross.

    if( BOUNDARY_DECREASE == 1 ) {
      if( linear != NULL ) {
        deriv_modified = TRUE;
        evaluate_quadratic_deriv_real( n_parameters, parameters,
                                       linear, square, n_cross_terms,
                                       cross_parms, cross_terms, derivative );
      }
    }

    ind = 0;
    for_less( surface, 0, deform->n_surfaces>0?1:0 ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_bound ) {
            bound = &deform->surfaces[surface].bound[i];

            if( bound->max_dist_weight > 0.0 ) {
                deriv_modified = TRUE;
                evaluate_boundary_search_fit_deriv(
                                              bound->image_weight_in,
                                              bound->image_weight_out,
                                              bound->max_dist_threshold,
                                              bound->max_dist_weight,
                                              bound->oversample,
                                              &boundary_flags[ind],
                                              &boundary_points[ind],
                                              0,
                       deform->surfaces[surface].surface.n_points,
                       deform->surfaces[surface].surface.n_points,
                       this_parms,
                       deform->surfaces[surface].surface.n_neighbours,
                       deform->surfaces[surface].surface.neighbours,
                       this_deriv,
                       // Added by June
                       bound->weights,
                       bound->max_weight_value,
                       deform->surfaces[surface].volume->adaptive_boundary_ratio);

                ind += deform->surfaces[surface].surface.n_points +
                       bound->oversample * deform->surfaces[surface].surface.n_edges +
                       ( bound->oversample * ( bound->oversample - 1 ) ) / 2 *
                       deform->surfaces[surface].surface.n_polygons;
            }
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.boundary_fit = compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    //////////////////////////////////////////////////////////////////
    // Added by June
    // constraint for WM surface intersection
    for_less( surface, 0, deform->n_surfaces>0?1:0 ) {
      this_parms = &parameters[start_parameter[surface]];
      this_deriv = &derivative[start_parameter[surface]];

      for_less( i, 0, deform->n_intersect_wm ) {
        deriv_modified = TRUE;
        evaluate_intersect_wm_fit_deriv(
                  deform->intersect_wm->objects,
                  deform->intersect_wm->weight,
                  this_parms,
                  deform->intersect_wm->wm_masked,
                  deform->intersect_wm->intersected,
                  &(deform->intersect_wm->n_intersected),
                  deform->surfaces[surface].surface.n_neighbours,
                  deform->surfaces[surface].surface.neighbours,
                  0,
                  deform->surfaces[surface].surface.n_points,
                  this_deriv );
      }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.surf_surf_fit += compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    //////////////////////////////////////////////////////////////////
    // Added by JUNE
    for_less( surface, 0, deform->n_surfaces>0?1:0 ) {
      this_parms = &parameters[start_parameter[surface]];
      this_deriv = &derivative[start_parameter[surface]];

      for_less( i, 0, deform->surfaces[surface].n_laplacian ) {
        deriv_modified = TRUE;
        evaluate_laplacian_fit_deriv(
                           deform->surfaces[surface].laplacian->weight,
                           deform->surfaces[surface].laplacian->volume,
                           deform->surfaces[surface].bound->volume,
                           deform->surfaces[surface].laplacian->from_value,
                           deform->surfaces[surface].laplacian->to_value,
                           deform->surfaces[surface].laplacian->type,
                           this_parms, 0,
                           deform->surfaces[surface].surface.n_points,
                           deform->surfaces[surface].surface.n_neighbours,
                           deform->surfaces[surface].surface.neighbours,
                           deform->surfaces[surface].laplacian->oversample,
                           this_deriv);
      }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.laplacian_fit = compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    //////////////////////////////////////////////////////////////////
    // Added by June
    for_less( surface, 0, deform->n_surfaces>0?1:0 ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_volume ) {
            deriv_modified = TRUE;
            bound = &deform->surfaces[surface].bound[i];
            anchor = &deform->surfaces[surface].anchors[i];

              // Added by June Sic Kim 9/12/2002
            evaluate_volume_fit_deriv(
                          deform->surfaces[surface].volume->weight,
                          deform->surfaces[surface].volume->max_weight,
                          deform->surfaces[surface].volume->direction,
                          bound->volume,
                          deform->surfaces[surface].laplacian->volume,
                          bound->threshold,
                          this_parms,
                          deform->surfaces[surface].surface.n_neighbours,
                          deform->surfaces[surface].surface.neighbours,
                          0,
                          deform->surfaces[surface].surface.n_points,
                          deform->surfaces[surface].surface.n_points,
                          this_deriv);
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.volume_fit = compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    /////////////////////////////////////////////////////////////////////

    for_less( surface, 0, deform->n_surfaces ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_curvature ) {
            deriv_modified = TRUE;
            curvature = &deform->surfaces[surface].curvature[i];

            evaluate_curvature_fit_deriv( curvature->curvature_weight,
                                          curvature->max_curvature_weight,
                                          curvature->min_curvature,
                                          curvature->max_curvature,
                                          curvature->n_points,
                                          this_parms,
                                          curvature->n_neighbours,
                                          curvature->neighbours,
                                          curvature->curvatures,
                                          this_deriv );
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.curvature_fit = compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    for_less( surface, 0, deform->n_surfaces ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_bend ) {
            deriv_modified = TRUE;
            bend = &deform->surfaces[surface].bend[i];

            evaluate_bend_fit_deriv( bend->bend_weight,
                                     bend->max_bend_weight,
                                     bend->min_bend,
                                     bend->max_bend,
                                     bend->n_points, this_parms,
                                     bend->n_neighbours,
                                     bend->neighbours,
                                     bend->model_x, bend->model_y,
                                     this_deriv );
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.bend_fit = compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    for_less( i, 0, deform->n_inter_surfaces ) {
        deriv_modified = TRUE;
        inter = &deform->inter_surfaces[i];

        evaluate_inter_surface_fit_deriv(
                   inter->weight,
                   inter->max_weight1,
                   inter->max_weight2,
                   inter->n_weight_steps,
                   inter->n_connections, inter->connections,
                   deform->surfaces[inter->surface_index1].surface.n_neighbours,
                   deform->surfaces[inter->surface_index1].surface.neighbours,
                   &parameters[start_parameter[inter->surface_index1]],
                   &parameters[start_parameter[inter->surface_index2]],
                   &derivative[start_parameter[inter->surface_index1]],
                   &derivative[start_parameter[inter->surface_index2]] );

        if( inter->oversample > 0 ) {
            evaluate_oversampled_inter_surface_fit_deriv(
                   inter->weight,
                   inter->max_weight1,
                   inter->max_weight2,
                   inter->n_weight_steps,
                   inter->n_connections, inter->connections,
                   inter->oversample,
                   deform->surfaces[inter->surface_index1].surface.n_points,
                   deform->surfaces[inter->surface_index1].surface.n_neighbours,
                   deform->surfaces[inter->surface_index1].surface.neighbours,
                   &parameters[start_parameter[inter->surface_index1]],
                   &parameters[start_parameter[inter->surface_index2]],
                   &derivative[start_parameter[inter->surface_index1]],
                   &derivative[start_parameter[inter->surface_index2]] );
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.inter_surface_fit = compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    if( NO_ANCHOR != 1 ) {
      for_less( surface, 0, deform->n_surfaces ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_anchors ) {
            deriv_modified = TRUE;
            anchor = &deform->surfaces[surface].anchors[i];

            evaluate_anchor_fit_deriv(0,
                       anchor->weight, anchor->max_dist_weight,
                       anchor->n_anchor_points, anchor->anchor_points,
                       this_parms, anchor->max_weight_value, this_deriv );
        }
      }

      eval.anchor_fit = compute_deriv_mag( n_parameters, derivative );

      if( deriv_modified ) {
        deriv_modified = FALSE;
        for_less( p, 0, n_parameters ) {
          full_deriv[p] += derivative[p];
          derivative[p] = 0.0;
        }
      }
    }

    for_less( surface, 0, deform->n_surfaces ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_weight_points ) {
            deriv_modified = TRUE;
            weight_point = &deform->surfaces[surface].weight_points[i];

            evaluate_weight_point_fit_deriv(
                       weight_point->weight, weight_point->max_dist_weight,
                       weight_point->n_weight_points,
                       weight_point->weight_points,
                       this_parms, this_deriv );
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.weight_point_fit = compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    for_less( surface, 0, deform->n_surfaces>0?1:0 ) {
        this_parms = &parameters[start_parameter[surface]];
        this_deriv = &derivative[start_parameter[surface]];

        for_less( i, 0, deform->surfaces[surface].n_self_intersects ) {
            deriv_modified = TRUE;
            self = &deform->surfaces[surface].self_intersects[i];

            evaluate_self_intersect_fit_deriv(
                       self->n_weights, self->weights, self->min_distances,
                       &si_lookup[surface][i],
                       self->use_tri_tri_dist,
                       self->square_flag,
                       deform->surfaces[surface].surface.triangles,
                       this_parms, this_deriv );
        }
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.self_intersect_fit = compute_deriv_mag( n_parameters, derivative );

      if( getenv( "PROJECT_DERIV" ) != NULL ) {
        if( sscanf( getenv( "PROJECT_DERIV" ), "%lf", &remove_ratio ) != 1 )
            remove_ratio = 1.0;
        remove_self_intersect_components( remove_ratio,
                                          n_parameters, full_deriv, derivative);
      } else {
        for_less( p, 0, n_parameters ) {
            full_deriv[p] += derivative[p];
        }
      }

      for_less( p, 0, n_parameters ) {
        derivative[p] = 0.0;
      }
    }

    for_less( i, 0, deform->n_surf_surfs ) {
        deriv_modified = TRUE;
        surf = &deform->surf_surfs[i];

        evaluate_surf_surf_fit_deriv(
                   surf->n_weights, surf->weights, surf->min_distances,
                   &ss_lookup[i],
                   deform->surfaces[surf->surface_index1].surface.n_points,
                   deform->surfaces[surf->surface_index1].surface.triangles,
                   &parameters[start_parameter[surf->surface_index1]],
                   &derivative[start_parameter[surf->surface_index1]],
                   deform->surfaces[surf->surface_index2].surface.n_points,
                   deform->surfaces[surf->surface_index2].surface.triangles,
                   &parameters[start_parameter[surf->surface_index2]],
                   &derivative[start_parameter[surf->surface_index2]] );
    }

    if( deriv_modified ) {
      deriv_modified = FALSE;
      eval.surf_surf_fit += compute_deriv_mag( n_parameters, derivative );

      for_less( p, 0, n_parameters ) {
        full_deriv[p] += derivative[p];
        derivative[p] = 0.0;
      }
    }

    FREE( derivative );

    if( fit_info != NULL )
        *fit_info = eval;

#ifdef DEBUG
    if( getenv( "DERIV_STEP" ) != NULL ) {
        int               p;
        Real              deriv_step, fit, fit1, fit2, save_param, used_delta;
        Real              deriv1, deriv2, min_deriv, max_deriv, tolerance;
        Real              diff, max_diff;
        fit_eval_struct   f_ignore;

        if( sscanf( getenv("DERIV_STEP"), "%lf %lf", &deriv_step,
            &tolerance ) != 2 )
            return;

        fit = evaluate_fit( deform, n_parameters, start_parameter,
                            parameters,
                            NULL, NULL,
                            0.0, linear, square, n_cross_terms,
                            cross_parms, cross_terms,
                            0.0, 0.0, si_lookup, &f_ignore );

        max_diff = deriv_step * FABS( fit );

        for_less( p, 0, n_parameters )
        {
            save_param = parameters[p];

            used_delta = 1.0e-4;

            do
            {
                used_delta /= 2.0;

                parameters[p] = save_param - used_delta;

                fit1 = evaluate_fit( deform, n_parameters, start_parameter,
                                     parameters, NULL, NULL,
                                     0.0, linear, square, n_cross_terms,
                                     cross_parms, cross_terms,
                                     0.0, 0.0, si_lookup, &f_ignore );

                parameters[p] = save_param + used_delta;

                fit2 = evaluate_fit( deform, n_parameters, start_parameter,
                                     parameters, NULL, NULL,
                                     0.0, linear, square, n_cross_terms,
                                     cross_parms, cross_terms,
                                     0.0, 0.0, si_lookup, &f_ignore );

                diff = FABS( fit2 - fit1 );
            }
            while( diff > max_diff );

            parameters[p] = save_param;

            deriv1 = (fit - fit1) / used_delta;
            deriv2 = (fit2 - fit) / used_delta;

            min_deriv = MIN( deriv1, deriv2 );
            max_deriv = MAX( deriv1, deriv2 );

            if( !numerically_close( min_deriv, max_deriv, tolerance ) ||
                (!numerically_close( min_deriv, full_deriv[p], tolerance ) &&
                 !numerically_close( max_deriv, full_deriv[p], tolerance )) )
            {
                print( "D %d: %g  %g  %g\n", p,
                       min_deriv, full_deriv[p], max_deriv );
            }
        }
    }
#endif
}

private  void  initialize_vector_interpolation(
    Vector   *v1,
    Vector   *v2,
    Vector   *axis,
    Real     *angle )
{
    Real     x, y;
    Vector   offset;

    CROSS_VECTORS( *axis, *v1, *v2 );

    if( null_Vector(axis) )
        *angle = 0.0;
    else
    {
        NORMALIZE_VECTOR( *axis, *axis );
        x = DOT_VECTORS( *v2, *v1 );
        SCALE_VECTOR( offset, *v1, x );
        SUB_VECTORS( offset, *v2, offset );
        y = MAGNITUDE( offset );

        *angle = 2.0 * PI - compute_clockwise_rotation( x, y );
        if( *angle == 2.0 * PI )
            *angle = 0.0;

        if( *angle < 0.0 || *angle > PI )
            handle_internal_error( "angle < 0.0 || 180.0" );
    }
}

private void  get_vector_interpolation(
    Vector  *v1,
    Vector  *axis,
    Real    angle,
    Real    alpha,
    Vector  *interp )
{
    Real       x, y, z;
    Transform  transform;

    if( angle == 0.0 )
        *interp = *v1;
    else
    {
        make_rotation_about_axis( axis, alpha * angle, &transform );

        transform_vector( &transform,
                          (Real) Vector_x(*v1),
                          (Real) Vector_y(*v1),
                          (Real) Vector_z(*v1), &x, &y, &z );

        fill_Vector( *interp, x, y, z );
    }

    NORMALIZE_VECTOR( *interp, *interp );
}


private  void  get_oversample_boundaries(
    Volume                      volume,
    voxel_coef_struct           *voxel_lookup,
    bitlist_3d_struct           *done_bits,
    bitlist_3d_struct           *surface_bits,
    Real                        isovalue,
    Normal_directions           normal_direction,
    Real                        max_outward,
    Real                        max_inward,
    Real                        distance_offset,
    BOOLEAN                     check_direction_flag,
    clip_struct                 *clip_search,
    int                         oversample,
    int                         n_neighbours[],
    int                         *neighbours[],
    Real                        parameters[],
    int                         p1,
    int                         p2,
    int                         p3,
    Point                       neigh_points[],
    Smallest_int                boundary_flags[],
    Point                       boundary_points[] )
{
    int              ind, n, neigh, w1, w2, p1_index, p2_index, p3_index;
    int              n_nodes_to_ignore, nodes_to_ignore[3];
    Real             dist, alpha1, alpha2, angle_p12, angle_p32;
    Point            search_point;
    Vector           p1_normal, p2_normal, p3_normal;
    Vector           search_normal, normal1, normal2;
    Vector           perp_p12, perp_p32, perp;
    Real             inward_search, outward_search, value;
    Real             angle, out_dist, in_dist;
    Point            point1, point2, point3, start1, start2;
    BOOLEAN          found;

    /*--- find p1 normal */

    for_less( n, 0, n_neighbours[p1] )
    {
        neigh = neighbours[p1][n];
        fill_Point( neigh_points[n],
                    parameters[IJ(neigh,0,3)],
                    parameters[IJ(neigh,1,3)],
                    parameters[IJ(neigh,2,3)] );
    }

    find_polygon_normal( n_neighbours[p1], neigh_points, &p1_normal );

    /*--- find p2 normal */

    for_less( n, 0, n_neighbours[p2] )
    {
        neigh = neighbours[p2][n];
        fill_Point( neigh_points[n],
                    parameters[IJ(neigh,0,3)],
                    parameters[IJ(neigh,1,3)],
                    parameters[IJ(neigh,2,3)] );
    }

    find_polygon_normal( n_neighbours[p2], neigh_points, &p2_normal );

    /*--- find p3 normal */

    for_less( n, 0, n_neighbours[p3] )
    {
        neigh = neighbours[p3][n];
        fill_Point( neigh_points[n],
                    parameters[IJ(neigh,0,3)],
                    parameters[IJ(neigh,1,3)],
                    parameters[IJ(neigh,2,3)] );
    }

    find_polygon_normal( n_neighbours[p3], neigh_points, &p3_normal );

    p1_index = IJ( p1, 0, 3 );
    p2_index = IJ( p2, 0, 3 );
    p3_index = IJ( p3, 0, 3 );

    fill_Point( point1,
                parameters[p1_index+0], parameters[p1_index+1],
                parameters[p1_index+2] );

    fill_Point( point2,
                parameters[p2_index+0], parameters[p2_index+1],
                parameters[p2_index+2] );

    fill_Point( point3,
                parameters[p3_index+0], parameters[p3_index+1],
                parameters[p3_index+2] );

    initialize_vector_interpolation( &p1_normal, &p2_normal,
                                     &perp_p12, &angle_p12 );
    initialize_vector_interpolation( &p3_normal, &p2_normal,
                                     &perp_p32, &angle_p32 );

    ind = 0;

    for_inclusive( w1, 0, oversample )
    {
        alpha1 = (Real) w1 / (Real) (oversample+1);

        INTERPOLATE_POINTS( start1, point1, point2, alpha1 );
        INTERPOLATE_POINTS( start2, point3, point2, alpha1 );

        get_vector_interpolation( &p1_normal, &perp_p12, angle_p12, alpha1,
                                  &normal1 );
        get_vector_interpolation( &p3_normal, &perp_p32, angle_p32, alpha1,
                                  &normal2 );

        initialize_vector_interpolation( &normal1, &normal2, &perp, &angle );

        for_inclusive( w2, 1, oversample - w1 + 1 )
        {
            if( w1 == 0 && w2 == oversample + 1 ||
                w2 == oversample - w1 + 1 && p2 < p3 )
            {
                continue;
            }

            alpha2 = (Real) w2 / (Real) (oversample - w1 +1);

            INTERPOLATE_POINTS( search_point, start1, start2, alpha2 );

            get_vector_interpolation( &normal1, &perp, angle, alpha2,
                                      &search_normal );

            outward_search = max_outward - distance_offset;
            inward_search = max_inward + distance_offset;

            if( check_direction_flag )
            {
                evaluate_volume_in_world( volume,
                                      RPoint_x(search_point),
                                      RPoint_y(search_point),
                                      RPoint_z(search_point),
                                      EVAL_LINEAR, FALSE,
                                      0.0, &value,
                                      NULL, NULL, NULL,
                                      NULL, NULL, NULL, NULL, NULL, NULL );

                if( normal_direction == TOWARDS_LOWER &&
                    value < isovalue ||
                    normal_direction == TOWARDS_HIGHER &&
                    value > isovalue )
                {
                    outward_search = 0.1;
                }
                else if( normal_direction == TOWARDS_HIGHER &&
                         value < isovalue ||
                         normal_direction == TOWARDS_LOWER &&
                         value > isovalue )
                {
                    inward_search = 0.1;
                }
            }

            out_dist = outward_search;
            in_dist = inward_search;

            if( clip_search != NULL )
            {
                n_nodes_to_ignore = 0;
                if( w2 != oversample - w1 + 1 )
                {
                    nodes_to_ignore[n_nodes_to_ignore] = p1;
                    ++n_nodes_to_ignore;
                }

                if( w1 != 0 )
                {
                    nodes_to_ignore[n_nodes_to_ignore] = p2;
                    ++n_nodes_to_ignore;
                }

                if( w2 != 0 )
                {
                    nodes_to_ignore[n_nodes_to_ignore] = p3;
                    ++n_nodes_to_ignore;
                }

                clip_search_line( clip_search,
                                  n_nodes_to_ignore, nodes_to_ignore,
                                  &search_point, &search_normal,
                                  out_dist, in_dist, &out_dist, &in_dist );
            }

            found = find_isosurface_boundary_in_direction( volume, voxel_lookup,
                                            done_bits, surface_bits,
                                            &search_point,
                                            &search_normal,
                                            out_dist,
                                            in_dist, 0,
                                            isovalue, normal_direction, &dist );

            if( found )
            {
                dist = dist + distance_offset;
                GET_POINT_ON_RAY( boundary_points[ind], search_point,
                                  search_normal, dist );
                if( dist >= 0.0 )
                    boundary_flags[ind] = (Smallest_int) BOUNDARY_IS_OUTSIDE;
                else
                    boundary_flags[ind] = (Smallest_int) BOUNDARY_IS_INSIDE;
            }
            else
                boundary_flags[ind] = (Smallest_int) BOUNDARY_NOT_FOUND;

            ++ind;
        }
    }
}

private  void  find_image_boundaries(
    Volume                      volume,
    voxel_coef_struct           *voxel_lookup,
    bitlist_3d_struct           *done_bits,
    bitlist_3d_struct           *surface_bits,
    Real                        isovalue,
    Normal_directions           normal_direction,
    Real                        image_weight_out,
    Real                        image_weight_in,
    Real                        max_outward,
    Real                        max_inward,
    Real                        distance_offset,
    BOOLEAN                     check_direction_flag,
    BOOLEAN                     normal_direction_only,
    BOOLEAN                     clip_to_surface,
    int                         oversample,
    int                         n_nodes,
    int                         n_neighbours[],
    int                         *neighbours[],
    int                         start_parameter,
    Real                        parameters[],
    Smallest_int                active_flags[],
    Smallest_int                boundary_flags[],
    Point                       boundary_points[],
    Real                        *constant,
    Real                        linear[],
    Real                        square[],
    int                         *n_cross_terms[],
    int                         **cross_parms[],
    Real                        **cross_terms[] )
{
    int                node, n, neigh, neigh1, neigh2, ind, p_index, pos;
    int                max_neighbours, n1_index, n2_index, w1, w2;
    Real               dist, lx1, ly1, lz1, lx2, ly2, lz2, cons, image_weight;
    Real               weight1, weight2, weight3, max_diff;
    Real               lx3, ly3, lz3, nx, ny, nz, pd, len, alpha1, alpha2;
    Real               sx1, sy1, sz1, sx2, sy2, sz2, sx3, sy3, sz3;
    Real               lx1y1, lx1z1, lx1x2, lx1y2, lx1z2, lx1x3, lx1y3, lx1z3;
    Real               ly1z1, ly1x2, ly1y2, ly1z2, ly1x3, ly1y3, ly1z3;
    Real               lz1x2, lz1y2, lz1z2, lz1x3, lz1y3, lz1z3;
    Real               lx2y2, lx2z2, lx2x3, lx2y3, lx2z3;
    Real               ly2z2, ly2x3, ly2y3, ly2z3;
    Real               lz2x3, lz2y3, lz2z3;
    Real               lx3y3, lx3z3;
    Real               ly3z3;
    Real               out_dist, in_dist;
    Real               bx, by, bz, inward_search, outward_search, value;
    Point              origin, bound;
    Smallest_int       *found_flags;
    Point              *neigh_points, *boundaries;
    Vector             normal;
    BOOLEAN            found, active;
    progress_struct    progress;
    clip_struct        *clip_search;

    initialize_progress_report( &progress, FALSE, n_nodes,
                                "Finding Image Boundaries" );

    // weights have been normalized for oversampling rate in surface_fit.c.
    max_diff = MAX( image_weight_out * max_outward * max_outward,
                    image_weight_in * max_inward * max_inward );

    // Normalize weight by local area around each vertex.
    // This means that mesh will move towards the boundary at
    // the same rate if size of triangles are different locally.

    Real * areas;
    ALLOC( areas, n_nodes );
    for_less( node, 0, n_nodes ) areas[node] = 0.0;

    for_less( node, 0, n_nodes ) {
      int ind0 = IJ(node,0,3);
      for_less( n, 0, n_neighbours[node] ) {
        int n1 = neighbours[node][n];
        int n2 = neighbours[node][(n+1)%n_neighbours[node]];
        if( !( node < n1 && node < n2 ) ) continue;
        int ind1 = IJ(n1,0,3);
        int ind2 = IJ(n2,0,3);
        Real e0 = sqrt( ( parameters[ind0+0] - parameters[ind1+0] ) *
                        ( parameters[ind0+0] - parameters[ind1+0] ) +
                        ( parameters[ind0+1] - parameters[ind1+1] ) *
                        ( parameters[ind0+1] - parameters[ind1+1] ) +
                        ( parameters[ind0+2] - parameters[ind1+2] ) *
                        ( parameters[ind0+2] - parameters[ind1+2] ) );
        Real e1 = sqrt( ( parameters[ind2+0] - parameters[ind1+0] ) *
                        ( parameters[ind2+0] - parameters[ind1+0] ) +
                        ( parameters[ind2+1] - parameters[ind1+1] ) *
                        ( parameters[ind2+1] - parameters[ind1+1] ) +
                        ( parameters[ind2+2] - parameters[ind1+2] ) *
                        ( parameters[ind2+2] - parameters[ind1+2] ) );
        Real e2 = sqrt( ( parameters[ind0+0] - parameters[ind2+0] ) *
                        ( parameters[ind0+0] - parameters[ind2+0] ) +
                        ( parameters[ind0+1] - parameters[ind2+1] ) *
                        ( parameters[ind0+1] - parameters[ind2+1] ) +
                        ( parameters[ind0+2] - parameters[ind2+2] ) *
                        ( parameters[ind0+2] - parameters[ind2+2] ) );
        Real s = 0.5 * ( e0 + e1 + e2 );
        Real local_area = sqrt( fabs( s * ( s - e0 ) * ( s - e1 ) * ( s - e2 ) ) +
                                1.0e-20 );
        areas[node] += local_area;
        areas[n1] += local_area;
        areas[n2] += local_area;
      }
    }

    Real area_avg = 0.0;
    for_less( node, 0, n_nodes ) {
      area_avg += areas[node];
    }
    area_avg /= (Real)n_nodes;
    for_less( node, 0, n_nodes ) {
      areas[node] /= area_avg;
    }

    max_neighbours = 0;
    for_less( node, 0, n_nodes )
        max_neighbours = MAX( max_neighbours, n_neighbours[node] );

    ALLOC( neigh_points, max_neighbours );

    if( clip_to_surface ) {
        clip_search = initialize_clip_search( n_nodes, n_neighbours, neighbours,
                                              parameters ); 
    } else
        clip_search = NULL; 

    cons = 0.0;

    for_less( node, 0, n_nodes ) {
        active = active_flags == NULL || active_flags[node];

        for_less( n, 0, n_neighbours[node] ) {
            neigh = neighbours[node][n];
            fill_Point( neigh_points[n],
                        parameters[IJ(neigh,0,3)],
			parameters[IJ(neigh,1,3)],
                        parameters[IJ(neigh,2,3)] );

            if( !active && active_flags[neigh] )
                active = TRUE;
        }

        if( !active ) continue;

        find_polygon_normal( n_neighbours[node], neigh_points, &normal );

        fill_Point( origin,
                    parameters[IJ(node,0,3)],
                    parameters[IJ(node,1,3)],
                    parameters[IJ(node,2,3)] );

        outward_search = max_outward - distance_offset;
        inward_search = max_inward + distance_offset;

        if( check_direction_flag ) {
            evaluate_volume_in_world( volume,
                                      parameters[IJ(node,0,3)],
                                      parameters[IJ(node,1,3)],
                                      parameters[IJ(node,2,3)],
                                      EVAL_LINEAR, FALSE, 0.0, &value,
                                      NULL, NULL, NULL,
                                      NULL, NULL, NULL, NULL, NULL, NULL );

            if( ( normal_direction == TOWARDS_LOWER && value < isovalue ) ||
                ( normal_direction == TOWARDS_HIGHER && value > isovalue ) ) {
// Claude: 0.1 is for 1mm voxel size. Maybe make this 10% of voxel size.
                outward_search = 0.1;
            } else if( ( normal_direction == TOWARDS_HIGHER && value < isovalue ) ||
                       ( normal_direction == TOWARDS_LOWER && value > isovalue ) ) {
                inward_search = 0.1;
            }
        }

        out_dist = outward_search;
        in_dist = inward_search;

        if( clip_to_surface ) {
            clip_search_line( clip_search, 1, &node, &origin, &normal,
                              out_dist, in_dist, &out_dist, &in_dist );
        }

        found = find_isosurface_boundary_in_direction( volume, voxel_lookup,
                                            done_bits, surface_bits,
                                            &origin, &normal,
                                            out_dist, in_dist, 0,
                                            isovalue, normal_direction, &dist );

        if( found ) {
            dist += distance_offset;

            GET_POINT_ON_RAY( bound, origin, normal, dist );

            if( dist >= 0.0 ) {
                image_weight = image_weight_out;
            } else {
                image_weight = image_weight_in;
            }
            image_weight *= areas[node];

            if( boundary_points != NULL ) {
                boundary_points[node] = bound;

                if( dist >= 0.0 ) {
                    boundary_flags[node] = BOUNDARY_IS_OUTSIDE;
                } else {
                    boundary_flags[node] = BOUNDARY_IS_INSIDE;
                }
            }

            p_index = IJ(node,0,3);
            bx = RPoint_x(bound);
            by = RPoint_y(bound);
            bz = RPoint_z(bound);

            // The fit function, in its simplest form, is:
            // dist( point, bound )^2 = (x-bx)^2 + (y-by)^2 + (z-bz)^2
            //                        = bx^2 + by^2 + bz^2        (cons)
            //                          -2*(x*bx + y*by + z*bz)   (linear)
            //                          + x^2 + y^2 + z^2         (quadratic)
            // The derivative can easily be evaluated from these factors.
            // There can also be cross term in the normal direction and with
            // oversampling (terms in x*y, x*z, y*z). The motivation for this
            // expansion is that boundary terms can be evaluated only every
            // few iterations, thus cheaper.

            if( normal_direction_only ) {

                // This is to go towards (bx,by,bz) always in the normal 
                // direction as defined on the first iteration when the 
                // boundary terms are evaluated (local direction unaffected
                // by stretch nor self-intersection).
                // The form of the fit function in the normal direction is:
                //   d   = ((x-xb)*n)n
                //   d^2 = ( (x-xb)*nx + (y-yb)*ny + (z-zb)*nz )^2
                //   (nx,ny,nz) = (x0-xb,y0-yb,z0-zb)/len
                // (It makes a bunch of messy cross terms).
                nx = parameters[p_index+0] - bx;
                ny = parameters[p_index+1] - by;
                nz = parameters[p_index+2] - bz;
                len = nx * nx + ny * ny + nz * nz;
                if( len == 0.0 ) {
                    nx = RPoint_x(normal);
                    ny = RPoint_y(normal);
                    nz = RPoint_z(normal);
                } else {
                    len = sqrt( len );
                    nx /= len;
                    ny /= len;
                    nz /= len;
                }

                pd = -(nx * bx + ny * by + nz * bz);

                cons += image_weight * pd * pd;
                linear[start_parameter+p_index+0] += image_weight * 2.0 * nx*pd;
                linear[start_parameter+p_index+1] += image_weight * 2.0 * ny*pd;
                linear[start_parameter+p_index+2] += image_weight * 2.0 * nz*pd;
                square[start_parameter+p_index+0] += image_weight * nx * nx;
                square[start_parameter+p_index+1] += image_weight * ny * ny;
                square[start_parameter+p_index+2] += image_weight * nz * nz;

                lx1y1 = image_weight * 2.0 * nx * ny;
                lx1z1 = image_weight * 2.0 * nx * nz;
                ly1z1 = image_weight * 2.0 * ny * nz;
                add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + p_index + 1, lx1y1, 5 );
                add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + p_index + 2, lx1z1, 5 );
                add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 1,
                                  start_parameter + p_index + 2, ly1z1, 5 );
            } else {
                // This is to point towards (bx,by,bz) from the current point.
                // At the first iteration, this is the same as in the normal
                // direction, but on subsequent iterations with the boundary
                // terms fixed, the direction may not be exactly in the 
                // normal direction (since other constraints such as stretch
                // and self-intersection may slide the node sideways).
		cons += image_weight * (bx * bx + by * by + bz * bz);
                linear[start_parameter+p_index+0] += image_weight * -2.0 * bx;
                linear[start_parameter+p_index+1] += image_weight * -2.0 * by;
                linear[start_parameter+p_index+2] += image_weight * -2.0 * bz;
                square[start_parameter+p_index+0] += image_weight;
                square[start_parameter+p_index+1] += image_weight;
                square[start_parameter+p_index+2] += image_weight;
            }
        } else {
            if( boundary_points != NULL ) {
                boundary_flags[node] = (Smallest_int) BOUNDARY_NOT_FOUND;
            }
            cons += max_diff;
        }

        update_progress_report( &progress, node+1 );
    }

    terminate_progress_report( &progress );

    ind = n_nodes;

    if( oversample > 0 ) {
        ALLOC( boundaries, (oversample+2) * (oversample+1) / 2 );
        ALLOC( found_flags, (oversample+2) * (oversample+1) / 2 );

        initialize_progress_report( &progress, FALSE, n_nodes,
                                    "Finding Oversample boundaries" );

        for_less( node, 0, n_nodes ) {
            active = active_flags == NULL || active_flags[node];

            for_less( n, 0, n_neighbours[node] ) {
                neigh = neighbours[node][n];
                fill_Point( neigh_points[n],
                            parameters[IJ(neigh,0,3)],
                            parameters[IJ(neigh,1,3)],
                            parameters[IJ(neigh,2,3)] );
            }

            find_polygon_normal( n_neighbours[node], neigh_points, &normal );

            for_less( n, 0, n_neighbours[node] ) {
                neigh1 = neighbours[node][n];
                neigh2 = neighbours[node][(n+1)%n_neighbours[node]];

                if( node > neigh1 || node > neigh2 )
                    continue;

                if( !active && !active_flags[neigh1] && !active_flags[neigh2] ) {
                    for_inclusive( w1, 0, oversample ) {
                        for_inclusive( w2, 1, oversample - w1 + 1 ) {
                            if( ( w1 == 0 && w2 == oversample + 1 ) ||
                                ( w2 == oversample - w1 + 1 && neigh1 < neigh2 ) ) {
                                continue;
                            }
                            ++ind;
                        }
                    }

                    continue;
                }

                p_index = IJ(node,0,3);
                n1_index = IJ(neigh1,0,3);
                n2_index = IJ(neigh2,0,3);

                get_oversample_boundaries( volume, voxel_lookup, done_bits,
                                     surface_bits, isovalue, normal_direction,
                                     max_outward,
                                     max_inward, distance_offset,
                                     check_direction_flag, clip_search,
                                     oversample,
                                     n_neighbours, neighbours,
                                     parameters,
                                     node, neigh1, neigh2, neigh_points,
                                     found_flags, boundaries );

                if( boundary_points != NULL ) {
                    pos = 0;
                    for_inclusive( w1, 0, oversample ) {
                        for_inclusive( w2, 1, oversample - w1 + 1 ) {
                            if( w1 == 0 && w2 == oversample + 1 ||
                                w2 == oversample - w1 + 1 && neigh1 < neigh2 ) {
                                continue;
                            }
                            boundary_points[ind+pos] = boundaries[pos];
                            boundary_flags[ind+pos] = found_flags[pos];
                            ++pos;
                        }
                    }
                }

                lx1 = 0.0;
                ly1 = 0.0;
                lz1 = 0.0;
                lx2 = 0.0;
                ly2 = 0.0;
                lz2 = 0.0;
                lx3 = 0.0;
                ly3 = 0.0;
                lz3 = 0.0;
                sx1 = 0.0;
                sx2 = 0.0;
                sx3 = 0.0;
                sy1 = 0.0;
                sy2 = 0.0;
                sy3 = 0.0;
                sz1 = 0.0;
                sz2 = 0.0;
                sz3 = 0.0;
                lx1y1 = 0.0;
                lx1z1 = 0.0;
                lx1x2 = 0.0;
                lx1y2 = 0.0;
                lx1z2 = 0.0;
                lx1x3 = 0.0;
                lx1y3 = 0.0;
                lx1z3 = 0.0;

                ly1z1 = 0.0;
                ly1x2 = 0.0;
                ly1y2 = 0.0;
                ly1z2 = 0.0;
                ly1x3 = 0.0;
                ly1y3 = 0.0;
                ly1z3 = 0.0;

                lz1x2 = 0.0;
                lz1y2 = 0.0;
                lz1z2 = 0.0;
                lz1x3 = 0.0;
                lz1y3 = 0.0;
                lz1z3 = 0.0;

                lx2y2 = 0.0;
                lx2z2 = 0.0;
                lx2x3 = 0.0;
                lx2y3 = 0.0;
                lx2z3 = 0.0;

                ly2z2 = 0.0;
                ly2x3 = 0.0;
                ly2y3 = 0.0;
                ly2z3 = 0.0;

                lz2x3 = 0.0;
                lz2y3 = 0.0;
                lz2z3 = 0.0;

                lx3y3 = 0.0;
                lx3z3 = 0.0;

                ly3z3 = 0.0;

                pos = 0;
                for_inclusive( w1, 0, oversample ) {
                    alpha1 = (Real) w1 / (Real) (oversample+1);
                    for_inclusive( w2, 1, oversample - w1 + 1 ) {
                        if( w1 == 0 && w2 == oversample + 1 ||
                            w2 == oversample - w1 + 1 && neigh1 < neigh2 ) {
                            continue;
                        }

                        alpha2 = (Real) w2 / (Real) (oversample-w1+1);

                        if( found_flags[pos] != BOUNDARY_NOT_FOUND ) {
                            if( found_flags[pos] == BOUNDARY_IS_OUTSIDE ) {
                                image_weight = image_weight_out;
                            } else {
                                image_weight = image_weight_in;
                            }

                            bx = RPoint_x(boundaries[pos]);
                            by = RPoint_y(boundaries[pos]);
                            bz = RPoint_z(boundaries[pos]);

                            weight1 = (1.0 - alpha1) * (1.0 - alpha2);
                            weight2 = alpha1;
                            weight3 = (1.0 - alpha1) * alpha2;

                            image_weight *= ( weight1 * areas[node] +
                                              weight2 * areas[neigh1] +
                                              weight3 * areas[neigh2] );

                            Real image_w1_weight = image_weight * weight1;
                            Real image_w2_weight = image_weight * weight2;
                            Real image_w3_weight = image_weight * weight3;

                            if( normal_direction_only ) {
                                nx = weight1 * parameters[p_index+0] +
                                     weight2 * parameters[n1_index+0] +
                                     weight3 * parameters[n2_index+0] - bx;
                                ny = weight1 * parameters[p_index+1] +
                                     weight2 * parameters[n1_index+1] +
                                     weight3 * parameters[n2_index+1] - by;
                                nz = weight1 * parameters[p_index+2] +
                                     weight2 * parameters[n1_index+2] +
                                     weight3 * parameters[n2_index+2] - bz;

                                len = nx * nx + ny * ny + nz * nz;
                                if( len == 0.0 ) {
                                    nx = RPoint_x(normal);
                                    ny = RPoint_y(normal);
                                    nz = RPoint_z(normal);
                                } else {
                                    len = sqrt( len );
                                    nx /= len;
                                    ny /= len;
                                    nz /= len;
                                }

                                pd = -(nx * bx + ny * by + nz * bz);

                                Real nx_nx = nx * nx;
                                Real nx_ny = nx * ny;
                                Real nx_nz = nx * nz;
                                Real ny_ny = ny * ny;
                                Real ny_nz = ny * nz;
                                Real nz_nz = nz * nz;

                                Real image_w1_w1_weight = image_w1_weight * weight1;
                                Real image_w1_w2_weight = image_w1_weight * weight2;
                                Real image_w1_w3_weight = image_w1_weight * weight3;
                                Real image_w2_w2_weight = image_w2_weight * weight2;
                                Real image_w2_w3_weight = image_w2_weight * weight3;
                                Real image_w3_w3_weight = image_w3_weight * weight3;

                                cons += image_weight * pd * pd;
                                lx1 += image_w1_weight * nx * pd;
                                ly1 += image_w1_weight * ny * pd;
                                lz1 += image_w1_weight * nz * pd;
                                lx2 += image_w2_weight * nx * pd;
                                ly2 += image_w2_weight * ny * pd;
                                lz2 += image_w2_weight * nz * pd;
                                lx3 += image_w3_weight * nx * pd;
                                ly3 += image_w3_weight * ny * pd;
                                lz3 += image_w3_weight * nz * pd;

                                sx1 += image_w1_w1_weight * nx_nx;
                                sy1 += image_w1_w1_weight * ny_ny;
                                sz1 += image_w1_w1_weight * nz_nz;
                                sx2 += image_w2_w2_weight * nx_nx;
                                sy2 += image_w2_w2_weight * ny_ny;
                                sz2 += image_w2_w2_weight * nz_nz;
                                sx3 += image_w3_w3_weight * nx_nx;
                                sy3 += image_w3_w3_weight * ny_ny;
                                sz3 += image_w3_w3_weight * nz_nz;

                                lx1y1 += image_w1_w1_weight * nx_ny;
                                lx1z1 += image_w1_w1_weight * nx_nz;
                                lx1x2 += image_w1_w2_weight * nx_nx;
                                lx1y2 += image_w1_w2_weight * nx_ny;
                                lx1z2 += image_w1_w2_weight * nx_nz;
                                lx1x3 += image_w1_w3_weight * nx_nx;
                                lx1y3 += image_w1_w3_weight * nx_ny;
                                lx1z3 += image_w1_w3_weight * nx_nz;

                                ly1z1 += image_w1_w1_weight * ny_nz;
                                ly1x2 += image_w1_w2_weight * nx_ny;
                                ly1y2 += image_w1_w2_weight * ny_ny;
                                ly1z2 += image_w1_w2_weight * ny_nz;
                                ly1x3 += image_w1_w3_weight * nx_ny;
                                ly1y3 += image_w1_w3_weight * ny_ny;
                                ly1z3 += image_w1_w3_weight * ny_nz;

                                lz1x2 += image_w1_w2_weight * nx_nz;
                                lz1y2 += image_w1_w2_weight * ny_nz;
                                lz1z2 += image_w1_w2_weight * nz_nz;
                                lz1x3 += image_w1_w3_weight * nx_nz;
                                lz1y3 += image_w1_w3_weight * ny_nz;
                                lz1z3 += image_w1_w3_weight * nz_nz;

                                lx2y2 += image_w2_w2_weight * nx_ny;
                                lx2z2 += image_w2_w2_weight * nx_nz;
                                lx2x3 += image_w2_w3_weight * nx_nx;
                                lx2y3 += image_w2_w3_weight * nx_ny;
                                lx2z3 += image_w2_w3_weight * nx_nz;

                                ly2z2 += image_w2_w2_weight * ny_nz;
                                ly2x3 += image_w2_w3_weight * nx_ny;
                                ly2y3 += image_w2_w3_weight * ny_ny;
                                ly2z3 += image_w2_w3_weight * ny_nz;

                                lz2x3 += image_w2_w3_weight * nx_nz;
                                lz2y3 += image_w2_w3_weight * ny_nz;
                                lz2z3 += image_w2_w3_weight * nz_nz;

                                lx3y3 += image_w3_w3_weight * nx_ny;
                                lx3z3 += image_w3_w3_weight * nx_nz;

                                ly3z3 += image_w3_w3_weight * ny_nz;
                            } else {
                                cons += image_weight *
                                        (bx * bx + by * by + bz * bz);

                                lx1 += -bx * image_w1_weight;
                                ly1 += -by * image_w1_weight;
                                lz1 += -bz * image_w1_weight;

                                lx2 += -bx * image_w2_weight;
                                ly2 += -by * image_w2_weight;
                                lz2 += -bz * image_w2_weight;

                                lx3 += -bx * image_w3_weight;
                                ly3 += -by * image_w3_weight;
                                lz3 += -bz * image_w3_weight;

                                sx1 += weight1 * image_w1_weight;
                                sx2 += weight2 * image_w2_weight;
                                sx3 += weight3 * image_w3_weight;

                                lx1x2 += weight1 * image_w2_weight;
                                lx1x3 += weight1 * image_w3_weight;
                                lx2x3 += weight2 * image_w3_weight;
                             }
                        } else {
                            cons += max_diff;
                        }

                        ++pos;
                    }
                }

                ind += pos;

                if( normal_direction_only ) {
                    linear[start_parameter+p_index+0] += lx1+lx1;
                    linear[start_parameter+p_index+1] += ly1+ly1;
                    linear[start_parameter+p_index+2] += lz1+lz1;
                    square[start_parameter+p_index+0] += sx1;
                    square[start_parameter+p_index+1] += sy1;
                    square[start_parameter+p_index+2] += sz1;

                    linear[start_parameter+n1_index+0] += lx2+lx2;
                    linear[start_parameter+n1_index+1] += ly2+ly2;
                    linear[start_parameter+n1_index+2] += lz2+lz2;
                    square[start_parameter+n1_index+0] += sx2;
                    square[start_parameter+n1_index+1] += sy2;
                    square[start_parameter+n1_index+2] += sz2;

                    linear[start_parameter+n2_index+0] += lx3+lx3;
                    linear[start_parameter+n2_index+1] += ly3+ly3;
                    linear[start_parameter+n2_index+2] += lz3+lz3;
                    square[start_parameter+n2_index+0] += sx3;
                    square[start_parameter+n2_index+1] += sy3;
                    square[start_parameter+n2_index+2] += sz3;

                    // multiplication by 2.
                    lx1y1 += lx1y1;
                    lx1z1 += lx1z1;
                    lx1x2 += lx1x2;
                    lx1y2 += lx1y2;
                    lx1z2 += lx1z2;
                    lx1x3 += lx1x3;
                    lx1y3 += lx1y3;
                    lx1z3 += lx1z3;

                    ly1z1 += ly1z1;
                    ly1x2 += ly1x2;
                    ly1y2 += ly1y2;
                    ly1z2 += ly1z2;
                    ly1x3 += ly1x3;
                    ly1y3 += ly1y3;
                    ly1z3 += ly1z3;

                    lz1x2 += lz1x2;
                    lz1y2 += lz1y2;
                    lz1z2 += lz1z2;
                    lz1x3 += lz1x3;
                    lz1y3 += lz1y3;
                    lz1z3 += lz1z3;

                    lx2y2 += lx2y2;
                    lx2z2 += lx2z2;
                    lx2x3 += lx2x3;
                    lx2y3 += lx2y3;
                    lx2z3 += lx2z3;

                    ly2z2 += ly2z2;
                    ly2x3 += ly2x3;
                    ly2y3 += ly2y3;
                    ly2z3 += ly2z3;

                    lz2x3 += lz2x3;
                    lz2y3 += lz2y3;
                    lz2z3 += lz2z3;

                    lx3y3 += lx3y3;
                    lx3z3 += lx3z3;

                    ly3z3 += ly3z3;

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + p_index + 1, lx1y1, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + p_index + 2, lx1z1, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + n1_index + 0, lx1x2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + n1_index + 1, lx1y2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + n1_index + 2, lx1z2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + n2_index + 0, lx1x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + n2_index + 1, lx1y3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 0,
                                  start_parameter + n2_index + 2, lx1z3, 5 );

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 1,
                                  start_parameter + p_index + 2, ly1z1, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 1,
                                  start_parameter + n1_index + 0, ly1x2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 1,
                                  start_parameter + n1_index + 1, ly1y2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 1,
                                  start_parameter + n1_index + 2, ly1z2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 1,
                                  start_parameter + n2_index + 0, ly1x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 1,
                                  start_parameter + n2_index + 1, ly1y3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 1,
                                  start_parameter + n2_index + 2, ly1z3, 5 );

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 2,
                                  start_parameter + n1_index + 0, lz1x2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 2,
                                  start_parameter + n1_index + 1, lz1y2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 2,
                                  start_parameter + n1_index + 2, lz1z2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 2,
                                  start_parameter + n2_index + 0, lz1x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 2,
                                  start_parameter + n2_index + 1, lz1y3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + p_index + 2,
                                  start_parameter + n2_index + 2, lz1z3, 5 );

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 0,
                                  start_parameter + n1_index + 1, lx2y2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 0,
                                  start_parameter + n1_index + 2, lx2z2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 0,
                                  start_parameter + n2_index + 0, lx2x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 0,
                                  start_parameter + n2_index + 1, lx2y3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 0,
                                  start_parameter + n2_index + 2, lx2z3, 5 );

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 1,
                                  start_parameter + n1_index + 2, ly2z2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 1,
                                  start_parameter + n2_index + 0, ly2x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 1,
                                  start_parameter + n2_index + 1, ly2y3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 1,
                                  start_parameter + n2_index + 2, ly2z3, 5 );

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 2,
                                  start_parameter + n2_index + 0, lz2x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 2,
                                  start_parameter + n2_index + 1, lz2y3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n1_index + 2,
                                  start_parameter + n2_index + 2, lz2z3, 5 );

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n2_index + 0,
                                  start_parameter + n2_index + 1, lx3y3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n2_index + 0,
                                  start_parameter + n2_index + 2, lx3z3, 5 );

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                  cross_parms, cross_terms,
                                  start_parameter + n2_index + 1,
                                  start_parameter + n2_index + 2, ly3z3, 5 );
                } else {
                    linear[start_parameter+p_index+0] += lx1+lx1;
                    linear[start_parameter+p_index+1] += ly1+ly1;
                    linear[start_parameter+p_index+2] += lz1+lz1;
                    square[start_parameter+p_index+0] += sx1;
                    square[start_parameter+p_index+1] += sx1;
                    square[start_parameter+p_index+2] += sx1;

                    linear[start_parameter+n1_index+0] += lx2+lx2;
                    linear[start_parameter+n1_index+1] += ly2+ly2;
                    linear[start_parameter+n1_index+2] += lz2+lz2;
                    square[start_parameter+n1_index+0] += sx2;
                    square[start_parameter+n1_index+1] += sx2;
                    square[start_parameter+n1_index+2] += sx2;

                    linear[start_parameter+n2_index+0] += lx3+lx3;
                    linear[start_parameter+n2_index+1] += ly3+ly3;
                    linear[start_parameter+n2_index+2] += lz3+lz3;
                    square[start_parameter+n2_index+0] += sx3;
                    square[start_parameter+n2_index+1] += sx3;
                    square[start_parameter+n2_index+2] += sx3;

                    lx1x2 += lx1x2;  // multiplication by 2.
                    lx1x3 += lx1x3;
                    lx2x3 += lx2x3;

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                                  cross_parms, cross_terms,
                                                  start_parameter + p_index + 0,
                                                  start_parameter + n1_index +0,
                                                  lx1x2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                                  cross_parms, cross_terms,
                                                  start_parameter + p_index + 1,
                                                  start_parameter + n1_index +1,
                                                  lx1x2, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                                  cross_parms, cross_terms,
                                                  start_parameter + p_index + 2,
                                                  start_parameter + n1_index +2,
                                                  lx1x2, 5 );

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                                  cross_parms, cross_terms,
                                                  start_parameter + p_index + 0,
                                                  start_parameter + n2_index +0,
                                                  lx1x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                                  cross_parms, cross_terms,
                                                  start_parameter + p_index + 1,
                                                  start_parameter + n2_index +1,
                                                  lx1x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                                  cross_parms, cross_terms,
                                                  start_parameter + p_index + 2,
                                                  start_parameter + n2_index +2,
                                                  lx1x3, 5 );

                    add_to_quadratic_cross_term_real( n_cross_terms,
                                                  cross_parms, cross_terms,
                                                  start_parameter + n1_index +0,
                                                  start_parameter + n2_index +0,
                                                  lx2x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                                  cross_parms, cross_terms,
                                                  start_parameter + n1_index +1,
                                                  start_parameter + n2_index +1,
                                                  lx2x3, 5 );
                    add_to_quadratic_cross_term_real( n_cross_terms,
                                                  cross_parms, cross_terms,
                                                  start_parameter + n1_index +2,
                                                  start_parameter + n2_index +2,
                                                  lx2x3, 5 );
                }
            }

            update_progress_report( &progress, node+1 );
        }

        terminate_progress_report( &progress );

        FREE( boundaries );
        FREE( found_flags );
    }

    FREE( neigh_points );
    FREE( areas );

    if( clip_search != NULL ) {
        delete_clip_search( clip_search );
    }

    *constant += cons;
}

public  void  find_boundary_points(
    Deform_struct               *deform,
    int                         n_parameters,
    int                         start_parameter[],
    Real                        parameters[],
    Smallest_int                active_flags[],
    Smallest_int                boundary_flags[],
    Point                       boundary_points[],
    Real                        *constant,
    Real                        linear[],
    Real                        square[],
    int                         *n_cross_terms[],
    int                         **cross_parms[],
    Real                        **cross_terms[] )
{
    int                         s, d, ind;
    surface_struct              *surface;
    surface_bound_struct        *bound;
    BOOLEAN                     save_points;

    zero_quadratic_real( n_parameters, constant, linear, square, *n_cross_terms,
                         *cross_parms, *cross_terms );

    ind = 0;
    for_less( s, 0, deform->n_surfaces )
    {
        surface = &deform->surfaces[s].surface;
        for_less( d, 0, deform->surfaces[s].n_bound )
        {
            bound = &deform->surfaces[s].bound[d];

            save_points = (bound->max_dist_weight > 0.0);

            find_image_boundaries(
                   bound->volume, bound->voxel_lookup,
                   &bound->done_bits, &bound->surface_bits,
                   bound->threshold,
                   bound->normal_direction,
                   bound->image_weight_out, bound->image_weight_in,
                   bound->max_outward, bound->max_inward,
                   bound->offset, bound->check_direction_flag,
                   bound->normal_direction_only,
                   bound->clip_to_surface,
                   bound->oversample,
                   surface->n_points,
                   surface->n_neighbours, surface->neighbours,
                   start_parameter[s],
                   &parameters[start_parameter[s]],
                   (active_flags == NULL) ? NULL :
                              &active_flags[start_parameter[s]/N_DIMENSIONS],
                   save_points ? &boundary_flags[ind] : NULL,
                   save_points ? &boundary_points[ind] : NULL,
                   constant, linear, square, n_cross_terms, cross_parms,
                   cross_terms );

            if( save_points ) {
                ind += surface->n_points +
                       bound->oversample * surface->n_edges +
                       bound->oversample * (bound->oversample-1) / 2 *
                       surface->n_polygons;
            }
        }
    }
}

public  void  print_fit_info(
    fit_eval_struct  *f )
{
    if( f->boundary_fit != 0.0 )
        print( " B:%.4g", f->boundary_fit );
    if( f->gradient_fit != 0.0 )
        print( " G:%.4g", f->gradient_fit );
    if( f->gw_gradient_fit != 0.0 )
        print( " GW-G:%.4g", f->gw_gradient_fit );
    if( f->value_fit != 0.0 )
        print( " V:%.4g", f->value_fit );
    if( f->stretch_fit != 0.0 )
        print( " S:%.4g", f->stretch_fit );
    if( f->curvature_fit != 0.0 )
        print( " C:%.4g", f->curvature_fit );
    if( f->bend_fit != 0.0 )
        print( " b:%.4g", f->bend_fit );

    if( f->self_intersect_fit != 0.0 )
        print( " I:%.4g", f->self_intersect_fit );
    if( f->closest_self_intersect == 0.0 )
        print( " #### self intersection ####\n" );
    else if( f->self_intersect_fit != 0.0 &&
             f->closest_self_intersect > 0.0 )
        print( " (%.4g)", f->closest_self_intersect );

    if( f->surf_surf_fit != 0.0 )
        print( " ss:%.4g", f->surf_surf_fit );
    if( f->closest_surf_surf == 0.0 )
        print( "#" );
    else if( f->surf_surf_fit != 0.0 &&
             f->closest_surf_surf > 0.0 )
        print( " (%.4g)", f->closest_surf_surf );

    if( f->inter_surface_fit != 0.0 )
        print( " 2:%.4g", f->inter_surface_fit );
    if( f->anchor_fit != 0.0 )
        print( " A:%.4g", f->anchor_fit );
    if( f->weight_point_fit != 0.0 )
        print( " w:%.4g", f->weight_point_fit );
    if( f->volume_fit != 0.0 )
      printf( " v:%.4g", f->volume_fit );
    if( f->laplacian_fit != 0.0 )
      printf( " L:%.4g", f->laplacian_fit );
}

public  void  print_fit_deriv_info(
    fit_eval_struct  *f )
{
    if( f->boundary_fit != 0.0 )
        print( " B:%.2g", f->boundary_fit );
    if( f->gradient_fit != 0.0 )
        print( " G:%.2g", f->gradient_fit );
    if( f->gw_gradient_fit != 0.0 )
        print( " GW-G:%.2g", f->gw_gradient_fit );
    if( f->value_fit != 0.0 )
        print( " V:%.2g", f->value_fit );
    if( f->stretch_fit != 0.0 )
        print( " S:%.2g", f->stretch_fit );
    if( f->curvature_fit != 0.0 )
        print( " C:%.2g", f->curvature_fit );
    if( f->bend_fit != 0.0 )
        print( " b:%.2g", f->bend_fit );

    if( f->self_intersect_fit != 0.0 )
        print( " I:%.2g", f->self_intersect_fit );
    if( f->closest_self_intersect == 0.0 )
        print( " #### self intersection ####\n" );
    else if( f->self_intersect_fit != 0.0 &&
             f->closest_self_intersect > 0.0 )
        print( " (%.4g)", f->closest_self_intersect );

    if( f->surf_surf_fit != 0.0 )
        print( " ss:%.2g", f->surf_surf_fit );
    if( f->closest_surf_surf == 0.0 )
        print( " #### surf surf intersection ####\n" );
    else if( f->surf_surf_fit != 0.0 &&
             f->closest_surf_surf > 0.0 )
        print( " (%.4g)", f->closest_surf_surf );

    if( f->inter_surface_fit != 0.0 )
        print( " 2:%.2g", f->inter_surface_fit );
    if( f->anchor_fit != 0.0 )
        print( " A:%.2g", f->anchor_fit );
    if( f->weight_point_fit != 0.0 )
        print( " w:%.2g", f->weight_point_fit );
    if( f->volume_fit != 0.0 )
      printf(" v:%.2g", f->volume_fit );
    if( f->laplacian_fit != 0.0 )
      printf(" L:%.2g", f->laplacian_fit );
}

