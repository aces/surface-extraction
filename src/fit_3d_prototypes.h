#ifndef  DEF_FIT_3D_PROTOTYPES
#define  DEF_FIT_3D_PROTOTYPES

public  Real  get_base_length(
    int    n_neighbours,
    Point  neighbours[],
    Point  *centroid );

public  void  get_neighbours_normal(
    int      n_points,
    Point    points[],
    Vector   *normal );

public  void  compute_centroid(
    int     n_points,
    Point   points[],
    Point   *centroid );

public  void  get_bending_xy(
    Real   p0[],
    Real   p1[],
    Real   p2[],
    Real   p3[],
    Real   *x,
    Real   *y );

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
    Real                          coefs[] );

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
    fit_eval_struct               *fit_info );

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
    fit_eval_struct               *fit_info );

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
    Real                        **cross_terms[] );

public  void  print_fit_info(
    fit_eval_struct  *f );

public  void  print_fit_deriv_info(
    fit_eval_struct  *f );

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
    Real                        *boundary_distance );

public  void  initialize_lookup_volume_coeficients(
    voxel_coef_struct  *lookup,
    Volume             volume );

public  void  lookup_volume_coeficients(
    voxel_coef_struct  *lookup,
    int                x,
    int                y,
    int                z,
    Real               coefs[] );

public  void  delete_lookup_volume_coeficients(
    voxel_coef_struct  *lookup );

public  conjugate   initialize_conjugate_gradient(
    int       n_parameters );

public  void   delete_conjugate_gradient(
    conjugate   con );

public  BOOLEAN  get_conjugate_unit_direction(
    conjugate   con,
    Real        derivative[],
    Real        unit_dir[],
    BOOLEAN     *changed );

public  Real  sq_triangle_triangle_dist(
    Real   a0[],
    Real   a1[],
    Real   a2[],
    Real   b0[],
    Real   b1[],
    Real   b2[],
    int    *which_case,
    Real   *sq_dist);

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
    Real   deriv_b2[] );

public  Real  sq_triangle_point_dist(
    Real   a0[],
    Real   a1[],
    Real   a2[],
    Real   point[] );

public  void  sq_triangle_point_dist_deriv(
    Real   a0[],
    Real   a1[],
    Real   a2[],
    Real   point[],
    Real   deriv_a0[],
    Real   deriv_a1[],
    Real   deriv_a2[],
    Real   deriv_point[] );

public  Real  recompute_surf_surfs(
    Deform_struct                 *deform,
    int                           start_parameter[],
    Real                          parameters[],
    Real                          max_movement,
    Real                          line_dir[],
    surf_surf_lookup_struct       *ss_lookup );

public  void  initialize_surf_surf_lookup(
    surf_surf_lookup_struct  *ss_lookup );

public  void  delete_surf_surf_lookup(
    surf_surf_lookup_struct  *ss_lookup );

public  int  get_n_surf_surf_candidate(
    surf_surf_lookup_struct  *ss_lookup );

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
    Real                     *dist_sq );

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
    surf_surf_deriv_struct        *deriv );

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
    surf_surf_deriv_struct        *deriv_info );

public  Real sq_triangle_triangle_dist_estimate(
    Real    a0[],
    Real    a1[],
    Real    a2[],
    Real    b0[],
    Real    b1[],
    Real    b2[],
    Real    search_distance_sq );

public  Real  recompute_self_intersects(
    Deform_struct                 *deform,
    int                           grid_size,
    int                           start_parameter[],
    Real                          parameters[],
    Real                          max_movement,
    Real                          line_dir[],
    self_intersect_lookup_struct  **si_lookup );

public  void  initialize_self_intersect_lookup(
    self_intersect_lookup_struct  *si_lookup );

public  void  delete_self_intersect_lookup(
    self_intersect_lookup_struct  *si_lookup );

public  int  get_n_self_intersect_candidate(
    self_intersect_lookup_struct  *si_lookup );

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
    Real                          *dist_sq );

public  void   create_self_intersect_deriv_info(
    self_intersect_lookup_struct  *si_lookup,
    BOOLEAN                       use_tri_tri_dist,
    int                           which,
    Real                          dist_sq,
    int                           n_points,
    Real                          parameters[],
    int                           n_neighbours[],
    int                           *neighbours[],
    self_intersect_deriv_struct   *deriv );

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
    self_intersect_deriv_struct   *deriv_info );

public  clip_struct  *initialize_clip_search(
    int    n_nodes,
    int    n_neighbours[],
    int    *neighbours[],
    Real   parameters[] );

public  void  clip_search_line(
    clip_struct  *clip,
    int          n_nodes_to_ignore,
    int          nodes_to_ignore[],
    Point        *origin,
    Vector       *normal,
    Real         outer_distance,
    Real         inner_distance,
    Real         *clipped_outer_distance,
    Real         *clipped_inner_distance );

public  void  delete_clip_search(
    clip_struct  *clip );

public  Real  get_max_movement_per_unit_line(
    int    n_parameters,
    Real   line_dir[] );
#endif
