#ifndef  DEF_FIT_PROTOTYPES
#define  DEF_FIT_PROTOTYPES

public   void   initialize_deform(
    Deform_struct   *deform );

public  void  get_surface_normal(
    Vector          *normal_2d,
    int             n_neighbours,
    Point           neighbours[],
    Vector          *normal );

public  Real  get_surface_curvature(
    Point    *point,
    Vector   *point_normal,
    int      n_neighbours,
    Point    neighbours[] );

public  Real  get_model_to_surface_scale(
    int                  n_points,
    Point                points[],
    model_info_struct    *model );

public  Status  add_surface_to_deform(
    Deform_struct               *deform,
    Volume                      volume,
    Volume                      label_volume,
    Surface_types               surface_type,
    boundary_definition_struct  *boundary,
    Real                        multi_scale_weight,
    int                         n_scales,
    multi_scale_struct          scales[],
    Real                        oversample_weight,
    int                         n_oversamples,
    Real                        max_outwards_distance,
    Real                        max_inwards_distance,
    Real                        data_weight,
    BOOLEAN                     threshold_per_vertex,
    Real                        threshold,
    object_struct               *threshold_object,
    Real                        curvature1,
    Real                        threshold1,
    Real                        curvature2,
    Real                        threshold2,
    Real                        max_diff,
    Normal_directions           normal_direction,
    object_struct               *object,
    int                         n_models,
    model_info_struct           models[],
    Deformation_model_types     model_types[],
    object_struct               *model_objects[] );

public  void  add_self_intersection_constraint(
    Deform_struct      *deform,
    int                surface_index,
    Real               min_distance1,
    Real               min_distance2,
    Real               weight );

public  void  add_surface_surface_constraint(
    Deform_struct      *deform,
    int                surface_index1,
    int                surface_index2,
    Real               min_distance1,
    Real               min_distance2,
    Real               weight );

public  void  add_dynamic_points_constraint(
    Deform_struct      *deform,
    int                surface_index1,
    int                surface_index2,
    Real               min_distance,
    Real               desired_distance,
    Real               max_distance,
    Real               weight );

public  void  add_list_of_points_constraint(
    Deform_struct      *deform,
    int                surface_index1,
    int                surface_index2,
    Real               min_distance_ratio,
    Real               max_distance_ratio,
    Real               weight,
    STRING             filename );

public  void  add_indiv_static_points_constraint(
    Deform_struct      *deform,
    int                surface_index,
    STRING             filename,
    Real               min_distance,
    Real               desired_distance,
    Real               max_distance,
    Real               weight );

public  void  add_static_points_constraint(
    Deform_struct      *deform,
    int                surface_index,
    STRING             filename,
    Real               min_distance,
    Real               desired_distance,
    Real               max_distance,
    Real               weight );

public  Status  get_deform_from_arguments(
    int            argc,
    char           *argv[],
    Deform_struct  *deform,
    STRING         **output_filenames );

public   void   delete_deform(
    Deform_struct   *deform,
    STRING          *output_filenames );

public  void  reset_recording_boundaries( void );

public  void  end_looking_up_boundaries( void );

public  void  start_looking_up_boundaries( void );

public  Real  get_linear_constraint_function(
    Real   value,
    Real   min_value,
    Real   max_value,
    Real   weight );

public  Real  get_linear_constraint_function(
    Real   value,
    Real   min_value,
    Real   max_value,
    Real   weight );

public  void  precompute_boundary_for_deform(
    Deform_struct          *deform,
    Point                  **surface_points,
    int                    surface,
    int                    point_index,
    Point                  *point,
    edge_info_struct       *boundary );

public  Real  evaluate_point_fit(
    Deform_struct          *deform,
    Point                  **points,
    uniform_subdiv_struct  subdivs[],
    int                    surface,
    int                    point,
    Vector                 *normal,
    BOOLEAN                edge_info_present,
    edge_info_struct       *edge_info );

public  Real  evaluate_fit(
    Deform_struct          *deform,
    Point                  **points,
    uniform_subdiv_struct  subdivs[] );

public  Real  evaluate_one_self_intersect(
    int                    surface,
    int                    point,
    surface_struct         surfaces[],
    uniform_subdiv_struct  subdivs[],
    Point                  **points,
    self_intersect_struct  *constraint );

public  Real  evaluate_all_self_intersect(
    surface_struct         surfaces[],
    uniform_subdiv_struct  subdivs[],
    Point                  **points,
    self_intersect_struct  *constraint );

public  Real  evaluate_one_surface_surface_intersect(
    int                     surface,
    int                     point,
    surface_struct          surfaces[],
    uniform_subdiv_struct   subdivs[],
    Point                   **points,
    surface_surface_struct  *constraint );

public  Real  evaluate_all_surface_surface_intersect(
    surface_struct          surfaces[],
    uniform_subdiv_struct   subdivs[],
    Point                   **points,
    surface_surface_struct  *constraint );

public  Real  evaluate_one_static_point_distance(
    int                     surface,
    int                     point,
    Point                   **points,
    static_points_struct    *constraint );

public  Real  evaluate_one_dynamic_point_distance(
    int                     surface,
    int                     point,
    Point                   **points,
    dynamic_points_struct   *constraint );

public  Real  evaluate_all_dynamic_point_distances(
    surface_struct          surfaces[],
    Point                   **points,
    dynamic_points_struct   *constraint );

public  Real  evaluate_all_static_point_distances(
    surface_struct          surfaces[],
    Point                   **points,
    static_points_struct    *constraint );

public  BOOLEAN  segments_closest_on_perpendicular(
    Point   *p1,
    Point   *p2,
    Point   *q1,
    Point   *q2,
    Real    *sq_dist );

public  Real  sq_distance_between_points(
    Point  *p1,
    Point  *p2 );

public  Real  point_segment_sq_distance(
    Point   *p,
    Point   *q1,
    Point   *q2 );

public  BOOLEAN    segments_within_distance(
    Point   *p1,
    Point   *p2,
    Point   *q1,
    Point   *q2,
    Real    min_distance,
    Real    *distance );

public  Real  sq_distance_between_polygons(
    polygons_struct  *p1,
    int              poly1,
    Point            points1[],
    polygons_struct  *p2,
    int              poly2,
    Point            points2[] );

public  BOOLEAN    polygons_within_distance(
    polygons_struct  *p1,
    int              poly1,
    Point            points1[],
    polygons_struct  *p2,
    int              poly2,
    Point            points2[],
    Real             min_distance,
    Real             *distance );

public  int  create_offsets(
    BOOLEAN     two_d,
    Grid_types  grid_type,
    Real        max_distance,
    int         grid_size,
    Real        step_ratio,
    BOOLEAN     one_d_flag,
    Vector      *offsets[] );

public  void  create_relaxed_models(
    Deform_struct       *deform,
    Real                relax_ratio,
    model_info_struct   ***models,
    constraint_struct   **constraints );

public  void  set_relaxed_models(
    Deform_struct       *deform,
    model_info_struct   **models,
    constraint_struct   constraints[] );

public  void  delete_relaxed_models(
    Deform_struct       *deform,
    model_info_struct   **models,
    constraint_struct   constraints[] );

public  Real  *save_deformation_weights(
    Deform_struct       *deform );

public  void  restore_deformation_weights(
    Deform_struct       *deform,
    Real                weights[],
    BOOLEAN             restore_data_flag,
    BOOLEAN             restore_stretch_flag,
    BOOLEAN             restore_curvature_flag,
    BOOLEAN             restore_constraints_flag );

public  void  zero_deformation_weights(
    Deform_struct       *deform );

public  int  get_points_to_evaluate(
    Deform_struct    *deform,
    int              n_points,
    int              s,
    int              points[],
    BOOLEAN          do_neighbours,
    int              *surfaces_to_evaluate[],
    int              *points_to_evaluate[],
    int              *models_to_evaluate[] );

public  Real  evaluate_set_of_points(
    Deform_struct          *deform,
    Point                  **points,
    uniform_subdiv_struct  subdivs[],
    int                    s,
    int                    n_points_to_evaluate,
    int                    surfaces_to_evaluate[],
    int                    points_to_evaluate[],
    int                    models_to_evaluate[],
    Vector                 *normal,
    BOOLEAN                use_cache,
    edge_info_struct       **static_cache,
    edge_info_struct       *this_point_cache,
    Real                   max_value_allowed );

public  int  get_n_multiscale_features(
    BOOLEAN   use_derivs[] );

public  void  delete_multi_scale(
    multi_scale_struct   *multi );

public  void  compute_multiscale_feature(
    multi_scale_struct   *multi,
    Point                *point,
    float                features[] );

public  Status  input_multi_scale(
    STRING              filename,
    int                 *n_scales,
    multi_scale_struct  *scales[] );

public  Real  get_multiscale_feature_diff(
    int                 n_scales,
    multi_scale_struct  scales[],
    float               feature1[],
    float               feature2[] );

public  Real  find_geodesic_distance(
    Volume                  volume,
    Real                    threshold,
    Real                    max_distance,
    Real                    voxel1[],
    Real                    voxel2[],
    Vector                  *search_normal );

public  void  initialize_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    Point                  *min_limits,
    Point                  *max_limits,
    int                    nx,
    int                    ny,
    int                    nz );

public  void  delete_uniform_subdiv(
    uniform_subdiv_struct  *subdiv );

public  void  add_point_to_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    Point                  *point,
    int                    point_index );

public  void  remove_point_from_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    Point                  *point,
    int                    point_index );

public  void  move_point_in_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    Point                  *prev_point,
    Point                  *new_point,
    int                    point_index );

public  void  create_uniform_subdiv(
    uniform_subdiv_struct  *subdiv,
    int                    n_points,
    Point                  points[],
    int                    n_cells );

public  int  get_subdiv_points_near_point(
    uniform_subdiv_struct  *subdiv,
    Point                  *point,
    Real                   distance,
    int                    points_list[] );
#endif
