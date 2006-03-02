#include  <fit_3d.h>
#include  <conjugate.h>

#define  DEFAULT_TOLERANCE           1.0e-4
#define  DEFAULT_FUNCTION_TOLERANCE  1.0e-4
#define  INITIAL_STEP           1.0e-4
#define  SEARCH_RATIO           5.0

private  Status  input_surface(
    STRING    filename,
    BOOLEAN   model_flag,
    int       *n_points,
    Point     *points[],
    int       *n_polygons,
    int       *n_neighbours[],
    int       **neighbours[] );

private  Status  input_surface_with_normals(
    STRING    filename,
    BOOLEAN   model_flag,
    int       *n_points,
    Point     *points[],
    Vector    *normals[],
    int       *n_polygons,
    int       *n_neighbours[],
    int       **neighbours[] );

private  Status  lookup_volume(
    STRING             filename,
    Volume             *volume,
    voxel_coef_struct  **voxel_lookup );

private  float   *create_model_lengths(
    int                n_points,
    int                n_neighbours[],
    int                *neighbours[],
    Point              model_points[],
    Real               stretch_scale,
    BOOLEAN            use_equal_lengths );

private  void   create_model_bends(
    int                n_points,
    int                n_neighbours[],
    int                *neighbours[],
    Point              model_points[],
    float              *model_x[],
    float              *model_y[] );

private  float   *create_curvatures(
    int                n_points,
    int                n_neighbours[],
    int                *neighbours[],
    Point              model_points[],
    Real               curvature_scale );

private  BOOLEAN   fit_polygons(
    Deform_struct      *deform,
    Real               tolerance,
    Real               function_tolerance,
    int                n_iters,
    int                recompute_every,
    Real               si_step,
    Real               max_step,
    Real               movement_threshold,
    int                n_movements,
    int                point_grid_size,
    int                n_deriv_smoothing_steps,
    Real               deriv_smoothing,
    BOOLEAN            print_deriv,
    BOOLEAN            print_closest,
    BOOLEAN            print_initial,
    STRING             node_dist_filename,
    BOOLEAN            timing_flag,
    Real               successful_movement_threshold );

private  int  count_edges(
    int    n_points,
    int    n_neighbours[],
    int    *neighbours[] );

private  void  delete_volume_lookup( void );
private  void  delete_surface_lookup( void );
private  void  delete_surface_lookup_model_points( void );

private  void  delete_surface_neighbours( void );

private  void  usage(
    char   executable_name[] )
{
    STRING  usage_format = "\
Usage:     %s  help not implemented yet \n\
                  \n\n";

    print_error( usage_format, executable_name );
}

static  int   n_surfaces_read = 0;
static  int  n_volumes_read = 0;
static  STRING log_filename = "";
static  int   n_recompute_intersect = 0;

private FILE *open_file_to_write_info(char *);
private int  write_file_info(FILE *, fit_eval_struct *);
private int  close_file_info(FILE *);

int  main(
    int    argc,
    char   *argv[] )
{
    STRING               arg, volume_filename, surface_direction;
    STRING               surf_surf_filename;
    STRING               model_filename, node_dist_filename;
    STRING               input_filename, output_filename;
    STRING               wm_obj_filename;
    STRING               *input_filenames, *output_filenames, filename;
    STRING               depth_filename;
    STRING               laplacian_filename, laplacian_gradient_filename;
    STRING               chamfer_filename;
    STRING               stretch_factor_name, curvature_factor_name;
    STRING               bend_factor_name, midpoint_filename;
    int                  n_bend, p1, p2, midpoint, n_weight_steps;
    int                  n_objects, s, i, n_stretch, n_curvature, point;
    int                  n_iters, recompute_every, ind, n_edges, j, best_index;
    int                  sizes[N_DIMENSIONS], n, neigh;
    int                  surface_index1, surface_index2;
    Real                 tolerance, min_offset, max_offset, x, y, z;
    Real                 min_offset1, min_offset2, max_offset1, max_offset2;
    Real                 function_tolerance, factor1, factor2;
    Real                 sum_surface_len, sum_model_len;
    Real                 max_step, si_step, factor, weight, min_distance;
    Real                 force_curvature, curvature_factor, stretch_factor;
    Real                 bend_factor;
    Real                 len, model_x, model_y, desired_distance;
    Real                 movement_threshold, successful_movement_threshold;
    Real                 value_differential_offset, value_differential_ratio;
    Real                 stretch_differential_offset, stretch_differential_ratio;
    Real                 deriv_smoothing;
    int                  n_deriv_smoothing_steps;
    int                  point_grid_size, p, n_movements;
    FILE                 *file, *file1;
    File_formats         format;
    object_struct        **object_list;
    polygons_struct      *polygons;
    Point                *points, **surface_points;
    Deform_struct        deform;
    one_surface_struct   surf;
    surface_struct       surface;
    surface_bound_struct bound;
    gradient_struct      gradient;
    surface_value_struct value;
    stretch_struct       stretch;
    bend_struct          bend;
    curvature_struct     curvature, *curv_ptr;
    inter_surface_struct inter;
    connection_struct    conn;
    anchor_struct        anchor;
    anchor_point_struct  a;
    weight_point_struct  weight_point;
    weight_struct        w;
    self_intersect_struct  *self;
    surf_surf_struct     *ss;
    intersect_wm_struct  intersect_wm;
    volume_info_struct   volume_s;
    laplacian_struct     laplacian;
    BOOLEAN              use_tri_tri_dist, use_equal_lengths;
    BOOLEAN              print_deriv, print_closest, print_initial;
    BOOLEAN              add_to_flag, timing_flag, did_all_iterations;
    STRING               two_or_three;
    int                  surface_total_number = 2;
    STRING               wm_masked_filename;

    n_surfaces_read = 0;
    n_volumes_read = 0;

    deform.n_surfaces = 0;
    deform.n_inter_surfaces = 0;
    deform.n_surf_surfs = 0;
    deform.n_intersect_wm = 0;

    timing_flag = FALSE;
    use_equal_lengths = FALSE;

    n_iters = 100;
    tolerance = DEFAULT_TOLERANCE;
    function_tolerance = DEFAULT_FUNCTION_TOLERANCE;
    max_step = -1.0;
    si_step = 0.5;
    point_grid_size = 30;
    recompute_every = 5;
    use_tri_tri_dist = TRUE;
    successful_movement_threshold = 0.0;

    input_filenames = NULL;
    output_filenames = NULL;

    print_deriv = getenv( "PRINT_DERIV" ) != NULL;
    print_closest = getenv( "PRINT_CLOSEST" ) != NULL;
    print_initial = getenv( "PRINT_INITIAL" ) != NULL;
    node_dist_filename = NULL;

    value_differential_ratio = 1.0e-8;
    value_differential_offset = 0.2;

    stretch_differential_ratio = 1.0e-8;
    stretch_differential_offset = 0.1;

    movement_threshold = -1.0;
    n_movements = 0;

    deriv_smoothing = 0.0;
    n_deriv_smoothing_steps = 0;

    initialize_argument_processing( argc, argv );

    while( get_string_argument( NULL, &arg ) )
    {
        if( equal_strings( arg, "-mode" ) )
        {
          if( !get_string_argument( NULL, &two_or_three ) )
          {
            print_error( "Error in -mode arguments.\n" );
            usage( argv[0] );
            return( 1 );
          }
          if( equal_strings( two_or_three, "TWO" ) || 
              equal_strings( two_or_three, "two" ) )
          {
            surface_total_number = 2;
          }
          else if( equal_strings( two_or_three, "THREE" ) || 
                equal_strings( two_or_three, "three" ) )
          {
            surface_total_number = 3;
          }
          // For the WM surface extraction
          else if( equal_strings( two_or_three, "WHITE" ) ||
                equal_strings( two_or_three, "white" ) )
          {
            surface_total_number = 1;
          }
          else
          {
            print_error( "Error in -mode arguments.\nDefine the number of surfaces (two or three)\n" );
            usage( argv[0] );
            return( 1 );
          }
        }
        else if( equal_strings( arg, "-surface" ) )
        {
          if( surface_total_number <= 2 ){
            if( !get_string_argument( NULL, &input_filename ) ||
                !get_string_argument( NULL, &output_filename ) )
            {
                print_error( "Error in -surface arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
          }
          else if( surface_total_number == 3 ){
            if( !get_string_argument( NULL, &input_filename ) ||
                !get_string_argument( NULL, &output_filename ) ||
                !get_string_argument( NULL, &surf_surf_filename ) )
            {
                print_error( "Error in -surface arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
          }
          else{
            print_error( "Error in -surface arguments.\n" );
            usage( argv[0] );
            return( 1 );
          }
            if( input_surface( input_filename, FALSE,
                               &surface.n_points, &surface.points,
                               &surface.n_polygons,
                               &surface.n_neighbours, &surface.neighbours )
                               != OK )
            {
                return( 1 );
            }

            surface.n_midpoints = 0;
            surf.surface = surface;
            surf.static_flag = FALSE;
            surf.n_bound = 0;
            surf.n_laplacian = 0;
            surf.n_volume = 0;
            surf.n_gradient = 0;
            surf.n_value = 0;
            surf.n_stretch = 0;
            surf.n_curvature = 0;
            surf.n_bend = 0;
            surf.n_anchors = 0;
            surf.n_weight_points = 0;
            surf.n_self_intersects = 0;

            ADD_ELEMENT_TO_ARRAY( deform.surfaces, deform.n_surfaces,
                                  surf, 1);
            --deform.n_surfaces;
            ADD_ELEMENT_TO_ARRAY( input_filenames, deform.n_surfaces,
                                  input_filename, 1 );
            --deform.n_surfaces;
            ADD_ELEMENT_TO_ARRAY( output_filenames, deform.n_surfaces,
                                  output_filename, 1 );

            if( surface_total_number == 3 )
            {
              if( input_surface( surf_surf_filename, FALSE,
                               &surface.n_points, &surface.points,
                               &surface.n_polygons,
                               &surface.n_neighbours, &surface.neighbours )
                  != OK )
              {
                return( 1 );
              }

              surface.n_midpoints = 0;
              surf.surface = surface;
              surf.static_flag = FALSE;
              surf.n_bound = 0;
              surf.n_laplacian = 0;
              surf.n_volume = 0;
              surf.n_gradient = 0;
              surf.n_value = 0;
              surf.n_stretch = 0;
              surf.n_curvature = 0;
              surf.n_bend = 0;
              surf.n_anchors = 0;
              surf.n_weight_points = 0;
              surf.n_self_intersects = 0;
              ADD_ELEMENT_TO_ARRAY( deform.surfaces, deform.n_surfaces,
                                    surf, 1);
            }
        }
	else if( equal_strings( arg, "-ban_wm" ) )
	{
	  if( !get_string_argument( NULL, &wm_obj_filename ) ||
          !get_string_argument( NULL, &wm_masked_filename ) ||
	      !get_real_argument( 0.0, &intersect_wm.weight) )
	  {
	    print_error( "Error in -ban_wm arguments.\n" );
	    usage( argv[0] );
	    return( 1 );
	  }

	  if( input_graphics_file( wm_obj_filename, &format,
          &intersect_wm.n_objects, &intersect_wm.objects ) != OK  || 
          intersect_wm.n_objects != 1 )
	  {
        print_error( "Error in -ban_wm.\n" );
	    return( 1 );
	  }

      if( lookup_volume( wm_masked_filename, &intersect_wm.wm_masked,
                         &intersect_wm.voxel_lookup ) != OK )
      {
        return( 1 );
      }

      SET_ARRAY_SIZE( intersect_wm.intersected, intersect_wm.n_intersected,
                      40960, DEFAULT_CHUNK_SIZE );
      intersect_wm.n_intersected = 0;
      intersect_wm.direction = (surface_total_number==1) ? -1 : 1;

      ADD_ELEMENT_TO_ARRAY( deform.intersect_wm, 
                            deform.n_intersect_wm, intersect_wm, 1 );

      //delete_object_list( n_objects, object_list );
	}
        else if( equal_strings( arg, "-static" ) )
        {
            deform.surfaces[deform.n_surfaces-1].static_flag = TRUE;
        }
        else if( equal_strings( arg, "-midpoints" ) )
        {
            if( deform.n_surfaces == 0 )
            {
                print_error( "No surface to specify midpoints for.\n" );
                return( 1 );
            }

            if( !get_string_argument( NULL, &midpoint_filename ) )
            {
                print_error( "Error in midpoint filename.\n" );
                return( 1 );
            }

            if( open_file( midpoint_filename, READ_FILE, ASCII_FORMAT, &file )
                != OK )
                return( 1 );

            while( input_int( file, &p1 ) == OK &&
                   input_int( file, &p2 ) == OK &&
                   input_int( file, &midpoint ) == OK )
            {
                ADD_ELEMENT_TO_ARRAY( deform.surfaces[deform.n_surfaces-1].
                                      surface.midpoints,
                                      deform.surfaces[deform.n_surfaces-1].
                                      surface.n_midpoints, p1,
                                      DEFAULT_CHUNK_SIZE );
                ADD_ELEMENT_TO_ARRAY( deform.surfaces[deform.n_surfaces-1].
                                      surface.midpoints,
                                      deform.surfaces[deform.n_surfaces-1].
                                      surface.n_midpoints, p2,
                                      DEFAULT_CHUNK_SIZE );
                ADD_ELEMENT_TO_ARRAY( deform.surfaces[deform.n_surfaces-1].
                                      surface.midpoints,
                                      deform.surfaces[deform.n_surfaces-1].
                                      surface.n_midpoints, midpoint,
                                      DEFAULT_CHUNK_SIZE );
            }

            close_file( file );
        }
        else if( equal_strings( arg, "-laplacian" ) )
        {
          if( !get_string_argument( NULL, &laplacian_filename ) ||
              !get_real_argument( 0.0, &laplacian.weight ) ||
              !get_real_argument( 0.0, &laplacian.from_value ) ||
              !get_real_argument( 0.0, &laplacian.to_value ) ||
              !get_real_argument( 0.0, &laplacian.deriv_factor ) ||
              !get_real_argument( 1.0, &laplacian.oversample ) )
          {
            print_error( "Error in -laplacian arguments.\n" );
            usage( argv[0] );
            return( 1 );
          }
          if( lookup_volume( laplacian_filename, &laplacian.volume,
                             &laplacian.voxel_lookup ) != OK )
          {
            return( 1 );
          }
          laplacian.direction = (surface_total_number==1) ? -1 : 1;
          ADD_ELEMENT_TO_ARRAY(
              deform.surfaces[0].laplacian,
              deform.surfaces[0].n_laplacian, laplacian, 1 );
        }
        else if( equal_strings( arg, "-volume" ) )
        {
          if( !get_real_argument( 0.0, &volume_s.weight ) ||
              !get_real_argument( 0.0, &volume_s.max_weight ) ||
              !get_real_argument( 0.0, &volume_s.adaptive_anchor_ratio ) ||
              !get_real_argument( 0.0, &volume_s.adaptive_boundary_ratio ) )
          {
            print_error( "Error in -volume arguments.\n" );
            usage( argv[0] );
            return( 1 );
          }
          volume_s.direction = (surface_total_number==1) ? -1 : 1;
          ADD_ELEMENT_TO_ARRAY( deform.surfaces[0/*deform.n_surfaces-1*/].
                                volume,
                                deform.surfaces[0/*deform.n_surfaces-1*/].
                                n_volume, volume_s, DEFAULT_CHUNK_SIZE );
        }
        else if( equal_strings( arg, "-boundary" ) )
        {
            if( !get_real_argument( 0.0, &bound.image_weight_out ) ||
                !get_real_argument( 0.0, &bound.image_weight_in ) ||
                !get_string_argument( NULL, &volume_filename ) ||
                !get_real_argument( 0.0, &bound.threshold ) ||
                !get_string_argument( NULL, &surface_direction ) ||
                !get_real_argument( 0.0, &bound.max_outward ) ||
                !get_real_argument( 0.0, &bound.max_inward ) ||
                !get_real_argument( 0.0, &bound.max_dist_weight ) ||
                !get_real_argument( 0.0, &bound.max_dist_threshold ) ||
                !get_int_argument( 0, &bound.oversample ) )
                //!get_string_argument( NULL, &depth_filename ) )
            {
                print_error( "Error in -boundary arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( lookup_volume( volume_filename, &bound.volume,
                               &bound.voxel_lookup ) != OK )
            {
              return( 1 );
            }

            if( surface_direction[0] == '-' )
                bound.normal_direction = TOWARDS_LOWER;
            else if( surface_direction[0] == '+' )
                bound.normal_direction = TOWARDS_HIGHER;
            else
                bound.normal_direction = ANY_DIRECTION;

            bound.offset = 0.0;
            bound.check_direction_flag = FALSE;
            bound.normal_direction_only = FALSE;
            bound.clip_to_surface = FALSE;

            get_volume_sizes( bound.volume, sizes );
            create_bitlist_3d( sizes[X], sizes[Y], sizes[Z], &bound.done_bits );
            create_bitlist_3d( sizes[X], sizes[Y], sizes[Z],
                               &bound.surface_bits );

            // Added by June Sic Kim at 6/12/2002
            /*if( open_file( depth_filename, READ_FILE, ASCII_FORMAT, &file ) != OK )
            {
                print_error( "Error opening file: %s\n", depth_filename );
                return( 1 );
            }
 
           while( input_real( file, &x ) == OK )
            {
              if( bound.max_weight_value < x )
              {
                bound.max_weight_value = x;
              }
                ADD_ELEMENT_TO_ARRAY( bound.weights,
                                      bound.n_weights,
                                      x, DEFAULT_CHUNK_SIZE );
            }

            (void) close_file( file );*/


            ADD_ELEMENT_TO_ARRAY(
                deform.surfaces[0/*deform.n_surfaces-1*/].bound,
                deform.surfaces[0/*deform.n_surfaces-1*/].n_bound, bound, 1 );
        }
        else if( equal_strings( arg, "-boundary_offset" ) )
        {
            if( deform.n_surfaces == 0 ||
                deform.surfaces[deform.n_surfaces-1].n_bound == 0 ||
                !get_real_argument( 0.0, &deform.surfaces[deform.n_surfaces-1].
                                    bound[deform.surfaces[deform.n_surfaces-1].
                                          n_bound-1].offset ) )
            {
                print_error( "Error in -boundary_offset.\n" );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-normal_only" ) )
        {
            if( deform.n_surfaces == 0 ||
                deform.surfaces[deform.n_surfaces-1].n_bound == 0 )
            {
                print_error( "Error in -normal_only.\n" );
                return( 1 );
            }
            deform.surfaces[deform.n_surfaces-1].
                                    bound[deform.surfaces[deform.n_surfaces-1].
                                      n_bound-1].normal_direction_only = TRUE;
        }
        else if( equal_strings( arg, "-clip_to_surface" ) )
        {
            if( deform.n_surfaces == 0 ||
                deform.surfaces[deform.n_surfaces-1].n_bound == 0 )
            {
                print_error( "Error in -clip_to_surface.\n" );
                return( 1 );
            }
            deform.surfaces[deform.n_surfaces-1].
                                    bound[deform.surfaces[deform.n_surfaces-1].
                                      n_bound-1].clip_to_surface = TRUE;
        }
        else if( equal_strings( arg, "-check_direction" ) )
        {
            if( deform.n_surfaces == 0 ||
                deform.surfaces[deform.n_surfaces-1].n_bound == 0 )
            {
                print_error( "Error in -check_direction.\n" );
                return( 1 );
            }

            deform.surfaces[deform.n_surfaces-1].
                        bound[deform.surfaces[deform.n_surfaces-1].n_bound-1].
                                check_direction_flag = TRUE;
        }
        else if( equal_strings( arg, "-gradient" ) )
        {
            if( !get_real_argument( 0.0, &gradient.image_weight ) ||
                !get_string_argument( NULL, &volume_filename ) ||
                !get_int_argument( 0, &gradient.continuity ) ||
                !get_real_argument( 0.0, &gradient.threshold ) ||
                !get_real_argument( 0.0, &gradient.min_diff ) ||
                !get_real_argument( 0.0, &gradient.max_diff ) ||
                !get_real_argument( 0.0, &factor ) ||
                !get_real_argument( 0.0, &value_differential_offset ) ||
                !get_real_argument( 0.0, &value_differential_ratio ) )
            {
                print_error( "Error in -gradient arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( gradient.image_weight > 0.0 )
                gradient.max_diff_weight = factor * gradient.image_weight;
            else
                gradient.max_diff_weight = factor;

            if( lookup_volume( volume_filename, &gradient.volume,
                               &gradient.voxel_lookup ) != OK )
            {
                return( 1 );
            }

            gradient.differential_offset = value_differential_offset;
            gradient.differential_ratio = value_differential_ratio;

            ADD_ELEMENT_TO_ARRAY(
                deform.surfaces[deform.n_surfaces-1].gradient,
                deform.surfaces[deform.n_surfaces-1].n_gradient, gradient, 1 );
        }
        else if( equal_strings( arg, "-value_differential" ) )
        {
            if( !get_real_argument( 0.0, &value_differential_offset ) ||
                !get_real_argument( 0.0, &value_differential_ratio ) )
            {
                print_error( "Error in -value_differential arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-value" ) )
        {
            if( !get_real_argument( 0.0, &value.image_weight ) ||
                !get_string_argument( NULL, &volume_filename ) ||
                !get_int_argument( 0, &value.continuity ) ||
                !get_real_argument( 0.0, &value.threshold ) ||
                !get_real_argument( 0.0, &value.min_diff ) ||
                !get_real_argument( 0.0, &value.max_diff ) ||
                !get_real_argument( 0.0, &factor ) ||
                !get_int_argument( 0, &value.oversample ) )
            {
                print_error( "Error in -value arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( value.image_weight > 0.0 )
                value.max_diff_weight = factor * value.image_weight;
            else
                value.max_diff_weight = factor;

            if( lookup_volume( volume_filename, &value.volume,
                               &value.voxel_lookup ) != OK )
            {
                return( 1 );
            }

            value.differential_offset = value_differential_offset;
            value.differential_ratio = value_differential_ratio;

            ADD_ELEMENT_TO_ARRAY(
                deform.surfaces[deform.n_surfaces-1].value,
                deform.surfaces[deform.n_surfaces-1].n_value, value, 1 );
        }
        else if( equal_strings( arg, "-equal_lengths" ) )
            use_equal_lengths = TRUE;
        else if( equal_strings( arg, "-noequal_lengths" ) )
            use_equal_lengths = FALSE;
        else if( equal_strings( arg, "-stretch_differential" ) )
        {
            if( !get_real_argument( 0.0, &stretch_differential_offset ) ||
                !get_real_argument( 0.0, &stretch_differential_ratio ) )
            {
                print_error( "Error in -stretch_differential arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-stretch" ) )
        {
            if( !get_real_argument( 0.0, &stretch.stretch_weight ) ||
                !get_string_argument( NULL, &model_filename ) ||
                !get_string_argument( NULL, &stretch_factor_name ) ||
                !get_real_argument( 0.0, &factor ) ||
                !get_real_argument( 0.0, &stretch.min_stretch ) ||
                !get_real_argument( 0.0, &stretch.max_stretch ) )
            {
                print_error( "Error in -stretch arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( stretch.stretch_weight <= 0.0 && factor <= 0.0 )
                continue;

            stretch.differential_offset = stretch_differential_offset;
            stretch.differential_ratio = stretch_differential_ratio;

            if( stretch_factor_name[0] == '+' )
            {
                if( sscanf( &stretch_factor_name[1], "%lf", &stretch_factor )
                    != 1 )
                {
                    print_error( "Error in -stretch arguments.\n" );
                    usage( argv[0] );
                    return( 1 );
                }
                add_to_flag = TRUE;
            }
            else
            {
                if( sscanf( stretch_factor_name, "%lf", &stretch_factor ) != 1)
                {
                    print_error( "Error in -stretch arguments.\n" );
                    usage( argv[0] );
                    return( 1 );
                }
                add_to_flag = FALSE;
            }

            if( stretch.stretch_weight > 0.0 )
                stretch.max_stretch_weight = factor * stretch.stretch_weight;
            else
                stretch.max_stretch_weight = factor;

#ifndef  USE_STRETCH_SQRT
            if( stretch.min_stretch > -1.0 )
                stretch.min_stretch = (1.0 + stretch.min_stretch) *
                                      (1.0 + stretch.min_stretch) - 1.0;
            else
                stretch.min_stretch = -1.0;

            stretch.max_stretch = (1.0 + stretch.max_stretch) *
                                  (1.0 + stretch.max_stretch) - 1.0;
#endif


            if( input_surface( model_filename, TRUE,
                               &stretch.n_points, &points, NULL,
                               &stretch.n_neighbours, &stretch.neighbours )
                               != OK )
            {
                return( 1 );
            }

            if( stretch.n_points >
                deform.surfaces[deform.n_surfaces-1].surface.n_points )
            {
              printf("stretch model=%d, surface=%d\n",stretch.n_points, 
                     deform.surfaces[deform.n_surfaces-1].surface.n_points);
                print_error( "Stretch model has more points than surface.\n" );
                return( 1 );
            }

            /*--- check for special case, where we want to assign a scale
                  automatically */

            if( stretch_factor < 0.0 )
            {
                sum_model_len = 0.0;
                sum_surface_len = 0.0;
                for_less( p, 0, stretch.n_points )
                {
                    for_less( n, 0, stretch.n_neighbours[p] )
                    {
                        neigh = stretch.neighbours[p][n];
                        sum_model_len += distance_between_points( &points[p],
                                                      &points[neigh] );
                        sum_surface_len += distance_between_points(
                                  &deform.surfaces[deform.n_surfaces-1].
                                   surface.points[p],
                                  &deform.surfaces[deform.n_surfaces-1].
                                   surface.points[neigh] );
                    }
                }

                stretch_factor = -stretch_factor * sum_surface_len /
                                 sum_model_len;
                print( "Surface %d:  Using model %s scaled by %g for stretch lengths\n",
                       deform.n_surfaces-1, model_filename, stretch_factor );
            }

            stretch.model_lengths = create_model_lengths( stretch.n_points,
                                                          stretch.n_neighbours,
                                                          stretch.neighbours,
                                                          points,
                                                          stretch_factor,
                                                          use_equal_lengths );

            if( add_to_flag )
            {
                n_stretch = deform.surfaces[deform.n_surfaces-1].n_stretch;
                if( deform.surfaces[deform.n_surfaces-1].n_stretch == 0 ||
                    stretch.n_points !=
                    deform.surfaces[deform.n_surfaces-1].stretch[n_stretch-1].
                         n_points )
                {
                    print_error( "Error in stretch argument.\n" );
                    return( 1 );
                }

                n_edges = count_edges( stretch.n_points, stretch.n_neighbours,
                                       stretch.neighbours );
                for_less( ind, 0, n_edges )
                {
                    deform.surfaces[deform.n_surfaces-1].stretch[n_stretch-1].
                         model_lengths[ind] += stretch.model_lengths[ind];
                }

                FREE( stretch.model_lengths );
            }
            else
            {
                ADD_ELEMENT_TO_ARRAY(
                      deform.surfaces[0/*deform.n_surfaces-1*/].stretch,
                      deform.surfaces[0/*deform.n_surfaces-1*/].n_stretch, stretch, 1);
            }
        }
        else if( equal_strings( arg, "-curvature" ) )
        {
            if( !get_real_argument( 0.0, &curvature.curvature_weight ) ||
                !get_string_argument( NULL, &model_filename ) ||
                !get_string_argument( NULL, &curvature_factor_name ) ||
                !get_real_argument( 0.0, &factor ) ||
                !get_real_argument( 0.0, &curvature.min_curvature ) ||
                !get_real_argument( 0.0, &curvature.max_curvature ) )
            {
                print_error( "Error in -curvature arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( curvature.curvature_weight <= 0.0 && factor <= 0.0 )
                continue;

            if( curvature_factor_name[0] == '+' )
            {
                if( sscanf( &curvature_factor_name[1], "%lf",
                             &curvature_factor ) != 1 )
                {
                    print_error( "Error in -curvature arguments.\n" );
                    usage( argv[0] );
                    return( 1 );
                }
                add_to_flag = TRUE;
            }
            else
            {
                if( sscanf( curvature_factor_name, "%lf", &curvature_factor )
                                                != 1 )
                {
                    print_error( "Error in -curvature arguments.\n" );
                    usage( argv[0] );
                    return( 1 );
                }
                add_to_flag = FALSE;
            }

            if( curvature.curvature_weight > 0.0 )
            {
                curvature.max_curvature_weight = factor *
                                   curvature.curvature_weight;
            }
            else
                curvature.max_curvature_weight = factor;

            if( input_surface( model_filename, TRUE,
                               &curvature.n_points, &points, NULL,
                               &curvature.n_neighbours, &curvature.neighbours )
                               != OK )
            {
                return( 1 );
            }

            if( curvature.n_points >
                deform.surfaces[deform.n_surfaces-1].surface.n_points )
            {
                print_error( "Curvature model has more points than surface.\n");
                return( 1 );
            }

            curvature.curvatures = create_curvatures( curvature.n_points,
                                                      curvature.n_neighbours,
                                                      curvature.neighbours,
                                                      points,
                                                      curvature_factor );
            if( add_to_flag )
            {
                n_curvature = deform.surfaces[deform.n_surfaces-1].n_curvature;
                if( n_curvature == 0 || curvature.n_points !=
                    deform.surfaces[deform.n_surfaces-1].curvature
                                       [n_curvature-1].  n_points )
                {
                    print_error( "Error in curvature argument.\n" );
                    return( 1 );
                }

                for_less( point, 0, curvature.n_points )
                {
                    deform.surfaces[deform.n_surfaces-1].
                         curvature[n_curvature-1].
                         curvatures[point] += curvature.curvatures[point];
                }

                FREE( curvature.curvatures );
            }
            else
            {
                ADD_ELEMENT_TO_ARRAY(
                    deform.surfaces[deform.n_surfaces-1].curvature,
                    deform.surfaces[deform.n_surfaces-1].n_curvature,
                    curvature, 1 );
            }
        }
        else if( equal_strings( arg, "-bend" ) )
        {
            if( !get_real_argument( 0.0, &bend.bend_weight ) ||
                !get_string_argument( NULL, &model_filename ) ||
                !get_string_argument( NULL, &bend_factor_name ) ||
                !get_real_argument( 0.0, &factor ) ||
                !get_real_argument( 0.0, &bend.min_bend ) ||
                !get_real_argument( 0.0, &bend.max_bend ) )
            {
                print_error( "Error in -bend arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( bend_factor_name[0] == '+' )
            {
                if( sscanf( &bend_factor_name[1], "%lf", &bend_factor ) != 1 )
                {
                    print_error( "Error in -bend arguments.\n" );
                    usage( argv[0] );
                    return( 1 );
                }
                add_to_flag = TRUE;
            }
            else
            {
                if( sscanf( bend_factor_name, "%lf", &bend_factor ) != 1)
                {
                    print_error( "Error in -bend arguments.\n" );
                    usage( argv[0] );
                    return( 1 );
                }
                add_to_flag = FALSE;
            }

            if( bend.bend_weight > 0.0 )
                bend.max_bend_weight = factor * bend.bend_weight;
            else
                bend.max_bend_weight = factor;

            if( input_surface( model_filename, TRUE,
                               &bend.n_points, &points, NULL,
                               &bend.n_neighbours, &bend.neighbours )
                               != OK )
            {
                return( 1 );
            }

            if( bend.n_points >
                deform.surfaces[deform.n_surfaces-1].surface.n_points )
            {
                print_error( "Bend model has more points than surface.\n" );
                return( 1 );
            }

            create_model_bends( bend.n_points, bend.n_neighbours,
                                bend.neighbours, points,
                                &bend.model_x, &bend.model_y );

            n_edges = count_edges( bend.n_points, bend.n_neighbours,
                                   bend.neighbours );
            for_less( ind, 0, n_edges )
            {
                if( bend_factor == 0.0 )
                {
                    bend.model_x[ind] = (float) -1.0;
                    bend.model_y[ind] = (float) 0.0;
                }
                else
                {
                    bend.model_x[ind] = (float) ((Real) bend.model_x[ind] *
                                                 bend_factor);
                    bend.model_y[ind] = (float) ((Real) bend.model_y[ind] *
                                                 bend_factor);
                }
            }

            if( add_to_flag )
            {
                n_bend = deform.surfaces[deform.n_surfaces-1].n_bend;
                if( deform.surfaces[deform.n_surfaces-1].n_bend == 0 ||
                    bend.n_points !=
                    deform.surfaces[deform.n_surfaces-1].bend[n_bend-1].
                         n_points )
                {
                    print_error( "Error in bend argument.\n" );
                    return( 1 );
                }

                for_less( ind, 0, n_edges )
                {
                    model_x = (Real) deform.surfaces[deform.n_surfaces-1].
                              bend[n_bend-1].model_x[ind] +
                              (Real) bend.model_x[ind];
                    model_y = (Real) deform.surfaces[deform.n_surfaces-1].
                              bend[n_bend-1].model_y[ind] +
                              (Real) bend.model_y[ind];

                    len = sqrt( model_x * model_x + model_y * model_y );

                    deform.surfaces[deform.n_surfaces-1].bend[n_bend-1].
                            model_x[ind] = (float) (model_x / len);
                    deform.surfaces[deform.n_surfaces-1].bend[n_bend-1].
                            model_y[ind] = (float) (model_y / len);
                }

                FREE( bend.model_x );
                FREE( bend.model_y );
            }
            else
            {
                ADD_ELEMENT_TO_ARRAY(
                    deform.surfaces[deform.n_surfaces-1].bend,
                    deform.surfaces[deform.n_surfaces-1].n_bend, bend, 1);
            }
        }
        else if( equal_strings( arg, "-force_curvature" ) )
        {
            if( !get_real_argument( 0.0, &force_curvature ) ||
                deform.n_surfaces == 0 ||
                deform.surfaces[deform.n_surfaces-1].n_curvature == 0 )
            {
                print_error( "Error in -force_curvature arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
            curv_ptr = &deform.surfaces[deform.n_surfaces-1].curvature[
                        deform.surfaces[deform.n_surfaces-1].n_curvature-1];
            for_less( i, 0, curv_ptr->n_points )
                curv_ptr->curvatures[i] = (float) force_curvature;
        }
        else if( equal_strings( arg, "-inter" ) )
        {
            if( !get_real_argument( 0.0, &inter.weight ) ||
                !get_int_argument( 0, &inter.surface_index1 ) ||
                !get_int_argument( 0, &inter.surface_index2 ) ||
                !get_string_argument( NULL, &filename ) ||
                !get_real_argument( 0.0, &factor1 ) ||
                !get_real_argument( 0.0, &factor2 ) ||
                !get_int_argument( 0, &n_weight_steps ) ||
                !get_real_argument( 0.0, &min_offset2 ) ||
                !get_real_argument( 0.0, &min_offset1 ) ||
                !get_real_argument( 0.0, &max_offset1 ) ||
                !get_real_argument( 0.0, &max_offset2 ) )
            {
                print_error( "Error in -inter arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            inter.oversample = 0;
            inter.n_weight_steps = n_weight_steps;

            if( inter.weight > 0.0 )
            {
                inter.max_weight1 = factor1 * inter.weight;
                inter.max_weight2 = factor2 * inter.weight;
            }
            else
            {
                inter.max_weight1 = factor1;
                inter.max_weight2 = factor2;
            }

            if( inter.surface_index1 < 0 ||
                inter.surface_index2 < 0 ||
                inter.surface_index1 >= deform.n_surfaces ||
                inter.surface_index2 >= deform.n_surfaces )
            {
                print_error( "Surface index out of range in -inter\n" );
                usage( argv[0] );
                return( 1 );
            }

            inter.n_connections = 0;

            if( open_file( filename, READ_FILE, ASCII_FORMAT, &file ) != OK )
            {
                print_error( "Error opening file: %s\n", filename );
                return( 1 );
            }

            while( input_int( file, &conn.surface1_point ) == OK )
            {
                if( input_int( file, &conn.surface2_point ) != OK ||
                    input_real( file, &conn.desired_distance ) != OK )
                {
                    print_error( "Error reading file: %s\n", filename );
                    return( 1 );
                }

                conn.min_distance1 = min_offset1;
                conn.min_distance2 = min_offset2;
                conn.max_distance1 = max_offset1;
                conn.max_distance2 = max_offset2;

                if( conn.surface1_point < 0 || conn.surface1_point >=
                    deform.surfaces[inter.surface_index1].surface.n_points ||
                    conn.surface2_point < 0 || conn.surface2_point >=
                    deform.surfaces[inter.surface_index2].surface.n_points )
                {
                    print_error( "Point index out of range in: %s\n", filename);
                    return( 1 );
                }

                ADD_ELEMENT_TO_ARRAY( inter.connections, inter.n_connections,
                                      conn, DEFAULT_CHUNK_SIZE );
            }

            (void) close_file( file );

            ADD_ELEMENT_TO_ARRAY( deform.inter_surfaces,
                                  deform.n_inter_surfaces, inter,
                                  1 );
        }
        else if( equal_strings( arg, "-two" ) )
        {
            if( !get_real_argument( 0.0, &inter.weight ) ||
                !get_int_argument( 0, &inter.surface_index1 ) ||
                !get_int_argument( 0, &inter.surface_index2 ) ||
                !get_real_argument( 0.0, &desired_distance ) ||
                !get_real_argument( 0.0, &factor1 ) ||
                !get_real_argument( 0.0, &factor2 ) ||
                !get_int_argument( 0, &inter.n_weight_steps ) ||
                !get_real_argument( 0.0, &min_offset2 ) ||
                !get_real_argument( 0.0, &min_offset1 ) ||
                !get_real_argument( 0.0, &max_offset1 ) ||
                !get_real_argument( 0.0, &max_offset2 ) ||
                !get_int_argument( 0, &inter.oversample ) )
            {
                print_error( "Error in -two arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( inter.weight > 0.0 )
            {
                inter.max_weight1 = factor1 * inter.weight;
                inter.max_weight2 = factor2 * inter.weight;
            }
            else
            {
                inter.max_weight1 = factor1;
                inter.max_weight2 = factor2;
            }

            if( inter.surface_index1 < 0 ||
                inter.surface_index2 < 0 ||
                inter.surface_index1 >= deform.n_surfaces ||
                inter.surface_index2 >= deform.n_surfaces )
            {
                print_error( "Surface index out of range in -inter\n" );
                usage( argv[0] );
                return( 1 );
            }

            inter.n_connections = MIN( deform.surfaces[inter.surface_index1].
                                          surface.n_points,
                                       deform.surfaces[inter.surface_index2].
                                          surface.n_points );

            if( deform.surfaces[inter.surface_index1].surface.n_points !=
                deform.surfaces[inter.surface_index2].surface.n_points &&
                inter.oversample > 0 )
            {
                print_error( "Warning, setting oversample to 0 on -two\n" );
                inter.oversample = 0;
            }

            ALLOC( inter.connections, inter.n_connections );

            for_less( p, 0, inter.n_connections )
            {
                inter.connections[p].surface1_point = p;
                inter.connections[p].surface2_point = p;
                inter.connections[p].desired_distance = desired_distance;
                inter.connections[p].min_distance1 = min_offset1;
                inter.connections[p].min_distance2 = min_offset2;
                inter.connections[p].max_distance1 = max_offset1;
                inter.connections[p].max_distance2 = max_offset2;
            }

            ADD_ELEMENT_TO_ARRAY( deform.inter_surfaces,
                                  deform.n_inter_surfaces, inter, 1 );
        }
        else if( equal_strings( arg, "-anchor" ) )
        {
            if( !get_real_argument( 0.0, &anchor.weight ) ||
                !get_string_argument( NULL, &filename ) ||
                !get_real_argument( 0.0, &factor ) ||
                !get_real_argument( 0.0, &min_offset ) ||
                !get_real_argument( 0.0, &max_offset ) )
                //!get_string_argument( NULL, &depth_filename ) )
            {
                print_error( "Error in -anchor arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( anchor.weight > 0.0 )
                anchor.max_dist_weight = factor * anchor.weight;
            else
                anchor.max_dist_weight = factor;

            anchor.n_anchor_points = 0;

            if( open_file( filename, READ_FILE, ASCII_FORMAT, &file ) != OK )
            {
                print_error( "Error opening file: %s\n", filename );
                return( 1 );
            }

            // Added by June Sic Kim at 5/12/2002
            /*if( open_file( depth_filename, READ_FILE, ASCII_FORMAT, &file1 ) != OK )
            {
                print_error( "Error opening file: %s\n", filename );
                return( 1 );
            }*/
 
           while( input_int( file, &a.surface_point ) == OK )
            {
                if( input_real( file, &x ) != OK ||
                    input_real( file, &y ) != OK ||
                    input_real( file, &z ) != OK ||
                    input_real( file, &a.desired_distance ) != OK )
                {
                    print_error( "Error reading file: %s\n", filename );
                    return( 1 );
                }

                fill_Point( a.anchor_point, x, y, z );

                a.min_distance = min_offset;
                a.max_distance = max_offset;

                if( a.surface_point < 0 || a.surface_point >=
                    deform.surfaces[deform.n_surfaces-1].surface.n_points )
                {
                    print_error( "Point index (%d) out of range in: %s\n",
                                 a.surface_point, filename);
                    return( 1 );
                }

                // Added by June Sic Kim at 6/12/2002
                /*if( input_real( file1, &x ) != OK )
                {
                  print_error( "Error reading file: %s\n", depth_filename);
                  return( 1 );
                }
                if( anchor.max_weight_value < x )
                {
                  anchor.max_weight_value = x;
                }
                a.weight = x;
                a.max_weight = x;*/

                ADD_ELEMENT_TO_ARRAY( anchor.anchor_points,
                                      anchor.n_anchor_points,
                                      a, DEFAULT_CHUNK_SIZE );
            }

            (void) close_file( file );
            //(void) close_file( file1 );

            ADD_ELEMENT_TO_ARRAY( deform.surfaces[0/*deform.n_surfaces-1*/].anchors,
                                  deform.surfaces[0/*deform.n_surfaces-1*/].n_anchors,
                                 anchor, 1 );
        }
        else if( equal_strings( arg, "-weight_point" ) )
        {
            if( !get_real_argument( 0.0, &weight_point.weight ) ||
                !get_string_argument( NULL, &filename ) ||
                !get_real_argument( 0.0, &factor ) ||
                !get_real_argument( 0.0, &min_offset ) ||
                !get_real_argument( 0.0, &max_offset ) )
            {
                print_error( "Error in -weight_point arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( weight_point.weight > 0.0 )
                weight_point.max_dist_weight = factor * weight_point.weight;
            else
                weight_point.max_dist_weight = factor;

            weight_point.n_weight_points = 0;

            if( open_file( filename, READ_FILE, ASCII_FORMAT, &file ) != OK )
            {
                print_error( "Error opening file: %s\n", filename );
                return( 1 );
            }

            while( input_int( file, &w.n_surface_points ) == OK )
            {
                ALLOC( w.surface_points, w.n_surface_points );
                ALLOC( w.surface_weights, w.n_surface_points );

                for_less( i, 0, w.n_surface_points )
                {
                    if( input_int( file, &w.surface_points[i] ) != OK ||
                        input_real( file, &w.surface_weights[i] ) != OK ||
                        w.surface_points[i] < 0 || w.surface_points[i] >=
                        deform.surfaces[deform.n_surfaces-1].surface.n_points )
                    {
                        print_error( "Error in surface points file.\n" );
                        return( 1 );
                    }
                }
                    
                if( input_real( file, &x ) != OK ||
                    input_real( file, &y ) != OK ||
                    input_real( file, &z ) != OK ||
                    input_real( file, &w.desired_distance ) != OK )
                {
                    print_error( "Error reading file: %s\n", filename );
                    return( 1 );
                }

                fill_Point( w.anchor_point, x, y, z );

                w.min_distance = min_offset;
                w.max_distance = max_offset;

                ADD_ELEMENT_TO_ARRAY( weight_point.weight_points,
                                      weight_point.n_weight_points,
                                      w, DEFAULT_CHUNK_SIZE );
            }

            (void) close_file( file );

            ADD_ELEMENT_TO_ARRAY( deform.surfaces[deform.n_surfaces-1].
                                  weight_points,
                          deform.surfaces[deform.n_surfaces-1].n_weight_points,
                          weight_point, 1 );
        }
        else if( equal_strings( arg, "-self_intersect" ) )
        {
            if( !get_real_argument( 0.0, &weight ) ||
                !get_real_argument( 0.0, &min_distance ) )
            {
                print_error( "Error in -self_intersect arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( deform.surfaces[0/*deform.n_surfaces-1*/].n_self_intersects == 0 )
            {
              ALLOC( deform.surfaces[0/*deform.n_surfaces-1*/].self_intersects, 1);
                deform.surfaces[0/*deform.n_surfaces-1*/].n_self_intersects = 1;

                deform.surfaces[0/*deform.n_surfaces-1*/].self_intersects[0].
                                     n_weights = 0;
                deform.surfaces[0/*deform.n_surfaces-1*/].self_intersects[0].
                                     use_tri_tri_dist = use_tri_tri_dist;
            }

            self = &deform.surfaces[0/*deform.n_surfaces-1*/].self_intersects[0];

            self->square_flag = TRUE;
            if( weight < 0.0 )
            {
                self->square_flag = FALSE;
                weight = -weight;
            }
            SET_ARRAY_SIZE( self->weights,
                            self->n_weights, self->n_weights + 1, 1 );
            SET_ARRAY_SIZE( self->min_distances,
                            self->n_weights, self->n_weights + 1, 1 );
            self->weights[self->n_weights] = weight;
            self->min_distances[self->n_weights] = min_distance;
            ++self->n_weights;

            for_less( i, 0, self->n_weights-1 )
            {
                Real  swap;

                best_index = i;
                for_less( j, i+1, self->n_weights )
                {
                    if( self->min_distances[j] < self->min_distances[i] )
                        best_index = j;
                }

                swap = self->min_distances[best_index];
                self->min_distances[best_index] = self->min_distances[i];
                self->min_distances[i] = swap;

                swap = self->weights[best_index];
                self->weights[best_index] = self->weights[i];
                self->weights[i] = swap;
            }
        }
        else if( equal_strings( arg, "-surf_surf" ) )
        {
            if( !get_int_argument( 0, &surface_index1 ) ||
                !get_int_argument( 0, &surface_index2 ) ||
                !get_real_argument( 0.0, &weight ) ||
                !get_real_argument( 0.0, &min_distance ) )
            {
                print_error( "Error in -surf_surf arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }

            if( deform.n_surf_surfs == 0 ||
                deform.surf_surfs[deform.n_surf_surfs-1].surface_index1 !=
                                                         surface_index1 ||
                deform.surf_surfs[deform.n_surf_surfs-1].surface_index2 !=
                                                         surface_index2 )
            {
                SET_ARRAY_SIZE( deform.surf_surfs, deform.n_surf_surfs,
                                deform.n_surf_surfs+1, DEFAULT_CHUNK_SIZE );

                deform.surf_surfs[deform.n_surf_surfs].surface_index1 =
                                                   surface_index1;
                deform.surf_surfs[deform.n_surf_surfs].surface_index2 =
                                                   surface_index2;
                deform.surf_surfs[deform.n_surf_surfs].n_weights = 0;

                ++deform.n_surf_surfs;
            }

            ss = &deform.surf_surfs[deform.n_surf_surfs-1];

            SET_ARRAY_SIZE( ss->weights,
                            ss->n_weights, ss->n_weights + 1, 1 );
            SET_ARRAY_SIZE( ss->min_distances,
                            ss->n_weights, ss->n_weights + 1, 1 );
            ss->weights[ss->n_weights] = weight;
            ss->min_distances[ss->n_weights] = min_distance;
            ++ss->n_weights;

        }
        else if( equal_strings( arg, "-fitting" ) )
        {
            if( !get_int_argument( 0, &n_iters ) ||
                !get_int_argument( 0, &recompute_every ) ||
                !get_real_argument( 0.0, &tolerance ) )
            {
                print_error( "Error in -fitting arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-deriv_smoothing" ) )
        {
            if( !get_int_argument( 0, &n_deriv_smoothing_steps ) ||
                !get_real_argument( 0.0, &deriv_smoothing ) )
            {
                print_error( "Error in -deriv_smoothing arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-ftol" ) )
        {
            if( !get_real_argument( 0.0, &function_tolerance ) )
            {
                print_error( "Error in -ftol arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-max_step" ) )
        {
            if( !get_real_argument( 0.0, &max_step ) )
            {
                print_error( "Error in -max_step arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-step" ) )
        {
            if( !get_real_argument( 0.0, &si_step ) )
            {
                print_error( "Error in -step arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-grid_size" ) )
        {
            if( !get_int_argument( 0, &point_grid_size ) )
            {
                print_error( "Error in -grid_size arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-stop" ) )
        {
            if( !get_real_argument( 0.0, &movement_threshold ) ||
                !get_int_argument( 0, &n_movements ) )
            {
                print_error( "Error in -stop arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-movement_success" ) )
        {
            if( !get_real_argument( 0.0, &successful_movement_threshold ) )
            {
                print_error( "Error in -movement_success arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-output_only" ) )
        {
            if( !get_string_argument( 0, &node_dist_filename ) )
            {
                print_error( "Error in -output_only arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else if( equal_strings( arg, "-tri" ) )
            use_tri_tri_dist = TRUE;
        else if( equal_strings( arg, "-point" ) )
            use_tri_tri_dist = FALSE;
        else if( equal_strings( arg, "-print_deriv" ) )
            print_deriv = TRUE;
        else if( equal_strings( arg, "-noprint_deriv" ) )
            print_deriv = FALSE;
        else if( equal_strings( arg, "-print_initial" ) )
            print_initial = TRUE;
        else if( equal_strings( arg, "-noprint_initial" ) )
            print_initial = FALSE;
        else if( equal_strings( arg, "-print_closest" ) )
            print_closest = TRUE;
        else if( equal_strings( arg, "-noprint_closest" ) )
            print_closest = FALSE;
        else if( equal_strings( arg, "-timing" ) )
            timing_flag = TRUE;
        else if( equal_strings( arg, "-log" ) )
        {
            if( !get_string_argument( 0, &log_filename ) )
            {
                print_error( "Error in -log arguments.\n" );
                usage( argv[0] );
                return( 1 );
            }
        }
        else
        {

            print_error( "Unrecognized argument: %s.\n", arg );
            usage( argv[0] );
            return( 1 );
        }
    }

    delete_surface_lookup_model_points();

    did_all_iterations = fit_polygons( &deform, tolerance, function_tolerance,
                  n_iters,
                  recompute_every,
                  si_step, max_step, movement_threshold, n_movements,
                  point_grid_size,
                  n_deriv_smoothing_steps, deriv_smoothing,
                  print_deriv, print_closest,
                  print_initial, node_dist_filename, timing_flag,
                  successful_movement_threshold );

    delete_volume_lookup();
    delete_surface_neighbours();

    ALLOC( surface_points, deform.n_surfaces );
    for_less( s, 0, deform.n_surfaces )
        surface_points[s] = deform.surfaces[s].surface.points;

    // modified by June at 14/08/2003
    for_less( s, 0, deform.n_surfaces>0?1:0 )
    {
        for_less( i, 0, deform.surfaces[s].n_bound )
        {
            delete_bitlist_3d( &deform.surfaces[s].bound[i].done_bits );
            delete_bitlist_3d( &deform.surfaces[s].bound[i].surface_bits );
        }

        if( deform.surfaces[s].n_bound > 0 )
            FREE( deform.surfaces[s].bound );

        if( deform.surfaces[s].n_value > 0 )
            FREE( deform.surfaces[s].value );

        for_less( i, 0, deform.surfaces[s].n_stretch )
        {
            FREE( deform.surfaces[s].stretch[i].model_lengths );
        }

        if( deform.surfaces[s].n_stretch > 0 )
            FREE( deform.surfaces[s].stretch );

        for_less( i, 0, deform.surfaces[s].n_curvature )
        {
            FREE( deform.surfaces[s].curvature[i].curvatures );
        }

        if( deform.surfaces[s].n_curvature > 0 )
            FREE( deform.surfaces[s].curvature );

        for_less( i, 0, deform.surfaces[s].n_bend )
        {
            FREE( deform.surfaces[s].bend[i].model_x );
            FREE( deform.surfaces[s].bend[i].model_y );
        }

        if( deform.surfaces[s].n_bend > 0 )
            FREE( deform.surfaces[s].bend );

        for_less( i, 0, deform.surfaces[s].n_anchors )
        {
            if( deform.surfaces[s].anchors[i].n_anchor_points > 0 )
                FREE( deform.surfaces[s].anchors[i].anchor_points );
        }

        if( deform.surfaces[s].n_anchors > 0 )
            FREE( deform.surfaces[s].anchors );

        for_less( i, 0, deform.surfaces[s].n_weight_points )
        {
            for_less( p, 0, deform.surfaces[s].weight_points[i].n_weight_points)
            {
                FREE( deform.surfaces[s].weight_points[i].
                      weight_points[p].surface_points );
                FREE( deform.surfaces[s].weight_points[i].
                      weight_points[p].surface_weights );
            }

            if( deform.surfaces[s].weight_points[i].n_weight_points > 0 )
                FREE( deform.surfaces[s].weight_points[i].weight_points );
        }

        if( deform.surfaces[s].n_weight_points > 0 )
            FREE( deform.surfaces[s].weight_points );

        for_less( i, 0, deform.surfaces[s].n_self_intersects )
        {
            if( deform.surfaces[s].self_intersects[i].n_weights > 0 )
            {
                FREE( deform.surfaces[s].self_intersects[i].weights );
                FREE( deform.surfaces[s].self_intersects[i].min_distances );
            }
        }

        if( deform.surfaces[s].n_self_intersects > 0 )
            FREE( deform.surfaces[s].self_intersects );

        if( deform.surfaces[s].surface.n_midpoints > 0 )
            FREE( deform.surfaces[s].surface.midpoints );
    }

    for_less( i, 0, deform.n_surf_surfs )
    {
        if( deform.surf_surfs[i].n_weights > 0 )
        {
            FREE( deform.surf_surfs[i].weights );
            FREE( deform.surf_surfs[i].min_distances );
        }
    }

    if( deform.n_surf_surfs > 0 )
    {
        FREE( deform.surf_surfs );
    }

    for_less( i, 0, deform.n_inter_surfaces )
    {
        if( deform.inter_surfaces[i].n_connections > 0 )
            FREE( deform.inter_surfaces[i].connections );
    }

    if( deform.n_inter_surfaces > 0 )
        FREE( deform.inter_surfaces );

    for_less( s, 0, deform.n_surfaces>0?1:0 )
    {
        if( deform.surfaces[s].static_flag )
            continue;

        if( input_graphics_file( input_filenames[s],
                                 &format, &n_objects, &object_list ) != OK ||
                                 n_objects != 1 ||
                                 get_object_type(object_list[0]) != POLYGONS )
            return( 1 );

        polygons = get_polygons_ptr( object_list[0] );
        FREE( polygons->points );
        polygons->points = surface_points[s];

        compute_polygon_normals( polygons );

        (void) output_graphics_file( output_filenames[s],
                                     format, n_objects, object_list );

        delete_object_list( n_objects, object_list );
    }

    if( deform.n_surfaces > 0 )
        FREE( deform.surfaces );

    delete_surface_lookup();

    if( deform.n_surfaces > 0 )
    {
        FREE( input_filenames );
        FREE( output_filenames );
        FREE( surface_points );
    }

    output_alloc_to_file( NULL );

    printf("This code is modified by June Sic Kim at 30/04/2004\n");

    return( did_all_iterations ? 0 : 1 );
}

/*------- surface storage for sharing */

static  int       n_neighbour_sets = 0;
static  struct    {
                      int     n_points;
                      int     *n_neighbours;
                      int     **neighbours;
                  }   *neighbour_sets;

private  void  get_surface_neighbours(
    object_struct  *object,
    int            *n_neighbours_return[],
    int            **neighbours_return[] )
{
    int              i, p, n, *n_neighbours, **neighbours;
    polygons_struct  *surface;

    surface = get_polygons_ptr(object);
    create_polygon_point_neighbours( surface, FALSE,
                                     &n_neighbours, &neighbours, NULL, NULL );

    for_less( i, 0, n_neighbour_sets )
    {
        if( surface->n_points != neighbour_sets[i].n_points )
            continue;

        for_less( p, 0, surface->n_points )
        {
            for_less( n, 0, n_neighbours[p] )
            {
                if( neighbours[p][n] != neighbour_sets[i].neighbours[p][n] )
                    break;
            }

            if( n < n_neighbours[p] )
                break;
        }

        if( p >= surface->n_points )
            break;
    }

    if( i >= n_neighbour_sets )
    {
        SET_ARRAY_SIZE( neighbour_sets, n_neighbour_sets, n_neighbour_sets+1,
                        1 );

        neighbour_sets[n_neighbour_sets].n_points = surface->n_points;
        neighbour_sets[n_neighbour_sets].n_neighbours = n_neighbours;
        neighbour_sets[n_neighbour_sets].neighbours = neighbours;
        ++n_neighbour_sets;
    }
    else
    {
        delete_polygon_point_neighbours( surface, n_neighbours, neighbours,
                                         NULL, NULL );
    }

    *n_neighbours_return = neighbour_sets[i].n_neighbours;
    *neighbours_return = neighbour_sets[i].neighbours;
}

private  void  delete_surface_neighbours( void )
{
    int              i;
    polygons_struct  tmp;

    for_less( i, 0, n_neighbour_sets )
    {
        tmp.n_points = neighbour_sets[i].n_points;
        delete_polygon_point_neighbours( &tmp,
                                         neighbour_sets[i].n_neighbours,
                                         neighbour_sets[i].neighbours,
                                         NULL, NULL );
    }

    if( n_neighbour_sets > 0 )
        FREE( neighbour_sets );
}

/*------- surface storage for sharing */

static  struct    {
                      STRING  filename;
                      int     n_points;
                      Point   *points;
  // Add one line for normals
                      Vector  *normals;
                      int     n_polygons;
                      int     *n_neighbours;
                      int     **neighbours;
                      BOOLEAN model_flag;
                  }   *surface_lookup;

private  Status  input_surface(
    STRING    filename,
    BOOLEAN   model_flag,
    int       *n_points,
    Point     *points[],
    int       *n_polygons,
    int       *n_neighbours[],
    int       **neighbours[] )
{
    int              i, n_objects;
    object_struct    **object_list;
    polygons_struct  *surface;
    File_formats     format;
    STRING           expanded;

    expanded = expand_filename( filename );

    for_less( i, 0, n_surfaces_read )
    {
        if( equal_strings( expanded, surface_lookup[i].filename ) )
            break;
    }

    if( i >= n_surfaces_read )
    {
        if( input_graphics_file( expanded, &format, &n_objects,
                                 &object_list ) != OK || n_objects != 1 ||
            get_object_type(object_list[0]) != POLYGONS )
        {
            print_error( "Error in %s\n", expanded );
            return( ERROR );
        }

        SET_ARRAY_SIZE( surface_lookup, n_surfaces_read, n_surfaces_read+1,
                        DEFAULT_CHUNK_SIZE );

        surface = get_polygons_ptr( object_list[0] );

        if( surface->n_items * 3 != surface->end_indices[surface->n_items-1] )
            print( "\n---- warning, non-triangulated surface will not work with self-intersection testing\n\n" );

        get_surface_neighbours( object_list[0], n_neighbours, neighbours );
        surface_lookup[n_surfaces_read].filename = expanded;
        surface_lookup[n_surfaces_read].n_points = surface->n_points;
        surface_lookup[n_surfaces_read].points = surface->points;
        surface_lookup[n_surfaces_read].n_polygons = surface->n_items;
        surface_lookup[n_surfaces_read].model_flag = model_flag;
        surface_lookup[n_surfaces_read].n_neighbours = *n_neighbours;
        surface_lookup[n_surfaces_read].neighbours = *neighbours;

        ALLOC( surface->points, 1 );
        delete_object_list( n_objects, object_list );
        ++n_surfaces_read;
    }
    else
        delete_string( expanded );

    if( !model_flag )
        surface_lookup[n_surfaces_read].model_flag = FALSE;

    *n_points = surface_lookup[i].n_points;
    *points = surface_lookup[i].points;

    if( n_polygons != NULL )
        *n_polygons = surface_lookup[i].n_polygons;

    *n_neighbours = surface_lookup[i].n_neighbours;
    *neighbours = surface_lookup[i].neighbours;

    return( OK );
}

private  Status  input_surface_with_normals(
    STRING    filename,
    BOOLEAN   model_flag,
    int       *n_points,
    Point     *points[],
    Vector    *normals[],
    int       *n_polygons,
    int       *n_neighbours[],
    int       **neighbours[] )
{
    int              i, n_objects;
    object_struct    **object_list;
    polygons_struct  *surface;
    File_formats     format;
    STRING           expanded;

    expanded = expand_filename( filename );

    for_less( i, 0, n_surfaces_read )
    {
        if( equal_strings( expanded, surface_lookup[i].filename ) )
            break;
    }

    if( i >= n_surfaces_read )
    {
        if( input_graphics_file( expanded, &format, &n_objects,
                                 &object_list ) != OK || n_objects != 1 ||
            get_object_type(object_list[0]) != POLYGONS )
        {
            print_error( "Error in %s\n", expanded );
            return( ERROR );
        }

        SET_ARRAY_SIZE( surface_lookup, n_surfaces_read, n_surfaces_read+1,
                        DEFAULT_CHUNK_SIZE );

        surface = get_polygons_ptr( object_list[0] );

        if( surface->n_items * 3 != surface->end_indices[surface->n_items-1] )
            print( "\n---- warning, non-triangulated surface will not work with self-intersection testing\n\n" );

        get_surface_neighbours( object_list[0], n_neighbours, neighbours );
        surface_lookup[n_surfaces_read].filename = expanded;
        surface_lookup[n_surfaces_read].n_points = surface->n_points;
        surface_lookup[n_surfaces_read].points = surface->points;
        surface_lookup[n_surfaces_read].normals = surface->normals;
        surface_lookup[n_surfaces_read].n_polygons = surface->n_items;
        surface_lookup[n_surfaces_read].model_flag = model_flag;
        surface_lookup[n_surfaces_read].n_neighbours = *n_neighbours;
        surface_lookup[n_surfaces_read].neighbours = *neighbours;

        ALLOC( surface->points, 1 );
        ALLOC( surface->normals, 1 );
        delete_object_list( n_objects, object_list );
        ++n_surfaces_read;
    }
    else
        delete_string( expanded );

    if( !model_flag )
        surface_lookup[n_surfaces_read].model_flag = FALSE;

    *n_points = surface_lookup[i].n_points;
    *points = surface_lookup[i].points;
    *normals = surface_lookup[i].normals;

    if( n_polygons != NULL )
        *n_polygons = surface_lookup[i].n_polygons;

    *n_neighbours = surface_lookup[i].n_neighbours;
    *neighbours = surface_lookup[i].neighbours;

    return( OK );
}

private  void  delete_surface_lookup_model_points( void )
{
    int       i;

    for_less( i, 0, n_surfaces_read )
    {
        if( surface_lookup[i].model_flag )
            FREE( surface_lookup[i].points );
    }
}

private  void  delete_surface_lookup( void )
{
    int  i;

    for_less( i, 0, n_surfaces_read )
        delete_string( surface_lookup[i].filename );

    if( n_surfaces_read > 0 )
        FREE( surface_lookup );
}

/*------- volume storage for sharing */

static  struct  {
                    STRING            filename;
                    Volume            volume;
                    voxel_coef_struct *voxel_lookup;
                } *volume_lookup;

private  Status  lookup_volume(
    STRING             filename,
    Volume             *volume,
    voxel_coef_struct  **voxel_lookup )
{
    int          i;
    Volume       in_volume;
    STRING       expanded;
    Real         high_limit, high, low;
    Data_types   which_type;
    nc_type      convert_to_type;

    expanded = expand_filename( filename );

    for_less( i, 0, n_volumes_read )
    {
        if( equal_strings( expanded, volume_lookup[i].filename ) )
            break;
    }

    if( i >= n_volumes_read )
    {
        if( input_volume_header_only( expanded, 3, XYZ_dimension_names,
                                      &in_volume, NULL ) != OK )
            return( ERROR );

        which_type = get_volume_data_type( in_volume );

        convert_to_type = NC_BYTE;
        high_limit = 0.0;

        if( which_type == UNSIGNED_BYTE )
            convert_to_type = NC_UNSPECIFIED;
        else if( which_type == SIGNED_BYTE ||
                 which_type == UNSIGNED_SHORT ||
                 which_type == SIGNED_SHORT ||
                 which_type == UNSIGNED_INT ||
                 which_type == SIGNED_INT )
        {
            get_volume_voxel_range( in_volume, &low, &high );
            if( ROUND( high - low ) <= 255 )
                high_limit = (Real) ROUND( high - low );
        }

        delete_volume( in_volume );

        if( input_volume( expanded, 3, XYZ_dimension_names,
                          convert_to_type, FALSE, 0.0, high_limit,
                          TRUE, &in_volume, NULL ) != OK )
        {
            return( ERROR );
        }

        SET_ARRAY_SIZE( volume_lookup, n_volumes_read, n_volumes_read+1,
                        DEFAULT_CHUNK_SIZE );

        volume_lookup[n_volumes_read].filename = expanded;
        volume_lookup[n_volumes_read].volume = in_volume;
        ALLOC( volume_lookup[n_volumes_read].voxel_lookup, 1 );
        initialize_lookup_volume_coeficients( volume_lookup[n_volumes_read].
                                              voxel_lookup, in_volume );

        ++n_volumes_read;
    }
    else
        delete_string( expanded );

    *volume = volume_lookup[i].volume;

    if( voxel_lookup != NULL )
        *voxel_lookup = volume_lookup[i].voxel_lookup;

    return( OK );
}

private  void  delete_volume_lookup( void )
{
    int       i;

    for_less( i, 0, n_volumes_read )
    {
        delete_string( volume_lookup[i].filename );
        delete_volume( volume_lookup[i].volume );
        FREE( volume_lookup[i].voxel_lookup );
    }

    if( n_volumes_read > 0 )
        FREE( volume_lookup );
}

private  int  count_edges(
    int    n_points,
    int    n_neighbours[],
    int    *neighbours[] )
{
    int   n_edges, point, n;

    n_edges = 0;
    for_less( point, 0, n_points )
    {
        for_less( n, 0, n_neighbours[point] )
        {
            if( THIS_IS_UNIQUE_EDGE( point, neighbours[point][n] ) )
            {
                ++n_edges;
            }
        }
    }

    return( n_edges );
}

private  float   *create_curvatures(
    int                n_points,
    int                n_neighbours[],
    int                *neighbours[],
    Point              model_points[],
    Real               curvature_scale )
{
    int                         point, n;
    int                         max_neighbours;
    Real                        normal_len;
    float                       *curvatures;
    Point                       centroid;
    Point                       *neigh_points;
    Vector                      normal, offset;

    max_neighbours = 0;
    for_less( point, 0, n_points )
        max_neighbours = MAX( max_neighbours, n_neighbours[point] );

    ALLOC( curvatures, n_points );
    ALLOC( neigh_points, max_neighbours );
    for_less( point, 0, n_points )
    {
        for_less( n, 0, n_neighbours[point] )
            neigh_points[n] = model_points[neighbours[point][n]];

        compute_centroid( n_neighbours[point], neigh_points, &centroid );
        get_neighbours_normal( n_neighbours[point], neigh_points, &normal );
        SUB_POINTS( offset, model_points[point], centroid );

#ifdef USE_CURV_SQRT
        normal_len = MAGNITUDE( normal ) * get_base_length(
                                n_neighbours[point], neigh_points, &centroid );
#else
        normal_len = DOT_VECTORS( normal, normal );
#endif

        curvatures[point] = (float) (curvature_scale *
                                     DOT_VECTORS( offset, normal ) /
                                     normal_len);
    }

    FREE( neigh_points );

    return( curvatures );
}

private  float   *create_model_lengths(
    int                n_points,
    int                n_neighbours[],
    int                *neighbours[],
    Point              model_points[],
    Real               stretch_scale,
    BOOLEAN            use_equal_lengths )
{
    int          ind, point, n, n_edges;
    Real         sum_length, desired_length;
    float        *model_lengths;

    n_edges = count_edges( n_points, n_neighbours, neighbours );;
    ALLOC( model_lengths, n_edges );

    if( use_equal_lengths )
    {
        sum_length = 0.0;
        for_less( point, 0, n_points )
        {
            for_less( n, 0, n_neighbours[point] )
            {
                if( THIS_IS_UNIQUE_EDGE( point, neighbours[point][n] ) )
                {
                    sum_length += distance_between_points(
                                     &model_points[point],
                                     &model_points[neighbours[point][n]] );
                }
            }
        }

        desired_length = sum_length / (Real) n_edges;
    }

    ind = 0;
    for_less( point, 0, n_points )
    {
        for_less( n, 0, n_neighbours[point] )
        {
            if( THIS_IS_UNIQUE_EDGE( point, neighbours[point][n] ) )
            {
#ifdef USE_STRETCH_SQRT
                if( use_equal_lengths )
                {
                    model_lengths[ind] = (float)
                                      (stretch_scale * desired_length);
                }
                else
                {
                    model_lengths[ind] = (float)
                           (stretch_scale * distance_between_points(
                                     &model_points[point],
                                     &model_points[neighbours[point][n]] ) );
                }
#else
                if( use_equal_lengths )
                {
                    model_lengths[ind] = (float)
                           (stretch_scale * stretch_scale *
                            sq_distance_between_points(
                                         &model_points[point],
                                     &model_points[neighbours[point][n]] ) );
                }
                else
                {
                    model_lengths[ind] = (float)
                                     (stretch_scale * stretch_scale *
                                      desired_length * desired_length);
                }
#endif
                ++ind;
            }
        }
    }

    return( model_lengths );
}

private  void   create_model_bends(
    int                n_points,
    int                n_neighbours[],
    int                *neighbours[],
    Point              model_points[],
    float              *model_x[],
    float              *model_y[] )
{
    int          ind, point, n, n_edges, nn, *neighs;
    Real         x, y;
    Real         p0[N_DIMENSIONS], p1[N_DIMENSIONS];
    Real         p2[N_DIMENSIONS], p3[N_DIMENSIONS];

    n_edges = count_edges( n_points, n_neighbours, neighbours );;
    ALLOC( *model_x, n_edges );
    ALLOC( *model_y, n_edges );

    ind = 0;
    for_less( point, 0, n_points )
    {
        p0[X] = RPoint_x( model_points[point] );
        p0[Y] = RPoint_y( model_points[point] );
        p0[Z] = RPoint_z( model_points[point] );

        nn = n_neighbours[point];
        neighs = neighbours[point];
        for_less( n, 0, nn )
        {
            if( THIS_IS_UNIQUE_EDGE( point, neighs[n] ) )
            {
                p1[X] = RPoint_x( model_points[neighs[(n-1+nn)%nn]] );
                p1[Y] = RPoint_y( model_points[neighs[(n-1+nn)%nn]] );
                p1[Z] = RPoint_z( model_points[neighs[(n-1+nn)%nn]] );

                p2[X] = RPoint_x( model_points[neighs[n]] );
                p2[Y] = RPoint_y( model_points[neighs[n]] );
                p2[Z] = RPoint_z( model_points[neighs[n]] );

                p3[X] = RPoint_x( model_points[neighs[(n+1)%nn]] );
                p3[Y] = RPoint_y( model_points[neighs[(n+1)%nn]] );
                p3[Z] = RPoint_z( model_points[neighs[(n+1)%nn]] );

                get_bending_xy( p0, p1, p2, p3, &x, &y );
                (*model_x)[ind] = (float) x;
                (*model_y)[ind] = (float) y;

                ++ind;
            }
        }
    }
}

typedef struct
{
    int                           min_interval;
    int                           max_interval;
    Real                          interval_size;
    Real                          max_movement;
    int                           point_grid_size;
    int                           *start_parameter;

    BOOLEAN                       self_intersect_present;
    self_intersect_lookup_struct  ***si_lookups;

    BOOLEAN                       surf_surf_present;
    surf_surf_lookup_struct       **ss_lookups;
} line_lookup_struct;

private  void  initialize_line_lookup(
    line_lookup_struct  *line_lookup,
    Deform_struct       *deform,
    Real                max_movement,
    Real                movement_per_unit_line,
    int                 point_grid_size,
    int                 start_parameter[] )
{
    int      s;
    BOOLEAN  self_intersect_present;

    self_intersect_present = FALSE;
    for_less( s, 0, deform->n_surfaces )
    {
        if( deform->surfaces[s].n_self_intersects > 0 )
            self_intersect_present = TRUE;
    }

    line_lookup->self_intersect_present = self_intersect_present;
    line_lookup->surf_surf_present = deform->n_surf_surfs > 0;

    if( !self_intersect_present && !line_lookup->surf_surf_present )
        return;

    line_lookup->max_movement = max_movement;
    line_lookup->min_interval = 0;
    line_lookup->max_interval = -1;
    line_lookup->point_grid_size = point_grid_size;
    line_lookup->start_parameter = start_parameter;
    line_lookup->si_lookups = NULL;
    line_lookup->ss_lookups = NULL;

    if( movement_per_unit_line == 0.0 )
        line_lookup->interval_size = 1.0e30;
    else
        line_lookup->interval_size = 2.0 * max_movement / movement_per_unit_line;
}

static  Real  closest_distance1 = 0.0;
static  Real  closest_distance2 = 0.0;

private  void   get_line_lookup(
    line_lookup_struct  *line_lookup,
    Deform_struct       *deform,
    int                 n_parameters,
    Real                parameters[],
    Real                line_pos,
    Real                line_dir[],
    Real                *relative_pos,
    self_intersect_lookup_struct  ***self_lookup,
    surf_surf_lookup_struct       **surf_surf_lookup )
{
    int                           interval, n_insert, size, i, s;
    Real                          centre, *line_point, closest_dist;
    self_intersect_lookup_struct  **si_lookup;
    surf_surf_lookup_struct       *ss_lookup;

    *self_lookup = NULL;
    *surf_surf_lookup = NULL;

    if( !line_lookup->self_intersect_present &&
        !line_lookup->surf_surf_present )
    {
        *relative_pos = 0.0;
        return;
    }

    interval = FLOOR( line_pos / line_lookup->interval_size + 0.5 );

    if( interval < line_lookup->min_interval ||
        interval > line_lookup->max_interval )
    {
        if( line_lookup->min_interval > line_lookup->max_interval )
        {
            line_lookup->min_interval = interval;
            line_lookup->max_interval = interval;
            if( line_lookup->self_intersect_present )
            {
                SET_ARRAY_SIZE( line_lookup->si_lookups, 0, 1,
                                DEFAULT_CHUNK_SIZE );
                line_lookup->si_lookups[0] = NULL;
            }
            if( line_lookup->surf_surf_present )
            {
                SET_ARRAY_SIZE( line_lookup->ss_lookups, 0, 1,
                                DEFAULT_CHUNK_SIZE );
                line_lookup->ss_lookups[0] = NULL;
            }
        }
        else if( interval < line_lookup->min_interval )
        {
            n_insert = line_lookup->min_interval - interval;
            size = line_lookup->max_interval - line_lookup->min_interval + 1;

            if( line_lookup->self_intersect_present )
            {
                SET_ARRAY_SIZE( line_lookup->si_lookups, size, size + n_insert,
                                DEFAULT_CHUNK_SIZE );
                for_down( i, size-1, 0 )
                    line_lookup->si_lookups[i+n_insert] =
                                       line_lookup->si_lookups[i];
                for_less( i, 0, n_insert )
                    line_lookup->si_lookups[i] = NULL;
            }

            if( line_lookup->surf_surf_present )
            {
                SET_ARRAY_SIZE( line_lookup->ss_lookups, size, size + n_insert,
                                DEFAULT_CHUNK_SIZE );
                for_down( i, size-1, 0 )
                    line_lookup->ss_lookups[i+n_insert] =
                                       line_lookup->ss_lookups[i];
                for_less( i, 0, n_insert )
                    line_lookup->ss_lookups[i] = NULL;
            }

            line_lookup->min_interval = interval;
        }
        else
        {
            n_insert = interval - line_lookup->max_interval;
            size = line_lookup->max_interval - line_lookup->min_interval + 1;

            if( line_lookup->self_intersect_present )
            {
                SET_ARRAY_SIZE( line_lookup->si_lookups, size, size + n_insert,
                                DEFAULT_CHUNK_SIZE );
                for_less( i, 0, n_insert )
                    line_lookup->si_lookups[size + i] = NULL;
            }

            if( line_lookup->surf_surf_present )
            {
                SET_ARRAY_SIZE( line_lookup->ss_lookups, size, size + n_insert,
                                DEFAULT_CHUNK_SIZE );
                for_less( i, 0, n_insert )
                    line_lookup->ss_lookups[size + i] = NULL;
            }

            line_lookup->max_interval = interval;
        }
    }

    if( line_lookup->self_intersect_present &&
        line_lookup->si_lookups[interval-line_lookup->min_interval] == NULL )
    {
        ALLOC( si_lookup, deform->n_surfaces );
        for_less( s, 0, deform->n_surfaces )
        {
            if( deform->surfaces[s].n_self_intersects > 0 )
            {
                ALLOC( si_lookup[s],
                       deform->surfaces[s].n_self_intersects );
            }
        }

        line_lookup->si_lookups[interval-line_lookup->min_interval] = si_lookup;

        centre = (Real) interval * line_lookup->interval_size;

        ALLOC( line_point, n_parameters );
        for_less( i, 0, n_parameters )
            line_point[i] = parameters[i] + centre * line_dir[i];

        n_recompute_intersect++;
        closest_dist = recompute_self_intersects( deform,
                                   line_lookup->point_grid_size,
                                   line_lookup->start_parameter,
                                   line_point, line_lookup->max_movement,
                                   line_dir, si_lookup );

        FREE( line_point );

        if( interval == 0 )
            closest_distance1 = closest_dist;
    }

    if( line_lookup->surf_surf_present &&
        line_lookup->ss_lookups[interval-line_lookup->min_interval] == NULL )
    {
        ALLOC( ss_lookup, deform->n_surf_surfs );

        line_lookup->ss_lookups[interval-line_lookup->min_interval] = ss_lookup;

        centre = (Real) interval * line_lookup->interval_size;

        ALLOC( line_point, n_parameters );
        for_less( i, 0, n_parameters )
            line_point[i] = parameters[i] + centre * line_dir[i];

        closest_dist = recompute_surf_surfs( deform,
                                      line_lookup->start_parameter,
                                      line_point, line_lookup->max_movement,
                                      line_dir, ss_lookup );

        if( interval == 0 )
            closest_distance2 = closest_dist;

        FREE( line_point );
    }

    *relative_pos = line_pos - (Real) interval * line_lookup->interval_size;

    if( line_lookup->self_intersect_present )
    {
        *self_lookup = line_lookup->si_lookups
                            [interval-line_lookup->min_interval];
    }

    if( line_lookup->surf_surf_present )
    {
        *surf_surf_lookup = line_lookup->ss_lookups
                            [interval-line_lookup->min_interval];
    }
}

private  void  narrow_line_lookup_range(
    line_lookup_struct  *line_lookup,
    Deform_struct       *deform,
    Real                t_min,
    Real                t_max )
{
    int                            i, s, m, interval_min, interval_max;
    self_intersect_lookup_struct   **si_lookup;
    surf_surf_lookup_struct        *ss_lookup;

    if( !line_lookup->self_intersect_present &&
        !line_lookup->surf_surf_present )
        return;

    interval_min = FLOOR( t_min / line_lookup->interval_size + 0.5 );
    interval_max = FLOOR( t_max / line_lookup->interval_size + 0.5 );

    for_less( i, line_lookup->min_interval, interval_min )
    {
        if( line_lookup->self_intersect_present &&
            line_lookup->si_lookups[i-line_lookup->min_interval] != NULL )
        {
            si_lookup = line_lookup->si_lookups[i-line_lookup->min_interval];

            for_less( s, 0, deform->n_surfaces )
            {
                for_less( m, 0, deform->surfaces[s].n_self_intersects )
                    delete_self_intersect_lookup( &si_lookup[s][m] );
                if( deform->surfaces[s].n_self_intersects > 0 )
                    FREE( si_lookup[s] );
            }

            FREE( si_lookup );
            line_lookup->si_lookups[i-line_lookup->min_interval] = NULL;
        }

        if( line_lookup->surf_surf_present &&
            line_lookup->ss_lookups[i-line_lookup->min_interval] != NULL )
        {
            ss_lookup = line_lookup->ss_lookups[i-line_lookup->min_interval];

            for_less( s, 0, deform->n_surf_surfs )
                delete_surf_surf_lookup( &ss_lookup[s] );

            FREE( ss_lookup );

            line_lookup->ss_lookups[i-line_lookup->min_interval] = NULL;
        }
    }

    for_inclusive( i, interval_max+1, line_lookup->max_interval )
    {
        if( line_lookup->self_intersect_present &&
            line_lookup->si_lookups[i-line_lookup->min_interval] != NULL )
        {
            si_lookup = line_lookup->si_lookups[i-line_lookup->min_interval];

            for_less( s, 0, deform->n_surfaces )
            {
                for_less( m, 0, deform->surfaces[s].n_self_intersects )
                    delete_self_intersect_lookup( &si_lookup[s][m] );
                if( deform->surfaces[s].n_self_intersects > 0 )
                    FREE( si_lookup[s] );
            }

            FREE( si_lookup );
            line_lookup->si_lookups[i-line_lookup->min_interval] = NULL;
        }

        if( line_lookup->surf_surf_present &&
            line_lookup->ss_lookups[i-line_lookup->min_interval] != NULL )
        {
            ss_lookup = line_lookup->ss_lookups[i-line_lookup->min_interval];

            for_less( s, 0, deform->n_surf_surfs )
                delete_surf_surf_lookup( &ss_lookup[s] );

            FREE( ss_lookup );
            line_lookup->ss_lookups[i-line_lookup->min_interval] = NULL;
        }
    }
}

private  void  delete_line_lookup(
    line_lookup_struct  *line_lookup,
    Deform_struct       *deform )
{
    int                            i, s, m;
    self_intersect_lookup_struct   **si_lookup;
    surf_surf_lookup_struct        *ss_lookup;

    if( !line_lookup->self_intersect_present &&
        !line_lookup->surf_surf_present )
        return;

    for_inclusive( i, line_lookup->min_interval, line_lookup->max_interval )
    {
        if( line_lookup->self_intersect_present &&
            line_lookup->si_lookups[i-line_lookup->min_interval] != NULL )
        {
            si_lookup = line_lookup->si_lookups[i-line_lookup->min_interval];

            for_less( s, 0, deform->n_surfaces )
            {
                for_less( m, 0, deform->surfaces[s].n_self_intersects )
                    delete_self_intersect_lookup( &si_lookup[s][m] );
                if( deform->surfaces[s].n_self_intersects > 0 )
                    FREE( si_lookup[s] );
            }

            FREE( si_lookup );
        }

        if( line_lookup->surf_surf_present &&
            line_lookup->ss_lookups[i-line_lookup->min_interval] != NULL )
        {
            ss_lookup = line_lookup->ss_lookups[i-line_lookup->min_interval];

            for_less( s, 0, deform->n_surf_surfs )
                delete_surf_surf_lookup( &ss_lookup[s] );

            FREE( ss_lookup );
        }
    }
        
    if( line_lookup->self_intersect_present &&
        line_lookup->min_interval <= line_lookup->max_interval )
        FREE( line_lookup->si_lookups );
        
    if( line_lookup->surf_surf_present &&
        line_lookup->min_interval <= line_lookup->max_interval )
        FREE( line_lookup->ss_lookups );
}

private  Real   evaluate_along_line(
    Deform_struct        *deform,        
    int                  n_parameters,
    int                  start_parameter[],
    Real                 parameters[],
    Smallest_int         active_flags[],
    Smallest_int         evaluate_flags[],
    Real                 line_dir[],
    Real                 buffer[],
    Real                 dist,
    Real                 boundary_t_coefs[],
    Smallest_int         boundary_flags[],
    Point                boundary_points[],
    Real                 max_value,
    line_lookup_struct   *line_lookup,
    fit_eval_struct      *fit_info  )
{
    int                           p;
    Real                          fit, offset;
    self_intersect_lookup_struct  **si_lookup;
    surf_surf_lookup_struct       *ss_lookup;

    for_less( p, 0, n_parameters ){
      buffer[p] = parameters[p] + dist * line_dir[p];
    }

    get_line_lookup( line_lookup, deform, n_parameters,
                     parameters, dist, line_dir, &offset,
                     &si_lookup, &ss_lookup );

    fit = evaluate_fit( deform, start_parameter, buffer,
                        active_flags, evaluate_flags,
                        boundary_t_coefs, dist,
                        boundary_flags, boundary_points,
                        max_value, FABS(offset), si_lookup,
                        ss_lookup, fit_info );

    return( fit );
}

#ifdef DEBUG
#define DEBUG
private  void  output_list(
    int       n_in_list,
    Real      t_list[],
    Real      f_list[] )
{
    char           str[1000];
    char           filename[1000];
    int            i, j, best, start, end;
    object_struct  **objects;
    static         int  count = 0;
    Real           min_t, max_t;
    Real           min_f, max_f;
    Point          swap;
    lines_struct   *lines;

    (void) sprintf( filename, "list_%d.obj", count );
    ++count;

    start = n_in_list * .2;
    end = n_in_list * .8;

    min_t = t_list[start];
    max_t = t_list[start];
    min_f = f_list[start];
    max_f = f_list[start];
    for_less( i, start, end )
    {
        min_t = MIN( min_t, t_list[i] );
        max_t = MAX( max_t, t_list[i] );
        min_f = MIN( min_f, f_list[i] );
        max_f = MAX( max_f, f_list[i] );
    }

    ALLOC( objects, 1 );

    objects[0] = create_object( LINES );
    lines = get_lines_ptr( objects[0] );
    initialize_lines_with_size( lines, WHITE, end - start, FALSE );

    for_less( i, start, end )
    {
        fill_Point( lines->points[i-start],
                    100.0 * (t_list[i] - min_t) / (max_t - min_t),
                    100.0 * (f_list[i] - min_f) / (max_f - min_f),
                    0.0 );
    }

    for_less( i, 0, lines->n_points-1 )
    {
        best = i;
        for_less( j, i+1, lines->n_points )
        {
            if( Point_x(lines->points[j]) < Point_x(lines->points[best]) )
                best = j;
        }
        swap = lines->points[best];
        lines->points[best] = lines->points[i];
        lines->points[i] = swap;
    }

    (void) output_graphics_file( filename, BINARY_FORMAT,
                                 1, objects );

    delete_object_list( 1, objects );
}
#endif

#define  GOLDEN_RATIO   0.618034

private  Real   minimize_along_line(
    Real                current_fit,
    Deform_struct       *deform,        
    int                 n_parameters,
    int                 start_parameter[],
    Real                parameters[],
    Smallest_int        active_flags[],
    Smallest_int        evaluate_flags[],
    Real                line_dir[],
    Smallest_int        boundary_flags[],
    Point               boundary_points[],
    Real                boundary_coefs[3],
    Real                tolerance,
    Real                function_tolerance,
    line_lookup_struct  *line_lookup,
    Real                *step_taken,
    BOOLEAN             *changed,
    fit_eval_struct     *current_fit_info )
{
    static           int n_evals = 0;
    int              p, n_init;
    Real             t0, t1, t2, f0, f1, f2, t_next, f_next, step, fsize;
    Real             *test_parameters, swap, bottom, new_size, old_size;
    static Real      initial_step_forward = INITIAL_STEP;
    static Real      search_ratio = SEARCH_RATIO;
    fit_eval_struct  t0_fit_info, fit_info, t2_fit_info;
    BOOLEAN          done, prev_failed;
#ifdef DEBUG
    int              n_in_list = 0;
    Real             *t_list;
    Real             *f_list;
#endif

    if( getenv("SEARCH_RATIO") == NULL ||
        sscanf( getenv( "SEARCH_RATIO" ), "%lf", &search_ratio ) != 1 )
    {
        search_ratio = SEARCH_RATIO;
    }

    n_init = n_evals;
    n_recompute_intersect = 0;

    *changed = FALSE;

    ALLOC( test_parameters, n_parameters );

    t0 = 0.0;
    f0 = current_fit;
    t0_fit_info = *current_fit_info;

    t1 = initial_step_forward;
    f1 = evaluate_along_line( deform, n_parameters, start_parameter,
                              parameters,
                              active_flags, evaluate_flags,
                              line_dir,
                              test_parameters, t1,
                              boundary_coefs, boundary_flags, boundary_points,
                              current_fit, line_lookup, current_fit_info );
    ++n_evals;

    step = initial_step_forward;
    if( f1 > f0 )
    {
        swap = t1;
        t1 = t0;
        t0 = swap;

        swap = f1;
        f1 = f0;
        f0 = swap;

        *current_fit_info = t0_fit_info;

        step = -step;
    }

    done = FALSE;
    while( !done )
    {
        t2 = t1 + step;
        step *= search_ratio;

        f2 = evaluate_along_line( deform, n_parameters, start_parameter,
                                  parameters,
                                  active_flags, evaluate_flags,
                                  line_dir,
                                  test_parameters, t2,
                                  boundary_coefs, boundary_flags,
                                  boundary_points,
                                  current_fit, line_lookup, &t2_fit_info );
        ++n_evals;

        if( f2 > f1 )
            done = TRUE;
        else
        {
            t0 = t1;
            f0 = f1;
            t1 = t2;
            f1 = f2;
            *current_fit_info = t2_fit_info;
        }
    }

    if( t2 < t0 )
    {
        swap = t2;
        t2 = t0;
        t0 = swap;

        swap = f2;
        f2 = f0;
        f0 = swap;
    }

#ifdef DEBUG
    ADD_ELEMENT_TO_ARRAY( t_list, n_in_list, t0, 100 );
    ADD_ELEMENT_TO_ARRAY( t_list, n_in_list, t1, 100 );
    ADD_ELEMENT_TO_ARRAY( t_list, n_in_list, t2, 100 );
    n_in_list -= 3;
    ADD_ELEMENT_TO_ARRAY( f_list, n_in_list, f0, 100 );
    ADD_ELEMENT_TO_ARRAY( f_list, n_in_list, f1, 100 );
    ADD_ELEMENT_TO_ARRAY( f_list, n_in_list, f2, 100 );
#endif

    n_init = n_evals - n_init;
    prev_failed = FALSE;

    do
    {
        bottom = 2.0 * (t2*f0-t0*f2-t1*f0+t1*f2-t2*f1+t0*f1);

        if( getenv( "USE_GOLDEN_ALWAYS" ) != NULL )
            bottom = 0.0;

        if( bottom != 0.0 )
        {
            t_next = (-t0*t0*f2+t0*t0*f1+t1*t1*f2-f0*t1*t1+f0*t2*t2-f1*t2*t2)/
                     bottom;
        }

        if( prev_failed || bottom == 0.0 ||
            t_next <= t0 ||
            t_next >= t2 ||
            t_next == t1 )
        {
            if( t1 - t0 > t2 - t1 )
                t_next = t0 + (t1 - t0) * GOLDEN_RATIO;
            else
                t_next = t2 + (t1 - t2) * GOLDEN_RATIO;
        }

        f_next = evaluate_along_line( deform, n_parameters, start_parameter,
                                      parameters,
                                      active_flags, evaluate_flags,
                                      line_dir, test_parameters, t_next,
                                      boundary_coefs, boundary_flags,
                                      boundary_points,
                                      current_fit, line_lookup, &fit_info );
        ++n_evals;

#ifdef DEBUG
        ADD_ELEMENT_TO_ARRAY( t_list, n_in_list, t_next, 100 );
        --n_in_list;
        ADD_ELEMENT_TO_ARRAY( f_list, n_in_list, f_next, 100 );
#endif

        old_size = t2 - t0;

        if( f_next < f1 && t_next < t1 ||
            f_next == f1 && t_next < t1 && t1 - t0 < t2 - t_next )
        {
            t2 = t1;
            f2 = f1;
            t1 = t_next;
            f1 = f_next;
            *current_fit_info = fit_info;
        }
        else if( f_next < f1 && t_next >= t1 ||
                 f_next == f1 && t_next >= t1 && t2 - t1 < t_next - t0 )
        {
            t0 = t1;
            f0 = f1;
            t1 = t_next;
            f1 = f_next;
            *current_fit_info = fit_info;
        }
        else if( f_next > f1 && t_next < t1 ||
                 f_next == f1 && t_next < t1 && t1 - t0 >= t2 - t_next )
        {
            t0 = t_next;
            f0 = f_next;
        }
        else if( f_next > f1 && t_next >= t1 ||
                 f_next == f1 && t_next >= t1 && t2 - t1 >= t_next - t0 )
        {
            t2 = t_next;
            f2 = f_next;
        }

        new_size = t2 - t0;
        if( (new_size / old_size) > .99 )
            prev_failed = TRUE;
        else
            prev_failed = FALSE;

        fsize = MAX( f0, f2 ) - f1;

        narrow_line_lookup_range( line_lookup, deform, t0, t2 );
    }
    while( (tolerance < 0.0 || (t2 - t0) > tolerance) &&
           (function_tolerance < 0.0 || fsize > function_tolerance) );

    if( getenv( "SHORTEN" ) != NULL )
    {
        Real  shorten;
        if( sscanf( getenv( "SHORTEN" ), "%lf", &shorten ) != 1 )
            shorten = 1.0;
        t1 *= shorten;

        f1 = evaluate_along_line( deform, n_parameters, start_parameter,
                                  parameters,
                                  active_flags, evaluate_flags,
                                  line_dir, test_parameters, t1,
                                  boundary_coefs,
                                  boundary_flags, boundary_points,
                                  -1.0, line_lookup, current_fit_info );
    }

    *step_taken = t1;

    initial_step_forward = FABS(t1) / search_ratio;
    if( initial_step_forward < 1.0e-6 )
        initial_step_forward = 1.0e-6;

    if( t1 != 0.0 )
    {
        *changed = TRUE;
        current_fit = f1;
        for_less( p, 0, n_parameters )
            parameters[p] += t1 * line_dir[p];
    }

#ifdef PRINT_DEBUG
{
    int  i;
    for_less( i, 0, n_parameters )
        if( line_dir[i] != 0.0 )
            break;

    print( "Dist[%d] at minimum: %g\n", i, FABS(t1*line_dir[i]) );
}
#endif

    FREE( test_parameters );

    if( getenv( "BRACKET_RATIO" ) != NULL )
        print( "N evals: %d (%d) [%d]    %g   (%.17g %.17g)\n", 
        n_evals, n_init, n_recompute_intersect, t1, t0, t2 );
 
#ifdef DEBUG
    output_list( n_in_list, t_list, f_list );

    FREE( t_list );
    FREE( f_list );
#endif

    return( current_fit );
}

private  int  find_active_points(
    int            n_points,
    int            n_neighbours[],
    int            *neighbours[],
    int            n_start_points,
    int            start_points[],
    int            expansion_level,
    Smallest_int   active_flags[],
    Smallest_int   evaluate_flags[] )
{
    int                 i, current, n_done, n, neigh, level;

    if( n_start_points < 1 || start_points[0] < 0 ||
        start_points[0] >= n_points )
    {
        for_less( i, 0, n_points )
        {
            active_flags[i] = TRUE;
            evaluate_flags[i] = TRUE;
        }
        return( n_points );
    }

    for_less( i, 0, n_points )
        evaluate_flags[i] = FALSE;

    for_less( i, 0, n_start_points )
        evaluate_flags[start_points[i]] = TRUE;

    for_less( level, 0, expansion_level+1 )
    {
        for_less( current, 0, n_points )
            active_flags[current] = evaluate_flags[current];

        for_less( current, 0, n_points )
        {
            if( active_flags[current] )
                continue;

            for_less( n, 0, n_neighbours[current] )
            {
                neigh = neighbours[current][n];
                if( active_flags[neigh] )
                    break;
            }

            if( n < n_neighbours[current] )
                evaluate_flags[current] = TRUE;
        }
    }

    n_done = 0;
    for_less( i, 0, n_points )
    {
        if( active_flags[i] )
            ++n_done;
    }

    return( n_done );
}

private  BOOLEAN   fit_polygons(
    Deform_struct      *deform,
    Real               tolerance,
    Real               function_tolerance,
    int                n_iters,
    int                recompute_every,
    Real               si_step,
    Real               max_movement_step,
    Real               movement_threshold,
    int                n_movements,
    int                point_grid_size,
    int                n_deriv_smoothing_steps,
    Real               deriv_smoothing,
    BOOLEAN            print_deriv,
    BOOLEAN            print_closest,
    BOOLEAN            print_initial,
    STRING             node_dist_filename,
    BOOLEAN            timing_flag,
    Real               successful_movement_threshold )
{
    int                         point, iter, n_parameters, n_edges, w;
    int                         start, n_points1, n_points2;
    int                         *start_parameter, s, i;
    int                         n_since_recompute, n_done;
    int                         active_point, n_active_points;
    int                         p1, oversample, n_eval_points;
    int                         new_fit_num=0;
    Smallest_int                *active_flags, *evaluate_flags;
    STRING                      active_string;
    Real                        fit, min_value, max_value, diff, *new_fits;
    Real                        step_taken;
    Real                        *derivative, *parameters, rms_movement;
    Real                        max_movement, movement;
    Real                        dx, dy, dz, boundary_coefs[3], max_step, delta;
    Real                        *prev_parms;
    conjugate                   conj;
    BOOLEAN                     position_changed, direction_changed;
    BOOLEAN                     bound_present, recomputed;
    Smallest_int                *boundary_flags;
    Point                       *boundary_points;
    BOOLEAN                     boundary_searching;
    int                         n_boundary_points;
    Real                        offset;
    Real                        max_move_per_unit_line;
    Real                        closest_self_intersect, active_threshold;
    int                         *start_points, active_level;
    fit_eval_struct             fit_info, deriv_info;
    self_intersect_lookup_struct  **si_lookup;
    surf_surf_lookup_struct     *ss_lookup;
    line_lookup_struct          line_lookup;
    int                         current_origin;
    int                         one_at_a_time_method, n_steps_before_change;
    int                         last_change;
    BOOLEAN                     one_at_a_time;
    int                         *n_cross_terms, **cross_parms;
    Real                        constant, *linear, *square, **cross_terms;
    Real                        start_time, current_time;
    int                         n_deriv_iters, n, neigh, deriv_iter;
    Real                        fraction, avg_x, avg_y, avg_z;
    FILE                        *parameter_log;

    ALLOC( start_parameter, deform->n_surfaces );
    ALLOC( new_fits, 400 );

    for_less( s, 0, 100 )
    {
        new_fits[s] = 0;
    }

    n_parameters = 0;
    for_less( s, 0, deform->n_surfaces )
    {
        start_parameter[s] = n_parameters;
        n_parameters += 3 * deform->surfaces[s].surface.n_points;
    }

    one_at_a_time = (getenv( "ONE_AT_A_TIME" ) != NULL);
    active_string = getenv( "ACTIVE_POINT" );

    if( one_at_a_time )
    {
        if( sscanf( getenv( "ONE_AT_A_TIME" ), "%d %d %d %lf",
                              &one_at_a_time_method,
                              &active_level,
                              &n_steps_before_change,
                              &active_threshold ) != 4 )
            handle_internal_error( "one at a time" );

        last_change = 0;

        if( one_at_a_time_method == 2 )
        {
            ALLOC( start_points, n_parameters / N_DIMENSIONS );
            n_active_points = 0;
        }

        print( "One at a time: %d %d %d %g\n",
                one_at_a_time_method, active_level, n_steps_before_change,
                active_threshold );

        ALLOC( active_flags, n_parameters / N_DIMENSIONS );
        ALLOC( evaluate_flags, n_parameters / N_DIMENSIONS );

        current_origin = -1;
    }
    else if( active_string != NULL )
    {
        int   ind;

        ALLOC( active_flags, n_parameters / N_DIMENSIONS );
        ALLOC( evaluate_flags, n_parameters / N_DIMENSIONS );

        for_less( s, 0, n_parameters / N_DIMENSIONS )
        {
            active_flags[s] = FALSE;
            evaluate_flags[s] = FALSE;
        }

        s = 0;
        ind = 0;

        while( s < deform->n_surfaces &&
               sscanf( active_string, "%d %d %n",
                       &active_point, &active_level, &ind ) == 2 )
        {
            n_done = find_active_points( deform->surfaces[s].surface.n_points,
                                deform->surfaces[s].surface.n_neighbours,
                                deform->surfaces[s].surface.neighbours,
                                1, &active_point, active_level,
                                &active_flags[start_parameter[s]/3],
                                &evaluate_flags[start_parameter[s]/3] );

            print( "Surface [%d] Active point: %d %d\n", s,
                   active_point, n_done );
            ++s;
            active_string += ind;
        }
    }
    else
    {
        active_flags = NULL;
        evaluate_flags = NULL;
    }

    ALLOC( parameters, n_parameters );

    for_less( s, 0, deform->n_surfaces )
    {
        for_less( point, 0, deform->surfaces[s].surface.n_points )
        {
            parameters[start_parameter[s]+IJ(point,0,3)] =
                   RPoint_x(deform->surfaces[s].surface.points[point]);
            parameters[start_parameter[s]+IJ(point,1,3)] =
                   RPoint_y(deform->surfaces[s].surface.points[point]);
            parameters[start_parameter[s]+IJ(point,2,3)] =
                   RPoint_z(deform->surfaces[s].surface.points[point]);
        }
    }

    closest_self_intersect = -1.0;
    for_less( s, 0, deform->n_surfaces )
    {
        n_edges = count_edges( deform->surfaces[s].surface.n_points,
                               deform->surfaces[s].surface.n_neighbours,
                               deform->surfaces[s].surface.neighbours );

        deform->surfaces[s].surface.n_edges = n_edges;

        for_less( i, 0, deform->surfaces[s].n_bound )
        {
            oversample = deform->surfaces[s].bound[i].oversample;

            deform->surfaces[s].bound[i].image_weight_out /= (Real)
                         (deform->surfaces[s].surface.n_points +
                          oversample * n_edges +
                          oversample * (oversample-1) / 2 *
                          deform->surfaces[s].surface.n_polygons );

            deform->surfaces[s].bound[i].image_weight_in /= (Real)
                         (deform->surfaces[s].surface.n_points +
                          oversample * n_edges +
                          oversample * (oversample-1) / 2 *
                          deform->surfaces[s].surface.n_polygons );
        }

        for_less( i, 0, deform->surfaces[s].n_value )
        {
            get_volume_real_range( deform->surfaces[s].value[i].volume,
                                   &min_value, &max_value );
            diff = (max_value - min_value) / 10.0;
            oversample = deform->surfaces[s].value[i].oversample;
            n_eval_points = deform->surfaces[s].surface.n_points +
                            deform->surfaces[s].surface.n_polygons *
                            ((oversample+3) * (oversample+2)/2 - 3) -
                            oversample * deform->surfaces[s].surface.n_edges;

            deform->surfaces[s].value[i].image_weight /=
                              (Real) n_eval_points * diff * diff;

            deform->surfaces[s].value[i].max_diff_weight /= (Real)
                                   deform->surfaces[s].surface.n_points *
                                   diff * diff;
        }

        for_less( i, 0, deform->surfaces[s].n_stretch )
        {
            n_edges = count_edges( deform->surfaces[s].stretch[i].n_points,
                                   deform->surfaces[s].stretch[i].n_neighbours,
                                   deform->surfaces[s].stretch[i].neighbours );

            deform->surfaces[s].stretch[i].stretch_weight /= (Real) n_edges;
            deform->surfaces[s].stretch[i].max_stretch_weight /= (Real) n_edges;
        }

        for_less( i, 0, deform->surfaces[s].n_curvature )
        {
            deform->surfaces[s].curvature[i].curvature_weight /=
                            (Real) deform->surfaces[s].curvature[i].n_points;
            deform->surfaces[s].curvature[i].max_curvature_weight /=
                            (Real) deform->surfaces[s].curvature[i].n_points;
        }

        for_less( i, 0, deform->surfaces[s].n_bend )
        {
            n_edges = count_edges( deform->surfaces[s].bend[i].n_points,
                                   deform->surfaces[s].bend[i].n_neighbours,
                                   deform->surfaces[s].bend[i].neighbours );

            deform->surfaces[s].bend[i].bend_weight /= (Real) n_edges;
            deform->surfaces[s].bend[i].max_bend_weight /= (Real) n_edges;
        }

        for_less( i, 0, deform->surfaces[s].n_self_intersects )
        {
            for_less( w, 0, deform->surfaces[s].self_intersects[i].n_weights )
            {
                deform->surfaces[s].self_intersects[i].weights[w] /=
                          (Real) deform->surfaces[s].surface.n_points;
                if( closest_self_intersect < 0.0 ||
                    deform->surfaces[s].self_intersects[i].min_distances[w] <
                    closest_self_intersect )
                {
                    closest_self_intersect =
                       deform->surfaces[s].self_intersects[i].min_distances[w];
                }
            }
        }
    }

    for_less( i, 0, deform->n_surf_surfs )
    {
        int  s1, s2;

        s1 = deform->surf_surfs[i].surface_index1;
        s2 = deform->surf_surfs[i].surface_index2;
        for_less( w, 0, deform->surf_surfs[i].n_weights )
        {
            deform->surf_surfs[i].weights[w] /=
                      (Real) (deform->surfaces[s1].surface.n_points +
                              deform->surfaces[s2].surface.n_points);

            if( closest_self_intersect < 0.0 ||
                deform->surf_surfs[i].min_distances[w] < closest_self_intersect)
            {
                closest_self_intersect = deform->surf_surfs[i].min_distances[w];
            }
        }
    }

    for_less( i, 0, deform->n_inter_surfaces )
    {
        n_points1 = deform->surfaces[deform->inter_surfaces[i].surface_index1].
                    surface.n_points;
        n_points2 = deform->surfaces[deform->inter_surfaces[i].surface_index2].
                    surface.n_points;

        oversample = deform->inter_surfaces[i].oversample;

        if( oversample > 0 )
        {
            n_eval_points = n_points1 +
                            deform->surfaces[deform->inter_surfaces[i].
                                 surface_index1].surface.n_polygons *
                            ((oversample+3) * (oversample+2)/2 - 3) -
                            oversample * deform->surfaces[
                      deform->inter_surfaces[i].surface_index1].surface.n_edges;

        }
        else
        {
            n_eval_points = n_points1 + n_points2;
        }

        deform->inter_surfaces[i].weight /= (Real) n_eval_points;
        deform->inter_surfaces[i].max_weight1 /= (Real) n_eval_points;
        deform->inter_surfaces[i].max_weight2 /= (Real) n_eval_points;
    }

    if( closest_self_intersect > 0.0 && max_movement_step < 0.0 )
        max_movement_step = closest_self_intersect;

    boundary_searching = FALSE;
    n_boundary_points = 0;
    for_less( s, 0, deform->n_surfaces )
    {
        for_less( i, 0, deform->surfaces[s].n_bound )
        {
            if( deform->surfaces[s].bound[i].max_dist_weight > 0.0 )
            {
                oversample = deform->surfaces[s].bound[i].oversample;

                n_boundary_points += deform->surfaces[s].surface.n_points +
                            oversample * deform->surfaces[s].surface.n_edges +
                            oversample * (oversample-1) / 2 *
                            deform->surfaces[s].surface.n_polygons;
            }
            boundary_searching = TRUE;
        }
    }

    if( n_boundary_points > 0 )
    {
        print( "Allocating %d boundary points\n", n_boundary_points );
        ALLOC( boundary_points, n_boundary_points );
        ALLOC( boundary_flags, n_boundary_points );
    }

    ALLOC( derivative, n_parameters );
    for_less( p1, 0, n_parameters )
        derivative[p1] = 0.0;

    if( boundary_searching )
    {
        initialize_quadratic_real( n_parameters, &constant,
                                   &linear, &square, &n_cross_terms,
                                   &cross_parms, &cross_terms );
    }
    else
    {
        linear = NULL;
        square = NULL;
        n_cross_terms = NULL;
        cross_parms = NULL;
        cross_terms = NULL;
    }

    conj = initialize_conjugate_gradient( n_parameters );

    bound_present = FALSE;
    for_less( s, 0, deform->n_surfaces )
    {
        if( deform->surfaces[s].n_bound > 0 )
            bound_present = TRUE;
    }

    if( !bound_present )
        recompute_every = 0;

    n_since_recompute = recompute_every;

    if( si_step <= 0.0 )
        si_step = 1.0;

    if( point_grid_size < 1 )
        point_grid_size = 1;

    initialize_line_lookup( &line_lookup, deform,
                            si_step, 0.0, point_grid_size, start_parameter );

    get_line_lookup( &line_lookup, deform, n_parameters,
                     parameters, 0.0,
                     parameters, &offset, &si_lookup, &ss_lookup );

    step_taken = 0.0;

    if( timing_flag )
        start_time = current_cpu_seconds();

    if( movement_threshold > 0.0 && n_movements > 0 )
    {
        ALLOC( prev_parms, n_parameters );
        for_less( i, 0, n_parameters )
            prev_parms[i] = parameters[i];
    }

    
    position_changed = TRUE;
    iter = 0;
    while( iter < n_iters )
    {
        ++iter;

        if( one_at_a_time )
        {
            if( iter == 1 || iter == n_iters )
            {
                current_origin = -1;
                last_change = n_steps_before_change;
            }
            else
            {
                ++last_change;
                if( last_change >= n_steps_before_change )
                {
                    last_change = 0;
                    n_since_recompute = recompute_every;

                    if( one_at_a_time_method == 1 )
                        current_origin = get_random_int( n_parameters / 3 );
                    else if( one_at_a_time_method == 0 )
                    {
                        ++current_origin;
                        if( current_origin >= n_parameters / 3 )
                            current_origin = -1;
                    }
                }
            }

            if( iter != 1 && iter != n_iters && one_at_a_time_method == 2 )
            {
                n_done = find_active_points(
                           deform->surfaces[0].surface.n_points,
                           deform->surfaces[0].surface.n_neighbours,
                           deform->surfaces[0].surface.neighbours,
                           n_active_points, start_points,
                           active_level, active_flags, evaluate_flags );

                print( "N active: %d %d (%d)\n", n_active_points, n_done,
                                                 start_points[0] );
            }
            else
            {
                n_done = find_active_points(
                           deform->surfaces[0].surface.n_points,
                           deform->surfaces[0].surface.n_neighbours,
                           deform->surfaces[0].surface.neighbours,
                           1, &current_origin, active_level, active_flags,
                           evaluate_flags );
                print( "Current: %d %d\n", current_origin, n_done );
            }
        }

        recomputed = FALSE;

        // Modified by June 15/08/2003
        if( recompute_every > 0 && n_since_recompute >= recompute_every  && 
            (deform->surfaces[0].bound->max_outward!=0 || deform->surfaces[0].bound->max_inward!=0 || deform->surfaces[0].bound->max_dist_weight!=0) )
        {
            n_since_recompute = 0;

            find_boundary_points( deform, n_parameters,
                                  start_parameter, parameters,
                                  (one_at_a_time && one_at_a_time_method == 2) ?
                                  NULL : evaluate_flags, boundary_flags,
                                  boundary_points,
                                  &constant, linear, square,
                                  &n_cross_terms, &cross_parms,
                                  &cross_terms );

            recomputed = TRUE;

            if( iter == 1 && boundary_searching )
            {
                realloc_quadratic_cross_terms_real( n_parameters, n_cross_terms,
                                               &cross_parms, &cross_terms );
            }

        }

        get_line_lookup( &line_lookup, deform, n_parameters,
                         parameters,
                         step_taken, derivative, &offset,
                         &si_lookup, &ss_lookup );

        compute_boundary_line_coefficients(
                            n_parameters, parameters,
                            constant, linear, square, n_cross_terms,
                            cross_parms, cross_terms,
                            derivative, boundary_coefs );

        fit = evaluate_fit( deform, start_parameter,
                            parameters,  
                            active_flags, evaluate_flags,
                            boundary_coefs, 0.0,
                            boundary_flags, boundary_points,
                            0.0, FABS(offset), si_lookup, ss_lookup,
                            &fit_info );

        if( iter == 1 || recomputed && print_initial )
        {
            print( "Initial:             %g ", fit );
            print_fit_info( &fit_info );
            print( "\n" );
            (void) flush_file( stdout );
        }

        evaluate_fit_deriv( deform, n_parameters, start_parameter,
                            parameters, boundary_flags, boundary_points,
                            linear, square, n_cross_terms, 
                            cross_parms, cross_terms,
                            si_lookup, ss_lookup, derivative, &deriv_info );

        if( getenv( "DERIV_SMOOTH" ) == NULL ||
            sscanf( getenv( "DERIV_SMOOTH" ), "%d %lf", &n_deriv_iters,
                                                        &fraction ) != 2 )
        {
            n_deriv_iters = n_deriv_smoothing_steps;
            fraction = deriv_smoothing;
        }

        /* Added by June */
        //////////////////////////////////////////////////////////////
/*        if( new_fit_num>=20 )
        {
            new_fit_num=0;
        }
        Real new_fit;
        new_fit = 0;
        for_less( s, 0, 19 )
        {
            new_fit += new_fits[s];
        }
        new_fit /= 20.0;
        if(fabs(fit-new_fit)<0.00001)
        {
            break;
        }
        new_fits[new_fit_num++] = fit;
*/		//////////////////////////////////////////////////////////////

        if( n_deriv_iters > 0 && fraction > 0.0 )
        {
            int     n_points, n3, n_neighs, *neighs, i3;
            Real    *x_deriv_ptr, *y_deriv_ptr, *z_deriv_ptr;
            Real    *x_new_deriv, *y_new_deriv, *z_new_deriv, one_minus_f;
            Real    scale;

            one_minus_f = 1.0 - fraction;
            for_less( deriv_iter, 0, n_deriv_iters )
            for_less( s, 0, deform->n_surfaces )
            {
                x_deriv_ptr = &derivative[start_parameter[s]+0];
                y_deriv_ptr = &derivative[start_parameter[s]+1];
                z_deriv_ptr = &derivative[start_parameter[s]+2];
                n_points = deform->surfaces[s].surface.n_points;
                n3 = 3 * n_points;
                ALLOC( x_new_deriv, n3 );
                y_new_deriv = &x_new_deriv[1];
                z_new_deriv = &x_new_deriv[2];
                for_less( i, 0, n_points )
                {
                    avg_x = 0.0;
                    avg_y = 0.0;
                    avg_z = 0.0;
                    neighs = deform->surfaces[s].surface.neighbours[i];
                    n_neighs = deform->surfaces[s].surface.n_neighbours[i];
                    for_less( n, 0, n_neighs)
                    {
                        neigh = 3*neighs[n];
                        avg_x += x_deriv_ptr[neigh];
                        avg_y += y_deriv_ptr[neigh];
                        avg_z += z_deriv_ptr[neigh];
                    }
                    scale = fraction / (Real) n_neighs;

                    i3 = 3 * i;
                    x_new_deriv[i3] = one_minus_f * x_deriv_ptr[i3] +
                                      scale * avg_x;
                    y_new_deriv[i3] = one_minus_f * y_deriv_ptr[i3] +
                                      scale * avg_y;
                    z_new_deriv[i3] = one_minus_f * z_deriv_ptr[i3] +
                                      scale * avg_z;
                }
                for_less( i, 0, n3 )
                    x_deriv_ptr[i] = x_new_deriv[i];
                FREE( x_new_deriv );
            }
        }

        if( active_flags != NULL )
        {
            for_less( i, 0, n_parameters )
            {
                if( !active_flags[i/3] )
                    derivative[i] = 0.0;
            }
        }

        for_less( s, 0, deform->n_surfaces )
        {
            if( deform->surfaces[s].static_flag )
            {
                for_less( i, 0, deform->surfaces[s].surface.n_points )
                {
                    derivative[start_parameter[s] + 3*i + 0] = 0.0;
                    derivative[start_parameter[s] + 3*i + 1] = 0.0;
                    derivative[start_parameter[s] + 3*i + 2] = 0.0;
                }
            }
        }

        if( one_at_a_time && one_at_a_time_method == 2 &&
            last_change >= n_steps_before_change - 1 )
        {
            Real   len, max_len;

            max_len = -1.0;
            for_less( i, 0, n_parameters / 3 )
            {
                len = derivative[IJ(i,0,3)] * derivative[IJ(i,0,3)] +
                      derivative[IJ(i,1,3)] * derivative[IJ(i,1,3)] +
                      derivative[IJ(i,2,3)] * derivative[IJ(i,2,3)];

                if( len > max_len )
                    max_len = len;
            }

            n_active_points = 0;
            max_len *= active_threshold * active_threshold;
            for_less( i, 0, n_parameters / 3 )
            {
                len = derivative[IJ(i,0,3)] * derivative[IJ(i,0,3)] +
                      derivative[IJ(i,1,3)] * derivative[IJ(i,1,3)] +
                      derivative[IJ(i,2,3)] * derivative[IJ(i,2,3)];

                if( len >= max_len )
                {
                    start_points[n_active_points] = i;
                    ++n_active_points;
                }
            }
        }

        ++n_since_recompute;

/*
        if( one_at_a_time )
        {
            for_less( i, 0, n_parameters )
            {
                derivative[i] *= -1.0;
            }
        }
        else */ if( !get_conjugate_unit_direction( conj, derivative, derivative,
                                           &direction_changed ) )
        {
            if( recompute_every <= 0 || n_since_recompute == 1 )
                break;
            n_since_recompute = recompute_every;
            continue;
        }

        if( active_flags != NULL )
        {
            for_less( i, 0, n_parameters )
            {
                if( !active_flags[i/3] )
                    derivative[i] = 0.0;
            }
        }

        max_move_per_unit_line = get_max_movement_per_unit_line( n_parameters,
                                                                 derivative );

        delete_line_lookup( &line_lookup, deform );
        initialize_line_lookup( &line_lookup, deform,
                                si_step, max_move_per_unit_line,
                                point_grid_size, start_parameter);

        compute_boundary_line_coefficients(
                            n_parameters, parameters,
                            constant, linear, square, n_cross_terms,
                            cross_parms, cross_terms,
                            derivative, boundary_coefs );

        fit = minimize_along_line( fit, deform, n_parameters,
                                   start_parameter,
                                   parameters,
                                   active_flags, evaluate_flags,
                                   derivative,
                                   boundary_flags, boundary_points,
                                   boundary_coefs,
                                   tolerance, function_tolerance,
                                   &line_lookup, &step_taken, &position_changed,
                                   &fit_info );

        if( print_deriv )
        {
            print( "                      Deriv: " );
            print_fit_deriv_info( &deriv_info );
        }

        if( print_closest && closest_self_intersect >= 0.0 &&
            (closest_distance1 >= 0.0 || closest_distance2 >= 0.0) )
        {
            print(    "  Closest: %.4g ", closest_distance1 );
            print(    "  Closest: %.4g\n", closest_distance2 );
        }
        else if( print_deriv )
        {
            print( "\n" );
        }

        if( timing_flag )
        {
            current_time = current_cpu_seconds() - start_time;
            print( "%.3g ", current_time );
        }

        print( "Iter %4d: %.6g \t", iter, fit );
        print_fit_info( &fit_info );
        print( "\n" );
        (void) flush_file( stdout );

/* writing file_info log file by JUNE */
        if( !equal_strings(log_filename,"") )
        {
          parameter_log = open_file_to_write_info(log_filename);
          write_file_info( parameter_log, &fit_info );
          close_file_info(parameter_log);
        }

        if( movement_threshold > 0.0 && n_movements > 0 &&
            (iter % n_movements) == 0 )
        {
            if( iter > 0 )
            {
                max_step = 0.0;
                for_less( i, 0, n_parameters )
                {
                    delta = FABS( parameters[i] - prev_parms[i] );
                    if( delta > max_step )
                        max_step = delta;
                }
            }

            if( max_step < movement_threshold )
            {
                print( "Maximum change in last %d iterations: %g\n",
                        n_movements, max_step );
                break;
            }

            for_less( i, 0, n_parameters )
                prev_parms[i] = parameters[i];
        }

        if( !position_changed && !direction_changed )
        {
            if( recompute_every <= 0 || n_since_recompute == 1 )
                break;
            n_since_recompute = recompute_every;
        }
    }


    if( movement_threshold > 0.0 && n_movements > 0 )
    {
        FREE( prev_parms );
    }

    delete_line_lookup( &line_lookup, deform );

    FREE( derivative );

    if( boundary_searching )
    {
        delete_quadratic_real( n_parameters, linear, square, n_cross_terms,
                               cross_parms, cross_terms );
    }

    delete_conjugate_gradient( conj );

    rms_movement = 0.0;
    max_movement = 0.0;

    for_less( s, 0, deform->n_surfaces )
    {
        start = start_parameter[s];
        for_less( point, 0, deform->surfaces[s].surface.n_points )
        {
            dx = RPoint_x( deform->surfaces[s].surface.points[point] ) -
                           parameters[start+IJ(point,0,3)];
            dy = RPoint_y( deform->surfaces[s].surface.points[point] ) -
                           parameters[start+IJ(point,1,3)];
            dz = RPoint_z( deform->surfaces[s].surface.points[point] ) -
                           parameters[start+IJ(point,2,3)];

            movement = dx * dx + dy * dy + dz * dz;
            max_movement = MAX( max_movement, movement );
            rms_movement += movement;
                          
            fill_Point( deform->surfaces[s].surface.points[point],
                        parameters[start+IJ(point,0,3)],
                        parameters[start+IJ(point,1,3)],
                        parameters[start+IJ(point,2,3)] )
        }
    }

    if( max_movement > 0.0 )
        max_movement = sqrt( max_movement );

    if( n_parameters == 0 )
        rms_movement = 0.0;
    else
        rms_movement = sqrt( rms_movement / (Real) (n_parameters / 3) );

    print( "Movement: %g    %g\n", rms_movement, max_movement );

    FREE( start_parameter );
    FREE( parameters );
    FREE( new_fits );

    if( active_flags != NULL )
    {
        FREE( active_flags );
        FREE( evaluate_flags );
    }

    if( n_boundary_points > 0 )
    {
        FREE( boundary_points );
        FREE( boundary_flags );
    }

    if( one_at_a_time && one_at_a_time_method == 2 )
        FREE( start_points );

    return( iter == n_iters && max_movement >= successful_movement_threshold );
}

/* style of log file (each constraint has 8 bytes)
   B   v   S   C   b   I   SS   2   A   L */
int write_file_info( FILE *log_file, fit_eval_struct *fit_info )
{
    char *fbuf = (char*)malloc(100);

    memset(fbuf,0,100);
    if( log_file != NULL )
    {
        sprintf(fbuf,
            "%11.4f\t%8.4f\t%10.4f\t%8.4f\t%8.4f\t%10.4f\t%8.4f\t%7.3f\t%10.4f\t%10.4f\n",
            fit_info->boundary_fit, fit_info->volume_fit, fit_info->stretch_fit,
            fit_info->curvature_fit, fit_info->bend_fit, 
            fit_info->self_intersect_fit, fit_info->surf_surf_fit,
            fit_info->inter_surface_fit, fit_info->anchor_fit, 
            fit_info->laplacian_fit);

        fwrite(fbuf, 1, 100, log_file);
    }else{
        return 0;
    }
    free(fbuf);

    return 90;
}

FILE *open_file_to_write_info( char *filename )
{
    FILE *handle = NULL;

    if((handle = fopen(filename, "a+")) == NULL)
    {
        print_error("Can't create log file\n");
        return NULL;
    }

    return handle;
}

int close_file_info( FILE *handle)
{
    if( fclose(handle) != 0 )
    {
        print_error("Warning: log file is not closed!\n");
    }
    return 0;
}
