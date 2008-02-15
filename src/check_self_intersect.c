
/* Check if a surface self-intersects. Report the minimal distance
   of self-intersection and, if desired, a texture file showing the
   closest distance to intersection at the vertices.
 */


#include  <fit_3d.h>

private  Status  input_surface( STRING    filename,
                                BOOLEAN   model_flag,
                                int       *n_points,
                                Point     *points[],
                                int       *n_polygons,
                                int       *polygons[],
                                int       *n_neighbours[],
                                int       **neighbours[] );
private  void  delete_surface_lookup( void );
private  void  delete_surface_neighbours( void );

private  void  usage( char   executable_name[] ) {
    STRING  usage_format = "\
Usage: %s surface.obj [distance.txt]\n\n\n";

    print_error( usage_format, executable_name );
}

static  int   n_surfaces_read = 0;

// ---------------------------------------------------------------------
// Initialize a new line lookup structure for the self-surface and/or
// surface-surface intersects.
//
private  void  initialize_line_lookup(
    line_lookup_struct  *line_lookup,
    Deform_struct       *deform,
    Real                max_movement,
    int                 point_grid_size,
    int                 start_parameter[] ) {

    int      s;
    BOOLEAN  self_intersect_present;

    self_intersect_present = FALSE;
    for_less( s, 0, deform->n_surfaces ) {
        if( deform->surfaces[s].n_self_intersects > 0 )
            self_intersect_present = TRUE;
    }

    line_lookup->self_intersect_present = self_intersect_present;
    line_lookup->surf_surf_present = deform->n_surf_surfs > 0;
    line_lookup->closest_dist = 0.0;
    line_lookup->max_movement = max_movement;
    line_lookup->point_grid_size = point_grid_size;
    line_lookup->start_parameter = start_parameter;
    line_lookup->si_lookups = NULL;
    line_lookup->ss_lookups = NULL;
}

// ---------------------------------------------------------------------
// Recompute the line lookup tables (deform->n_surfaces of them) 
// for the self-surface and/or surface-surface intersects. The
// lookup tables (pairs of close triangles) are evaluated at the
// current surface position xyz.
//
private  void   get_line_lookup(
    line_lookup_struct  *line_lookup,
    Deform_struct       *deform,
    int                 n_parameters,
    Real                parameters[],
    Real                line_dir[] ) {

    int                 i, s;
    Real                closest_dist;

    line_lookup->closest_dist = line_lookup->max_movement;

    if( line_lookup->self_intersect_present ) {
        if( line_lookup->si_lookups == NULL ) {
          ALLOC( line_lookup->si_lookups, deform->n_surfaces );
          for_less( s, 0, deform->n_surfaces ) {
            if( deform->surfaces[s].n_self_intersects > 0 ) {
              ALLOC( line_lookup->si_lookups[s],
                     deform->surfaces[s].n_self_intersects );
              for( i = 0; i < deform->surfaces[s].n_self_intersects; i++ ) {
                // We pre-allocate memory for the lookup table, of
                // size 6 * n_polygons. This seems to be a safe upper
                // bound, but it can be made larger. It assumes 6
                // pairs per triangle, on average.
                line_lookup->si_lookups[s][i].n_pairs = 0;
                line_lookup->si_lookups[s][i].n_pairs_alloc = 
                    6*deform->surfaces[s].surface.n_polygons;
                ALLOC( line_lookup->si_lookups[s][i].p1s,
                       line_lookup->si_lookups[s][i].n_pairs_alloc );
                ALLOC( line_lookup->si_lookups[s][i].p2s,
                       line_lookup->si_lookups[s][i].n_pairs_alloc );
                ALLOC( line_lookup->si_lookups[s][i].cases,
                       line_lookup->si_lookups[s][i].n_pairs_alloc );
                ALLOC( line_lookup->si_lookups[s][i].min_line_dists,
                       line_lookup->si_lookups[s][i].n_pairs_alloc );
              }
            }
          }
        }
        closest_dist = recompute_self_intersects( deform,
                                     line_lookup->start_parameter,
                                     parameters, line_lookup->max_movement,
                                     line_dir, line_lookup->si_lookups );
        if( closest_dist < line_lookup->closest_dist ) {
          line_lookup->closest_dist = closest_dist;
        }
    }
}

// ---------------------------------------------------------------------
// Free the data structures for the line lookup table after we
// are all done.
//
private  void  delete_line_lookup( line_lookup_struct  *line_lookup,
                                   Deform_struct       *deform ) {

    int                            s, m;
    self_intersect_lookup_struct   **si_lookup;

    if( !line_lookup->self_intersect_present &&
        !line_lookup->surf_surf_present ) {
        return;
    }

    if( line_lookup->self_intersect_present &&
        line_lookup->si_lookups != NULL ) {
        si_lookup = line_lookup->si_lookups;

        for_less( s, 0, deform->n_surfaces ) {
            for_less( m, 0, deform->surfaces[s].n_self_intersects ) {
                delete_self_intersect_lookup( &si_lookup[s][m] );
            }
            if( deform->surfaces[s].n_self_intersects > 0 ) {
                FREE( si_lookup[s] );
            }
        }

        FREE( si_lookup );
        line_lookup->si_lookups = NULL;
    }
}


int main( int argc, char * argv[] ) {

    STRING               input_filename = NULL;
    STRING               output_filename = NULL;
    int                  p1, point;
    Deform_struct        deform;
    one_surface_struct   surf;
    surface_struct       surface;

    n_surfaces_read = 0;

    deform.n_surfaces = 0;
    deform.n_inter_surfaces = 0;
    deform.n_surf_surfs = 0;
    deform.n_intersect_wm = 0;

    initialize_argument_processing( argc, argv );

    if( argc >= 2 ) {
      if( !get_string_argument( NULL, &input_filename ) ) {
        print_error( "Error in -surface arguments.\n" );
        usage( argv[0] );
        return( 1 );
      }
      if( argc >= 3 ) {
        if( !get_string_argument( NULL, &output_filename ) ) {
          print_error( "Error in -surface arguments.\n" );
          usage( argv[0] );
          return( 1 );
        }
      }
    } else {
      usage( argv[0] );
      return( 1 );
    }

    if( input_surface( input_filename, FALSE,
                       &surface.n_points, &surface.points,
                       &surface.n_polygons, &surface.triangles,
                       &surface.n_neighbours, &surface.neighbours ) != OK ) {
      return( 1 );
    }
    // if( input_filename ) FREE( input_filename );

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
    surf.n_self_intersects = 1;

    ADD_ELEMENT_TO_ARRAY( deform.surfaces, deform.n_surfaces, surf, 1);

    ALLOC( deform.surfaces[0].self_intersects, 1);
    deform.surfaces[0].self_intersects[0].n_weights = 0;

    int * start_parameter, s;
    ALLOC( start_parameter, deform.n_surfaces );

    int n_parameters = 0;
    for_less( s, 0, deform.n_surfaces ) {
      start_parameter[s] = n_parameters;
      n_parameters += N_DIMENSIONS * deform.surfaces[s].surface.n_points;
    }

    line_lookup_struct   line_lookup;
    Real  si_step = 0.5;
    int   point_grid_size = 30;
    initialize_line_lookup( &line_lookup, &deform, si_step,
                            point_grid_size, start_parameter );

    Real   * derivative, * parameters;
    ALLOC( parameters, n_parameters );

    for_less( point, 0, deform.surfaces[0].surface.n_points ) {
        parameters[start_parameter[0]+IJ(point,0,3)] =
               RPoint_x(deform.surfaces[0].surface.points[point]);
        parameters[start_parameter[0]+IJ(point,1,3)] =
               RPoint_y(deform.surfaces[0].surface.points[point]);
        parameters[start_parameter[0]+IJ(point,2,3)] =
               RPoint_z(deform.surfaces[0].surface.points[point]);
    }

    ALLOC( derivative, n_parameters );
    for_less( p1, 0, n_parameters ) {
        derivative[p1] = 1.0;
    }

    get_line_lookup( &line_lookup, &deform, n_parameters,
                     parameters, derivative );
    FREE( derivative );
    FREE( parameters );
    FREE( start_parameter );

    self_intersect_lookup_struct * si_lookup = &line_lookup.si_lookups[0][0];
    int n_candidates = get_n_self_intersect_candidate( si_lookup );

    int bad_count = 0; 
    for( s = 0; s < n_candidates; s++ ) {
      Real dist = sqrt( si_lookup->min_line_dists[s] );
      if( dist <= 1.0e-08 ) bad_count++;
    }   

    printf( "Self-intersection distance = %f\n", line_lookup.closest_dist );
    printf( "Number of self_intersecting triangles = %d\n", bad_count );

    if( output_filename ) {
      Real * nodal_distance = NULL;
      ALLOC( nodal_distance, deform.surfaces[0].surface.n_points );
      for( s = 0; s < deform.surfaces[0].surface.n_points; s++ ) {
        nodal_distance[s] = si_step;
      }

      int * triangles = deform.surfaces[0].surface.triangles;
      for( s = 0; s < n_candidates; s++ ) {
        int poly1 = si_lookup->p1s[s];
        int poly2 = si_lookup->p2s[s];
        Real dist = sqrt( si_lookup->min_line_dists[s] );
        nodal_distance[triangles[3*poly1]] = MIN( nodal_distance[triangles[3*poly1]], dist );
        nodal_distance[triangles[3*poly1+1]] = MIN( nodal_distance[triangles[3*poly1+1]], dist );
        nodal_distance[triangles[3*poly1+2]] = MIN( nodal_distance[triangles[3*poly1+2]], dist );
        nodal_distance[triangles[3*poly2]] = MIN( nodal_distance[triangles[3*poly2]], dist );
        nodal_distance[triangles[3*poly2+1]] = MIN( nodal_distance[triangles[3*poly2+1]], dist );
        nodal_distance[triangles[3*poly2+2]] = MIN( nodal_distance[triangles[3*poly2+2]], dist );
      }   

      FILE * fp = fopen( output_filename, "wt" );
      for( s = 0; s < deform.surfaces[0].surface.n_points; s++ ) {
        fprintf( fp, "%f\n", nodal_distance[s] );
      }
      fclose( fp );
      FREE( nodal_distance );
    }

    delete_surface_neighbours();
    delete_line_lookup( &line_lookup, &deform );
    FREE( deform.surfaces[0].self_intersects );
    FREE( deform.surfaces );

    delete_surface_lookup();

    return( 0 );
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
                      Vector  *normals;
                      int     n_polygons;
                      int     *polygons;
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
    int       *polygons[],
    int       *n_neighbours[],
    int       **neighbours[] )
{
    int              i, n_objects;
    object_struct    **object_list;
    polygons_struct  *surface;
    File_formats     format;
    STRING           expanded;

    expanded = expand_filename( filename );

    for_less( i, 0, n_surfaces_read ) {
        if( equal_strings( expanded, surface_lookup[i].filename ) )
            break;
    }

    if( i >= n_surfaces_read ) {
        if( input_graphics_file( expanded, &format, &n_objects,
                                 &object_list ) != OK || n_objects != 1 ||
            get_object_type(object_list[0]) != POLYGONS ) {
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
        surface_lookup[n_surfaces_read].polygons = surface->indices;
        surface_lookup[n_surfaces_read].model_flag = model_flag;
        surface_lookup[n_surfaces_read].n_neighbours = *n_neighbours;
        surface_lookup[n_surfaces_read].neighbours = *neighbours;

        // dummy points, indices so as not to free the real ones in delete_object_list.
        ALLOC( surface->points, 1 );
        ALLOC( surface->indices, 1 );
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

    if( polygons != NULL )
        *polygons = surface_lookup[i].polygons;

    *n_neighbours = surface_lookup[i].n_neighbours;
    *neighbours = surface_lookup[i].neighbours;

    return( OK );
}

private  void  delete_surface_lookup( void ) {
    int  i;

    for_less( i, 0, n_surfaces_read ) {
        delete_string( surface_lookup[i].filename );
        FREE( surface_lookup[i].points );
        FREE( surface_lookup[i].polygons );
    }

    if( n_surfaces_read > 0 )
        FREE( surface_lookup );
}
