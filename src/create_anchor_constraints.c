#include  <internal_volume_io.h>
#include  <bicpl.h>
#include  <fit.h>

int  main(
    int    argc,
    char   *argv[] )
{
    STRING           input_filename, output_filename;
    int              n_objects, point;
    int              n_points;
    Point            *points;
    Real             desired_distance;
    File_formats     format;
    object_struct    **object_list;
    FILE             *file;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ) ||
        !get_real_argument( 0.0, &desired_distance ) )
    {
        print_error(
             "Usage: %s input.obj output.constraint desired_dist\n", argv[0] );
        return( 1 );
    }

    if( input_graphics_file( input_filename, &format, &n_objects,
                             &object_list ) != OK || n_objects != 1 )
    {
        print( "File %s must contain exactly one object\n",
               input_filename );
        return( 1 );
    }
    
    n_points = get_object_points( object_list[0], &points );

    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
        return( 1 );

    for_less( point, 0, n_points )
    {
        if( output_int( file, point ) != OK ||
            output_real( file, (Real) Point_x(points[point]) ) != OK ||
            output_real( file, (Real) Point_y(points[point]) ) != OK ||
            output_real( file, (Real) Point_z(points[point]) ) != OK ||
            output_real( file, desired_distance ) != OK ||
            output_newline( file ) )
        {
            return( 1 );
        }
    }

    (void) close_file( file );

    return( 0 );
}
