
/*
   SURFACE_ANGLES

   Find the angle between the W/G surfaces and compare to the
   surface normal. How much distortion is there during the 
   gray surface expansion? The output should be blurred after
   using depth_potential.

   surface_angles white.obj mid.obj gray.obj angles.txt

   Values: white.obj = white surface
           mid.obj = mid surface
           gray.obj = gray surface
           angles.txt = output angles

   By: Claude Lepage, February 2015.

   COPYRIGHT: McConnell Brain Imaging Center, 
              Montreal Neurological Institute,
              Department of Psychology,
              McGill University, Montreal, Quebec, Canada. 
  
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
*/

#include <stdio.h>
#include <volume_io.h>
#include <bicpl.h>

// Prototypes of functions in this file.

static void usage( char * );

// Main program.

int main( int argc, char * argv[] ) {

  int      i, j, k;

  int n_objects;
  object_struct  ** white_object;
  object_struct  ** mid_object;
  object_struct  ** gray_object;

  polygons_struct * white_surface;
  polygons_struct * mid_surface;
  polygons_struct * gray_surface;
  File_formats      format;
  STRING            expanded;

  // Parse the command line arguments for the file names.

  if( argc != 5 ) {
    usage( argv[0] );
    return( 1 );
  }

  char * angles_txt = argv[4];
  if( angles_txt == NULL ) {
    usage( argv[0] );
    return( 1 );
  }

  // Read the white surface.
  expanded = expand_filename( argv[1] );
  int err = input_graphics_file( expanded, &format, &n_objects,
                                 &white_object );
  delete_string( expanded );

  if( err != OK || n_objects != 1) {
    print_error( "Error reading file %s\n", argv[1] );
    return( 1 );
  }

  // Read the mid surface.
  expanded = expand_filename( argv[2] );
  err = input_graphics_file( expanded, &format, &n_objects, &mid_object );
  delete_string( expanded );

  if( err != OK || n_objects != 1) {
    print_error( "Error reading file %s\n", argv[2] );
    return( 1 );
  }

  // Read the gray surface.
  expanded = expand_filename( argv[3] );
  err = input_graphics_file( expanded, &format, &n_objects, &gray_object );
  delete_string( expanded );

  if( err != OK || n_objects != 1) {
    print_error( "Error reading file %s\n", argv[3] );
    return( 1 );
  }

  if( get_object_type(white_object[0]) != POLYGONS ||
      get_object_type(mid_object[0]) != POLYGONS ||
      get_object_type(gray_object[0]) != POLYGONS ) {
    print_error( "Error in contents of files %s %s %s\n", 
                 argv[1], argv[2], argv[3] );
    return( 1 );
  }

  white_surface = get_polygons_ptr( white_object[0] );
  mid_surface = get_polygons_ptr( mid_object[0] );
  gray_surface = get_polygons_ptr( gray_object[0] );

  if( white_surface->n_points != gray_surface->n_points ||
      white_surface->n_points != mid_surface->n_points ||
      white_surface->n_items != gray_surface->n_items ||
      white_surface->n_items != mid_surface->n_items ) {
    print_error( "Mismatched number of points/polygons in white, mid and gray surfaces.\n" );
    return( 1 );
  }

  // For each vertex, compute the ray between W/G surfaces and
  // compare its orientation angle relative to the normal vector
  // of the mid surface.

  FILE   * fp = fopen( angles_txt, "w+t" );

  for( j = 0; j < white_surface->n_points; j++ ) {

    Real vecx = gray_surface->points[j].coords[0] -
                white_surface->points[j].coords[0];
    Real vecy = gray_surface->points[j].coords[1] -
                white_surface->points[j].coords[1];
    Real vecz = gray_surface->points[j].coords[2] -
                white_surface->points[j].coords[2];

    Real norm = sqrt( vecx * vecx + vecy * vecy + vecz * vecz );
    Real cos_angle = ( vecx * mid_surface->normals[j].coords[0] +
                       vecy * mid_surface->normals[j].coords[1] +
                       vecz * mid_surface->normals[j].coords[2] ) / norm;
    // Real angle = acos( cos_angle ) * 180.0 / 3.14159265;
#if 0
if( cos_angle < 0.0 ) {
  printf( "Vertex j = %d\n", j );
  printf( "  w at %g %g %g\n", white_surface->points[j].coords[0],
          white_surface->points[j].coords[1], 
          white_surface->points[j].coords[2] );
  printf( "  m at %g %g %g\n", mid_surface->points[j].coords[0],
          mid_surface->points[j].coords[1], 
          mid_surface->points[j].coords[2] );
  printf( "  g at %g %g %g\n", gray_surface->points[j].coords[0],
          gray_surface->points[j].coords[1], 
          gray_surface->points[j].coords[2] );
  printf( "  m at %g %g %g to %g %g %g\n", 
          mid_surface->normals[j].coords[0],
          mid_surface->normals[j].coords[1], 
          mid_surface->normals[j].coords[2],
      mid_surface->points[j].coords[0]+0.5*mid_surface->normals[j].coords[0],
      mid_surface->points[j].coords[1]+0.5*mid_surface->normals[j].coords[1],
      mid_surface->points[j].coords[2]+0.5*mid_surface->normals[j].coords[2] );



  printf( "  norm = %g  cos_angle = %g\n", norm, cos_angle );
}
#endif
    fprintf( fp, "%f\n", 1.0/cos_angle );
  }
  fclose( fp );

  // Free the memory after usage.

  delete_object_list( 1, white_object );
  delete_object_list( 1, mid_object );
  delete_object_list( 1, gray_object );

  return 0;
}


// -------------------------------------------------------------------
// Help message on how to use this module.
//
static void usage( char * executable_name ) {

  STRING  usage_format = "\
Usage: %s white.obj mid.obj gray.obj angles.txt\n\
Values: white.obj = white input object file\n\
        mid.obj = mid input object file\n\
        gray.obj = gray input object file\n\
        angles.txt = output .txt file for vertex-wise distoration angles\n\n";

  print_error( usage_format, executable_name );
}

