#!xPERLx -w

# Part of the ASP cortex extraction suite. For a description of the
# algorithms used here see MacDonald et. al. "Automated 3-D Extraction
# of Inner and Outer Surfaces of Cerebral Cortex from MRI" in
# NeuroImage, 12, 340-356, 2000

# extract_white_surface takes a masked classified volume as an input,
# and writes out the various iterations of the deformation. The exact
# iterations are determined by the schedule included in this file.

# Author: David MacDonald
# Last Modified: Sep 2001, Jason Lerch <jason@bic.mni.mcgill.ca>
#    * removed any hardcoded paths
#    * uses MNI::DataDir
#    * soften the stretch weight constraint on last iteration for slightly 
#      better fitting.

require "xINCDIRx/deform_utils.pl";
use MNI::DataDir;

    $volume = shift;
    $output_file2 = shift;
    $isovalue2 = shift;
    $initial_model2 = shift;

    $start_n = shift;
    $end_n = shift;
    $dont_copy = shift;

    if( ! defined($isovalue2) )
    {
        die "Usage: $0  input.mnc output_prefix  isovalue [model]  [start] [end] [dont copy]\n";
    }

#--------- initialize volumes

#---------

    $output2 = &get_prefix( $output_file2 );

    if( !defined( $initial_model2 ) )
    {
      $model_data_dir = MNI::DataDir::dir('ASP');
      MNI::DataDir::check_data($model_data_dir, [qw(white_model_320.obj)]);
      $initial_model2 = "${model_data_dir}/white_model_320.obj";
    }

    $fit = "new_fit_3d ";
    $refine = "refine_mesh";

    $self_dist2 = 0.01;
    $n_selfs = 9;
    $self_factor = 1.0;

    $stop_threshold = 3e-2;
    $stop_iters = 10;

    $filter = 0;
    $n_per = 5;
    $tolerance = 1.0e-2;
    $tolerance = 1.0e-10;
    $f_tolerance = 1.0e-2;
    $f_tolerance = 1.0e-10;

    $iters_scale = 1.0;
    $break_scale = 1.0;
    $oo_scale = 0.5;
    $iters_override = 0;

    $stretch_scale = 1;
    $curvature_scale = 0;

    $wo = 100;

    @schedule = (

#size    sw    cw  n_itr  inc  offo offi  si over   sw   self
#-----  ----  ---- -----  ---- ---- ----  -- ----  ----  ----
   320, 1e3,   0,   1000,  50, $wo,  1, .5,  40,  1e0,  1,
  1280, 1e3,   0,   1000,  50, $wo,  1, .5,  20,  1e0,  1,
  1280, 1e2,   0,   2000,  50, $wo,  1,  1,  20,  1e0,  5,
  5120, 1e2,   0,   2000,  50, $wo,  1,  1,  10,  1e0,  3,
  5120, 1e1,   0,   2000,  50, $wo,  1,  1,  10,  1e0,  3,
 20480, 1e1,   0,   5000,  50, $wo,  1,  1,   5,  1e0,  1.5,
 81920, 2,     0,   1000,  50, $wo,  1,  1,   2,  1e0,  .75,
  );


    $sched_size =  11;

    $st = $schedule[0];

    if( ! defined($start_n) )
        { $start_n = 0; }

    if( ! defined($end_n) )
        { $end_n = @schedule / $sched_size - 1; }

#------ loop over each schedule

    if( !$dont_copy && $start_n <= 0 ) {

        if( $output_file2 ne "-" )
        {


            system_call( "set_object_colour $initial_model2 ${output2}_${st}.obj white" );
            system_call( "set_object_colour $initial_model2 ${output2}_m_${st}.obj white" );

        }
    }

    $prev_n = $st;

    for( $i = 0;  $i < @schedule;  $i += $sched_size )
    {
        $step = $i / $sched_size;

        #--- get the 14 components of the deformation schedule entry

        ( $size, $sw, $cw, $n_iters, $iter_inc, $offo, $offi,
          $si_step, $oversample, $self_weight, $self_dist ) =
                     @schedule[$i..$i+$sched_size-1];

        $sw *= $stretch_scale;
        $cw *= $curvature_scale;
        $oversample *= $oo_scale;

        if( $iters_override > 0 )
            { $n_iters = $iters_override; }
        else
            { $n_iters = int( $n_iters * $iters_scale ); }

        $self_weight *= $self_factor;

        $self2 = get_self_intersect( $self_weight, $n_selfs, $self_dist,
                                     $self_dist2 );

        if( $step < $start_n || $step > $end_n )
            { $prev_n = $size;  next; }

        if( $size != $prev_n && (!$dont_copy || $step > $start_n) )
        {
            if( $size < 320 )
            {
                system_call( "$refine ${output2}_${prev_n}.obj " .
                             " ${output2}_m_${prev_n}.obj" .
                             " ${output2}_m_${size}.obj $size " );
                system_call( "$refine ${output2}_${prev_n}.obj " .
                             " ${output2}_${prev_n}.obj" .
                             " ${output2}_${size}.obj $size " );
            }
            else
            {
                system_call( "subdivide_polygons " .
                             " ${output2}_m_${prev_n}.obj" .
                             " ${output2}_m_${size}.obj $size " );
                system_call( "subdivide_polygons " .
                             " ${output2}_${prev_n}.obj" .
                             " ${output2}_${size}.obj $size " );
            }
        }

        $prev_n = $size;

        #--- if the schedule size is greater than the current number of
        #--- polygons in the deforming surface, subdivide the deforming surface

        print( "Fitting $size-sized polygons at filter ${filter}, " .
               "max $n_iters iters.\n" );

        $b2 = " -boundary $offo $offi $volume " .
                           " $isovalue2 - 3 50 0 0 $oversample ";

        $iter_inc *= $break_scale;
        if( $iter_inc <= 0 )  { $iter_inc = $n_iters; }

        $n_failures = 0;

        for( $iter = 0;  $iter < $n_iters;  $iter += $iter_inc )
        {
            system( "echo Step ${size}: $iter / $n_iters    $sw $cw" );

            if( $size < 320 )
            {
                system_call( "$refine  ${output2}_${size}.obj " .
                             " ${output2}_m_${size}.obj" .
                             " ${output2}_m_${size}.obj $size " );
                system_call( "$refine  ${output2}_${size}.obj " .
                             " ${output2}_${size}.obj" .
                             " ${output2}_${size}.obj $size " );
            }


            $ni = $n_iters - $iter;
            if( $ni > $iter_inc )  { $ni = $iter_inc; }

            $sw2 = $sw;

            $surf2_info = " -surface ${output2}_$size ${output2}_$size " .
              " -stretch $sw2 ${output2}_m_$size -.9 0 0 0".
              " $b2 ".
              " $self2 ";

            $command = "$fit $surf2_info ".
                       " -print_deriv " .
                       " -step $si_step " .
                       " -fitting $ni $n_per $tolerance " .
                       " -ftol $f_tolerance " .
                       " -stop $stop_threshold $stop_iters ";

            $ret = system_call( "$command", 1 );

            system_call( "measure_surface_area ${output2}_$size" );

            if( $ret == 1 )
            {
                ++$n_failures;

                if( $n_failures == 2 )
                    { last; }
            }
            else
            {
                $n_failures = 0;
            }
        }
    }

    print( "Surface extraction finished.\n" );

    &clean_up();
