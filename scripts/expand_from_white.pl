#!xPERLx -w

require "xINCDIRx/deform_utils.pl";
use MNI::DataDir;

#    &print_out_script_info();

    $volume = shift;
    $white_surface = shift;
    $output = shift;
    $isovalue = shift;

    $laplacian_file = shift;
    $start_n = shift;
    $end_n = shift;
    $dont_copy = shift;

    if( ! defined($isovalue) )
    {
        die "Usage: $0  input.mnc white.obj output1_pref [start] [end] [dont copy]\n";
    }

#---------

    $fit = "new_fit_3d ";

    $oversample_reference = 20480;
    $n_polygons = `print_n_polygons $white_surface`;
    chop( $n_polygons );

    $model_data_dir = MNI::DataDir::dir('ASP');
    MNI::DataDir::check_data($model_data_dir, ["ellipsoid_${n_polygons}.obj.gz"]);

    $model = $white_surface;
    $model = "${model_data_dir}/ellipsoid_${n_polygons}.obj.gz";

    $self_dist2 = 0.01;
    $n_selfs = 9;
    $self_factor = 1.0;

    $stop_threshold = 3e-2;
    $stop_iters = 10;

    $value_differential_offset = 0.5;

    $filter = 0;
    $n_per = 1;
    $tolerance = 1.0e-2;
    $tolerance = 1.0e-10;
    $f_tolerance = 1.0e-2;
    $f_tolerance = 1.0e-10;

    $iters_scale = 1.0;
    $break_scale = 1.0;
    $oo_scale = 1.0 * sqrt( $oversample_reference / $n_polygons );
    $value_oo_scale = 1.0;
    $iters_override = 0;

    $stretch_scale = 1;
    $curvature_scale = 0;

    @schedule = (

#size   sw  n_it  inc offo offi  si over   sw   self
#----- ---- ----  --- ---  ----  -- ----  ----  ----
 1e1,   0, 1000, 100,  1,   1,   .5,  1,  1e1,  .25, 4, -1, -5, 2, 0,1e-5,0,0,1.0, 2.0, 1e-3, 1e3,
 1e1,   0, 1500,  50,  1,   1,   .5,  1,  1e1,  .25, 4, -1, -5, 2, 0,1e-5,0,0,1.0, 2.0, 1e-3, 1e4,
 1e1,   0,  500,  50,  1,   1,   .5,  1,  1e1,  .25, 4, -1, -5, 2, 0,1e-5,0,0,1.0, 2.0, 1e-3, 1e5,
  );


    $sched_size =  22;
 
    if( ! defined($start_n) )
        { $start_n = 0; }

    if( ! defined($end_n) )
        { $end_n = @schedule / $sched_size - 1; }

#------ loop over each schedule

    if( !$dont_copy && $start_n <= 0 ) {
        system_call( "set_object_colour $white_surface ${output} white" );
    }

    $once = 0;

    for( $i = 0;  $i < @schedule;  $i += $sched_size )
    {
        $step = $i / $sched_size;

        #--- get the components of the deformation schedule entry

        ( $sw, $cw, $n_iters, $iter_inc, $offo, $offi,
          $si_step, $oversample, $self_weight, $self_dist, $desired_dist,
          $anchor_weight, $min_dist, $max_dist,
          $v_weight, $v_max_weight,
          $adaptive_anchor_ratio, $adaptive_boundary_ratio,
          $laplacian_sampling, $laplacian_factor, $laplacian_weight,
          $surf_surf_weight ) =
                     @schedule[$i..$i+$sched_size-1];

        $sw *= $stretch_scale;
        $cw *= $curvature_scale;
        $oversample *= $oo_scale;
        $value_oversample = $oversample * $value_oo_scale;

        if( $iters_override > 0 )
            { $n_iters = $iters_override; }
        else
            { $n_iters = int( $n_iters * $iters_scale ); }

        $self_weight *= $self_factor;

        $self2 = get_self_intersect( $self_weight, $n_selfs, $self_dist,
                                     $self_dist2 );

        if( $step < $start_n || $step > $end_n )
            { next; }

        #--- if the schedule size is greater than the current number of
        #--- polygons in the deforming surface, subdivide the deforming surface

        print( "Fitting polygons at filter ${filter}, " .
               "max $n_iters iters.\n" );

        $b2 = " -boundary $offo $offi $volume " .
                           " $isovalue - 20 20 0 0 $oversample " .
              " -value_differential $value_differential_offset 1e-9" .
              " -value 0 $masked_volume 0 2.0 -10 0 1e4 $value_oversample";

        $b2 = " -boundary $offo $offi $volume " .
                           " $isovalue - 20 2 0 0 $oversample ";

        $iter_inc *= $break_scale;
        if( $iter_inc <= 0 )  { $iter_inc = $n_iters; }

        $surf2_info = " -mode three".
          " -surface ${output} ${output} $white_surface" .
          " -stretch $sw $input -1.0 0 0 0".
          " $b2 ".
          " $self2 ";

        if( $anchor_weight > 0 )
        {
            $tmp_anchor = "/tmp/anchor_${$}.txt";

            register_tmp_files( $tmp_anchor );

            system_call( "create_anchor_constraints " .
                         " $white_surface $tmp_anchor $desired_dist"  );
            $anchor = "-anchor $anchor_weight $tmp_anchor ".
                      "1e4 $min_dist $max_dist";
        }
        else
        {
            $anchor = "";
        }

        $n_failures = 0;

        for( $iter = 0;  $iter < $n_iters;  $iter += $iter_inc )
        {
            system( "echo Step: $iter / $n_iters    $sw $cw" );

            if( $once > 0 ){
              $surf2_info = " -mode three".
                " -surface ${output} ${output} ${white_surface}" .
                " -equal_lengths".
                " -stretch $sw ${output} -1.0 0 0 0".
                " $b2 ".
                " $self2 ";
	
              $ni = $n_iters - $iter;
              if( $ni > $iter_inc )  { $ni = $iter_inc; }
              $command = "$fit $surf2_info ".
                " $anchor " .
                " -volume $v_weight $v_max_weight $adaptive_anchor_ratio $adaptive_boundary_ratio " .
                " -laplacian $laplacian_file $laplacian_weight 0 20 $laplacian_factor $laplacian_sampling ".
                " -surf_surf 0 1 $surf_surf_weight 0.1" .
                " -print_deriv " .
                " -step $si_step " .
                " -fitting $ni $n_per $tolerance " .
                " -ftol $f_tolerance " .
                " -stop $stop_threshold $stop_iters ";
            }
            else
            {
              $ni = $n_iters - $iter;
              if( $ni > $iter_inc )  { $ni = $iter_inc; }

              $command = "$fit $surf2_info ".
                       " $anchor " .
                       " -volume $v_weight $v_max_weight $adaptive_anchor_ratio $adaptive_boundary_ratio " .
                       " -laplacian $laplacian_file $laplacian_weight 0 20 $laplacian_factor $laplacian_sampling ".
                       #" -surf_surf 0 1 $surf_surf_weight 0.1" .
                       " -print_deriv " .
                       " -step $si_step " .
                       " -fitting $ni $n_per $tolerance " .
                       " -ftol $f_tolerance " .
                       " -stop $stop_threshold $stop_iters ";
            }
            $once = 1;

            $ret = system_call( "$command", 1 );

            system_call( "measure_surface_area $output" );

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

        if( $anchor_weight > 0 )
        {
            unlink( $tmp_anchor );
        }
    }

    print( "Surface extraction finished.\n" );

    clean_up();
