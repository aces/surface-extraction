#!xPERLx -w

use strict;

require "deform_utils.pl";
use Getopt::Tabular;
use MNI::Startup;
use MNI::FileUtilities;
use MNI::Spawn;
use MNI::DataDir;
use MNI::FileUtilities qw(check_output_dirs);

MNI::Spawn::RegisterPrograms
  ( [qw/  rm
     mincmorph
     minccalc
     subdivide_polygons
     new_fit_3d_1
     make_gradient_volume
     set_object_colour/ ] )
  or exit 1;

my ($i, $step, $size, $sw, $cw, $n_iters, $iter_inc);
my ($si_step, $oversample, $self_weight, $self_dist);
my ($laplacian_sampling, $laplacian_factor, $laplacian_weight);
my ($self2, $surf2_info, $n_failures, $iter, $ni, $command, $ret, $b2);
my ($laplace_info);

# --- set the help & usage strings ---
my $help = <<HELP;
[refine]   : give a value >0 to generate the gray surface with 327680 triangles

HELP

my $usage = <<USAGE;
usage: $ProgramName cls.mnc white.obj gray.obj field.mnc [refine]
       $ProgramName -help to list options
USAGE

Getopt::Tabular::SetHelp( $help, $usage );

# --- process options ---
my( $masking_surface, $masked_input, $output_prefix );
my @options = 
  ( @DefaultArgs,     # from MNI::Startup
  );

GetOptions( \@options, \@ARGV ) 
  or exit 1;
die "$usage\n" unless @ARGV >= 4;

    my $cls = shift;
    my $white_surface = shift;
    my $input = shift;
    my $laplacian_file = shift;

    my $refine = shift;

    my $logfile = shift;
    my $start_n = shift;
    my $end_n = shift;
    my $dont_copy = shift;
my $output = $input;
    if( ! defined($laplacian_file) )
    {
        die "Usage: $0  input.mnc white.obj input_gray.obj output1_pref field.mnc [log_file] [start] [end] [dont copy]\n";
    }

#---------
    check_output_dirs($TmpDir);

    my $fit = "new_fit_3d ";

    my $oversample_reference = 20480;
    my $n_polygons = `print_n_polygons $white_surface`;
    chop( $n_polygons );

    my $model_data_dir = MNI::DataDir::dir('CLASP');
    MNI::DataDir::check_data($model_data_dir, ["ellipsoid_${n_polygons}.obj.gz"]);

    my $model = $white_surface;
    $model = "${model_data_dir}/ellipsoid_${n_polygons}.obj.gz";

    my $self_dist2 = 0.01;
    my $n_selfs = 9;
    my $self_factor = 1.0;

    my $stop_threshold = 3e-2;
    my $stop_iters = 10;

    my $filter = 0;
    my $n_per = 1;
    my $tolerance = 1.0e-2;
    $tolerance = 1.0e-10;
    my $f_tolerance = 1.0e-2;
    $f_tolerance = 1.0e-10;

    my $iters_scale = 1.0;
    my $break_scale = 1.0;
    my $oo_scale = 1.0 * sqrt( $oversample_reference / $n_polygons );
    my $iters_override = 0;

    my $stretch_scale = 1;
    my $curvature_scale = 0;

    if( ! defined($refine) )
        { $refine = 0; }

    my @schedule;
    if( $refine <= 81920 )
    {
      @schedule = (
#size   sw        n_it  inc  si over   sw   self   l_s  l_d   l_w
#----- ----       ----  ---  -- ----  ----  ----   ---  ---  ----
 81920, 1e1,   0,  200,  50,  .5,  1,  1e1,  .25, , 1.0, 0.0, 5e-6, 
 81920, 1e1,   0,  200,  50,  .5,  1,  1e1,  .25, , 1.0, 0.0, 1e-4, 
 81920, 1e1,   0,  300,  50,  .5,  1,  1e1,  .25, , 1.0, 0.0, 5e-4, 
 81920, 1e1,   0,  800,  50,  .5,  1,  1e1,  .25, , 1.0, 1.0, 1e-2,
  );
    }
    else
    {
      @schedule = (
#size   sw        n_it  inc  si over   sw   self   l_s  l_d   l_w
#----- ----       ----  ---  -- ----  ----  ----   ---  ---  ----
327680, 1e1,   0,  200,  50,  .5,  1,  1e1,  .25, , 1.0, 0.0, 5e-6, 
327680, 1e1,   0,  200,  50,  .5,  1,  1e1,  .25, , 1.0, 0.0, 1e-5, 
327680, 1e1,   0,  300,  50,  .5,  1,  1e1,  .25, , 1.0, 0.0, 1e-4, 
327680, 1e1,   0,  800,  50,  .5,  1,  1e1,  .25, , 1.0, 1.0, 1e-3,
327680, 1e1,   0,  500,  50,  .5,  1,  1e1,  .25, , 1.0, 1.0, 2e-3,
  );
    }


    my $sched_size =  12;

    my $log = "";
    if( defined($logfile) )
        { $log = " -log $logfile"; }

    if( ! defined($start_n) )
        { $start_n = 0; }

    if( ! defined($end_n) )
        { $end_n = @schedule / $sched_size - 1; }

#------ loop over each schedule

    if( !$dont_copy && $start_n <= 0 ) {
        Spawn(["set_object_colour", $white_surface, $output, "white"]);
    }

    my $prev_n = $schedule[0];
    my $movie_num = 1;
    for( $i = 0;  $i < @schedule;  $i += $sched_size )
    {
        $step = $i / $sched_size;

        #--- get the components of the deformation schedule entry

        ( $size, $sw, $cw, $n_iters, $iter_inc,
          $si_step, $oversample, $self_weight, $self_dist, 
          $laplacian_sampling, $laplacian_factor, $laplacian_weight ) =
                     @schedule[$i..$i+$sched_size-1];

        $sw *= $stretch_scale;
        $cw *= $curvature_scale;
        $oversample *= $oo_scale;
        #$value_oversample = $oversample * $value_oo_scale;

        if( $iters_override > 0 )
            { $n_iters = $iters_override; }
        else
            { $n_iters = int( $n_iters * $iters_scale ); }

        $self_weight *= $self_factor;

        $self2 = get_self_intersect( $self_weight, $n_selfs, $self_dist,
                                     $self_dist2 );

        #--- if the schedule size is greater than the current number of
        #--- polygons in the deforming surface, subdivide the deforming surface
        if( $step < $start_n || $step > $end_n )
            { $prev_n = $size;  next; }
        if( $size != $prev_n && (!$dont_copy || $step > $start_n) )
        {
          Spawn(["subdivide_polygons", $white_surface, "${TmpDir}/white_${size}.obj", $size]);
          Spawn(["subdivide_polygons", $output, "${output}_${size}.obj", 
                  $size]);
          $input = $output = "${output}_${size}.obj";
          $white_surface = "${TmpDir}/white_${size}.obj";
        }
        $prev_n = $size;

        print( "Fitting polygons at filter ${filter}, " .
               "max $n_iters iters.\n" );

        $b2 = " -boundary 1 1 $cls " .
          " 1.5 - 0 0 0 0 $oversample";
        $laplace_info = " -laplacian $laplacian_file $laplacian_weight 0 10 $laplacian_factor $laplacian_sampling ";

        $iter_inc *= $break_scale;
        if( $iter_inc <= 0 )  { $iter_inc = $n_iters; }

        $surf2_info = " -surface ${input} ${output} ${white_surface}" .
          " -equal_lengths".
          " -stretch $sw ${input} -1.0 0 0 0".
          " $b2 ".
          " $self2 ".
          " $laplace_info ";

        $n_failures = 0;

        for( $iter = 0;  $iter < $n_iters;  $iter += $iter_inc )
        {
            system( "echo Step: $iter / $n_iters    $sw $cw" );
	
            $ni = $n_iters - $iter;
            if( $ni > $iter_inc )  { $ni = $iter_inc; }

            $command = "$fit -mode three $surf2_info ".
                       " -print_deriv " .
                       " -step $si_step " .
                       " -fitting $ni $n_per $tolerance " .
                       " -ftol $f_tolerance " .
                       " -stop $stop_threshold $stop_iters ".
                       " $log ";

            $ret = system_call( "$command", 1 );

            system_call( "measure_surface_area $output" );
            #$string_num = sprintf "%04s", ${movie_num};
            $movie_num = ${movie_num} + 1;
            #system_call( "cp -f ${output} movie/gray_${string_num}.obj" );

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

    clean_up();

