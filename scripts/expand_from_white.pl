#!xPERLx -w

# takes as input a white matter surface and a classified volume, then
# expands out from the white matter surface to find the grey matter
# surface.

# Author: David MacDonald
# Modified: Jason Lerch <jason@bic.mni.mcgill.ca>

require "xINCDIRx/deform_utils.pl";

use MNI::DataDir;
use MNI::Startup;
use Getopt::Tabular;

#    &print_out_script_info();

my ($volume, $white_surface, $output);
my ($start_n, $end_n, $dont_copy, $isovalue);
my ($schedule_file, $model_data_dir);

# defaults:
$schedule_file = "xSYSCONFx/expand_from_white.cfg";
$model_data_dir = MNI::DataDir::dir('ASP');

$isovalue = 1.5;

# get command line arguments and set defaults
sub Initialise {

    my @options = 
	( @DefaultArgs, # from MNI::Startup
	  [ "-start_n", "int", 1, \$start_n,
	    "Schedule number to start iterating at." ],
	  [ "-end_n", "int", 1, \$end_n,
	    "Schedule number to end iterating at." ],
	  [ "-schedule", "string", 1, \$schedule_file,
	    "Schedule to use. [Default: $schedule_file]" ],
	  [ "-isovalue", "float", 1, \$isovalue,
	    "Isovalue of volume to fit. [Default: $isovalue]" ],
	  [ "-modeldir", "string", 1, \$model_data_dir,
	    "Location of model data. [Default: $model_data_dir]" ],
	  );
    
    GetOptions( \@options, \@ARGV ) or die "\n";

    $volume = shift;
    $white_surface = shift;
    $output = shift;
    $dont_copy = shift;

    MNI::DataDir::check_data($model_data_dir, 
			     ["ellipsoid_${n_polygons}.obj.gz"]);

    $model = $white_surface;
    $model = "${model_data_dir}/ellipsoid_${n_polygons}.obj.gz";


#      if( ! defined($isovalue) )
#      {
#          die "Usage: $0  input.mnc white.obj output1_pref [start] [end] [dont copy]\n";
#      }
}

#---------
    $fit = "new_fit_3d ";

    $oversample_reference = 20480;
    $n_polygons = `print_n_polygons $white_surface`;
    chop( $n_polygons );

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

my @schedule;

open SCHEDULE, $schedule_file 
    or die "ERROR: could not open schedule file $schedule_file\n";
while (<SCHEDULE>) {
    my $line = $_;
    chomp $line;
    unless ($line ~= /^\#/) { # ignore comment lines
	my @fields = split /,/, $line;
	foreach my $field (@fields) {
	    $field =~ s/\s//g;
	    push @schedule, $field;
	}
    }
}

print "SCHEDULE: @schedule\n";

    $sched_size =  14;
 
    if( ! defined($start_n) )
        { $start_n = 0; }

    if( ! defined($end_n) )
        { $end_n = @schedule / $sched_size - 1; }

#------ loop over each schedule

    if( !$dont_copy && $start_n <= 0 ) {
        system_call( "set_object_colour $white_surface ${output} white" );
    }

    for( $i = 0;  $i < @schedule;  $i += $sched_size )
    {
        $step = $i / $sched_size;

        #--- get the components of the deformation schedule entry

        ( $sw, $cw, $n_iters, $iter_inc, $offo, $offi,
          $si_step, $oversample, $self_weight, $self_dist, $desired_dist,
          $anchor_weight, $min_dist, $max_dist ) =
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

        $surf2_info = " -surface ${output} ${output} " .
          " -stretch $sw $model -.9 0 0 0".
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

            $ni = $n_iters - $iter;
            if( $ni > $iter_inc )  { $ni = $iter_inc; }

            $command = "$fit $surf2_info ".
                       " $anchor " .
                       " -print_deriv " .
                       " -step $si_step " .
                       " -fitting $ni $n_per $tolerance " .
                       " -ftol $f_tolerance " .
                       " -stop $stop_threshold $stop_iters ";

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
