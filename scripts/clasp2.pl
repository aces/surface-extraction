#!xPERLx -w
#
# Based on ~david/Surface_deformation/How_to_extract_cortical_surfaces.txt

use strict;

use Getopt::Tabular;
use MNI::Startup;
use MNI::FileUtilities;
use MNI::Spawn;

# Several of the sub-programs will look for static files using
# MNI::DataDir.  Set the directory here if the user does not have
# their own.
#
$ENV{MNI_DATAPATH} = 'xDATADIRx'
  unless defined $ENV{MNI_DATAPATH};


MNI::Spawn::RegisterPrograms
  ( [qw/  rm
     mv
     minccalc
     mincresample
     dilate_volume
     cortical_surface
     discretize_pve
     pve_script
     pve_curvature
     skel
     mask_cortical_white_matter
     surface_mask2 mask_volume
     extract_white_surface
     classify_correct
     calibrate_white
     make_asp_grid
     expand_from_white/ ] )
  or exit 1;


# --- set the help & usage strings ---
my $help = <<HELP;
The help string goes here.

HELP

my $usage = <<USAGE;
usage: $ProgramName [options] classified.mnc final.mnc [mask.obj]
       $ProgramName -help to list options
USAGE

Getopt::Tabular::SetHelp( $help, $usage );

# --- process options ---
my( $masking_surface, $masked_input, $output_prefix );
my( $laplace, $cls_correct, $skelCSF, $white );
my( $Mask );
my @options = 
  ( @DefaultArgs,     # from MNI::Startup
    ['-out', 'string', 1, \$output_prefix,
     'prefix for all output files'],
  );

GetOptions( \@options, \@ARGV ) 
  or exit 1;
die "$usage\n" unless @ARGV == 2 or @ARGV == 3;

# $input_volume should be a CLASSIFIED MINC volume
my $input_volume = shift;
my $t1_volume = shift;
my $objMask = shift;


# Masking surface and cortical-matter masked version of input.
# These are temporary files.
MNI::FileUtilities::check_output_path("${TmpDir}")
  or exit 1;
$masking_surface = "${TmpDir}input-mask.obj"
  unless defined $masking_surface;
$masked_input = "${TmpDir}input-masked"
  unless defined $masked_input;

# Output filename prefix
if ( !defined $output_prefix ) {
    my( $dir, $base, $ext ) = MNI::PathUtilities::split_path($input_volume);
    $output_prefix = $base;
}

# Generate mask image
my $Filled = "${TmpDir}/Filled.mnc";
if ( !defined $objMask ) {
    $Mask = "${TmpDir}/cortical_mask.mnc";
    Spawn(["cortical_surface", $input_volume, "${output_prefix}_cortex.obj"]);
    Spawn(["minccalc", "-expression", 'out=1;', $input_volume, $Filled]);
    Spawn(["surface_mask2", "-binary_mask", $Filled, "${output_prefix}_cortex.obj", $Mask]);
    Spawn(["dilate_volume", $Mask, "${TmpDir}/cortical_mask_tmp.mnc", "1", "26", "3"]);
    Spawn(["mv", "-f", "${TmpDir}/cortical_mask_tmp.mnc", $Mask]);
    Spawn(["rm", "-f", $Filled]);
}
else{
    Spawn(["minccalc", "-expression", 'out=1;', $input_volume, $Filled]);
    Spawn(["surface_mask2", "-binary_mask", $Filled, $objMask, "${TmpDir}/cortical_mask.mnc"]);
    Spawn(["rm", "-f", $Filled]);
    $Mask = "${TmpDir}/cortical_mask.mnc";
    Spawn(["dilate_volume", $Mask, "${TmpDir}/cortical_mask_tmp.mnc", "1", "26", "3"]);
    Spawn(["mv", "-f", "${TmpDir}/cortical_mask_tmp.mnc", $Mask]);
}
Spawn(["mincresample", "-clobber", "-like", $input_volume, $Mask, "${TmpDir}/mask.mnc"]);
Spawn(["mv", "-f", "${TmpDir}/mask.mnc", $Mask]);
Spawn(["minccalc", "-expression", 'if(A[1]>0){out=A[0];}else{out=0;}', $input_volume, $Mask, "${TmpDir}/cls_masked.mnc"]);
$input_volume = "${TmpDir}/cls_masked.mnc";


#---------------------------------------------------------------------------
#  Step 1:   2 hours    Correcting a cls volume to classify putamen as WM
#                       Partial volume estimation
#---------------------------------------------------------------------------

# Partial volume estimation using pve2
$cls_correct = "${TmpDir}/cls_correct.mnc";
$skelCSF = "${TmpDir}skel_CSF.mnc";
Spawn(["pve_curvature", $t1_volume, $input_volume, $Mask, "${TmpDir}/curve"]);
Spawn(["pve_script", "-curve", "${TmpDir}/curve_cg.mnc", "-subcortical", $t1_volume, "${TmpDir}/pve", "-mask", $Mask, "-image", $input_volume]);
Spawn(["skel", "${TmpDir}/pve_csf.mnc", $skelCSF]);
Spawn(["discretize_pve", "${TmpDir}/pve_csf.mnc", "${TmpDir}/pve_wm.mnc", "${TmpDir}/pve_gm.mnc", "${TmpDir}/pve_sc.mnc", $cls_correct]);
#$input_volume = $cls_correct;


#---------------------------------------------------------------------------
#  Step 2:   1/2 hour   Creating a mask surface to chop of non-cortical white
#                       matter
#---------------------------------------------------------------------------

# Deforms from model file `white_matter_mask.obj' using new_fit_3d
Spawn( "mask_cortical_white_matter $cls_correct  $masking_surface  2.5");


#---------------------------------------------------------------------------
#  Step 3:   2 minutes  Masking the classified volume with the masking_surface
#---------------------------------------------------------------------------

# conglomerate routine
Spawn("surface_mask2  $cls_correct  $masking_surface  $masked_input");
#Spawn("rm -f $masking_surface");

# Set all voxels with value in range (1,2.5) to zero
# conglomerate routine
Spawn("mask_volume    $masked_input $masked_input $masked_input 1 2.5 0");


#---------------------------------------------------------------------------
#  Step 4:   20 hours   Shrink wrapping a sphere to the MASKED white matter
#                       matter,  (note that the second argument below is a
#                                 prefix)
#                       creating white_surface_{320,1280,5120,20480,81920}.obj
#                       you can delete the model files that are created:
#                            white_surface_m_{320,1280,5120,20480,81920}.obj
#---------------------------------------------------------------------------

# Deforms from model file `white_model_320.obj' using new_fit_3d
Spawn("extract_white_surface  $masked_input   ${output_prefix}_white  2.5");
#Spawn("rm -f $masked_input");

#---------------------------------------------------------------------------
#  Step 5:   x hours    Calibrate white_surface wite a gradient field
#---------------------------------------------------------------------------

$white = "${output_prefix}_white.obj";
Spawn(["calibrate_white", $t1_volume, $input_volume, $skelCSF, "${output_prefix}_white_81920.obj", $white]);


#---------------------------------------------------------------------------
#  Step 6:   1 hour     Create a Laplacian field from the WM surface to the
#                       outer boundary of gray matter
#---------------------------------------------------------------------------

$laplace = "${TmpDir}laplace.mnc";
Spawn(["make_asp_grid", $skelCSF, $white, $cls_correct, $laplace]);


#---------------------------------------------------------------------------
#  Step 7:   20 hours   Expand a copy of the white surface out to the gray
#                       boundary.  Note that now we can delete the temporary
#                       masked white matter ($masked_input), and use the
#                       original classified volume for this step.
#---------------------------------------------------------------------------

# Model file is ellipsoid_${n_polygons}.obj.gz
# (n_polygons will be 81920)
#Spawn("expand_from_white  $input_volume ${output_prefix}_white_81920.obj"
#      . " ${output_prefix}_gray_81920.obj  1.5");
#Spawn(["expand_from_white", $input_volume, $white, "${output_prefix}_gray_81920.obj", "1.5", $laplace]);
Spawn(["expand_from_white", $cls_correct, $white, "${output_prefix}_gray_81920.obj", $laplace]);
