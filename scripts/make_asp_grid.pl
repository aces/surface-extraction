#!xPERLx -w

# creates a grid from the classified image. This grid can then
# be used to initialise the solving of Laplace's equation.
#
# Authors: June-sic Kim <luck3d@bic.mni.mcgill.ca> and
#          Jason Lerch <jason@bic.mni.mcgill.ca>
#
# Sep 2003

use strict;
use MNI::Startup;
use Getopt::Tabular;
use MNI::Spawn;
use MNI::DataDir;
use MNI::FileUtilities qw(check_output_dirs);

# ===== Global Variables =====
my ($usage, $help);
my ($skelCSF, $cls, $wmSurface, $grid, $laplace, $output);
my ($expression);
my ($objMask, $unmasked);
my ($clsMasked, $wmLine, $wmMask, $wmMask2, $filledImage, $rslCSF);

# ===== Defaults =====
$unmasked = 0;

# ===== Argument Processing =====

$usage = "$ProgramName [options] skeletonized_csf.mnc wm_surface.obj classified.mnc output.mnc\n";
$help = "Help still to be written";

my @leftOverArgs;
my @argTbl = 
    (
     @DefaultArgs,
     ["Masking options", "section"],
     ["-obj_mask", "string", 1, \$objMask,
      "Mask to be used to mask the classified file with. By default the ".
      "classified file is already masked."],
     ["-unmasked", "boolean", undef, \$unmasked,
      "Create a mask and apply it to the classified image. By default the " .
      "classified file is already masked."],
     );
GetOptions(\@argTbl, \@ARGV, \@leftOverArgs) or die "\n";

$skelCSF = shift @leftOverArgs or die $usage;
$wmSurface = shift @leftOverArgs or die $usage;
$cls = shift @leftOverArgs or die $usage;
$output = shift @leftOverArgs or die $usage;

# register the programs
RegisterPrograms(["minccalc", "mincresample", "scan_object_to_volume",
                  "surface_mask2", "cortical_surface", "laplacian_thickness"]);

if ($Clobber) {
    AddDefaultArgs("minccalc", ["-clobber"]);
}

# create necessary tmp directory
check_output_dirs($TmpDir);

# ===== Main program execution =====
if ($unmasked) {
    # compute a mask
    $objMask = "${TmpDir}/cortex.obj";
    Spawn(["cortical_surface", $cls, $objMask, 1.5]);
}
if ($objMask) {
    # apply the mask
    $clsMasked = "${TmpDir}/masked_cls.mnc";
    Spawn(["surface_mask2", $cls, $objMask, $clsMasked]);
}
else {
    # the default assumption that the input classified file is already masked
    $clsMasked = $cls;
}

# compute the WM lines
# first, create a filled image
$filledImage = "${TmpDir}/filled.mnc";
Spawn(["minccalc", "-expression", 'out = 0;', $clsMasked, $filledImage]);

# now add the lines to the filled image
$wmLine = "${TmpDir}/wm_lines.mnc";
Spawn(["scan_object_to_volume", $filledImage, $wmSurface, $wmLine]);

# create a binary white matter mask
$wmMask = "${TmpDir}/wm_mask.mnc";
$wmMask2 = "${TmpDir}/wm_mask2.mnc";
Spawn(["surface_mask2", "-binary_mask", $cls, $wmSurface, $wmMask2]);
Spawn(["mincresample", "-like", ,$cls, $wmMask2, $wmMask]);

# resample the CSF skel map to be like the classified map
$rslCSF = "${TmpDir}/csf_rsl.mnc";
Spawn(["mincresample", "-nearest_neighbour", "-like", 
       $clsMasked, $skelCSF, $rslCSF]);

# create the grid itself
$expression = 'if(A[2]>0 && A[3]==0){out=0;}else if(A[0]>0 && A[3]==0){out=10;}else if(A[1]==0){out=10;}else{out=5;}';
$grid = "${TmpDir}/grid.mnc";
Spawn(["minccalc", "-expression", $expression, $rslCSF, $clsMasked,
       $wmMask, $wmLine,$grid]);

# create the laplacian field
$laplace = "${TmpDir}/laplace.mnc";
Spawn(["laplacian_thickness", "-potential_only", "-volume-double", "-gradients-double", "-from_grid", $grid, "-convergence", "-0.1", "-max_iterations", "500", $output]);

# correct the field
Spawn(["mincresample", "-like", $cls, $output, $laplace]);
$expression = 'if(A[0]==0){out=-0.1;}else{out=A[1];}';
Spawn(["minccalc", "-clobber", "-expression", $expression, $grid, $laplace, $output]);
