#!xPERLx -w

# correct a classified volume: putamen is classified as WM.
# partial volume estimation: generate a skeletonized volume of CSF.
#
# Authors: June-sic Kim <luck3d@bic.mni.mcgill.ca>
#
# May 2004

use strict;
use MNI::Startup;
use Getopt::Tabular;
use MNI::Spawn;
use MNI::DataDir;
use MNI::FileUtilities qw(check_output_dirs);

# ===== Global Variables =====
my ($usage, $help);
my ($skelCSF, $cls, $cls_correct, $final);
my ($objMask, $unmasked);
my ($clsMasked, $filledImage);
my ($blurred, $curvature, $cortical_mask, $subcortical, $pve_classify);

# ===== Defaults =====
$unmasked = 0;

# ===== Argument Processing =====

$usage = "$ProgramName [options] classified.mnc final.mnc classified_correct.mnc skeletonized_csf.mnc\n";
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
     ["-subcortical", "string", 1, \$subcortical,
      "Mask to be used to correct classified file (putamen)."],
     );
GetOptions(\@argTbl, \@ARGV, \@leftOverArgs) or die "\n";

$cls = shift @leftOverArgs or die $usage;
$final = shift @leftOverArgs or die $usage;
$cls_correct = shift @leftOverArgs or die $usage;
$skelCSF = shift @leftOverArgs or die $usage;

# register the programs
RegisterPrograms(["minccalc", "mincresample", "pve2", "discretize_pve.pl",
                  "skel", "mincblur", "surface_mask2", "cortical_surface",
                  "make-curvature-volume"]);

if ($Clobber) {
    AddDefaultArgs("minccalc", ["-clobber"]);
}

# create necessary tmp directory
check_output_dirs($TmpDir);
$cortical_mask = "${TmpDir}/mask.mnc";
$filledImage = "${TmpDir}/filled_white.mnc";
Spawn(["minccalc", "-expression", 'out=1;', $cls, $filledImage]);

# ===== Main program execution =====
if ($unmasked) {
    # compute a mask
    $objMask = "${TmpDir}cortex.obj";
    Spawn(["cortical_surface", $cls, $objMask]);
    Spawn(["surface_mask2", "-binary_mask", $filledImage, $objMask, $cortical_mask]);
}
if ($objMask) {
    # apply the mask
    $clsMasked = "${TmpDir}/masked_cls.mnc";
    Spawn(["surface_mask2", $cls, $objMask, $clsMasked]);
    Spawn(["surface_mask2", "-binary_mask", $filledImage, $objMask, $cortical_mask]);
}
else {
    # the default assumption that the input classified file is already masked
    $clsMasked = $cls;
    Spawn(["minccalc", "-expression", 'if(A[0]>0){out=1;}else{out=0;}', $clsMasked, $cortical_mask]);
}
if ($subcortical) {
}
else {
    $subcortical = "/home/bic/vsingh/minc_files/ipt2.mnc";
}

# correct the classified image (putamen-->WM)
# first, create a filled image
$filledImage = "${TmpDir}/filled.mnc";
Spawn(["minccalc", "-expression", 'out = 0;', $clsMasked, $filledImage]);

# create the blurred image
$blurred = "${TmpDir}/fwhm_4";
Spawn(["mincblur", "-fwhm", "4", $final, $blurred]);
$blurred = "${TmpDir}/fwhm_4_blur.mnc";

# create the curvature map
$curvature = "${TmpDir}/fwhm_4";
Spawn(["make-curvature-volume", $blurred, $curvature]);
$curvature = "${TmpDir}/fwhm_4_k1.mnc";

# run partial volume estimation
$pve_classify = "${TmpDir}/pve_cls";
Spawn(["pve2", "-image", "$clsMasked", "-curve", $curvature, "-subcortical", $subcortical, "-mask", $cortical_mask, $final, $pve_classify]);

# skeletonize CSF image
Spawn(["skel", "${pve_classify}_csf.mnc", $skelCSF]);

# discretize PVE volume
Spawn(["discretize_pve.pl", "${pve_classify}_csf.mnc", "${pve_classify}_gm.mnc", "${pve_classify}_wm.mnc", "${pve_classify}_sc.mnc", $cls_correct]);
