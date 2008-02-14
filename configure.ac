# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

# The arguments are package name, and package version.
AC_INIT([surface-extraction],[2.3.2], [Andrew Janke <a.janke@gmail.com>])
AC_CONFIG_SRCDIR([src/surface_fit.c])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

mni_REQUIRE_BICPL

AC_PATH_PROGS(PERL, [perl], no)

# Finish up by writing output files.
AC_OUTPUT(
scripts/utils.pl
scripts/deform_utils.pl
scripts/asp
scripts/calibrate_white
scripts/clasp
scripts/classify_correct
scripts/cortical_thickness
scripts/expand_from_white
scripts/extract_white_surface
scripts/make_asp_grid
scripts/mask_cortical_white_matter
scripts/pve_defrag
scripts/surface_qc
Makefile
src/Makefile
epm-header
)

