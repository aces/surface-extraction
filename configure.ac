dnl Process this file with autoconf to produce a configure script.

AC_INIT
AC_CONFIG_SRCDIR([src/new_fit_3d.c])
AM_INIT_AUTOMAKE(ASP, 1.0)
AM_CONFIG_HEADER(config.h)

smr_WITH_BUILD_PATH

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PROG_CC
AC_PROG_CXX

mni_REQUIRE_BICPL

AC_PATH_PROGS(PERL, perl5 perl)

AC_CONFIG_FILES([Makefile src/Makefile scripts/Makefile models/Makefile])
AC_OUTPUT

AC_MSG_RESULT([====================])
AC_MSG_RESULT([Once the program has built make sure that the model files are accessible through the MNI_DATAPATH variable])