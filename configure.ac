dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Herwig++],[SVN],[herwig@projects.hepforge.org],[Herwig++])
AC_CONFIG_SRCDIR([UnderlyingEvent/UA5Handler.h])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_HEADERS([Config/config.h])

AC_CANONICAL_HOST

case "${host}" in
  *-darwin[[0156]].*) 
    AC_MSG_ERROR([Herwig++ requires OS X 10.3 or later])
    ;;
  *-darwin*) 
    if test "x$MACOSX_DEPLOYMENT_TARGET" != "x10.3"; then 
      AC_MSG_ERROR(
      [Please add 'MACOSX_DEPLOYMENT_TARGET=10.3' to the configure line.]) 
    fi 
    ;;
esac

dnl === disable debug symbols by default =====
if test "x$CXXFLAGS" = "x"; then
   CXXFLAGS=-O2
fi
if test "x$CFLAGS" = "x"; then
   CFLAGS=-O2
fi
dnl Looptools manual requires no optimization
if test "x$FFLAGS" = "x"; then
   FFLAGS=-O0
fi
dnl ==========================================

AC_LANG(C++)

AM_INIT_AUTOMAKE([1.9])
dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LN_S

## needed for looptools
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_PROG_CC
##

## can be removed eventually
AC_DISABLE_STATIC

AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

dnl ####################################
dnl ####################################

AC_PATH_PROG(PERL, perl)

HERWIG_CHECK_GSL

HERWIG_CHECK_THEPEG

HERWIG_LOOPTOOLS

HERWIG_PDF_PATH

HERWIG_CHECK_HEPMC

HERWIG_CHECK_KTJET

HERWIG_CHECK_FASTJET

HERWIG_CHECK_ROOT

HERWIG_CHECK_EVTGEN

HERWIG_VERSIONSTRING

HERWIG_COMPILERFLAGS

AM_CONDITIONAL(NEED_APPLE_FIXES,
		[test "xx${host/darwin/foundit}xx" != "xx${host}xx"])

AC_CONFIG_FILES([UnderlyingEvent/Makefile
		Models/Makefile
		Models/StandardModel/Makefile
		Models/RSModel/Makefile
		Models/General/Makefile
		Models/Susy/Makefile
		Models/Susy/NMSSM/Makefile
		Models/UED/Makefile
		Interfaces/Makefile
		Decay/Makefile
		Decay/FormFactors/Makefile
		Decay/Tau/Makefile
		Decay/Baryon/Makefile
		Decay/VectorMeson/Makefile
		Decay/Perturbative/Makefile
                Decay/ScalarMeson/Makefile
                Decay/TensorMeson/Makefile
                Decay/WeakCurrents/Makefile
		Decay/Partonic/Makefile
		Decay/General/Makefile	
		Decay/Radiation/Makefile
                Doc/refman.conf
		PDT/Makefile
		PDF/Makefile
                MatrixElement/Makefile
		Shower/SplittingFunctions/Makefile
		Shower/Couplings/Makefile
		Shower/Default/MECorrections/Makefile
		Shower/Default/Makefile
		Shower/Base/Makefile
		Shower/CKKW/Clustering/Makefile
		Shower/CKKW/Reweighting/Makefile
		Shower/CKKW/Vertices/Makefile
		Shower/CKKW/Processes/Makefile
		Shower/CKKW/Qtilde/Makefile
		Shower/CKKW/Makefile
                Shower/Makefile
                Utilities/Makefile
		Utilities/Histogram2/Makefile
                Hadronization/Makefile
		lib/Makefile
		include/Makefile
		src/Makefile
		src/defaults/Makefile
		src/herwig-config
		Doc/Makefile
		Doc/HerwigDefaults.in
		Looptools/Makefile
		Analysis/Makefile
		Analysis/BSM/Makefile
		Analysis2/Makefile
                Analysis2/Data/Makefile
		src/Makefile-UserModules
		src/defaults/Analysis.in
		Makefile])

AC_CONFIG_FILES([Doc/fixinterfaces.pl],[chmod +x Doc/fixinterfaces.pl])

AC_OUTPUT
