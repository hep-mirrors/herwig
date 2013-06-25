dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([Herwig++],[devel],[herwig@projects.hepforge.org],[Herwig++])
AC_CONFIG_SRCDIR([Utilities/HerwigStrategy.cc])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([Config/config.h])
dnl AC_PRESERVE_HELP_ORDER
AC_CANONICAL_HOST

case "${host}" in
  *-darwin[[0156]].*) 
    AC_MSG_ERROR([Herwig++ requires OS X 10.3 or later])
    ;;
  *-darwin7.*) 
    if test "x$MACOSX_DEPLOYMENT_TARGET" != "x10.3"; then 
      AC_MSG_ERROR(
      [Please add 'MACOSX_DEPLOYMENT_TARGET=10.3' to the configure line.]) 
    fi 
    ;;
esac

dnl === disable debug symbols by default =====
if test "x$CXXFLAGS" = "x"; then
   CXXFLAGS=-O3
fi
if test "x$CFLAGS" = "x"; then
   CFLAGS=-O3
fi
dnl Looptools manual requires optimization off
if test "x$FCFLAGS" = "x"; then
   FCFLAGS=-O0
fi
dnl ==========================================

AC_LANG([C++])

AM_INIT_AUTOMAKE([1.9 gnu dist-bzip2 -Wall])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

dnl Checks for programs.
AC_PROG_CXX([g++])
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LN_S

dnl modified search order
AC_PROG_FC([gfortran g95 g77]) 
dnl xlf95 f95 fort ifort ifc efc pgf95 lf95 ftn xlf90 f90 pgf90 pghpf epcf90 xlf f77 frt pgf77 cf77 fort77 fl32 af77])
AC_LANG_PUSH([Fortran])
AC_MSG_CHECKING([if the Fortran compiler ($FC) works])
AC_COMPILE_IFELSE(
   	AC_LANG_PROGRAM([],[      print *[,]"Hello"]),
	[AC_MSG_RESULT([yes])],
	[AC_MSG_RESULT([no])
	 AC_MSG_ERROR([A Fortran compiler is required to build Herwig++.])
	]
)
AC_LANG_POP([Fortran])


LT_PREREQ([2.2])
LT_INIT([disable-static dlopen pic-only])

dnl ####################################
dnl ####################################

dnl for Doc/fixinterfaces.pl
AC_PATH_PROG(PERL, perl)

HERWIG_CHECK_GSL

HERWIG_CHECK_THEPEG

BOOST_REQUIRE
BOOST_FIND_HEADER([boost/array.hpp])
BOOST_FIND_HEADER([boost/numeric/ublas/io.hpp])
BOOST_FIND_HEADER([boost/numeric/ublas/matrix.hpp])
BOOST_FIND_HEADER([boost/numeric/ublas/matrix_proxy.hpp])
BOOST_FIND_HEADER([boost/numeric/ublas/matrix_sparse.hpp])
BOOST_FIND_HEADER([boost/numeric/ublas/symmetric.hpp])
BOOST_FIND_HEADER([boost/numeric/ublas/vector.hpp])
BOOST_FIND_HEADER([boost/operators.hpp])
BOOST_FIND_HEADER([boost/progress.hpp])
BOOST_FIND_HEADER([boost/scoped_array.hpp])
BOOST_FIND_HEADER([boost/scoped_ptr.hpp])
BOOST_FIND_HEADER([boost/utility.hpp])

HERWIG_COMPILERFLAGS

HERWIG_LOOPTOOLS

HERWIG_PDF_PATH

FASTJET_CHECK_FASTJET

HERWIG_VERSIONSTRING

HERWIG_CHECK_ABS_BUG

HERWIG_ENABLE_MODELS

HERWIG_ENABLE_DIPOLE

SHARED_FLAG=-shared
AM_CONDITIONAL(NEED_APPLE_FIXES,
		[test "xx${host/darwin/foundit}xx" != "xx${host}xx"])
if test "xx${host/darwin/foundit}xx" != "xx${host}xx"; then
   APPLE_DSO_FLAGS=-Wl,-undefined,dynamic_lookup
   SHARED_FLAG=-bundle
fi
AC_SUBST([APPLE_DSO_FLAGS])
AC_SUBST([SHARED_FLAG])

AC_CONFIG_FILES([UnderlyingEvent/Makefile
		Models/Makefile
		Models/StandardModel/Makefile
		Models/RSModel/Makefile
		Models/General/Makefile
		Models/Susy/Makefile
		Models/Susy/NMSSM/Makefile
		Models/Susy/RPV/Makefile
		Models/UED/Makefile
		Models/LH/Makefile
		Models/LHTP/Makefile
		Models/Transplanckian/Makefile
		Models/Leptoquarks/Makefile
		Models/Zprime/Makefile
		Models/TTbAsymm/Makefile
		Models/ADD/Makefile
		Models/Sextet/Makefile
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
                Doc/refman.h
		PDT/Makefile
		PDF/Makefile
                MatrixElement/Makefile
                MatrixElement/General/Makefile
                MatrixElement/Lepton/Makefile
                MatrixElement/Hadron/Makefile
                MatrixElement/DIS/Makefile
                MatrixElement/Powheg/Makefile
                MatrixElement/Gamma/Makefile
                MatrixElement/Matchbox/Makefile
                MatrixElement/Matchbox/Base/Makefile
                MatrixElement/Matchbox/Utility/Makefile
                MatrixElement/Matchbox/Phasespace/Makefile
                MatrixElement/Matchbox/Dipoles/Makefile
                MatrixElement/Matchbox/InsertionOperators/Makefile
                MatrixElement/Matchbox/Matching/Makefile
                MatrixElement/Matchbox/Builtin/Makefile
                MatrixElement/Matchbox/Builtin/Amplitudes/Makefile
                MatrixElement/Matchbox/External/Makefile
                MatrixElement/Matchbox/External/BLHAGeneric/Makefile
                MatrixElement/Matchbox/Tests/Makefile
		Exsample2/Makefile
		Shower/SplittingFunctions/Makefile
		Shower/Default/Makefile
		Shower/Base/Makefile
                Shower/Makefile
		DipoleShower/Makefile
		DipoleShower/Base/Makefile
		DipoleShower/Kernels/Makefile
		DipoleShower/Kinematics/Makefile
		DipoleShower/Utility/Makefile
		DipoleShower/AlphaS/Makefile
                Utilities/Makefile
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
		src/Makefile-UserModules
		src/defaults/Analysis.in
		Contrib/Makefile
		Contrib/make_makefiles.sh
		Tests/Makefile
		Makefile])

AC_CONFIG_FILES([Doc/fixinterfaces.pl],[chmod +x Doc/fixinterfaces.pl])

HERWIG_OVERVIEW

AC_CONFIG_COMMANDS([summary],[cat config.herwig])

AC_OUTPUT
