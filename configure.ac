dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([Herwig],[7.0.4],[herwig@projects.hepforge.org],[Herwig])
AC_CONFIG_SRCDIR([Utilities/HerwigStrategy.cc])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([Config/config.h])
dnl AC_PRESERVE_HELP_ORDER
AC_CANONICAL_HOST

dnl === disable debug symbols by default =====
if test "x$CXXFLAGS" = "x"; then
   CXXFLAGS=-O3
fi
if test "x$CFLAGS" = "x"; then
   CFLAGS=-O3
fi

AC_LANG([C++])

AM_INIT_AUTOMAKE([1.11 subdir-objects gnu dist-bzip2 no-dist-gzip -Wall])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

dnl Checks for C++ compiler. Handle C++11 flags.
AC_PROG_CXX([g++])

dnl this block can go once we decide to have C++11 always on
AC_MSG_CHECKING([whether to include C++11 flag for testing])
AC_ARG_ENABLE(stdcxx11,
        AC_HELP_STRING([--enable-stdcxx11],
        [turn on C++11 flag (only for testing, do not use in production!)]),
        [],
        [enable_stdcxx11=no]
        )


if test "x$enable_stdcxx11" = "xyes"; then
AC_MSG_RESULT([yes])
dnl remove the wrapper if block and "optional" once we decide to have C++11 always on
AX_CXX_COMPILE_STDCXX([11],[noext],[optional])
if test "x$HAVE_CXX11" != "x1"; then
AC_MSG_ERROR([compiler does not recognize requested c++11 option])
fi
else
AC_MSG_RESULT([no])
fi

dnl check for POSIX
AC_CHECK_HEADER([unistd.h],[],
      [AC_MSG_ERROR([Herwig needs "unistd.h". Non-POSIX systems are not supported.])])

dnl Checks for programs.
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
	 AC_MSG_ERROR([A Fortran compiler is required to build Herwig.])
	]
)
AC_LANG_POP([Fortran])


LT_PREREQ([2.2.6])
LT_INIT([disable-static dlopen pic-only])

dnl ####################################
dnl ####################################

dnl for Doc/fixinterfaces.pl
AC_PATH_PROG(PERL, perl)

dnl for Models/Feynrules
AM_PATH_PYTHON([2.6],, [:])
AM_CONDITIONAL([HAVE_PYTHON], [test "x$PYTHON" != "x:"])

HERWIG_CHECK_GSL

HERWIG_CHECK_THEPEG

BOOST_REQUIRE([1.41])
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
BOOST_FILESYSTEM([mt])
BOOST_TEST()

HERWIG_CHECK_VBFNLO

HERWIG_CHECK_NJET

HERWIG_CHECK_GOSAM

HERWIG_CHECK_GOSAM_CONTRIB

HERWIG_CHECK_OPENLOOPS

HERWIG_CHECK_MADGRAPH

HERWIG_COMPILERFLAGS

HERWIG_LOOPTOOLS

HERWIG_PDF_PATH

FASTJET_CHECK_FASTJET

HERWIG_CHECK_ABS_BUG

HERWIG_ENABLE_MODELS

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
		Models/Feynrules/Makefile
		Models/Feynrules/python/Makefile-FR
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
                MatrixElement/Reweighters/Makefile
                MatrixElement/Matchbox/Makefile
                MatrixElement/Matchbox/Base/Makefile
                MatrixElement/Matchbox/Utility/Makefile
                MatrixElement/Matchbox/Phasespace/Makefile
                MatrixElement/Matchbox/Dipoles/Makefile
                MatrixElement/Matchbox/InsertionOperators/Makefile
                MatrixElement/Matchbox/Matching/Makefile
                MatrixElement/Matchbox/Cuts/Makefile
                MatrixElement/Matchbox/Scales/Makefile
                MatrixElement/Matchbox/Scales/MatchboxScale.cc
                MatrixElement/Matchbox/ColorFull/Makefile
                MatrixElement/Matchbox/CVolver/Makefile
                MatrixElement/Matchbox/Builtin/Makefile
                MatrixElement/Matchbox/Builtin/Amplitudes/Makefile
                MatrixElement/Matchbox/Tests/Makefile
                MatrixElement/Matchbox/External/Makefile
                MatrixElement/Matchbox/External/BLHAGeneric/Makefile
                MatrixElement/Matchbox/External/VBFNLO/Makefile
                MatrixElement/Matchbox/External/NJet/Makefile
		MatrixElement/Matchbox/External/GoSam/Makefile
                MatrixElement/Matchbox/External/OpenLoops/Makefile
                MatrixElement/Matchbox/External/MadGraph/Makefile
                MatrixElement/Matchbox/External/MadGraph/mg2herwig.py
		Sampling/Makefile
		Sampling/CellGrids/Makefile
                Shower/Makefile
		Shower/Matching/Makefile
		DipoleShower/Makefile
		DipoleShower/Base/Makefile
		DipoleShower/Kernels/Makefile
		DipoleShower/Kinematics/Makefile
		DipoleShower/Utility/Makefile
		DipoleShower/AlphaS/Makefile
                Utilities/Makefile
                Utilities/XML/Makefile
                Utilities/Statistics/Makefile
                Hadronization/Makefile
		lib/Makefile
		include/Makefile
		src/Makefile
		src/defaults/Makefile
		src/snippets/Makefile
		src/Matchbox/Makefile
		src/herwig-config
		Doc/Makefile
		Doc/HerwigDefaults.in
		Looptools/Makefile
		Analysis/Makefile
		API/Makefile
		src/Makefile-UserModules
		src/defaults/Analysis.in
		src/defaults/MatchboxDefaults.in
		src/defaults/setup.gosam.in
		src/Matchbox/LO-DefaultShower.in
		src/Matchbox/LO-DipoleShower.in
		src/Matchbox/MCatLO-DefaultShower.in
		src/Matchbox/MCatLO-DipoleShower.in
		src/Matchbox/LO-NoShower.in
		src/Matchbox/MCatNLO-DefaultShower.in
		src/Matchbox/MCatNLO-DipoleShower.in
		src/Matchbox/NLO-NoShower.in
		src/Matchbox/Powheg-DefaultShower.in
		src/Matchbox/Powheg-DipoleShower.in
		Contrib/Makefile
		Contrib/make_makefiles.sh
		Tests/Makefile
		Makefile])

AC_CONFIG_LINKS([Doc/BSMlibs.in:Doc/BSMlibs.in])
AC_CONFIG_FILES([Doc/fixinterfaces.pl],[chmod +x Doc/fixinterfaces.pl])

HERWIG_OVERVIEW

AC_CONFIG_COMMANDS([summary],[cat config.herwig])

AC_OUTPUT
