dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([Herwig],[7.2.1],[herwig@projects.hepforge.org],[Herwig])
AC_CONFIG_SRCDIR([Utilities/HerwigStrategy.cc])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([Config/config.h])
dnl AC_PRESERVE_HELP_ORDER
AC_CANONICAL_HOST

dnl === disable debug symbols by default =====
if test "x$CXXFLAGS" = "x"; then
   CXXFLAGS="-O2 -DBOOST_UBLAS_NDEBUG"
fi
if test "x$CFLAGS" = "x"; then
   CFLAGS=-O2
fi

AC_LANG([C++])

AM_INIT_AUTOMAKE([1.11 subdir-objects gnu dist-bzip2 no-dist-gzip -Wall -Wno-portability])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

dnl Checks for C++ compiler. Handle C++11 flags.
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX([11],[noext],[mandatory])

dnl check for POSIX
AC_CHECK_HEADER([unistd.h],[],
      [AC_MSG_ERROR([Herwig needs "unistd.h". Non-POSIX systems are not supported.])])
AC_CHECK_HEADER([sys/stat.h],[],
      [AC_MSG_ERROR([Herwig needs "sys/stat.h". Non-POSIX systems are not supported.])])

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
AC_FC_WRAPPERS


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
BOOST_FIND_HEADER([boost/numeric/ublas/io.hpp])
dnl Boost 1.64 is missing a required header to make these work
dnl we just assume they're there if io.hpp has been found OK above
dnl BOOST_FIND_HEADER([boost/numeric/ublas/matrix.hpp])
dnl BOOST_FIND_HEADER([boost/numeric/ublas/matrix_proxy.hpp])
dnl BOOST_FIND_HEADER([boost/numeric/ublas/matrix_sparse.hpp])
dnl BOOST_FIND_HEADER([boost/numeric/ublas/symmetric.hpp])
dnl BOOST_FIND_HEADER([boost/numeric/ublas/vector.hpp])
BOOST_FIND_HEADER([boost/operators.hpp])
BOOST_TEST()

HERWIG_CHECK_VBFNLO

HERWIG_CHECK_NJET

HERWIG_CHECK_GOSAM

HERWIG_CHECK_GOSAM_CONTRIB

HERWIG_CHECK_OPENLOOPS

HERWIG_CHECK_MADGRAPH

HERWIG_CHECK_EVTGEN
HERWIG_CHECK_PYTHIA

HERWIG_COMPILERFLAGS

HERWIG_LOOPTOOLS

FASTJET_CHECK_FASTJET

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
		Decay/EvtGen/Makefile
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
		MatrixElement/Matchbox/External/MadGraph/mg2herwig
		MatrixElement/FxFx/Makefile
		Sampling/Makefile
		Sampling/CellGrids/Makefile
		Shower/Makefile
		Shower/QTilde/Makefile
		Shower/QTilde/Matching/Makefile
		Shower/Dipole/Makefile
		Shower/Dipole/Base/Makefile
		Shower/Dipole/Kernels/Makefile
		Shower/Dipole/Kinematics/Makefile
		Shower/Dipole/Utility/Makefile
		Shower/Dipole/AlphaS/Makefile
		Shower/Dipole/SpinCorrelations/Makefile
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
		src/defaults/Decays.in
		src/defaults/decayers.in
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
		src/Merging/Makefile
		Shower/Dipole/Merging/Makefile
		src/defaults/MatchboxMergingDefaults.in
		Contrib/Makefile
		Contrib/make_makefiles.sh
		Tests/Makefile
		Makefile])

AC_CONFIG_LINKS([Doc/BSMlibs.in:Doc/BSMlibs.in])
AC_CONFIG_FILES([Doc/fixinterfaces.pl],[chmod +x Doc/fixinterfaces.pl])
AC_CONFIG_HEADERS([PDF/SaSPhotonPDF.cc])
HERWIG_OVERVIEW

AC_CONFIG_COMMANDS([summary],[cat config.herwig])

AC_OUTPUT
