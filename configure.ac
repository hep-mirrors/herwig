dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Herwig++],[2.0-alpha3],[http://www.hep.phy.cam.ac.uk/theory/Herwig++/],[Herwig++])
AC_CONFIG_AUX_DIR([Config])

AC_CONFIG_SRCDIR([UnderlyingEvent/UA5.h])
AC_CONFIG_HEADERS([Config/config.h])

AC_CANONICAL_HOST

case "${host}" in
  *-darwin[[0156]].*) 
    AC_MSG_ERROR([Herwig++ requires OS X 10.3 or later])
    ;;
  *-darwin*)
    if test "x$MACOSX_DEPLOYMENT_TARGET" != "x10.3"; then
      AC_MSG_ERROR([Please export the MACOSX_DEPLOYMENT_TARGET variable, set to 10.3])
    fi
esac

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

AC_CHECK_CLHEP
AC_CHECK_THEPEG

AC_CHECK_KTJET
AC_CHECK_AMEGIC

## temporary workaround as long as looptools isn't actually necessary
case "{$host}" in
     *-darwin*) DO_LOOPTOOLS= ;;
     *)		DO_LOOPTOOLS=Looptools ;;
esac
AC_SUBST(DO_LOOPTOOLS)

dnl ###########################
dnl ###########################

AM_CPPFLAGS="-I\$(top_builddir)/include $THEPEGINCLUDE $CLHEPINCLUDE"
AM_CXXFLAGS="-ansi -pedantic -Wall"

# version info should really be specified one-by-one for each lib!
AM_LDFLAGS="-version-info 1:0:0"

# don't know how else to do this. see bug #14
if test -n "$GCC"; then
   case "${host}" in
   *-linux*)	AM_LDFLAGS="$AM_LDFLAGS -Wl,--disable-new-dtags" ;;
   *)		;;
   esac
fi

AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_LDFLAGS)

AC_CONFIG_FILES([UnderlyingEvent/Makefile
		Models/Makefile
		Models/StandardModel/Makefile
		Models/RSModel/Makefile
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
		PDT/Makefile
		PDF/Makefile
                MatrixElement/Makefile
		Shower/SplittingFunctions/Makefile
                Shower/Makefile
                Utilities/Makefile
                Hadronization/Makefile
		Helicity/Makefile
                Helicity/Correlations/Makefile
		Helicity/Vertex/Makefile
                Helicity/Vertex/RSModel/Makefile
                Helicity/Vertex/Scalar/Makefile
                Helicity/Vertex/StandardModel/Makefile
                Helicity/Vertex/Tensor/Makefile
                Helicity/Vertex/Vector/Makefile
                Helicity/WaveFunction/Makefile
		lib/Makefile
		include/Makefile
		src/Makefile
		Doc/Makefile
		Looptools/Makefile
		Analysis/Makefile
		Analysis/Makefile-shared
		Makefile])

AC_CONFIG_FILES([Doc/fixinterfaces.pl],[chmod +x Doc/fixinterfaces.pl])

AC_OUTPUT
