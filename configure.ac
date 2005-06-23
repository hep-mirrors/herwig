dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Herwig++],[0.9.5],[http://www.hep.phy.cam.ac.uk/theory/Herwig++/],[Herwig++])
AC_CONFIG_AUX_DIR([Config])

AC_LANG(C++)

AC_CONFIG_SRCDIR([UnderlyingEvent/UA5.h])
AC_CONFIG_HEADERS([Config/config.h])

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

AC_CHECK_CLHEP
AC_CHECK_THEPEG

AC_CHECK_KTJET
AC_CHECK_AMEGIC

dnl ###########################
dnl ###########################

AM_CPPFLAGS="-I\$(top_builddir)/include $THEPEGINCLUDE $CLHEPINCLUDE"
AM_CXXFLAGS="-ansi -pedantic -Wall"
# should really be specified one-by-one for each lib!
## AM_LDFLAGS="-version-info 1:0:0 -Wl,--enable-new-dtags"
AM_LDFLAGS="-version-info 1:0:0"
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
		PDT/Makefile
		PDF/Makefile
                MatrixElement/Makefile
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
		Looptools/Makefile
		Makefile])
AC_OUTPUT
