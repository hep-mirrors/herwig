dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([Herwig++],[0.9.5],[http://www.hep.phy.cam.ac.uk/theory/Herwig++/],[Herwig++])
AC_CONFIG_AUX_DIR([Config])

AC_LANG(C++)

AC_CONFIG_SRCDIR([UnderlyingEvent/UA5.h])
AC_CONFIG_HEADERS([Config/config.h])

AM_INIT_AUTOMAKE
dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET #???
AC_PROG_LN_S
AC_DISABLE_STATIC

## needed for looptools
AC_PROG_F77
## Looptools seems to need this if you want cache.c
AC_PROG_CC

AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

dnl ####################################
dnl ####################################

AC_MSG_CHECKING([CLHEPPATH is set])
if test -z "$CLHEPPATH"; then
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([Set CLHEPPATH to point to the path of CLHEP before running configure.])
else
  AC_MSG_RESULT("$CLHEPPATH")
fi

AC_MSG_CHECKING([CLHEPLIB is])
CLHEPLIB=`$CLHEPPATH/bin/clhep-config --libs | awk '{print $2}'`
CLHEPLDFLAGS=`$CLHEPPATH/bin/clhep-config --libs | awk '{print $1}'`
AC_MSG_RESULT("$CLHEPLIB")

AC_MSG_CHECKING([CLHEPINCLUDE is])
CLHEPINCLUDE=`$CLHEPPATH/bin/clhep-config --include`
AC_MSG_RESULT("$CLHEPINCLUDE")

dnl ###############################

AC_MSG_CHECKING([THEPEGPATH is set])
if test -z "$THEPEGPATH"; then
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([Set THEPEGPATH to point to the path of ThePEG before running configure.])
else
  AC_MSG_RESULT("$THEPEGPATH")
fi

AC_MSG_CHECKING([THEPEGLIB is])
THEPEGLIB="-lThePEG"
THEPEGLDFLAGS="-L$THEPEGPATH/lib/ThePEG"
AC_MSG_RESULT("$THEPEGLIB")

AC_MSG_CHECKING([THEPEGINCLUDE is])
THEPEGINCLUDE=-I$THEPEGPATH/include
AC_MSG_RESULT("$THEPEGINCLUDE")

dnl ###############################
dnl ###############################

dnl Now lets see if the libraries work properly
oldLIBS="$LIBS"
oldLDFLAGS="$LDFLAGS"
oldCPPFLAGS="$CPPFLAGS"
LIBS="$LIBS $CLHEPLIB"
LDFLAGS="$LDFLAGS $CLHEPLDFLAGS"
CPPFLAGS="$CPPFLAGS $CLHEPINCLUDE"

dnl check CLHEP first
AC_MSG_CHECKING([that CLHEP works])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <CLHEP/Random/Random.h>]], [[using namespace CLHEP; HepRandom r; r.flat();]])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no) 
AC_MSG_ERROR(CLHEP must be installed to continue.)
])

dnl check ThePEG now
LIBS="$oldLIBS $THEPEGLIB $CLHEPLIB -ldl"
LDFLAGS="$oldLDFLAGS $THEPEGLDFLAGS $CLHEPLDFLAGS"
CPPFLAGS="$oldCPPFLAGS $CLHEPINCLUDE $THEPEGINCLUDE"
AC_MSG_CHECKING([that ThePEG works])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <ThePEG/Repository/Repository.h>
]], [[breakThePEG();
]])],[AC_MSG_RESULT(yes)
],[AC_MSG_RESULT(no)
AC_MSG_ERROR(ThePEG must be installed to continue.)])

LIBS="$oldLIBS"
LDFLAGS="$oldLDFLAGS"
CPPFLAGS="$oldCPPFLAGS"

dnl ###########################
dnl ###########################

AC_MSG_CHECKING([KTJETPATH is])
if test -z "$KTJETPATH"; then
  AC_MSG_RESULT([*** No KtJet path set... won't build KtJet interface ***])
else
  AC_MSG_RESULT("$KTJETPATH")

  AC_MSG_CHECKING([KTJETLIBS is])
  if test -z "$KTJETLIBS"; then
    KTJETLIBS="-L$KTJETPATH/lib -lKtEvent"
  fi
  AC_MSG_RESULT("$KTJETLIBS")
fi

AM_CONDITIONAL(WANT_LIBKTJET,[test ! -z "$KTJETPATH"])

AC_MSG_CHECKING([AMEGICPATH is])
if test -z "$AMEGICPATH"; then	
  AC_MSG_RESULT([*** No AMEGIC path set... won't build AMEGIC interface ***])
else
  AMEGICLIBS="-lAmegic -lAmplitude -lAmplTools -lISR -lModel -lPhasespace -lProcess -lString -lVector -lZfunctions"
  AC_MSG_RESULT("$AMEGICPATH")
fi

AM_CONDITIONAL(WANT_LIBAMEGIC,[test ! -z "$AMEGICPATH"])


dnl ##########################
dnl ##########################

AM_CPPFLAGS="-I\$(top_builddir)/include $THEPEGINCLUDE $CLHEPINCLUDE"
AM_CXXFLAGS="-ansi -pedantic -Wall"
# should really be specified one-by-one for each lib!
## AM_LDFLAGS="-version-info 1:0:0 -Wl,--enable-new-dtags"
AM_LDFLAGS="-version-info 1:0:0"
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_LDFLAGS)

dnl AC_SUBST(CLHEPPATH)
AC_SUBST(CLHEPLIB)
AC_SUBST(CLHEPLDFLAGS)
dnl AC_SUBST(CLHEPINCLUDE)

dnl AC_SUBST(THEPEGPATH)
AC_SUBST(THEPEGLIB)
AC_SUBST(THEPEGLDFLAGS)
dnl AC_SUBST(THEPEGINCLUDE)

AC_SUBST(KTJETPATH)
AC_SUBST(KTJETLIBS)

AC_SUBST(AMEGICPATH)
AC_SUBST(AMEGICLIBS)

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
