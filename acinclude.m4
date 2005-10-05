AC_DEFUN([AC_CHECK_CLHEP],
[
AC_MSG_CHECKING([for CLHEPPATH])
if test -z "$CLHEPPATH"; then
  AC_MSG_RESULT([none])
  AC_MSG_ERROR([CLHEPPATH not set])
fi
AC_MSG_RESULT([$CLHEPPATH])

AC_MSG_CHECKING([for CLHEPLIB])
if test -z "$CLHEPLIB"; then
  for filename in $CLHEPPATH/lib/libCLHEP-?.?.?.?.{so,dylib} $CLHEPPATH/lib/libCLHEP.{so,dylib}; do
    if test -e $filename; then
	CLHEPLIB=`basename $filename | sed -e 's/^lib/-l/' -e 's/\.\(so\|dylib\)$//'`
    fi
  done
  if test -z "$CLHEPLIB"; then
    AC_MSG_RESULT([none])
    AC_MSG_ERROR([Cannot find libCLHEP at $CLHEPPATH/lib.])
  fi
fi
CLHEPLDFLAGS=-L$CLHEPPATH/lib
AC_MSG_RESULT([$CLHEPLIB])

AC_MSG_CHECKING([for CLHEPINCLUDE]) 
if test -z "$CLHEPINCLUDE"; then
  CLHEPINCLUDE=-I$CLHEPPATH/include
fi
AC_MSG_RESULT([$CLHEPINCLUDE])

# Now lets see if the libraries work properly
oldLIBS="$LIBS"
oldLDFLAGS="$LDFLAGS"
oldCPPFLAGS="$CPPFLAGS"
LIBS="$LIBS $CLHEPLIB"
LDFLAGS="$LDFLAGS $CLHEPLDFLAGS"
CPPFLAGS="$CPPFLAGS $CLHEPINCLUDE"

# check CLHEP first
AC_MSG_CHECKING([that CLHEP works])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <CLHEP/Random/Random.h>]],[[using namespace CLHEP; HepRandom r; r.flat();]])],[AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no]) 
AC_MSG_ERROR([CLHEP must be installed to continue.])
])

LIBS="$oldLIBS"
LDFLAGS="$oldLDFLAGS"
CPPFLAGS="$oldCPPFLAGS"

AC_SUBST(CLHEPLIB)
AC_SUBST(CLHEPLDFLAGS)
])


AC_DEFUN([AC_CHECK_THEPEG],
[
AC_MSG_CHECKING([for THEPEGPATH])
if test -z "$THEPEGPATH"; then
  AC_MSG_RESULT([none])
  AC_MSG_ERROR([THEPEGPATH not set])
fi
AC_MSG_RESULT([$THEPEGPATH])

THEPEGLIB="-lThePEG"
THEPEGLDFLAGS="-L$THEPEGPATH/lib/ThePEG"

AC_MSG_CHECKING([for THEPEGINCLUDE])
if test -z "$THEPEGINCLUDE"; then
  THEPEGINCLUDE=-I$THEPEGPATH/include
fi
AC_MSG_RESULT([$THEPEGINCLUDE])

dnl ###############################
dnl ###############################

dnl Now lets see if the libraries work properly
oldLIBS="$LIBS"
oldLDFLAGS="$LDFLAGS"
oldCPPFLAGS="$CPPFLAGS"
LIBS="$oldLIBS $THEPEGLIB $CLHEPLIB -ldl"
LDFLAGS="$oldLDFLAGS $THEPEGLDFLAGS $CLHEPLDFLAGS"
CPPFLAGS="$oldCPPFLAGS $CLHEPINCLUDE $THEPEGINCLUDE"
AC_MSG_CHECKING([that ThePEG works])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <ThePEG/Repository/Repository.h>
]], [[breakThePEG();
]])],[AC_MSG_RESULT([yes])
],[AC_MSG_RESULT([no])
AC_MSG_ERROR([ThePEG must be installed to continue.])])

LIBS="$oldLIBS"
LDFLAGS="$oldLDFLAGS"
CPPFLAGS="$oldCPPFLAGS"

AC_SUBST(THEPEGLIB)
AC_SUBST(THEPEGLDFLAGS)
])


AC_DEFUN([AC_CHECK_KTJET],[AC_MSG_CHECKING([KTJETPATH is])
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

AC_SUBST(KTJETPATH)
AC_SUBST(KTJETLIBS)

])

AC_DEFUN([AC_CHECK_AMEGIC],
[
AC_MSG_CHECKING([AMEGICPATH is])
if test -z "$AMEGICPATH"; then	
  AC_MSG_RESULT([*** No AMEGIC path set... won't build AMEGIC interface ***])
else
  AMEGICLIBS="-lAmegic -lAmplitude -lAmplTools -lISR -lModel -lPhasespace -lProcess -lString -lVector -lZfunctions"
  AC_MSG_RESULT("$AMEGICPATH")
fi

AM_CONDITIONAL(WANT_LIBAMEGIC,[test ! -z "$AMEGICPATH"])

AC_SUBST(AMEGICPATH)
AC_SUBST(AMEGICLIBS)

])

AC_DEFUN([AC_CHECK_EVTGEN],
[
AC_MSG_CHECKING([EVTGENPATH is])
if test -z "$EVTGENPATH"; then
  AC_MSG_RESULT([*** No EvtGen path set... won't build EvtGen interface ***])
else
  EVTGENLIBS="-lHwEvtGen -lEvtGenBase -lEvtGenModels"
  AC_MSG_RESULT("EVTGENPATH")
fi

AM_CONDITIONAL(WANT_LIBEVTGEN,[test ! -z "$EVTGENPATH"])

AC_SUBST(EVTGENPATH)
AC_SUBST(EVTGENLIBS)
])
