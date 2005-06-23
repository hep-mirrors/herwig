AC_DEFUN([AC_CHECK_CLHEP],
[
AC_MSG_CHECKING([CLHEPPATH is set])
if test -z "$CLHEPPATH"; then
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([Set CLHEPPATH to point to the path of CLHEP before running configure.])
else
  AC_MSG_RESULT("$CLHEPPATH")
fi


AC_MSG_CHECKING([CLHEPLIB is])
CLHEPLIB=`$CLHEPPATH/bin/clhep-config --libs | cut -d' ' -f2`
CLHEPLDFLAGS=`$CLHEPPATH/bin/clhep-config --libs | cut -d' ' -f1`
AC_MSG_RESULT("$CLHEPLIB")

AC_MSG_CHECKING([CLHEPINCLUDE is]) 
CLHEPINCLUDE=`$CLHEPPATH/bin/clhep-config --include` 
AC_MSG_RESULT("$CLHEPINCLUDE")

# Now lets see if the libraries work properly
oldLIBS="$LIBS"
oldLDFLAGS="$LDFLAGS"
oldCPPFLAGS="$CPPFLAGS"
LIBS="$LIBS $CLHEPLIB"
LDFLAGS="$LDFLAGS $CLHEPLDFLAGS"
CPPFLAGS="$CPPFLAGS $CLHEPINCLUDE"

# check CLHEP first
AC_MSG_CHECKING([that CLHEP works])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <CLHEP/Random/Random.h>]], [[using namespace CLHEP; HepRandom r; r.flat();]])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no) 
AC_MSG_ERROR(CLHEP must be installed to continue.)
])

LIBS="$oldLIBS"
LDFLAGS="$oldLDFLAGS"
CPPFLAGS="$oldCPPFLAGS"

AC_SUBST(CLHEPLIB)
AC_SUBST(CLHEPLDFLAGS)
])


AC_DEFUN([AC_CHECK_THEPEG],
[
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
