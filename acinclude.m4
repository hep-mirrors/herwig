dnl ##### CLHEP #####
AC_DEFUN([HERWIG_CHECK_CLHEP],
[
AC_MSG_CHECKING([for CLHEP])
AC_ARG_WITH(clhep,
        AC_HELP_STRING([--with-clhep=DIR],[location of CLHEP installation]),
        [],
	[with_clhep=no])
AC_MSG_RESULT([$with_clhep])
CLHEPINCLUDE=""
CLHEPPATH=""
CLHEPLIB=""
if test "x$with_clhep" != "xno"; then
	CLHEPPATH=$with_clhep

	AC_MSG_CHECKING([for CLHEPLIB])
	if test -z "$CLHEPLIB"; then
	  for filename in $CLHEPPATH/lib/libCLHEP-?.?.?.?.{so,dylib} $CLHEPPATH/lib/libCLHEP.{so,dylib}
	  do
		if test -e $filename; then
		   CLHEPLIB=`basename $filename | sed -e 's/^lib/-l/' -e 's/\.\(so\|dylib\)$//'`
		fi
	  done
	  if test -z "$CLHEPLIB"; then
	      AC_MSG_RESULT([none])
	      AC_MSG_ERROR([Cannot find libCLHEP at $CLHEPPATH/lib.])
	  fi
	fi
	CLHEPLDFLAGS="-L$CLHEPPATH/lib -R$CLHEPPATH/lib"
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
	AC_MSG_ERROR([CLHEP not OK. See 'config.log' for details.])
	])
	
	LIBS="$oldLIBS"
	LDFLAGS="$oldLDFLAGS"
	CPPFLAGS="$oldCPPFLAGS"

	AC_SUBST(CLHEPLIB)
	AC_SUBST(CLHEPLDFLAGS)
	AC_SUBST(CLHEPINCLUDE)
fi
])


dnl ##### HEPMC #####
AC_DEFUN([HERWIG_CHECK_HEPMC],
[
AC_REQUIRE([HERWIG_CHECK_CLHEP])
AC_MSG_CHECKING([for HepMC location])
HEPMCINCLUDE=""
HEPMCLIBS="-lHepMC"
hepmclinkname=HepMC
AC_ARG_WITH(hepmc,
        AC_HELP_STRING([--with-hepmc=DIR],[Location of HepMC installation @<:@default=CLHEP@:>@]),
        [],
	[if test -z $CLHEPINCLUDE; then with_hepmc=no; else with_hepmc="CLHEP"; fi])

if test "x$with_hepmc" = "xno"; then
	AC_MSG_RESULT([HepMC support disabled.])
else
	if test "$with_hepmc" = "CLHEP"; then
		if test -f "${CLHEPINCLUDE#-I}/CLHEP/HepMC/GenEvent.h"; then
			AC_MSG_RESULT([part of CLHEP])
			HEPMCINCLUDE=$CLHEPINCLUDE/CLHEP
			HEPMCLIBS="$CLHEPLDFLAGS $CLHEPLIB"
			hepmclinkname=CLHEP
		else
			AC_MSG_RESULT([not found in CLHEP, use '--with-hepmc=' explicitly.])
			with_hepmc=no
		fi
	else
		AC_MSG_RESULT([$with_hepmc])
		HEPMCINCLUDE=-I$with_hepmc/include
		HEPMCLIBS="-L$with_hepmc/lib -R$with_hepmc/lib -lHepMC"
	fi
fi

if test "x$with_hepmc" != "xno"; then
	# Now lets see if the libraries work properly
	oldLIBS="$LIBS"
	oldLDFLAGS="$LDFLAGS"
	oldCPPFLAGS="$CPPFLAGS"
	LIBS="$LIBS $CLHEPLIB $HEPMCLIBS"
	LDFLAGS="$LDFLAGS $CLHEPLDFLAGS"
	CPPFLAGS="$CPPFLAGS $CLHEPINCLUDE $HEPMCINCLUDE"

	# check HepMC
	AC_MSG_CHECKING([that HepMC works])
	AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <HepMC/GenEvent.h>
]],[[HepMC::GenEvent();]])],[AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no]) 
	AC_MSG_ERROR([Use '--with-hepmc=' to set a path or use '--without-hepmc'.])
	])
	
	LIBS="$oldLIBS"
	LDFLAGS="$oldLDFLAGS"
	CPPFLAGS="$oldCPPFLAGS"
fi


AM_CONDITIONAL(HAVE_HEPMC,[test "x$with_hepmc" != "xno"])
AC_SUBST(HEPMCINCLUDE)
AC_SUBST(HEPMCLIBS)
AC_CONFIG_LINKS([Config/HepMCHelper.h:Config/HepMCHelper_$hepmclinkname.h])
])

dnl ##### THEPEG #####
AC_DEFUN([HERWIG_CHECK_THEPEG],
[
AC_MSG_CHECKING([for libThePEG in])
AC_ARG_WITH(thepeg,
        AC_HELP_STRING([--with-thepeg=DIR],[location of ThePEG installation]),
        [],
	[with_thepeg="${prefix}"])
AC_MSG_RESULT([$with_thepeg])

if test "x$with_thepeg" = "xno"; then
	AC_MSG_ERROR([Cannot build Herwig++ without ThePEG. Please set --with-thepeg.])
fi

THEPEGLDFLAGS="-L${with_thepeg}/lib/ThePEG"
THEPEGPATH="${with_thepeg}"

oldldflags="$LDFLAGS"
oldlibs="$LIBS"

LDFLAGS=$THEPEGLDFLAGS
AC_CHECK_LIB([ThePEG],[debugThePEG],[],
	[AC_MSG_ERROR([No ThePEG libraries in $THEPEGLDFLAGS. Please set --with-thepeg.])])

AC_SUBST([THEPEGLIB],[-lThePEG])
AC_SUBST(THEPEGLDFLAGS)
AC_SUBST(THEPEGPATH)

LIBS="$oldlibs"
LDFLAGS="$oldldflags"

AC_MSG_CHECKING([for ThePEG headers in])
AC_ARG_WITH([thepeg-headers],
        AC_HELP_STRING([--with-thepeg-headers=DIR],[location of ThePEG include directory]),
        [],
	[with_thepeg_headers="${with_thepeg}/include"])
AC_MSG_RESULT([$with_thepeg_headers])

if test "x$with_thepeg_headers" = "xno"; then
	AC_MSG_ERROR([Cannot build Herwig++ without ThePEG headers. Please set --with-thepeg-headers.])
fi

THEPEGINCLUDE="-I$with_thepeg_headers"

oldcppflags="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $THEPEGINCLUDE"
AC_CHECK_HEADER([ThePEG/Config/ThePEG.h],[],
	[AC_MSG_ERROR([No ThePEG headers in $with_thepeg_headers. Please set --with-thepeg-headers.])])
CPPFLAGS="$oldcppflags"

AC_SUBST(THEPEGINCLUDE)
])

dnl ##### KTJET #####
AC_DEFUN([HERWIG_CHECK_KTJET],[
AC_REQUIRE([HERWIG_CHECK_CLHEP])

KTJETPATH=""
KTJETLIBS=""
KTJETINCLUDE=""
LOAD_KTJET=""

AC_MSG_CHECKING([for KtJet])

AC_ARG_WITH(ktjet,
        AC_HELP_STRING([--with-ktjet=DIR],[location of KtJet installation]),
        [],
	[with_ktjet=no])


if test "x$with_ktjet" = "xno"; then
	AC_MSG_RESULT([not required])	
else
	if test -z "$CLHEPINCLUDE"; then
		AC_MSG_RESULT([need CLHEP])
		AC_MSG_ERROR([KtJet needs CLHEP headers. Please set --with-clhep.])
	fi

	AC_MSG_RESULT([required])
	KTJETPATH="$with_ktjet"

	AC_MSG_CHECKING([KtJet library name is])
	if test -f "$KTJETPATH/lib/libKtJet.a"; then
		ktjetname=KtJet
		AC_MSG_RESULT([KtJet])
	elif test -f "$KTJETPATH/lib/libKtEvent.a"; then
		ktjetname=KtEvent
		AC_MSG_RESULT([KtEvent])
	else
		AC_MSG_RESULT([?])
		AC_MSG_ERROR([No KtJet library found in $KTJETPATH/lib.])
	fi

	ktjetrpath=""
	if test -e "$KTJETPATH/lib/lib$ktjetname.so"; then
		ktjetrpath="-R$KTJETPATH/lib"
	fi

	oldlibs=$LIBS
	oldcxxflags=$CXXFLAGS
	LIBS=""
	CXXFLAGS="-L$KTJETPATH/lib $ktjetrpath -l$ktjetname $CLHEPLDFLAGS"
	AC_CHECK_LIB([$ktjetname],[abort],
		     [],
		     [
			AC_MSG_ERROR([lib$ktjename not working. See 'config.log' for details.])
		     ],
		     [$CLHEPLIB])   
	KTJETLIBS="$CXXFLAGS $LIBS"
	LIBS=$oldlibs
	CXXFLAGS=$oldcxxflags

	AC_SUBST(KTJETLIBS)
	
	AC_MSG_CHECKING([KtJet headers])
	if test -f "$KTJETPATH/include/KtJet/KtJet.h"; then
		KTJETINCLUDE="-I$KTJETPATH/include"
	elif test -f "$KTJETPATH/KtJet/KtJet.h"; then
		KTJETINCLUDE="-I$KTJETPATH"
	else
		AC_MSG_RESULT([not found.])
		AC_MSG_ERROR([No KtJet headers. Please set KTJETINCLUDE explicitly.])
	fi
	AC_MSG_RESULT([$KTJETINCLUDE])
	KTJETINCLUDE="$KTJETINCLUDE $CLHEPINCLUDE"

	oldcppflags="$CPPFLAGS"
	CPPFLAGS="$CPPFLAGS $KTJETINCLUDE"
	AC_CHECK_HEADER([KtJet/KtJet.h],[],
	[AC_MSG_ERROR([Problem with KtJet headers in $KTJETINCLUDE.])])
	CPPFLAGS="$oldcppflags"

	AC_SUBST(KTJETINCLUDE)
	
	LOAD_KTJET="read KtJetAnalysis.in"
	AC_SUBST(LOAD_KTJET)
fi

AM_CONDITIONAL(WANT_LIBKTJET,[test ! -z "$KTJETPATH"])
])

dnl ##### LOOPTOOLS #####
AC_DEFUN([HERWIG_LOOPTOOLS],
[
AC_MSG_CHECKING([whether to build Looptools dependent parts])
AC_ARG_ENABLE(looptools,
        AC_HELP_STRING([--disable-looptools],[turn off Looptools-dependent parts]),
        [],
        [enable_looptools=yes]
        )
AC_MSG_RESULT([$enable_looptools])
AM_CONDITIONAL(WANT_LOOPTOOLS,[test "x$enable_looptools" = "xyes"])
])

dnl ##### PDF PATH #####
AC_DEFUN([HERWIG_PDF_PATH],
[
AC_MSG_CHECKING([which Herwig++ PDF data to use])
AC_ARG_WITH(pdf,
        AC_HELP_STRING([--with-pdf=DIR],[installation path of Herwig++PDF data tarball]),
        [],
        [with_pdf=${prefix}]
        )
HERWIG_PDF_DEFAULT=${with_pdf}/share/Herwig++PDF/mrst/1998/lo05a.dat

if test -f "${HERWIG_PDF_DEFAULT}"; then
	AC_MSG_RESULT([$with_pdf])
	localPDFneeded=false
else
	AC_MSG_RESULT([Using built-in PDF data set. For other data sets, set --with-pdf.])
	HERWIG_PDF_DEFAULT=../PDF/mrst/1998/lo05a.dat
	localPDFneeded=true
fi
AM_CONDITIONAL(WANT_LOCAL_PDF,[test "x$localPDFneeded" = "xtrue"])
AC_SUBST(HERWIG_PDF_DEFAULT)
])

dnl ##### EVTGEN #####
AC_DEFUN([HERWIG_CHECK_EVTGEN],
[
AC_MSG_CHECKING([for EVTGEN])
AC_ARG_WITH(evtgen,
        AC_HELP_STRING([--with-evtgen=DIR],[installation path of EvtGen]),
        [],
        [with_evtgen=no])

EVTGENPATH=
LOAD_EVTGEN=""

if test "x$with_evtgen" = "xno"; then
	AC_MSG_RESULT([not required])
else
	AC_MSG_RESULT([$with_evtgen])

	oldLIBS="$LIBS"
	tmpcxxflags=$CXXFLAGS

	# Now lets see if the libraries work properly
	

	CXXFLAGS="$CXXFLAGS -L${with_evtgen}/lib"

	AC_CHECK_LIB([evtgenlhc],[abort],
		[],
		[
			AC_MSG_ERROR([libevtgenlhc could not be found at ${with_evtgen}/lib])
		])
	LIBS="$oldLIBS"
	CXXFLAGS=$tmpcxxflags
	EVTGENPATH=$with_evtgen
	LOAD_EVTGEN="library HwEvtGen.so"
fi
AM_CONDITIONAL(WANT_EVTGEN,[test "x$with_evtgen" != "xno"])
AC_SUBST(EVTGENPATH)
AC_SUBST(LOAD_EVTGEN)
])]


dnl ###### GSL ######
AC_DEFUN([HERWIG_CHECK_GSL],
[
AC_MSG_CHECKING([for gsl location])
GSLINCLUDE=""
GSLLIBS=""

AC_ARG_WITH(gsl,
        AC_HELP_STRING([--with-gsl=DIR],[location of gsl installation @<:@default=system libs@:>@]),
        [],
	[with_gsl=no])

if test "x$with_gsl" = "xno"; then
	AC_MSG_RESULT([in system libraries])
	oldlibs="$LIBS"
	AC_CHECK_LIB(m,main)
	AC_CHECK_LIB(gslcblas,main)
	AC_CHECK_LIB(gsl,main,[],
			[
			AC_MSG_ERROR([Cannot find libgsl. Please install the GNU scientific library.])
			]
		     )
	GSLLIBS="$LIBS"
	LIBS=$oldlibs
else
	if test -e "$with_gsl/lib/libgsl.a" -a -d "$with_gsl/include/gsl"; then
		AC_MSG_RESULT([found in $with_gsl])
		GSLLIBS="-L$with_gsl/lib -R$with_gsl/lib -lgslcblas -lgsl"
		GSLINCLUDE="-I$with_gsl/include"
	else
		AC_MSG_RESULT([not found])
		AC_MSG_ERROR([Can't find $with_gsl/lib/libgsl.a or the headers in $with_gsl/include])
	fi
fi

dnl AM_CONDITIONAL(HAVE_GSL,[test "x$with_hepmc" != "xno"])
AC_SUBST(GSLINCLUDE)
AC_SUBST(GSLLIBS)
])

AC_DEFUN([HERWIG_VERSIONSTRING],
[
if test -d $srcdir/.svn; then
	AC_CHECK_PROG(have_svnversion,[svnversion],[yes],[no])
fi
AM_CONDITIONAL(USE_SVNVERSION,[test "x$have_svnversion" = "xyes"])
])
