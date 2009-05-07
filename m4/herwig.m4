# check for gcc bug http://gcc.gnu.org/bugzilla/show_bug.cgi?id=34130
AC_DEFUN([HERWIG_CHECK_ABS_BUG],
[
AC_REQUIRE([HERWIG_COMPILERFLAGS])
if test "$GCC" = "yes"; then
AC_MSG_CHECKING([for gcc abs bug])
AC_RUN_IFELSE(
	AC_LANG_PROGRAM(
		[ int foo (int i) { return -2 * __builtin_abs(i - 2); }	],
		[ if ( foo(1) != -2 || foo(3) != -2 ) return 1; ]
	),
	[ AC_MSG_RESULT([not found. Compiler is ok.]) ],
	[
	AC_MSG_RESULT([found. Builtin abs() is buggy.])
	AC_MSG_CHECKING([if -fno-builtin-abs works])
	oldcxxflags=$CXXFLAGS
	CXXFLAGS="$CXXFLAGS -fno-builtin-abs"
	AC_RUN_IFELSE(
		AC_LANG_PROGRAM(
			[
			#include <cstdlib>
			int foo (int i) { return -2 * std::abs(i - 2); }
			],
			[
			if (foo(1) != -2 || foo(3) != -2) return 1; 
			]
		),
		[
		AC_MSG_RESULT([yes. Setting -fno-builtin-abs.])
		AM_CXXFLAGS="$AM_CXXFLAGS -fno-builtin-abs"
		AM_CFLAGS="$AM_CFLAGS -fno-builtin-abs"
		],
		[
		AC_MSG_RESULT([no. Setting -fno-builtin.])
		AC_MSG_WARN([
*****************************************************************************
For this version of gcc, -fno-builtin-abs alone did not work to avoid the 
gcc abs() bug. Instead, all gcc builtin functions are now disabled.
Update gcc if possible.
*****************************************************************************])
		AM_CXXFLAGS="$AM_CXXFLAGS -fno-builtin"
		AM_CFLAGS="$AM_CFLAGS -fno-builtin"
		]
	)
	CXXFLAGS=$oldcxxflags
	]
)
fi
])

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
		   CLHEPLIB=`basename $filename | sed -e 's/^lib/-l/' -e 's/\.so//' -e 's/\.dylib//'`
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
AC_MSG_CHECKING([for HepMC location])
HEPMCINCLUDE=""
CREATE_HEPMC="#create"
HEPMCLIBS="-lHepMC"

AC_ARG_WITH(hepmc,
        AC_HELP_STRING([--with-hepmc=DIR],[Location of HepMC installation @<:@default=system libs@:>@]),
        [],
	[with_hepmc=system])

if test "x$with_hepmc" = "xno"; then
	AC_MSG_RESULT([HepMC support disabled.])
elif test "x$with_hepmc" = "xsystem"; then
        AC_MSG_RESULT([in system libraries])
	oldlibs="$LIBS"
	AC_CHECK_LIB(HepMC,main,
		[],
		[with_hepmc=no
		 AC_MSG_WARN([
HepMC not found in system libraries])
		])
	HEPMCLIBS="$LIBS"
	LIBS=$oldlibs
else
	AC_MSG_RESULT([$with_hepmc])
	HEPMCINCLUDE=-I$with_hepmc/include
	HEPMCLIBS="-L$with_hepmc/lib -R$with_hepmc/lib -lHepMC"
fi

if test "x$with_hepmc" != "xno"; then
	# Now lets see if the libraries work properly
	oldLIBS="$LIBS"
	oldLDFLAGS="$LDFLAGS"
	oldCPPFLAGS="$CPPFLAGS"
	LIBS="$LIBS $HEPMCLIBS"
	LDFLAGS="$LDFLAGS"
	CPPFLAGS="$CPPFLAGS $HEPMCINCLUDE"

	# check HepMC
	AC_MSG_CHECKING([that HepMC works])
	AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <HepMC/GenEvent.h>
]],[[HepMC::GenEvent();]])],[AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no]) 
	AC_MSG_ERROR([Use '--with-hepmc=' to set a path or use '--without-hepmc'.])
	])

	AC_CHECK_HEADERS([HepMC/PdfInfo.h],[],[AC_MSG_ERROR([Need HepMC with PdfInfo support.])],[#include <algorithm>])
	AC_CHECK_HEADERS([HepMC/IO_GenEvent.h])
	AC_CHECK_HEADERS([HepMC/IO_ExtendedAscii.h])

	LIBS="$oldLIBS"
	LDFLAGS="$oldLDFLAGS"
	CPPFLAGS="$oldCPPFLAGS"
	CREATE_HEPMC="create"
fi

AM_CONDITIONAL(HAVE_HEPMC,[test "x$with_hepmc" != "xno"])
AC_SUBST(HEPMCINCLUDE)
AC_SUBST(HEPMCLIBS)
AC_SUBST(CREATE_HEPMC)
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

THEPEGLDFLAGS="-L${with_thepeg}/lib/ThePEG $GSLLIBS"
THEPEGPATH="${with_thepeg}"

oldldflags="$LDFLAGS"
oldlibs="$LIBS"

LDFLAGS="$LDFLAGS $THEPEGLDFLAGS"
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
CREATE_KTJET="#create"

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
	CREATE_KTJET="create"
	AC_SUBST(LOAD_KTJET)
	AC_SUBST(CREATE_KTJET)
fi

AM_CONDITIONAL(WANT_LIBKTJET,[test ! -z "$KTJETPATH"])
])

dnl ##### FastJet #####
AC_DEFUN([HERWIG_CHECK_FASTJET],[

FASTJETPATH=""
FASTJETLIBS=""
FASTJETINCLUDE=""
LOAD_FASTJET=""

AC_MSG_CHECKING([for FastJet])

AC_ARG_WITH(fastjet,
        AC_HELP_STRING([--with-fastjet=DIR],[location of FastJet installation]),
        [],
	[with_fastjet=no])

if test "x$with_fastjet" = "xno"; then
	AC_MSG_RESULT([not required])	
else

	AC_MSG_RESULT([required])
	FASTJETPATH="$with_fastjet"

	# check for libraries

	AC_MSG_CHECKING([for FastJet library])

	if test -z "$FASTJETLIB" ; then
		FASTJETLIB="$FASTJETPATH/lib"
	fi

	if test -e "$FASTJETLIB/libfastjet.a"; then
	  	AC_MSG_RESULT([found static library])

		AC_MSG_CHECKING([for FastJet dynamic library])

		if test ! -e "$FASTJETLIB/libfastjet.so"; then
			AC_MSG_RESULT([not found, please create libfastjet.so from $FASTJETLIB/libfastjet.a])
			AC_MSG_ERROR([No FastJet library found in $FASTJETLIB.])
		fi
		AC_MSG_RESULT([$FASTJETLIB/libfastjet.so])
	else
		AC_MSG_RESULT([?])
		AC_MSG_ERROR([No FastJet library found in $FASTJETLIB.])
	fi

	FASTJETLIBS="-L$FASTJETLIB -lfastjet"

	AC_SUBST(FASTJETLIBS)

	# check for headers

	AC_MSG_CHECKING([for FastJet headers])

	if test -z "$FASTJETINCLUDE" ; then
		FASTJETINCLUDE="$FASTJETPATH/include"
	fi

	if test -f "$FASTJETINCLUDE/fastjet/ClusterSequence.hh"; then
		FASTJETINCLUDE="-I$FASTJETINCLUDE"
	else
		AC_MSG_RESULT([not found.])
		AC_MSG_ERROR([No FastJet headers.])
	fi

	AC_MSG_RESULT([$FASTJETINCLUDE])

	AC_SUBST(FASTJETINCLUDE)

fi

AM_CONDITIONAL(WANT_LIBFASTJET,[test ! -z "$FASTJETPATH"])
])

dnl ##### LOOPTOOLS #####
AC_DEFUN([HERWIG_LOOPTOOLS],
[
AC_REQUIRE([AC_PROG_FC])
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([HERWIG_COMPILERFLAGS])

AC_MSG_CHECKING([if Looptools build works])
enable_looptools=yes

if test "x$GCC" = "xyes"; then
   case "${host}" in
      x86_64-*)
	AM_FCFLAGS="$AM_FCFLAGS -fdefault-integer-8"
      	;;
   esac

   AC_LANG_PUSH([Fortran])
   	oldFCFLAGS="$FCFLAGS"
   	FCFLAGS="$AM_FCFLAGS"
   	AC_COMPILE_IFELSE(
	   	AC_LANG_PROGRAM([],[      print *[,]"Hello"]),
		[],
		[AC_MSG_RESULT([no])
 		 AC_MSG_ERROR([needs gfortran on 64bit machines])]
	)
	FCFLAGS="$oldFCFLAGS"
  AC_LANG_POP([Fortran])
fi
AC_MSG_RESULT([$enable_looptools])

AC_SUBST([F77],[$FC])
AC_SUBST([FFLAGS],[$FCFLAGS])
AC_SUBST([AM_FFLAGS],[$AM_FCFLAGS])
AC_SUBST([FLIBS],[$FCLIBS])
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
HERWIG_PDF_DEFAULT=${with_pdf}/share/Herwig++PDF/mrst/2001/lo2002.dat
HERWIG_PDF_NLO=${with_pdf}/share/Herwig++PDF/mrst/2002/mrst2002nlo.dat

if test -f "${HERWIG_PDF_DEFAULT}"; then
	AC_MSG_RESULT([$with_pdf])
	localPDFneeded=false
else
	AC_MSG_RESULT([Using built-in PDF data set. For other data sets, set --with-pdf.])
	HERWIG_PDF_DEFAULT=PDF/mrst/2001/lo2002.dat
	HERWIG_PDF_NLO=PDF/mrst/2002/mrst2002nlo.dat
	localPDFneeded=true
fi
AM_CONDITIONAL(WANT_LOCAL_PDF,[test "x$localPDFneeded" = "xtrue"])
AC_SUBST(HERWIG_PDF_DEFAULT)
AC_SUBST(HERWIG_PDF_NLO)
])

dnl ###### GSL ######
AC_DEFUN([HERWIG_CHECK_GSL],
[
AC_MSG_CHECKING([for gsl location])
GSLINCLUDE=""
GSLLIBS=""

AC_ARG_WITH(gsl,
        AC_HELP_STRING([--with-gsl=DIR],[location of gsl installation @<:@default=system libs@:>@]),
        [],
	[with_gsl=system])

if test "x$with_gsl" = "xno"; then
AC_MSG_ERROR([libgsl is required. Please install the GNU scientific library and header files.])
fi

if test "x$with_gsl" = "xsystem"; then
	AC_MSG_RESULT([in system libraries])
	oldlibs="$LIBS"
	AC_CHECK_LIB(m,main)
	AC_CHECK_LIB(gslcblas,main)
	AC_CHECK_LIB(gsl,main,[],
			[
			AC_MSG_ERROR([Cannot find libgsl. Please install the GNU scientific library and header files or use --with-gsl=.])
			]
		     )
	GSLLIBS="$LIBS"
	LIBS=$oldlibs
else
	if test "`uname -m`" = "x86_64" -a -e "$with_gsl/lib64/libgsl.a" -a -d "$with_gsl/include/gsl"; then
		AC_MSG_RESULT([found in $with_gsl])
		GSLLIBS="-L$with_gsl/lib64 -R$with_gsl/lib64 -lgslcblas -lgsl"
		GSLINCLUDE="-I$with_gsl/include"
	elif test -e "$with_gsl/lib/libgsl.a" -a -d "$with_gsl/include/gsl"; then
		AC_MSG_RESULT([found in $with_gsl])
		GSLLIBS="-L$with_gsl/lib -R$with_gsl/lib -lgslcblas -lgsl"
		GSLINCLUDE="-I$with_gsl/include"
	else
		AC_MSG_RESULT([not found])
		AC_MSG_ERROR([Can't find $with_gsl/lib/libgsl.a or the headers in $with_gsl/include])
	fi
fi

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

dnl ##### COMPILERFLAGS #####
AC_DEFUN([HERWIG_COMPILERFLAGS],
[
AC_REQUIRE([HERWIG_CHECK_THEPEG])

AM_CPPFLAGS="-I\$(top_builddir)/include $THEPEGINCLUDE \$(GSLINCLUDE)"

AC_MSG_CHECKING([for debugging mode])
AC_ARG_ENABLE(debug,
        AC_HELP_STRING([--enable-debug],[debug mode, use --enable-debug=slow for additional options that slow down the run.]),
        [],
        [enable_debug=no]
        )
AC_MSG_RESULT([$enable_debug])

if test "x$enable_debug" = "xno"; then
	AM_CPPFLAGS="$AM_CPPFLAGS -DNDEBUG"
else
	debugflags="-g"
fi

dnl -Wfloat-equal -fvisibility-inlines-hidden -Wctor-dtor-privacy -Weffc++
if test -n "$GCC"; then
	warnflags="-ansi -pedantic -Wall -W"

	if test "x$enable_debug" = "xslow"; then
		debugflags="$debugflags -fno-inline"
		AM_CPPFLAGS="$AM_CPPFLAGS -D_GLIBCXX_DEBUG"
	fi
fi

AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CFLAGS,  ["$warnflags $debugflags"])
AC_SUBST(AM_CXXFLAGS,["$warnflags $debugflags"])
AC_SUBST(AM_FCFLAGS,  ["$debugflags"])
AC_SUBST(AM_LDFLAGS)
])

AC_DEFUN([HERWIG_ENABLE_MODELS],
[
AC_MSG_CHECKING([for BSM models to include])

AC_ARG_ENABLE(models,
        AC_HELP_STRING([--enable-models=LIST],[Comma-separated list of BSM models to enable. Options are (mssm ued rs) or --disable-models to turn them all off.]),
        [],
        [enable_models=all]
        )
if test "x$enable_models" = "xyes" -o "x$enable_models" = "xall"; then
   all=yes
fi
AC_MSG_RESULT([$enable_models])

if test ! "$all"; then
   oldIFS="$IFS"
   IFS=","
   for i in $enable_models; do
       declare $i=yes
   done
   IFS="$oldIFS"
fi

AC_SUBST([CREATE_BSM_ANALYSIS],["# create"])
if test "$mssm" -a "$ued"; then
   CREATE_BSM_ANALYSIS="create"
fi

AM_CONDITIONAL(WANT_MSSM,[test "$mssm" -o "$all"])
AM_CONDITIONAL(WANT_UED,[test "$ued" -o "$all"])
AM_CONDITIONAL(WANT_RS,[test "$rs" -o "$all"])
])

AC_DEFUN([HERWIG_OVERVIEW],
[
echo    "*****************************************************"
echo    "*** $PACKAGE_STRING configuration summary"
echo    "***"
echo 	"*** BSM models:		$enable_models"
echo 	"*** Herwig debug mode:	$enable_debug"
echo    "***"
echo    "*** GSL:		$with_gsl"
echo    "***"
echo    "*** ThePEG:		$with_thepeg"
echo    "*** ThePEG headers:	$with_thepeg_headers"
echo    "***"
echo    "*** CLHEP:		$with_clhep"
echo    "*** HepMC:		$with_hepmc"
echo    "***"
echo    "*** KtJet:		$with_ktjet"
dnl echo    "*** FastJet:		$with_fastjet"
echo    "*****************************************************"
])
