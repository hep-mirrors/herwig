# check for gcc bug http://gcc.gnu.org/bugzilla/show_bug.cgi?id=34130
AC_DEFUN([HERWIG_CHECK_ABS_BUG],
[
AC_REQUIRE([HERWIG_COMPILERFLAGS])
if test "$GCC" = "yes"; then
AC_MSG_CHECKING([for gcc abs bug])
AC_RUN_IFELSE([
	AC_LANG_PROGRAM(
		[[ int foo (int i) { return -2 * __builtin_abs(i - 2); } ]],
		[[ if ( foo(1) != -2 || foo(3) != -2 ) return 1; ]]
	)],
	[ AC_MSG_RESULT([not found. Compiler is ok.]) ],
	[
	AC_MSG_RESULT([found. Builtin abs() is buggy.])
	AC_MSG_CHECKING([if -fno-builtin-abs works])
	oldcxxflags=$CXXFLAGS
	CXXFLAGS="$CXXFLAGS -fno-builtin-abs"
	AC_RUN_IFELSE([
		AC_LANG_PROGRAM(
			[[
			#include <cstdlib>
			int foo (int i) { return -2 * std::abs(i - 2); }
			]],
			[[
			if (foo(1) != -2 || foo(3) != -2) return 1; 
			]]
		)],
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

dnl ##### THEPEG #####
AC_DEFUN([HERWIG_CHECK_THEPEG],
[
defaultlocation="${prefix}"
test "x$defaultlocation" = xNONE && defaultlocation="${ac_default_prefix}"
AC_MSG_CHECKING([for libThePEG in])
AC_ARG_WITH(thepeg,
        AC_HELP_STRING([--with-thepeg=DIR],[location of ThePEG installation]),
        [],
	[with_thepeg="${defaultlocation}"])
AC_MSG_RESULT([$with_thepeg])

if test "x$with_thepeg" = "xno"; then
	AC_MSG_ERROR([Cannot build Herwig++ without ThePEG. Please set --with-thepeg.])
fi

THEPEGLDFLAGS="-L${with_thepeg}/lib/ThePEG"
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

AC_MSG_CHECKING([for HepMCAnalysis.so in ThePEG])


if test -x "$THEPEGPATH/lib/ThePEG/HepMCAnalysis.so" ; then
     	CREATE_HEPMC="create"
	AC_MSG_RESULT([found])
else
	CREATE_HEPMC="# create"
	AC_MSG_RESULT([not found])
fi

AC_SUBST([CREATE_HEPMC])
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
      x86_64-*|*-darwin1*)
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
HERWIG_PDF_PREFIX=${with_pdf}/share/Herwig++PDF

if test -f "${HERWIG_PDF_PREFIX}/mrst/2008/mrstMCal.dat"; then
	AC_MSG_RESULT([$with_pdf])
	localPDFneeded=false
else
	AC_MSG_RESULT([Using built-in PDF data set. For other data sets, set --with-pdf.])
	HERWIG_PDF_PREFIX=PDF
	localPDFneeded=true
fi
HERWIG_PDF_DEFAULT=${HERWIG_PDF_PREFIX}/mrst/2008/mrstMCal.dat
HERWIG_PDF_NLO=${HERWIG_PDF_PREFIX}/mrst/2002/mrst2002nlo.dat
HERWIG_PDF_POMERON=${HERWIG_PDF_PREFIX}/diffraction/

AM_CONDITIONAL(WANT_LOCAL_PDF,[test "x$localPDFneeded" = "xtrue"])
AC_SUBST(HERWIG_PDF_DEFAULT)
AC_SUBST(HERWIG_PDF_NLO)
AC_SUBST(HERWIG_PDF_POMERON)
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

dnl do an actual capability check on ld instead of this workaround
case "${host}" in
  *-darwin*) 
     ;;
  *)
     AM_LDFLAGS="-Wl,--enable-new-dtags"
     ;;
esac

AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CFLAGS,  ["$warnflags $debugflags"])
AC_SUBST(AM_CXXFLAGS,["$warnflags $debugflags"])
AC_SUBST(AM_FCFLAGS,  ["$debugflags"])
AC_SUBST(AM_LDFLAGS)
])

AC_DEFUN([HERWIG_ENABLE_MODELS],
[
AC_MSG_CHECKING([for BSM models to include])

LOAD_RS=""
LOAD_SUSY=""
LOAD_NMSSM=""
LOAD_TRP=""
LOAD_UED=""
LOAD_ADD=""
LOAD_SEXTET=""

AC_ARG_ENABLE(models,
        AC_HELP_STRING([--enable-models=LIST],[Comma-separated list of BSM models to enable. Options are (mssm nmssm ued rs trp add leptoquarks sextet) or --disable-models to turn them all off.]),
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

if test "$nmssm" -o "$all"; then
   LOAD_NMSSM="library HwNMSSM.so"
   mssm=yes
fi
AC_SUBST(LOAD_NMSSM)

if test "$rs" -o "$all" ; then
   LOAD_RS="library HwRSModel.so"
fi
AC_SUBST(LOAD_RS)

if test "$mssm" -o "$all"; then
   LOAD_SUSY="library HwSusy.so"
fi
AC_SUBST(LOAD_SUSY)

if test "$trp" -o "$all"; then
   LOAD_TRP="library HwTransplanck.so"
fi
AC_SUBST(LOAD_TRP)

if test "$ued" -o "$all"; then
   LOAD_UED="library HwUED.so"
fi
AC_SUBST(LOAD_UED)

if test "$add" -o "$all"; then
   LOAD_ADD="library HwADDModel.so"
fi
AC_SUBST(LOAD_ADD)

if test "$leptoquarks" -o "$all"; then
   LOAD_LEPTOQUARKS="library HwLeptoquarkModel.so"
fi
AC_SUBST(LOAD_LEPTOQUARKS)

if test "$sextet" -o "$all"; then
   LOAD_SEXTET="library HwSextetModel.so"
fi
AC_SUBST(LOAD_SEXTET)

AM_CONDITIONAL(WANT_MSSM,[test "$mssm" -o "$all"])
AM_CONDITIONAL(WANT_NMSSM,[test "$nmssm" -o "$all"])
AM_CONDITIONAL(WANT_UED,[test "$ued" -o "$all"])
AM_CONDITIONAL(WANT_RS,[test "$rs" -o "$all"])
AM_CONDITIONAL(WANT_Leptoquark,[test "$leptoquarks" -o "$all"])
AM_CONDITIONAL(WANT_TRP,[test "$trp" -o "$all"])
AM_CONDITIONAL(WANT_ADD,[test "$add" -o "$all"])
AM_CONDITIONAL(WANT_SEXTET,[test "$sextet" -o "$all"])
])

AC_DEFUN([HERWIG_OVERVIEW],
[
FCSTRING=`$FC --version | head -1`
CXXSTRING=`$CXX --version | head -1`
CCSTRING=`$CC --version | head -1`
cat << _HW_EOF_ > config.herwig
*****************************************************
*** $PACKAGE_STRING configuration summary
*** Please include this information in bug reports!
***--------------------------------------------------
*** Prefix:		$prefix
***
*** BSM models:		$enable_models
*** Herwig debug mode:	$enable_debug
***
*** GSL:		$with_gsl
***
*** ThePEG:		$with_thepeg
*** ThePEG headers:	$with_thepeg_headers
***
*** Fastjet:		${fjconfig}
***
*** Host:		$host
*** CXX:		$CXXSTRING
*** FC:			$FCSTRING
*** CC:			$CCSTRING
*****************************************************
_HW_EOF_
])
