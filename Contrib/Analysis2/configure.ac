dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([Analysis2],[SVN],[sp@particle.uni-karlsruhe.de],[Analysis2])
AC_CONFIG_SRCDIR([Histogram2/Histogram2.h])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

case "${host}" in
  *-darwin[[0156]].*) 
    AC_MSG_ERROR([Analysis2 requires OS X 10.3 or later])
    ;;
  *-darwin7.*) 
    if test "x$MACOSX_DEPLOYMENT_TARGET" != "x10.3"; then 
      AC_MSG_ERROR(
      [Please add 'MACOSX_DEPLOYMENT_TARGET=10.3' to the configure line.]) 
    fi 
    ;;
esac

dnl === disable debug symbols by default =====
if test "x$CXXFLAGS" = "x"; then
   CXXFLAGS=-O3
fi
if test "x$CFLAGS" = "x"; then
   CFLAGS=-O3
fi
dnl ==========================================

AC_LANG([C++])

AM_INIT_AUTOMAKE([1.9 gnu dist-bzip2 -Wall])
dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LN_S

AC_PROG_LIBTOOL

dnl LT_PREREQ([2.2])
dnl LT_INIT([disable-static dlopen pic-only])

ANALYSIS2_CHECK_GSL

ANALYSIS2_CHECK_THEPEG

ANALYSIS2_CHECK_CLHEP
ANALYSIS2_CHECK_KTJET
ANALYSIS2_CHECK_FASTJET

ANALYSIS2_COMPILERFLAGS

AC_CONFIG_FILES([Histogram2/Makefile
         	 Analysis2/Makefile
		 Analysis2/Data/Makefile
      		 src/Makefile
      		 include/Makefile
		 Makefile])

ANALYSIS2_OVERVIEW

AC_OUTPUT
