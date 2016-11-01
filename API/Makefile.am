pkglib_LTLIBRARIES = libHerwigAPI.la


libHerwigAPI_la_SOURCES = \
RunDirectories.h RunDirectories.cc HerwigAPI.h HerwigAPI.cc \
HerwigUI.h

libHerwigAPI_la_LDFLAGS = $(AM_LDFLAGS) -version-info 1:0:0
libHerwigAPI_la_CPPFLAGS = $(AM_CPPFLAGS) \
-DHERWIG_PKGDATADIR="\"$(pkgdatadir)\"" \
-DHERWIG_PKGLIBDIR="\"$(pkglibdir)\"" \
-DTHEPEG_PKGLIBDIR="\"$(THEPEGLIBPATH)\""
