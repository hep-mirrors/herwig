pkglib_LTLIBRARIES = HwMatchboxOpenLoops.la
HwMatchboxOpenLoops_la_SOURCES = \
OpenLoopsAmplitude.h OpenLoopsAmplitude.cc
HwMatchboxOpenLoops_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 12:0:0

HwMatchboxOpenLoops_la_CPPFLAGS = $(AM_CPPFLAGS) \
-DOPENLOOPSLIBS="\"$(OPENLOOPSLIBS)\"" \
-DOPENLOOPSPREFIX="\"$(OPENLOOPSPREFIX)\""
