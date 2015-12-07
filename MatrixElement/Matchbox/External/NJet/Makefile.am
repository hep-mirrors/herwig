pkglib_LTLIBRARIES = HwMatchboxNJet.la
HwMatchboxNJet_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 12:0:0

HwMatchboxNJet_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(NJETINCLUDEPATH) \
-DNJET_PREFIX="\"$(NJETPREFIX)\"" \
-DNJET_LIBS="\"$(NJETLIBPATH)\""

HwMatchboxNJet_la_SOURCES = \
NJetsAmplitude.h NJetsAmplitude.cc
