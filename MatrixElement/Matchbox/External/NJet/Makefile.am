pkglib_LTLIBRARIES = HwMatchboxNJet.la
HwMatchboxNJet_la_LDFLAGS = -module -version-info 11:1:0

HwMatchboxNJet_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(NJETINCLUDEPATH) \
-DNJET_PREFIX="\"$(NJETPREFIX)\""

HwMatchboxNJet_la_SOURCES = \
NJetsAmplitude.h NJetsAmplitude.cc
