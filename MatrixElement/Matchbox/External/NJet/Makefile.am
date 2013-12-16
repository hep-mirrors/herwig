pkglib_LTLIBRARIES = HwMatchboxNJet.la
HwMatchboxNJet_la_LDFLAGS = -module -version-info 11:1:0

HwMatchboxNJet_la_LIBADD = \
-L$(NJETLIBPATH) -lnjet2

HwMatchboxNJet_la_CPPFLAGS = $(AM_CPPFLAGS)
HwMatchboxNJet_la_CPPFLAGS += -I$(NJETINCLUDEPATH) -I$(NJETINCLUDEPATH)

HwMatchboxNJet_la_SOURCES = \
NJetsAmplitude.h NJetsAmplitude.cc

