SUBDIRS = \
Amplitudes \
PP2Jets

pkglib_LTLIBRARIES = HwMatchboxNLOJet.la
HwMatchboxNLOJet_la_LDFLAGS = -module -version-info 11:1:0

HwMatchboxNLOJet_la_LIBADD = \
-L$(NLOJETLIBPATH) -lnlo++ \
PP2Jets/libHwNLOJetPP2Jets.la \
Amplitudes/libHwNLOJetAmplitudes.la

HwMatchboxNLOJet_la_CPPFLAGS = $(AM_CPPFLAGS)
HwMatchboxNLOJet_la_CPPFLAGS += -I$(NLOJETINCLUDEPATH) -I$(NLOJETINCLUDEPATH)/nlo++

HwMatchboxNLOJet_la_SOURCES = \
NLOJetRandomWrapper.h \
NLOJetPhasespace.h NLOJetPhasespace.tcc NLOJetPhasespace.cc \
NLOJetAmplitude.h NLOJetAmplitude.tcc NLOJetAmplitude.cc \
NLOJetMEBase.h NLOJetMEBase.tcc NLOJetMEBase.cc
