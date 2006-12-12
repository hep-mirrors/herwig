pkglib_LTLIBRARIES = 

if WANT_LIBKTJET
AM_CPPFLAGS+= $(KTJETINCLUDE)
pkglib_LTLIBRARIES += HwKtJet.la
endif

if WANT_EVTGEN
AM_CPPFLAGS+= -I$(EVTGENPATH)/include
pkglib_LTLIBRARIES += HwEvtGen.la
endif

HwKtJet_la_SOURCES = KtJetInterface.cc KtJetInterface.h
HwKtJet_la_LIBADD = -L$(KTJETPATH)/lib -lKtEvent
HwKtJet_la_LDFLAGS = -module

HwEvtGen_la_SOURCES = \
EvtGen.cc EvtGenDecayer.cc EvtGenDecayer.fh \
EvtGenDecayer.h EvtGenDecayer.icc EvtGen.fh EvtGen.h \
EvtGen.icc EvtGenRandom.h
HwEvtGen_la_LIBADD = -L$(EVTGENPATH)/lib -lEvtGen
HwEvtGen_la_LDFLAGS = -module
HwEvtGen_la_CXXFLAGS = $(AM_CXXFLAGS) -Wno-non-virtual-dtor -Wno-unused
