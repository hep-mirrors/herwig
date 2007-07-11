pkglib_LTLIBRARIES = 

if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwKtJet.la
endif

if WANT_EVTGEN
AM_CPPFLAGS+= -I$(EVTGENPATH)
pkglib_LTLIBRARIES += HwEvtGen.la
endif

HwKtJet_la_SOURCES = KtJetInterface.cc KtJetInterface.h
HwKtJet_la_LIBADD = $(KTJETLIBS)
HwKtJet_la_LDFLAGS = -module -version-info 2:0:0
HwKtJet_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE)

HwEvtGen_la_SOURCES = \
EvtGen.cc EvtGen.icc EvtGen.fh EvtGen.h \
EvtGenRandom.h \
EvtGenDecayer.cc EvtGenDecayer.fh \
EvtGenDecayer.h EvtGenDecayer.icc 

HwEvtGen_la_LIBADD = -L$(EVTGENPATH)/lib  -levtgenlhc
HwEvtGen_la_LDFLAGS = -R$(EVTGENPATH)/lib -module -version-info 1:0:0
HwEvtGen_la_CXXFLAGS = $(AM_CXXFLAGS) -Wno-non-virtual-dtor -Wno-unused
