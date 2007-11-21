pkglib_LTLIBRARIES = 

if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwKtJet.la
endif

HwKtJet_la_SOURCES = KtJetInterface.cc KtJetInterface.h
HwKtJet_la_LIBADD = $(KTJETLIBS) -L$(CLHEPPATH)  $(CLHEPLIB)
HwKtJet_la_LDFLAGS = -module -version-info 3:0:0
HwKtJet_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)

if WANT_LIBROOT
pkglib_LTLIBRARIES += HwRoot.la
endif

HwRoot_la_SOURCES = RootInterface.cc RootInterface.h
HwRoot_la_LIBADD =  $(ROOTLIBS)
HwRoot_la_LDFLAGS = -R$(ROOTLIBPATH) -module -version-info 3:0:0
HwRoot_la_CPPFLAGS = $(AM_CPPFLAGS) $(ROOTINCLUDE)

if WANT_EVTGEN
pkglib_LTLIBRARIES += HwEvtGen.la
endif

HwEvtGen_la_SOURCES = \
EvtGen.cc EvtGen.icc EvtGen.fh EvtGen.h \
EvtGenRandom.h \
EvtGenDecayer.cc EvtGenDecayer.fh \
EvtGenDecayer.h EvtGenDecayer.icc 
HwEvtGen_la_LIBADD = -L$(EVTGENPATH)/lib  -levtgenlhc
HwEvtGen_la_LDFLAGS = -R$(EVTGENPATH)/lib -module -version-info 2:0:0
HwEvtGen_la_CXXFLAGS = $(AM_CXXFLAGS) -Wno-non-virtual-dtor -Wno-unused
HwEvtGen_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(EVTGENPATH)

pkglib_LTLIBRARIES += HwPOWHEG.la
HwPOWHEG_la_LDFLAGS = -R$(EVTGENPATH)/lib -module 

HwPOWHEG_la_SOURCES =\
POWHEGReader.cc POWHEGReader.fh \
POWHEGReader.h POWHEGReader.icc 
