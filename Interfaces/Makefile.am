pkglib_LTLIBRARIES = 

if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwKtJet.la
endif

HwKtJet_la_SOURCES = KtJetInterface.cc KtJetInterface.h
HwKtJet_la_LIBADD = $(KTJETLIBS) $(CLHEPLIB)
HwKtJet_la_LDFLAGS = -module -version-info 3:0:0
HwKtJet_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)

pkglib_LTLIBRARIES += HwPOWHEG.la
HwPOWHEG_la_LDFLAGS = -R$(EVTGENPATH)/lib -module 

HwPOWHEG_la_SOURCES =\
POWHEGReader.cc POWHEGReader.fh \
POWHEGReader.h POWHEGReader.icc 
