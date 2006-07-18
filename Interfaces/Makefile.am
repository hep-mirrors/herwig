pkglib_LTLIBRARIES = 

if WANT_LIBAMEGIC
pkglib_LTLIBRARIES += HwAmegic.la 
endif

if WANT_LIBKTJET
AM_CPPFLAGS+= $(KTJETINCLUDE)
pkglib_LTLIBRARIES += HwKtJet.la
endif

HwAmegic_la_SOURCES = AmegicInterface.cc  AmegicInterface.icc \
		AmegicInterface.h
HwAmegic_la_LDFLAGS = -module

HwKtJet_la_SOURCES = KtJetInterface.cc KtJetInterface.h
HwKtJet_la_LIBADD = -L$(KTJETPATH)/lib -lKtEvent
HwKtJet_la_LDFLAGS = -module
