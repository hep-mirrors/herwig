pkglib_LTLIBRARIES = 

if WANT_LIBKTJET
AM_CPPFLAGS+= $(KTJETINCLUDE)
pkglib_LTLIBRARIES += HwKtJet.la
endif

HwKtJet_la_SOURCES = KtJetInterface.cc KtJetInterface.h
HwKtJet_la_LIBADD = -L$(KTJETPATH)/lib -lKtEvent
HwKtJet_la_LDFLAGS = -module
