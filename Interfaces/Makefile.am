pkglib_LTLIBRARIES = 

if WANT_LIBKTJET
pkglib_LTLIBRARIES += HwKtJet.la
endif

HwKtJet_la_SOURCES = KtJetInterface.cc KtJetInterface.h
HwKtJet_la_LIBADD = $(KTJETLIBS) -L$(CLHEPPATH)  $(CLHEPLIB)
HwKtJet_la_LDFLAGS = -module -version-info 3:0:0
HwKtJet_la_CPPFLAGS = $(AM_CPPFLAGS) $(KTJETINCLUDE) -I$(KTJETPATH)
