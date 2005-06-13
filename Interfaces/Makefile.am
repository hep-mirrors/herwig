pkglib_LTLIBRARIES = 

if WANT_LIBAMEGIC
pkglib_LTLIBRARIES += libHwAmegic.la 
endif

if WANT_LIBKTJET
pkglib_LTLIBRARIES += libHwKtJet.la
endif

libHwAmegic_la_SOURCES = AmegicInterface.cc  AmegicInterface.icc \
		AmegicInterface.h
libHwKtJet_la_SOURCES = KtJetInterface.cc KtJetInterface.h
