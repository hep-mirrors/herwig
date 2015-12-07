if WANT_BSM
pkglib_LTLIBRARIES = HwLHTPModel.la
endif
HwLHTPModel_la_SOURCES = \
LHTPModel.fh     LHTPModel.cc      LHTPModel.h \
LHTPWWHVertex.cc  LHTPWWHVertex.h \
LHTPFFGVertex.cc  LHTPFFGVertex.h \
LHTPFFPVertex.cc  LHTPFFPVertex.h \
LHTPFFWVertex.cc  LHTPFFWVertex.h \
LHTPFFZVertex.cc  LHTPFFZVertex.h \
LHTPWWWVertex.cc  LHTPWWWVertex.h \
LHTPHHHVertex.cc  LHTPHHHVertex.h \
LHTPWHHVertex.cc  LHTPWHHVertex.h \
LHTPFFHVertex.cc  LHTPFFHVertex.h
HwLHTPModel_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 4:0:0
HwLHTPModel_la_LIBADD= $(GSLLIBS)
HwLHTPModel_la_CPPFLAGS= $(AM_CPPFLAGS) $(GSLINCLUDE)
