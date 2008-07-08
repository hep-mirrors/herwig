if WANT_LHTP
pkglib_LTLIBRARIES = HwLHTPModel.la
endif
HwLHTPModel_la_SOURCES = \
LHTPModel.fh     LHTPModel.cc      LHTPModel.h \
LHTPWWHVertex.fh LHTPWWHVertex.cc  LHTPWWHVertex.h \
LHTPFFGVertex.fh LHTPFFGVertex.cc  LHTPFFGVertex.h \
LHTPFFPVertex.fh LHTPFFPVertex.cc  LHTPFFPVertex.h \
LHTPFFWVertex.fh LHTPFFWVertex.cc  LHTPFFWVertex.h \
LHTPFFZVertex.fh LHTPFFZVertex.cc  LHTPFFZVertex.h \
LHTPWWWVertex.fh LHTPWWWVertex.cc  LHTPWWWVertex.h
HwLHTPModel_la_LDFLAGS = -module -version-info 3:0:0
HwLHTPModel_la_LIBADD= $(GSLLIBS)
HwLHTPModel_la_CPPFLAGS= $(AM_CPPFLAGS) $(GSLINCLUDE)
