if WANT_LHTP
pkglib_LTLIBRARIES = HwLHTPModel.la
endif
HwLHTPModel_la_SOURCES = \
LHTPModel.fh     LHTPModel.cc      LHTPModel.h      LHTPModel.icc  \
LHTPWWHVertex.fh LHTPWWHVertex.cc  LHTPWWHVertex.h  LHTPWWHVertex.icc \
LHTPFFGVertex.fh LHTPFFGVertex.cc  LHTPFFGVertex.h  LHTPFFGVertex.icc \
LHTPFFPVertex.fh LHTPFFPVertex.cc  LHTPFFPVertex.h  LHTPFFPVertex.icc \
LHTPFFWVertex.fh LHTPFFWVertex.cc  LHTPFFWVertex.h  LHTPFFWVertex.icc  \
LHTPFFZVertex.fh LHTPFFZVertex.cc  LHTPFFZVertex.h  LHTPFFZVertex.icc 
HwLHTPModel_la_LDFLAGS = -module -version-info 3:0:0
HwLHTPModel_la_LIBADD= $(GSLLIBS)
HwLHTPModel_la_CPPFLAGS= $(AM_CPPFLAGS) $(GSLINCLUDE)
