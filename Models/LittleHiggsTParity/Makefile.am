pkglib_LTLIBRARIES = HwLHTP.la
HwLHTP_la_SOURCES = \
LHTPModel.fh     LHTPModel.cc      LHTPModel.h      LHTPModel.icc  \
LHTPWWHVertex.fh LHTPWWHVertex.cc  LHTPWWHVertex.h  LHTPWWHVertex.icc \
LHTPFFGVertex.fh LHTPFFGVertex.cc  LHTPFFGVertex.h  LHTPFFGVertex.icc \
LHTPFFPVertex.fh LHTPFFPVertex.cc  LHTPFFPVertex.h  LHTPFFPVertex.icc \
LHTPFFWVertex.fh LHTPFFWVertex.cc  LHTPFFWVertex.h  LHTPFFWVertex.icc  \
LHTPFFZVertex.fh LHTPFFZVertex.cc  LHTPFFZVertex.h  LHTPFFZVertex.icc 
HwLHTP_la_LDFLAGS = -module -version-info 3:0:0
HwLHTP_la_LIBADD= $(GSLLIBS)
HwLHTP_la_CPPFLAGS= $(AM_CPPFLAGS) $(GSLINCLUDE)
