if WANT_Anomalous
pkglib_LTLIBRARIES = HwAnomalousCouplings.la
endif
HwAnomalousCouplings_la_SOURCES = \
AnomalousWWWVertex.h AnomalousWWWVertex.cc
HwAnomalousCouplings_la_LDFLAGS = -module -version-info 1:0:0
HwAnomalousCouplings_la_LIBADD= $(GSLLIBS)
HwAnomalousCouplings_la_CPPFLAGS= $(AM_CPPFLAGS) $(GSLINCLUDE)
