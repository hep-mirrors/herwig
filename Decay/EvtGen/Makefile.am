pkglib_LTLIBRARIES = HwEvtGenInterface.la

HwEvtGenInterface_la_SOURCES = \
EvtGenInterface.cc EvtGenInterface.h EvtGenInterface.fh \
EvtGenRandom.h \
EvtGenDecayer.cc EvtGenDecayer.h

HwEvtGenInterface_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 2:0:0
HwEvtGenInterface_la_LIBADD = $(EVTGENLIBS) $(PYTHIA8LIB)
HwEvtGenInterface_la_CPPFLAGS = $(AM_CPPFLAGS) $(EVTGENINCLUDE) \
-DEVTGEN_PREFIX="\"$(EVTGENPREFIX)\"" \
-DPYTHIA8DATA="\"$(PYTHIA8DATA)\""
