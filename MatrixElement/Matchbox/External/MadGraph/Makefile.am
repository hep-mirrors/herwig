pkglib_LTLIBRARIES = HwMatchboxMadGraph.la
HwMatchboxMadGraph_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 3:1:0

HwMatchboxMadGraph_la_SOURCES = \
MadGraphAmplitude.h MadGraphAmplitude.cc

madgraphdir = ${pkgdatadir}/MadGraphInterface

INPUTFILES = \
InterfaceMadGraph.f.in \
param_card.dat.in \
param_card_heft.dat.in

dist_madgraph_DATA = $(INPUTFILES)

HwMatchboxMadGraph_la_CPPFLAGS = $(AM_CPPFLAGS) \
-DHERWIG_BINDIR="\"$(bindir)\"" \
-DHERWIG_INCLUDEDIR="\"$(includedir)\"" \
-DHERWIG_PKGDATADIR="\"$(pkgdatadir)\"" \
-DMADGRAPH_PREFIX="\"$(MADGRAPHPREFIX)\""

EXTRA_DIST = mg2herwig.py

install-exec-local:
	$(install_sh_SCRIPT) $(builddir)/mg2herwig.py $(DESTDIR)$(bindir)/mg2herwig

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/mg2herwig
