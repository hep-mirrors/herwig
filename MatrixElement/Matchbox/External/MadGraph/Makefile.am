pkglib_LTLIBRARIES = HwMatchboxMadGraph.la
HwMatchboxMadGraph_la_LDFLAGS = -module -version-info 1:0:0

HwMatchboxMadGraph_la_SOURCES = \
MadGraphAmplitude.h MadGraphAmplitude.cc

madgraphdir = ${pkgdatadir}/MadGraphInterface

INPUTFILES = \
InterfaceMadGraph.f.in \
param_card.dat.in

dist_madgraph_DATA = $(INPUTFILES)

install-exec-local:
	$(install_sh_SCRIPT) $(builddir)/mg2Matchbox.py $(DESTDIR)$(bindir)/mg2Matchbox.py

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/mg2Matchbox.py

