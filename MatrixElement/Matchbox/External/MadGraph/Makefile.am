EXTRA_DIST = \
MadGraphAmplitude.h \
MadGraphAmplitude.cc \
mg2herwig.py


madgraphdir = ${pkgdatadir}/MadGraphInterface

INPUTFILES = \
InterfaceMadGraph.f.in \
param_card.dat.in \
param_card_heft.dat.in

dist_madgraph_DATA = $(INPUTFILES)



install-exec-local:
	$(install_sh_SCRIPT) $(builddir)/mg2herwig.py $(DESTDIR)$(bindir)/mg2herwig

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/mg2herwig
