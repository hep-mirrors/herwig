pkglib_LTLIBRARIES = HwMatchboxMadGraphfortran.la
HwMatchboxMadGraphfortran_la_LDFLAGS = -module -version-info 1:0:0

HwMatchboxMadGraphfortran_la_SOURCES = \
MadGraphAmplitudefortran.h MadGraphAmplitudefortran.cc

madgraphfortrandir = ${pkgdatadir}/MadGraphInterfacefortran

INPUTFILES = \
InterfaceMadGraphfortran.f.in \
param_card.dat.in

dist_madgraphfortran_DATA = $(INPUTFILES)

install-exec-local:
	$(install_sh_SCRIPT) $(builddir)/mg2Matchboxfortran.py $(DESTDIR)$(bindir)/mg2Matchboxfortran.py

uninstall-local:
	rm -f  $(DESTDIR)$(bindir)/mg2Matchboxfortran.py

