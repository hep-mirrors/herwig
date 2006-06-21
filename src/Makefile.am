bin_PROGRAMS = Herwig++
Herwig___SOURCES = Herwig++.cc
Herwig___LDFLAGS = -export-dynamic $(THEPEGLDFLAGS) $(CLHEPLDFLAGS)
Herwig___LDADD = $(top_builddir)/lib/libHerwig.la $(THEPEGLIB) $(CLHEPLIB) -ldl

dist_pkgdata_DATA = \
Decays.in LHC.in TVT.in Hadronization.in \
HerwigDefaults.in   Shower.in \
HerwigParticles.in  SMWZDecays.in \
Hw64Decays.in       ThePEGParticles.in \
LEPbb.in            UnderlyingEvent.in \
LEP.in Herwig++.cc Herwig.Makefile  \
SMWZDecays.in Model.in ILC.in \
StandardModelVertices.in Analysis.in

pkgdata_DATA = PDF.in

CLEANFILES = HerwigDefaults.rpo LEP.run LEP.log LEP.out \
  LHC.run LHC.log LHC.out TVT.run TVT.log TVT.out \
  ILC.run ILC.log ILC.out

myENV= HERWIG_USER_MODULES=$(top_builddir)/lib

## checking targets ##

check-local: check-LHC check-LEP check-TVT check-ILC


check-LEP: 
	$(myENV) ./Herwig++ init
	$(myENV) ./Herwig++ read LEP.in
	$(myENV) ./Herwig++ run LEP.run -N100

check-ILC: 
	$(myENV) ./Herwig++ init
	$(myENV) ./Herwig++ read ILC.in
	$(myENV) ./Herwig++ run ILC.run -N100

check-LHC: 
	$(myENV) ./Herwig++ init
	$(myENV) ./Herwig++ read LHC.in
	$(myENV) ./Herwig++ run LHC.run -N100

check-TVT: 
	$(myENV) ./Herwig++ init
	$(myENV) ./Herwig++ read TVT.in
	$(myENV) ./Herwig++ run TVT.run -N100

if WANT_LOCAL_PDF
install-data-hook:
	rm -f $(DESTDIR)$(pkgdatadir)/PDF.in
	sed -e s:../PDF/mrst:$(DESTDIR)$(pkgdatadir)/PDF/mrst:g \
             PDF.in > $(DESTDIR)$(pkgdatadir)/PDF.in

endif
