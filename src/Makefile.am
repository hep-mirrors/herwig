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
SMWZDecays.in Model.in \
StandardModelVertices.in Analysis.in

CLEANFILES = HerwigDefaults.rpo LEP.run LEP.log LEP.out \
  LHC.run LHC.log LHC.out TVT.run TVT.log TVT.out

myENV= HERWIG_USER_MODULES=$(top_builddir)/lib

## checking targets ##

check-local: check-LHC check-LEP check-TVT


check-LEP: 
	$(myENV) ./Herwig++ init
	$(myENV) ./Herwig++ read LEP.in
	$(myENV) ./Herwig++ run LEP.run -N100

check-LHC: 
	$(myENV) ./Herwig++ init
	$(myENV) ./Herwig++ read LHC.in
	$(myENV) ./Herwig++ run LHC.run -N100

check-TVT: 
	$(myENV) ./Herwig++ init
	$(myENV) ./Herwig++ read TVT.in
	$(myENV) ./Herwig++ run TVT.run -N100

######################

install-data-hook:
	rm -f $(DESTDIR)$(pkgdatadir)/HerwigDefaults.in
	sed -e s:../PDF/mrst:$(DESTDIR)$(pkgdatadir)/PDF/mrst:g \
          $(srcdir)/HerwigDefaults.in > $(DESTDIR)$(pkgdatadir)/HerwigDefaults.in
