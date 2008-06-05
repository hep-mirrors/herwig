BUILT_SOURCES = done-all-links

defaultsdir = ${pkgdatadir}/defaults

INPUTFILES = Analysis.in Analysis2.in LEPAnalysis2.in\
baryon_decays.in baryons.in boson_decays.in \
bosons.in Cuts.in decayers.in \
Decays.in DiffractiveParticles.in diquarks.in \
Hadronization.in HerwigDefaults.in HerwigCleanup.in \
lepton_decays.in leptons.in \
masses.in MatrixElements.in meson_decays.in mesons.in Model.in \
Particles.in QEDRadiation.in quark_decays.in quarks.in RSModelVertices.in \
Shower.in StandardModelVertices.in UnderlyingEvent.in widths.in Partons.in

dist_defaults_DATA = $(INPUTFILES)

EXTRA_DIST = PDF.in.in

defaults_DATA = PDF.in
CLEANFILES = PDF.in done-all-links

## For an explanation of this magic, see autoconf book 4.7.2
if WANT_LOCAL_PDF
edit    = sed -e "s,@HERWIG_PDF_DEFAULT\@,`cd $(top_srcdir) && pwd`/$(HERWIG_PDF_DEFAULT)," 
editNLO = sed -e "s,@HERWIG_PDF_NLO\@,`cd $(top_srcdir) && pwd`/$(HERWIG_PDF_NLO),"
installname = $(DESTDIR)$(pkgdatadir)/$(HERWIG_PDF_DEFAULT)
installnameNLO = $(DESTDIR)$(pkgdatadir)/$(HERWIG_PDF_NLO)
install-data-hook:
	rm -f $(DESTDIR)$(defaultsdir)/PDF.in
	sed -e 's,@HERWIG_PDF_DEFAULT\@,$(installname),' \
             $(srcdir)/PDF.in.in > $(srcdir)/PDF.in.tmp
	sed -e 's,@HERWIG_PDF_NLO\@,$(installnameNLO),' \
             $(srcdir)/PDF.in.tmp > $(DESTDIR)$(defaultsdir)/PDF.in
	rm -f $(srcdir)/PDF.in.tmp
else
edit    = sed -e 's,@HERWIG_PDF_DEFAULT\@,$(HERWIG_PDF_DEFAULT),'
editNLO = sed -e 's,@HERWIG_PDF_NLO\@,$(HERWIG_PDF_NLO),'
endif

PDF.in: Makefile $(srcdir)/PDF.in.in
	@echo "Updating PDF.in"
	@rm -f PDF.in PDF.in.tmp
	@$(edit) $(srcdir)/PDF.in.in > PDF.in.tmp
	@$(editNLO) PDF.in.tmp > $(srcdir)/PDF.in
	@rm PDF.in.tmp

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
