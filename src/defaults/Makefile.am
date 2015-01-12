BUILT_SOURCES = done-all-links

defaultsdir = ${pkgdatadir}/defaults

INPUTFILES = Analysis.in \
baryon_decays.in baryons.in boson_decays.in \
bosons.in Cuts.in decayers.in \
Decays.in DiffractiveParticles.in diquarks.in \
Hadronization.in HerwigDefaults.in HerwigCleanup.in \
lepton_decays.in leptons.in \
masses.in MatrixElements.in meson_decays.in mesons.in Model.in BSM.in \
Particles.in QEDRadiation.in quark_decays.in quarks.in \
setup.gosam.in Shower.in StandardModelVertices.in UnderlyingEvent.in widths.in Partons.in \
MatchboxDefaults.in DipoleShowerDefaults.in MatchboxLoopInduced.in

dist_defaults_DATA = $(INPUTFILES)

EXTRA_DIST = PDF.in.in

defaults_DATA = PDF.in
CLEANFILES = PDF.in done-all-links

## For an explanation of this magic, see autoconf book 4.7.2
if WANT_LOCAL_PDF
edit    = sed -e "s,@HERWIG_PDF_DEFAULT\@,`cd $(top_srcdir) && pwd`/$(HERWIG_PDF_DEFAULT)," -e "s,@HERWIG_PDF_NLO\@,`cd $(top_srcdir) && pwd`/$(HERWIG_PDF_NLO)," -e "s,@HERWIG_PDF_POMERON\@,`cd $(top_srcdir) && pwd`/$(HERWIG_PDF_POMERON),"
installname = $(DESTDIR)$(pkgdatadir)/$(HERWIG_PDF_DEFAULT)
installnameNLO = $(DESTDIR)$(pkgdatadir)/$(HERWIG_PDF_NLO)
installnamePOMERON = $(DESTDIR)$(pkgdatadir)/$(HERWIG_PDF_POMERON)
install-data-hook:
	rm -f $(DESTDIR)$(defaultsdir)/PDF.in
	sed -e 's,@HERWIG_PDF_DEFAULT\@,$(installname),' \
	    -e 's,@HERWIG_PDF_NLO\@,$(installnameNLO),' \
	    -e 's,@HERWIG_PDF_POMERON\@,$(installnamePOMERON),' \
             $(srcdir)/PDF.in.in > $(DESTDIR)$(defaultsdir)/PDF.in
else
edit    = sed -e 's,@HERWIG_PDF_DEFAULT\@,$(HERWIG_PDF_DEFAULT),' -e 's,@HERWIG_PDF_NLO\@,$(HERWIG_PDF_NLO),' -e 's,@HERWIG_PDF_POMERON\@,$(HERWIG_PDF_POMERON),'
endif

PDF.in: Makefile $(srcdir)/PDF.in.in
	@echo "Updating PDF.in"
	@rm -f PDF.in PDF.in.tmp
	@$(edit) $(srcdir)/PDF.in.in > PDF.in.tmp
	@mv PDF.in.tmp PDF.in

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
