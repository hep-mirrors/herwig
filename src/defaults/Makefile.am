BUILT_SOURCES = done-all-links

defaultsdir = ${pkgdatadir}/defaults

INPUTFILES = Analysis.in \
baryon_decays.in baryons.in boson_decays.in \
bosons.in Cuts.in decayers.in \
Decays.in DiffractiveParticles.in diquarks.in \
Hadronization.in HerwigDefaults.in \
lepton_decays.in leptons.in \
masses.in MatrixElements.in meson_decays.in mesons.in Model.in BSM.in \
Particles.in QEDRadiation.in quark_decays.in quarks.in \
setup.gosam.in Shower.in StandardModelVertices.in UnderlyingEvent.in widths.in Partons.in \
MatchboxDefaults.in DipoleShowerDefaults.in MatchboxLoopInduced.in Samplers.in \
EvtGenDecayer.in EvtGenBDecays.in HerwigBDecays.in MatchboxMergingDefaults.in 

dist_defaults_DATA = $(INPUTFILES)

EXTRA_DIST = PDF.in.in

defaults_DATA = PDF.in
CLEANFILES = PDF.in done-all-links

## For an explanation of this magic, see autoconf book 4.7.2

edit    = sed -e "s,@HERWIG_PDF_POMERON\@,`cd $(top_srcdir) && pwd`/PDF/diffraction/,"

## may still be broken and need $(DESTDIR)
installnamePOMERON = $(pkgdatadir)/PDF/diffraction/

install-data-hook:
	rm -f $(DESTDIR)$(defaultsdir)/PDF.in
	sed -e 's,@HERWIG_PDF_POMERON\@,$(installnamePOMERON),' \
             $(srcdir)/PDF.in.in > $(DESTDIR)$(defaultsdir)/PDF.in

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
