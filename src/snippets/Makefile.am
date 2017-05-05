BUILT_SOURCES = done-all-links

snippetsdir = ${pkgdatadir}/snippets

INPUTFILES = \
DipoleShowerFiveFlavours.in \
DipoleShowerFourFlavours.in \
Particles-SetLonglivedParticlesStable.in \
PDF-CT10.in \
PDF-NNPDF30NLO.in \
Matchbox.in \
MonacoSampler.in \
DipoleMerging.in \
EECollider.in EPCollider.in PPCollider.in \
MB.in \
MB-DipoleShower.in \
SoftModel.in \
Diffraction.in

dist_snippets_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
