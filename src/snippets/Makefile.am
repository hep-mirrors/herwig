BUILT_SOURCES = done-all-links

snippetsdir = ${pkgdatadir}/snippets

INPUTFILES = \
CellGridSampler.in \
Diffraction.in \
DipoleMerging.in \
DipoleShowerFiveFlavours.in \
DipoleShowerFourFlavours.in \
EECollider.in \
EPCollider.in \
HepMCFixedOrder.in \
HepMC.in \
Matchbox.in \
MB-DipoleShower.in \
MB.in \
MonacoSampler.in \
Particles-SetLonglivedParticlesStable.in \
PDF-CT10.in \
PDF-NNPDF30NLO.in \
PPCollider.in \
RivetFixedOrder.in \
Rivet.in \
SoftModel.in \
SoftTune.in \
YFS.in \
BaryonicReconnection.in \
Dipole_AutoTune_prel.in

dist_snippets_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
