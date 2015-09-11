BUILT_SOURCES = done-all-links

snippetsdir = ${pkgdatadir}/snippets

INPUTFILES = \
DipoleShowerParameters.in \
DipoleShowerFiveFlavours.in \
DipoleShowerFourFlavours.in \
Particles-SetLonglivedParticlesStable.in \
PDF-CT10.in \
PDF-NNPDF30NLO.in \
Tune-UE-EE-5-CTEQ6L1.in \
Tune-UE-EE-5-MRST.in

dist_snippets_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
