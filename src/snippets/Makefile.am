BUILT_SOURCES = done-all-links

snippetsdir = ${pkgdatadir}/snippets

INPUTFILES = \
Analysis-HepMCOutput.in \
Analysis-RivetOutput.in \
CM-13TeV.in \
CM-14TeV.in \
CM-7TeV.in \
CM-8TeV.in \
DipoleShowerParameters.in \
General-IncreaseVerbosity.in \
LHC-StandardFinalization.in \
LHC-StandardInitialization.in \
Particles-SetLonglivedParticlesStable.in \
PDF-CT10.in \
PDF-NNPDF30NLO.in \
Shower-OptionsForMCatNLOMatching.in \
Shower-OptionsForPowhegMatching.in \
Shower-SwitchOffDecay.in \
Shower-SwitchOffHadronization.in \
Shower-SwitchOffPartonShower.in \
Shower-SwitchOffUE.in \
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
