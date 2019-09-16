BUILT_SOURCES = done-all-links

Matchboxdir = ${pkgdatadir}/Matchbox

INPUTFILES = \
CMEC.in \
DefaultEEJets.in \
DefaultEPJets.in \
DefaultPPJets.in \
DiagonalCKM.in \
FiveFlavourNoBMassScheme.in \
FiveFlavourScheme.in \
FourFlavourScheme.in \
GoSam-GoSam.in \
HiggsEffective.in \
HJets.in \
IdentifiedBs.in \
InclusiveDurhamJets.in \
IncreaseVerbosity.in \
KrkNLO-DipoleShower.in \
LO-DefaultShower.in \
LO-DipoleShower.in \
LO.in \
LO-NoShower.in \
MadGraph-GoSam.in \
MadGraph-MadGraph.in \
MadGraph-NJet.in \
MadGraph-OpenLoops.in \
MCatLO-DefaultShower.in \
MCatLO-DipoleShower.in \
MCatNLO-DefaultShower.in \
MCatNLO-Dipole-HardAlphaSTune.in \
MCatNLO-DipoleShower.in \
MMHT2014.in \
CT14.in \
MuDown.in \
MuQDown.in \
MuQUp.in \
MuUp.in \
NJet-NJet.in \
NLO-NoShower.in \
NonDiagonalCKM.in \
OnShellHProduction.in \
OnShellTopProduction.in \
OnShellWProduction.in \
OnShellZProduction.in \
OpenLoops-OpenLoops.in \
Powheg-DefaultShower.in \
Powheg-DipoleShower.in \
Powheg.in \
PQCDLevel.in \
StandardModelLike.in \
ShowerBenchmarks.in \
VBFDiagramsOnly.in \
VBFNLO.in \
VBFNLOPhasespace.in

dist_Matchbox_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
