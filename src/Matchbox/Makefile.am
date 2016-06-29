BUILT_SOURCES = done-all-links

Matchboxdir = ${pkgdatadir}/Matchbox

INPUTFILES = \
DefaultEEJets.in \
DefaultEPJets.in \
DefaultPPJets.in \
EECollider.in \
EPCollider.in \
FiveFlavourScheme.in \
FiveFlavourNoBMassScheme.in \
FourFlavourScheme.in \
GenericCollider.in \
HiggsEffective.in \
HJets.in \
IdentifiedBs.in \
InclusiveDurhamJets.in \
IncreaseVerbosity.in \
LO-DefaultShower.in \
LO-DipoleShower.in \
MCatLO-DefaultShower.in \
MCatLO-DipoleShower.in \
LO.in \
LO-NoShower.in \
MadGraph-MadGraph.in \
GoSam-GoSam.in \
OpenLoops-OpenLoops.in \
MadGraph-GoSam.in \
MadGraph-NJet.in \
NJet-NJet.in \
MadGraph-OpenLoops.in \
MCatNLO-DefaultShower.in \
MCatNLO-DipoleShower.in \
MuDown.in \
MuQDown.in \
MuQUp.in \
MuUp.in \
NLO-NoShower.in \
OnShellHProduction.in \
OnShellWProduction.in \
OnShellZProduction.in \
OnShellTopProduction.in \
Powheg-DefaultShower.in \
Powheg-DipoleShower.in \
Powheg.in \
PPCollider.in \
PQCDLevel.in \
VBFNLO.in \
VBFNLOPhasespace.in \
VBFDiagramsOnly.in \
MMHT2014.in \
MCatNLO-Default-HardAlphaSTune.in \
MCatNLO-Default-ShowerAlphaSTune.in \
Powheg-Default-HardAlphaSTune.in \
Powheg-Default-ShowerAlphaSTune.in \
MCatNLO-Dipole-HardAlphaSTune.in \
MCatNLO-Dipole-ShowerAlphaSTune.in \
StandardModelLike.in \
DiagonalCKM.in \
NonDiagonalCKM.in

dist_Matchbox_DATA = $(INPUTFILES)

CLEANFILES = done-all-links

done-all-links: $(INPUTFILES)
	@echo "Linking input files"
	@for i in $(INPUTFILES); do \
	if test -f $(srcdir)/$$i -a ! -e $$i; then \
	$(LN_S) -f $(srcdir)/$$i; fi; done
	@touch done-all-links
