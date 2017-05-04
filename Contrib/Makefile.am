EXTRA_DIST = \
AcerDetInterface \
AlpGen \
AnomalousHVV \
DecayAnalysis \
FxFx \
HiggsPair \
HiggsPairOL \
LeptonME \
MultiWeight \
PGSInterface \
RadiativeZPrime \
ShowerVeto \
TauAnalysis \
VBFTest

dist-hook:
	rm -rf `find $(distdir) -name '.svn' -or -name 'Makefile' -or -name '.hg'`

all:  
	bash make_makefiles.sh
