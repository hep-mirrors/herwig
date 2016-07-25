EXTRA_DIST = \
AcerDetInterface \
AlpGen \
Analysis2 \
AnomalousHVV \
DecayAnalysis \
FxFx \
HiggsPair \
HiggsPairOL \
LeptonME \
PGSInterface \
RadiativeZPrime \
TauAnalysis \
MultiWeight \
ShowerVeto

dist-hook:
	rm -rf `find $(distdir) -name '.svn' -or -name 'Makefile' -or -name '.hg'`

all:  
	bash make_makefiles.sh
