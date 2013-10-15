EXTRA_DIST = \
AcerDetInterface \
AlpGen \
Analysis2 \
AnomalousHVV \
DecayAnalysis \
HiggsPair \
LeptonME \
MCPWNLO \
PGSInterface \
RadiativeZPrime \
TauAnalysis \
VBFTest

dist-hook:
	rm -rf `find $(distdir) -name '.svn' -or -name 'Makefile'`

all:  
	bash make_makefiles.sh
