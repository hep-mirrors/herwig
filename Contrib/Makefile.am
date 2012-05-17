EXTRA_DIST = \
AcerDetInterface \
Analysis2 \
AnomalousHVV \
DecayAnalysis \
LeptonME \
PGSInterface \
RadiativeZPrime \
TauAnalysis \
AlpGen


dist-hook:
	rm -rf `find $(distdir) -name '.svn' -or -name 'Makefile'`

all:  
	bash make_makefiles.sh
