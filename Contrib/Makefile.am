EXTRA_DIST = Analysis2 AnomalousHVV DecayAnalysis LeptonME \
  RadiativeZPrime RootInterface TauAnalysis

dist-hook:
	rm -rf `find $(distdir) -name '.svn' -or -name 'Makefile'`

all:  
	bash make_makefiles.sh
