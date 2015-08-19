EXTRA_DIST = \
AcerDetInterface \
AlpGen \
Analysis2 \
AnomalousHVV \
DecayAnalysis \
HiggsPair \
HiggsPairOL \
LeptonME \
PGSInterface \
RadiativeZPrime \
TauAnalysis \
HJets++ \
MultiWeight

dist-hook:
	rm -rf `find $(distdir) -name '.svn' -or -name 'Makefile' -or -name '.hg'`
        chmod u+w $(distdir)/HJets++ $(distdir)/HJets++/autom4te.cache
	cd $(distdir)/HJets++; autoreconf -vi

all:  
	bash make_makefiles.sh
