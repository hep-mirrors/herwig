AUTOMAKE_OPTIONS = -Wno-portability
AM_LDFLAGS += -module -avoid-version -rpath /dummy/path/not/used

EXTRA_DIST = Inputs python Rivet 

dist-hook:
	rm -rf $(distdir)/Inputs/.svn
	rm -rf $(distdir)/python/.svn
	rm -rf $(distdir)/Rivet/.svn

EXTRA_LTLIBRARIES = LeptonTest.la GammaTest.la HadronTest.la DISTest.la

if WANT_LIBFASTJET
EXTRA_LTLIBRARIES += HadronJetTest.la LeptonJetTest.la
HadronJetTest_la_SOURCES = \
Hadron/VHTest.h Hadron/VHTest.cc\
Hadron/VTest.h Hadron/VTest.cc\
Hadron/HTest.h Hadron/HTest.cc
HadronJetTest_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE) \
-I$(FASTJETPATH)
HadronJetTest_la_LIBADD = $(FASTJETLIBS) 
LeptonJetTest_la_SOURCES = \
Lepton/TopDecay.h Lepton/TopDecay.cc
LeptonJetTest_la_CPPFLAGS = $(AM_CPPFLAGS) $(FASTJETINCLUDE) \
-I$(FASTJETPATH)
LeptonJetTest_la_LIBADD = $(FASTJETLIBS) 
endif

LeptonTest_la_SOURCES = \
Lepton/VVTest.h Lepton/VVTest.cc \
Lepton/VBFTest.h Lepton/VBFTest.cc \
Lepton/VHTest.h Lepton/VHTest.cc \
Lepton/FermionTest.h Lepton/FermionTest.cc

GammaTest_la_SOURCES = \
Gamma/GammaMETest.h  Gamma/GammaMETest.cc \
Gamma/GammaPMETest.h Gamma/GammaPMETest.cc

DISTest_la_SOURCES = \
DIS/DISTest.h  DIS/DISTest.cc

HadronTest_la_SOURCES = \
Hadron/HadronVVTest.h  Hadron/HadronVVTest.cc\
Hadron/HadronVBFTest.h  Hadron/HadronVBFTest.cc\
Hadron/WHTest.h  Hadron/WHTest.cc\
Hadron/ZHTest.h  Hadron/ZHTest.cc\
Hadron/VGammaTest.h  Hadron/VGammaTest.cc\
Hadron/ZJetTest.h  Hadron/ZJetTest.cc\
Hadron/WJetTest.h  Hadron/WJetTest.cc\
Hadron/QQHTest.h  Hadron/QQHTest.cc


REPO = $(top_builddir)/src/HerwigDefaults.rpo
HERWIG = $(top_builddir)/src/Herwig++
HWREAD = $(HERWIG) read --repo $(REPO) -L $(builddir)/.libs
HWRUN = $(HERWIG) run

tests : tests-LEP tests-DIS tests-LHC tests-Gamma

if WANT_LIBFASTJET
tests-LEP : test-LEP-VV test-LEP-VH test-LEP-VBF test-LEP-BB test-LEP-Quarks test-LEP-Leptons \
            test-LEP-default test-LEP-Powheg test-LEP-TopDecay
else
tests-LEP : test-LEP-VV test-LEP-VH test-LEP-VBF test-LEP-BB test-LEP-Quarks test-LEP-Leptons
endif

tests-DIS : test-DIS-Charged test-DIS-Neutral

if WANT_LIBFASTJET
tests-LHC : test-LHC-WW test-LHC-WZ test-LHC-ZZ test-LHC-ZGamma test-LHC-WGamma \
	    test-LHC-ZH test-LHC-WH test-LHC-ZJet test-LHC-WJet test-LHC-Z test-LHC-W test-LHC-ZZVBF test-LHC-VBF \
	    test-LHC-WWVBF test-LHC-bbH test-LHC-ttH test-LHC-GammaGamma test-LHC-GammaJet test-LHC-Higgs \
	    test-LHC-HiggsJet test-LHC-QCDFast test-LHC-QCD test-LHC-Top test-LHC-Bottom \
            test-LHC-WHJet test-LHC-ZHJet test-LHC-HJet test-LHC-ZShower test-LHC-WShower\
            test-LHC-WHJet-Powheg test-LHC-ZHJet-Powheg test-LHC-HJet-Powheg \
            test-LHC-ZShower-Powheg test-LHC-WShower-Powheg
else
tests-LHC : test-LHC-WW test-LHC-WZ test-LHC-ZZ test-LHC-ZGamma test-LHC-WGamma \
	    test-LHC-ZH test-LHC-WH test-LHC-ZJet test-LHC-WJet test-LHC-Z test-LHC-W test-LHC-ZZVBF test-LHC-VBF \
	    test-LHC-WWVBF test-LHC-bbH test-LHC-ttH test-LHC-GammaGamma test-LHC-GammaJet test-LHC-Higgs \
	    test-LHC-HiggsJet test-LHC-QCDFast test-LHC-QCD test-LHC-Top
endif

tests-Gamma : test-Gamma-FF test-Gamma-WW test-Gamma-P

if WANT_LIBFASTJET 
test-LEP-% : Inputs/LEP-%.in LeptonTest.la LeptonJetTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}
else
test-LEP-% : Inputs/LEP-%.in LeptonTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}
endif

Rivet-LHC-Matchbox-% : Rivet/LHC-Matchbox-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet/LEP-Powheg-%.in : Rivet/LEP-Powheg.in
	infile=$(subst LEP-Powheg,,$(notdir $(subst .in,,$@))); \
	PROC=''; \
	if [ $${infile} == '-10' ]; then PROC='set /Herwig/MatrixElements/PowhegMEee2gZ2qq:MaximumFlavour 4'; fi; \
	sed s"!RUNNAME!LEP-Powheg$${infile}!" < Rivet/LEP-Powheg.in | sed s"!PARAMETERS!LEP/LEP$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/LEP-Matchbox-%.in : Rivet/LEP-Matchbox.in
	infile=$(subst LEP-Matchbox,,$(notdir $(subst .in,,$@))); \
	PROC="do Factory:Process e- e+ -> j j"; \
	if [ $${infile} == '-10' ]; then PROC="do Factory:Process e- e+ -> u ubar\ndo Factory:Process e- e+ -> d dbar\ndo Factory:Process e- e+ -> c cbar\ndo Factory:Process e- e+ -> s sbar"; fi; \
	sed s"!RUNNAME!LEP-Matchbox$${infile}!" < Rivet/LEP-Matchbox.in | sed s"!PARAMETERS!LEP/LEP$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/LEP-Dipole-%.in : Rivet/LEP-Dipole.in
	infile=$(subst LEP-Dipole,,$(notdir $(subst .in,,$@))); \
	PROC="do Factory:Process e- e+ -> j j"; \
	if [ $${infile} == '-10' ]; then PROC="do Factory:Process e- e+ -> u ubar\ndo Factory:Process e- e+ -> d dbar\ndo Factory:Process e- e+ -> c cbar\ndo Factory:Process e- e+ -> s sbar"; fi; \
	sed s"!RUNNAME!LEP-Dipole$${infile}!" < Rivet/LEP-Dipole.in | sed s"!PARAMETERS!LEP/LEP$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/LEP-%.in : Rivet/LEP.in
	infile=$(subst LEP,,$(notdir $(subst .in,,$@))); \
	PROC=''; \
	if [ $${infile} == '-10' ]; then PROC='set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4'; fi; \
	sed s"!RUNNAME!LEP$${infile}!" < Rivet/LEP.in | sed s"!PARAMETERS!LEP/LEP$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/DIS-Dipole-%.in : Rivet/DIS-Dipole.in
	infile=$(subst DIS-Dipole,,$(notdir $(subst .in,,$@))); \
	PROC="do /Herwig/MatrixElements/Dipole/DISFactory:Process e+ p e+ j"; \
	if [ $${infile} ~= 'e-' ]; then PROC="do /Herwig/MatrixElements/Dipole/DISFactory:Process e- p e- j"; fi; \
	sed s"!RUNNAME!DIS-Dipole$${infile}!" < Rivet/DIS-Dipole.in | sed s"!PARAMETERS!DIS/DIS$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/DIS-Matchbox-%.in : Rivet/DIS-Matchbox.in
	infile=$(subst DIS-Matchbox,,$(notdir $(subst .in,,$@))); \
	PROC="do /Herwig/MatrixElements/Matchbox/DISFactory:Process e+ p e+ j"; \
	if [ $${infile} ~= 'e-' ]; then PROC="do /Herwig/MatrixElements/Matchbox/DISFactory:Process e- p e- j"; fi; \
	sed s"!RUNNAME!DIS-Matchbox$${infile}!" < Rivet/DIS-Matchbox.in | sed s"!PARAMETERS!DIS/DIS$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/DIS-Powheg-%.in : Rivet/DIS-Powheg.in
	infile=$(subst DIS-Powheg,,$(notdir $(subst .in,,$@))); \
	PROC=''; \
	if [ $${infile} ~= '-NoME' ]; then PROC='set /Herwig/Shower/Evolver:MECorrMode 0'; fi; \
	sed s"!RUNNAME!DIS-Powheg$${infile}!" < Rivet/DIS-Powheg.in | sed s"!PARAMETERS!DIS/DIS$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/DIS-%.in : Rivet/DIS.in
	infile=$(subst DIS,,$(notdir $(subst .in,,$@))); \
	PROC=''; \
	if [ $${infile} ~= '-NoME' ]; then PROC='set /Herwig/Shower/Evolver:MECorrMode 0'; fi; \
	sed s"!RUNNAME!DIS$${infile}!" < Rivet/DIS.in | sed s"!PARAMETERS!DIS/DIS$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/BFactory-Dipole-%.in : Rivet/LEP-Dipole.in
	infile=$(subst BFactory-Dipole,,$(notdir $(subst .in,,$@))); \
	PROC="do Factory:Process e- e+ -> u ubar\ndo Factory:Process e- e+ -> d dbar\ndo Factory:Process e- e+ -> c cbar\ndo Factory:Process e- e+ -> s sbar"; \
	if [ $${infile} == '-10.58' ]; then PROC='do Factory:Process e- e+ -> u ubar\ndo Factory:Process e- e+ -> d dbar\ndo Factory:Process e- e+ -> c cbar\ndo Factory:Process e- e+ -> s sbar\ninsert /Herwig/Generators/EventGenerator:EventHandler:SubProcessHandlers 0 /Herwig/MatrixElements/SimpleEE\ncreate Herwig::MEee2VectorMeson /Herwig/MatrixElements/MEUpsilon HwMELepton.so\nset /Herwig/MatrixElements/MEUpsilon:VectorMeson /Herwig/Particles/Upsilon(4S)\nset /Herwig/MatrixElements/MEUpsilon:Coupling 0.0004151809\ninsert /Herwig/MatrixElements/SimpleEE:MatrixElements 0 /Herwig/MatrixElements/MEUpsilon'; elif [ $${infile} == '-Upsilon' ]; then PROC='erase /Herwig/MatrixElements/SimpleEE:MatrixElements 0'; fi; \
	sed s"!RUNNAME!BFactory-Dipole$${infile}!" < Rivet/LEP-Dipole.in | sed s"!PARAMETERS!BFactory/BFactory$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/BFactory-Matchbox-%.in : Rivet/LEP-Matchbox.in
	infile=$(subst BFactory-Matchbox,,$(notdir $(subst .in,,$@))); \
	PROC="do Factory:Process e- e+ -> u ubar\ndo Factory:Process e- e+ -> d dbar\ndo Factory:Process e- e+ -> c cbar\ndo Factory:Process e- e+ -> s sbar"; \
	if [ $${infile} == '-10.58' ]; then PROC='do Factory:Process e- e+ -> u ubar\ndo Factory:Process e- e+ -> d dbar\ndo Factory:Process e- e+ -> c cbar\ndo Factory:Process e- e+ -> s sbar\ninsert /Herwig/Generators/EventGenerator:EventHandler:SubProcessHandlers 0 /Herwig/MatrixElements/SimpleEE\ncreate Herwig::MEee2VectorMeson /Herwig/MatrixElements/MEUpsilon HwMELepton.so\nset /Herwig/MatrixElements/MEUpsilon:VectorMeson /Herwig/Particles/Upsilon(4S)\nset /Herwig/MatrixElements/MEUpsilon:Coupling 0.0004151809\ninsert /Herwig/MatrixElements/SimpleEE:MatrixElements 0 /Herwig/MatrixElements/MEUpsilon'; elif [ $${infile} == '-Upsilon' ]; then PROC='erase /Herwig/MatrixElements/SimpleEE:MatrixElements 0'; fi; \
	sed s"!RUNNAME!BFactory-Matchbox$${infile}!" < Rivet/LEP-Matchbox.in | sed s"!PARAMETERS!BFactory/BFactory$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/BFactory-Powheg-%.in : Rivet/LEP-Powheg.in
	infile=$(subst BFactory-Powheg,,$(notdir $(subst .in,,$@))); \
	PROC='set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4'; \
	if [ $${infile} == '-10.58' ]; then PROC='set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4\ncreate Herwig::MEee2VectorMeson /Herwig/MatrixElements/MEUpsilon HwMELepton.so\nset /Herwig/MatrixElements/MEUpsilon:VectorMeson /Herwig/Particles/Upsilon(4S)\nset /Herwig/MatrixElements/MEUpsilon:Coupling 0.0004151809\ninsert /Herwig/MatrixElements/SimpleEE:MatrixElements 0 /Herwig/MatrixElements/MEUpsilon'; elif [ $${infile} == '-Upsilon' ]; then PROC='erase /Herwig/MatrixElements/SimpleEE:MatrixElements 0'; fi; \
	sed s"!RUNNAME!BFactory-Powheg$${infile}!" < Rivet/LEP-Powheg.in | sed s"!PARAMETERS!BFactory/BFactory$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet/BFactory-%.in : Rivet/LEP.in
	infile=$(subst BFactory,,$(notdir $(subst .in,,$@))); \
	PROC='set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4'; \
	if [ $${infile} == '-10.58' ]; then PROC='set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4\ncreate Herwig::MEee2VectorMeson /Herwig/MatrixElements/MEUpsilon HwMELepton.so\nset /Herwig/MatrixElements/MEUpsilon:VectorMeson /Herwig/Particles/Upsilon(4S)\nset /Herwig/MatrixElements/MEUpsilon:Coupling 0.0004151809\ninsert /Herwig/MatrixElements/SimpleEE:MatrixElements 0 /Herwig/MatrixElements/MEUpsilon'; fi; \
	sed s"!RUNNAME!BFactory$${infile}!" < Rivet/LEP.in | sed s"!PARAMETERS!BFactory/BFactory$${infile}!" | sed s"!PROCESS!$${PROC}!" > $@

Rivet-LEP-Matchbox-% : Rivet/LEP-Matchbox-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet-LEP-Dipole-% : Rivet/LEP-Dipole-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet-BFactory-Matchbox-% : Rivet/BFactory-Matchbox-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet-LEP-% : Rivet/LEP-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

Rivet-BFactory-% : Rivet/BFactory-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

Rivet-TVT-% : Rivet/TVT-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

Rivet-DIS-% : Rivet/DIS-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}


Rivet-LHC-% : Rivet/LHC-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

Rivet-Star-% : Rivet/Star-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

Rivet-SppS-% : Rivet/SppS-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

Rivet-ISR-% : Rivet/ISR-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

Rivet-LEP: $(shell echo Rivet-LEP{,-Powheg,-Matchbox,-Dipole}-{10,14,22,35,44,91,130,133,136,161,172,177,183,189,192,196,197,200,202,206,91-nopi})
	rm -rf Rivet-LEP
	python/merge-LEP LEP
	python/merge-LEP LEP-Powheg
	python/merge-LEP LEP-Matchbox
	python/merge-LEP LEP-Dipole
	rivet-mkhtml -o Rivet-LEP LEP.yoda:Hw++ LEP-Powheg.yoda:Hw++-Powheg LEP-Matchbox.yoda:Hw++-Matchbox LEP-Dipole.yoda:Hw++-Dipole

Rivet-BFactory: $(shell echo Rivet-BFactory{,-Powheg,-Matchbox,-Dipole}-{10.52,10.52-sym,10.54,10.45,10.58}) \
                Rivet-BFactory-Upsilon Rivet-BFactory-Upsilon2 Rivet-BFactory-Upsilon4 \
                Rivet-BFactory-Tau
	rm -rf Rivet-BFactory
	python/merge-BFactory BFactory
	python/merge-BFactory BFactory-Powheg
	python/merge-BFactory BFactory-Matchbox
	python/merge-BFactory BFactory-Dipole
	rivet-mkhtml -o Rivet-BFactory BFactory.yoda:Hw++ BFactory-Powheg.yoda:Hw++-Powheg BFactory-Matchbox.yoda:Hw++-Matchbox BFactory-Dipole.yoda:Hw++-Dipole

Rivet-DIS: $(shell echo Rivet-DIS{,-NoME,-Powheg,-Matchbox,-Dipole}-{e--LowQ2,e+-LowQ2,e+-HighQ2})
	rm -rf Rivet-DIS
	python/merge-DIS DIS 
	python/merge-DIS DIS-Powheg
	python/merge-DIS DIS-NoME
	python/merge-DIS DIS-Matchbox
	python/merge-DIS DIS-Dipole
	rivet-mkhtml -o Rivet-DIS DIS.yoda:Hw++ DIS-Powheg.yoda:Hw++-Powheg DIS-NoME.yoda:Hw++-NoME DIS-Matchbox.yoda:Hw++-Matchbox DIS-Dipole.yoda:Hw++-Dipole

Rivet-TVT-WZ: Rivet-TVT-Run-I-Z     Rivet-TVT-Powheg-Run-I-Z     Rivet-TVT-Matchbox-Run-I-Z     \
              Rivet-TVT-Run-I-W     Rivet-TVT-Powheg-Run-I-W     Rivet-TVT-Matchbox-Run-I-W     \
              Rivet-TVT-Run-I-WZ    Rivet-TVT-Powheg-Run-I-WZ    Rivet-TVT-Matchbox-Run-I-WZ    \
              Rivet-TVT-Run-II-Z-e  Rivet-TVT-Powheg-Run-II-Z-e  Rivet-TVT-Matchbox-Run-II-Z-e  \
              Rivet-TVT-Run-II-Z-mu Rivet-TVT-Powheg-Run-II-Z-mu Rivet-TVT-Matchbox-Run-II-Z-mu \
              Rivet-TVT-Run-II-W    Rivet-TVT-Powheg-Run-II-W    Rivet-TVT-Matchbox-Run-II-W
	rm -rf Rivet-TVT-WZ
	python/merge-TVT-EW TVT-Run-II-W.yoda TVT-Run-II-Z-{e,mu}.yoda\
                            TVT-Run-I-{W,Z,WZ}.yoda -o TVT-WZ.yoda
	python/merge-TVT-EW TVT-Powheg-Run-II-W.yoda   TVT-Powheg-Run-II-Z-{e,mu}.yoda\
                            TVT-Powheg-Run-I-{W,Z,WZ}.yoda -o TVT-Powheg-WZ.yoda
	python/merge-TVT-EW TVT-Matchbox-Run-II-W.yoda   TVT-Matchbox-Run-II-Z-{e,mu}.yoda\
                            TVT-Matchbox-Run-I-{W,Z,WZ}.yoda -o TVT-Matchbox-WZ.yoda
	rivet-mkhtml -o Rivet-TVT-WZ TVT-WZ.yoda:Hw++ TVT-Powheg-WZ.yoda:Hw++-Powheg TVT-Matchbox-WZ.yoda:Hw++-Matchbox

Rivet-TVT-Photon: Rivet-TVT-Run-II-DiPhoton Rivet-TVT-Run-II-PromptPhoton
#                  Rivet-TVT-Run-I-PromptPhoton 
	rm -rf Rivet-TVT-Photon 
	yodamerge TVT-Run-II-DiPhoton.yoda TVT-Run-II-PromptPhoton.yoda\
                           -o TVT-Photon.yoda
	rivet-mkhtml -o Rivet-TVT-Photon TVT-Photon.yoda:Hw++


Rivet-TVT-Jets: Rivet-TVT-Run-II-Jets-1  Rivet-TVT-Run-II-Jets-2 \
                Rivet-TVT-Run-II-Jets-3  Rivet-TVT-Run-II-Jets-4 \
                Rivet-TVT-Run-II-Jets-5  Rivet-TVT-Run-II-Jets-6 \
                Rivet-TVT-Run-II-Jets-7  Rivet-TVT-Run-II-Jets-8 \
                Rivet-TVT-Run-II-Jets-9  Rivet-TVT-Run-II-Jets-10\
                Rivet-TVT-Run-II-Jets-11 Rivet-TVT-Run-II-UE \
                Rivet-TVT-Run-I-Jets-1   Rivet-TVT-Run-I-Jets-2 \
                Rivet-TVT-Run-I-Jets-3   Rivet-TVT-Run-I-Jets-4 \
                Rivet-TVT-Run-I-Jets-5   Rivet-TVT-Run-I-Jets-6 \
                Rivet-TVT-Run-I-Jets-7   Rivet-TVT-Run-I-Jets-8\
                Rivet-TVT-Run-I-UE       Rivet-TVT-Run-II-Jets-0\
                Rivet-TVT-630-UE         Rivet-TVT-630-Jets-1 \
                Rivet-TVT-630-Jets-2     Rivet-TVT-630-Jets-3 \
                Rivet-TVT-300-UE         Rivet-TVT-300-Jets-1 \
                Rivet-TVT-900-UE         Rivet-TVT-900-Jets-1
	python/merge-TVT-Energy TVT
	rivet-merge-CDF_2012_NOTE10874 TVT-300-Energy.yoda TVT-900-Energy.yoda TVT-1960-Energy.yoda
	flat2yoda RatioPlots.dat -o TVT-RatioPlots.yoda
	rm -rf Rivet-TVT-Jets
	python/merge-TVT-Jets TVT
	rivet-mkhtml -o Rivet-TVT-Jets TVT-Jets.yoda:Hw++

Rivet-LHC-Jets: Rivet-LHC-7-Jets-0  Rivet-LHC-7-Jets-1  Rivet-LHC-7-Jets-2 \
		Rivet-LHC-7-Jets-3  Rivet-LHC-7-Jets-4 \
		Rivet-LHC-7-Jets-5  Rivet-LHC-7-Jets-6 \
		Rivet-LHC-7-Jets-7  Rivet-LHC-7-Jets-8 \
		Rivet-LHC-7-Jets-9  Rivet-LHC-7-Jets-10 \
		Rivet-LHC-7-Jets-11 Rivet-LHC-7-Jets-12 \
		Rivet-LHC-7-Jets-13 Rivet-LHC-7-UE Rivet-LHC-8-UE \
		Rivet-LHC-2360-UE Rivet-LHC-900-UE Rivet-LHC-2760-UE \
		Rivet-LHC-7-Heavy-0 \
		Rivet-LHC-7-Heavy-1 Rivet-LHC-7-Heavy-2 \
		Rivet-LHC-7-Heavy-3 Rivet-LHC-7-Heavy-4 \
		Rivet-LHC-7-Heavy-5 Rivet-LHC-7-Heavy-6 \
		Rivet-LHC-7-Heavy-7 Rivet-LHC-7-Heavy-8 \
		Rivet-LHC-7-Top-L   Rivet-LHC-7-Top-SL \
		Rivet-LHC-7-UE-Long Rivet-LHC-900-UE-Long \
		Rivet-LHC-7-Jets-All-1 Rivet-LHC-7-Jets-All-2 \
		Rivet-LHC-7-Jets-All-3 Rivet-LHC-7-Jets-All-4 \
		Rivet-LHC-7-Jets-All-5
	rm -rf Rivet-LHC-Jets
	python/merge-LHC-Jets
	rivet-mkhtml -o Rivet-LHC-Jets LHC-Jets.yoda:Hw++

Rivet-Star: Rivet-Star-UE     Rivet-Star-Jets-1 \
            Rivet-Star-Jets-2 Rivet-Star-Jets-3 \
            Rivet-Star-Jets-4
	rm -rf Rivet-Star
	python/merge-Star Star
	rivet-mkhtml -o Rivet-Star Star.yoda

Rivet-SppS: Rivet-ISR-44-UE   Rivet-ISR-63-UE   \
	    Rivet-SppS-63-UE  Rivet-ISR-30-UE   \
	    Rivet-ISR-53-UE   Rivet-SppS-200-UE \
	    Rivet-SppS-500-UE Rivet-SppS-900-UE \
	    Rivet-SppS-546-UE Rivet-SppS-53-UE
	rm -rf Rivet-SppS
	python/merge-SppS SppS
	rivet-mkhtml -o Rivet-SppS SppS.yoda

Rivet-LHC-EW: Rivet-LHC-W-e    Rivet-LHC-Powheg-W-e    Rivet-LHC-Matchbox-W-e    \
              Rivet-LHC-W-mu   Rivet-LHC-Powheg-W-mu   Rivet-LHC-Matchbox-W-mu   \
              Rivet-LHC-Z-e    Rivet-LHC-Powheg-Z-e    Rivet-LHC-Matchbox-Z-e    \
              Rivet-LHC-Z-mu   Rivet-LHC-Powheg-Z-mu   Rivet-LHC-Matchbox-Z-mu   \
              Rivet-LHC-Z-LowMass-e    Rivet-LHC-Powheg-Z-LowMass-e    Rivet-LHC-Matchbox-Z-LowMass-e    \
              Rivet-LHC-Z-MedMass-e    Rivet-LHC-Powheg-Z-MedMass-e    Rivet-LHC-Matchbox-Z-MedMass-e    \
              Rivet-LHC-Z-LowMass-mu   Rivet-LHC-Powheg-Z-LowMass-mu   Rivet-LHC-Matchbox-Z-LowMass-mu   \
              Rivet-LHC-WZ     Rivet-LHC-Powheg-WZ     Rivet-LHC-Matchbox-WZ     \
              Rivet-LHC-WW-emu Rivet-LHC-Powheg-WW-emu Rivet-LHC-Matchbox-WW-emu \
              Rivet-LHC-WW-ll  Rivet-LHC-Powheg-WW-ll  Rivet-LHC-Matchbox-WW-ll  \
              Rivet-LHC-ZZ-ll  Rivet-LHC-Powheg-ZZ-ll  Rivet-LHC-Matchbox-ZZ-ll  \
              Rivet-LHC-ZZ-lv  Rivet-LHC-Powheg-ZZ-lv  Rivet-LHC-Matchbox-ZZ-lv  \
              Rivet-LHC-Matchbox-7-W-Jet-1-e Rivet-LHC-Matchbox-7-W-Jet-2-e \
              Rivet-LHC-Matchbox-7-W-Jet-3-e Rivet-LHC-Matchbox-7-Z-Jet-1-e \
              Rivet-LHC-Matchbox-7-Z-Jet-2-e Rivet-LHC-Matchbox-7-Z-Jet-3-e \
              Rivet-LHC-7-W-Jet-1-e Rivet-LHC-7-W-Jet-2-e \
              Rivet-LHC-7-W-Jet-3-e Rivet-LHC-7-Z-Jet-1-e \
              Rivet-LHC-7-Z-Jet-2-e Rivet-LHC-7-Z-Jet-3-e \
              Rivet-LHC-W-Z-e Rivet-LHC-W-Z-mu \
              Rivet-LHC-Powheg-W-Z-e Rivet-LHC-Powheg-W-Z-mu \
              Rivet-LHC-Matchbox-W-Z-e Rivet-LHC-Matchbox-W-Z-mu \
              Rivet-LHC-Matchbox-Z-b Rivet-LHC-Matchbox-Z-bb Rivet-LHC-Matchbox-W-b 
	rm -rf Rivet-LHC-EW;
	python/merge-LHC-EW LHC-{W-e,W-mu,Z-e,Z-mu,Z-LowMass-e,Z-MedMass-e,Z-LowMass-mu,W-Z-e,W-Z-mu,WW-emu,WW-ll,WZ,ZZ-ll,ZZ-lv}.yoda LHC-7-{W,Z}-Jet-{1,2,3}-e.yoda -o LHC-EW.yoda;
	python/merge-LHC-EW LHC-Matchbox-{W-e,W-mu,Z-e,Z-mu,Z-LowMass-e,Z-MedMass-e,Z-LowMass-mu,W-Z-e,W-Z-mu,WW-emu,WW-ll,WZ,ZZ-ll,ZZ-lv}.yoda LHC-Matchbox-7-{W,Z}-Jet-{1,2,3}-e.yoda -o LHC-Matchbox-EW.yoda;
	python/merge-LHC-EW LHC-Powheg-{W-e,W-mu,Z-e,Z-mu,Z-LowMass-e,Z-MedMass-e,Z-LowMass-mu,W-Z-e,W-Z-mu,WW-emu,WW-ll,WZ,ZZ-ll,ZZ-lv}.yoda -o LHC-Powheg-EW.yoda;
	rivet-mkhtml -o Rivet-LHC-EW LHC-EW.yoda:Hw++ LHC-Powheg-EW.yoda:Hw++-Powheg LHC-Matchbox-EW.yoda:Hw++-Matchbox LHC-Matchbox-Z-b.yoda:Hw++-Matchbox-Zb LHC-Matchbox-Z-bb.yoda:Hw++-Matchbox-Zbb LHC-Matchbox-W-b.yoda:Hw++-Matchbox-W-bb;

Rivet-LHC-Photon: Rivet-LHC-7-PromptPhoton-1  Rivet-LHC-7-PromptPhoton-2 Rivet-LHC-7-PromptPhoton-3 \
                  Rivet-LHC-7-PromptPhoton-4 Rivet-LHC-7-DiPhoton Rivet-LHC-7-GammaGamma
	rm -rf Rivet-LHC-Photon  
	python/merge-LHC-Photon -o LHC-Photon.yoda
	rivet-mkhtml -o Rivet-LHC-Photon LHC-Photon.yoda:Hw++

Rivet-LHC-Higgs: Rivet-LHC-Powheg-ggH Rivet-LHC-ggH Rivet-LHC-ggHJet \
		 Rivet-LHC-Powheg-VBF Rivet-LHC-VBF Rivet-LHC-WH Rivet-LHC-ZH \
		 Rivet-LHC-Powheg-WH Rivet-LHC-Powheg-ZH
	rm -rf Rivet-LHC-Higgs
	rivet-mkhtml -o Rivet-LHC-Higgs LHC-Powheg-ggH.yoda:gg-Powheg LHC-ggH.yoda:gg LHC-ggHJet.yoda:HJet \
                        LHC-Powheg-VBF.yoda:VBF-Powheg LHC-VBF.yoda:VBF LHC-WH.yoda:WH LHC-ZH.yoda:ZH \
                        LHC-Powheg-WH.yoda:WH-Powheg LHC-Powheg-ZH.yoda:ZH-Powheg

tests-Rivet : Rivet-LEP Rivet-BFactory Rivet-DIS Rivet-TVT-WZ Rivet-TVT-Photon Rivet-TVT-Jets Rivet-LHC-Jets Rivet-Star Rivet-SppS Rivet-LHC-EW Rivet-LHC-Photon


test-Gamma-% : Inputs/Gamma-%.in GammaTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

test-DIS-% : Inputs/DIS-%.in DISTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

if WANT_LIBFASTJET 
test-LHC-% : Inputs/LHC-%.in HadronTest.la GammaTest.la HadronJetTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}
else
test-LHC-% : Inputs/LHC-%.in HadronTest.la GammaTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}
endif

clean-local:
	rm -f *.out *.log *.tex *.top *.run *.dump *.mult *.Bmult *.yoda
