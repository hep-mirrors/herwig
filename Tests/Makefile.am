AM_LDFLAGS += -module -avoid-version -rpath /dummy/path/not/used

EXTRA_DIST = Inputs python Rivet 

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
HERWIG = $(top_builddir)/src/Herwig
HWREAD = $(HERWIG) read --repo $(REPO) -L $(builddir)/.libs -i $(top_builddir)/src
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
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib -i ../$(top_builddir)/src ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet-TVT-Matchbox-% : Rivet/TVT-Matchbox-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib -i ../$(top_builddir)/src ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet-TVT-Dipole-% : Rivet/TVT-Dipole-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib -i ../$(top_builddir)/src ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet-LHC-Dipole-% : Rivet/LHC-Dipole-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib -i ../$(top_builddir)/src ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet/LEP-%.in :
	python/make_input_files.py $(notdir $(subst .in,,$@))

Rivet/DIS-%.in : 
	python/make_input_files.py $(notdir $(subst .in,,$@))

Rivet/BFactory-%.in:
	python/make_input_files.py $(notdir $(subst .in,,$@))

Rivet/TVT-%.in:
	python/make_input_files.py $(notdir $(subst .in,,$@))

Rivet/LHC-%.in:
	python/make_input_files.py $(notdir $(subst .in,,$@))

Rivet/Star-%.in:
	python/make_input_files.py $(notdir $(subst .in,,$@))

Rivet/SppS-%.in:
	python/make_input_files.py $(notdir $(subst .in,,$@))

Rivet/ISR-%.in:
	python/make_input_files.py $(notdir $(subst .in,,$@))

Rivet-LEP-Matchbox-% : Rivet/LEP-Matchbox-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib -i ../$(top_builddir)/src ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet-LEP-Dipole-% : Rivet/LEP-Dipole-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib -i ../$(top_builddir)/src ../$<; \
	../$(HERWIG) run $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}; \
	mv $(notdir $(subst .in,.yoda,$<)) ..; \
	cd ..

Rivet-BFactory-Matchbox-% : Rivet/BFactory-Matchbox-%.in
	if [ ! -d Rivet-$(notdir $(subst .in,,$<))  ]; then mkdir Rivet-$(notdir $(subst .in,,$<)); fi;
	cd Rivet-$(notdir $(subst .in,,$<)); echo `pwd`; \
	../$(HERWIG) read --repo ../$(REPO) -L ../$(top_builddir)/lib -i ../$(top_builddir)/src ../$<; \
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

Rivet-inputfiles: $(shell echo Rivet/LEP{,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{9.4,12,13,17,27.6,29,30.2,30.7,30.75,30,31.3,34.8,43.6,50,52,55,56,57,60.8,60,61.4,10,12.8,22,26.8,35,44,48.0,91,93.0,130,133,136,161,172,177,183,189,192,196,197,200,202,206,91-nopi}.in) \
                  $(shell echo Rivet/LEP{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Powheg,-Matchbox-Powheg}-14.in) \
	          $(shell echo Rivet/LEP{,-Dipole}-{10.5,11.96,12.8,13.96,16.86,21.84,26.8,28.48,35.44,48.0,97.0}-gg.in) \
                  $(shell echo Rivet/BFactory{,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg}-{10.52,10.52-sym,10.54,10.45}.in) \
                  $(shell echo Rivet/BFactory{,-Dipole}-{Upsilon,Upsilon2,Upsilon4,Tau,10.58-res}.in) \
                  $(shell echo Rivet/DIS{,-NoME,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{e--LowQ2,e+-LowQ2,e+-HighQ2}.in) \
                  $(shell echo Rivet/TVT{,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{Run-I-Z,Run-I-W,Run-I-WZ,Run-II-Z-e,Run-II-Z-{,LowMass-,HighMass-}mu,Run-II-W}.in) \
	          $(shell echo Rivet/TVT{,-Dipole}-Run-II-{DiPhoton-GammaGamma,DiPhoton-GammaJet,PromptPhoton}.in) \
	          $(shell echo Rivet/TVT-Powheg-Run-II-{DiPhoton-GammaGamma,DiPhoton-GammaJet}.in) \
                  $(shell echo Rivet/TVT{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{Run-II-Jets-{0..11},Run-I-Jets-{1..8}}.in ) \
	          $(shell echo Rivet/TVT{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{630-Jets-{1..3},300-Jets-1,900-Jets-1}.in ) \
                  $(shell echo Rivet/TVT{,-Dipole}-{Run-I,Run-II,300,630,900}-UE.in) \
                  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-DiJets-{1..7}-{A,B,C}.in ) \
                  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{7,8,13}-Jets-{0..10}.in ) \
	          $(shell echo Rivet/LHC{,-Dipole}-{900,2360,2760,7,8,13}-UE.in ) \
	          $(shell echo Rivet/LHC{,-Dipole}-{900,7,13}-UE-Long.in ) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-Charm-{1..5}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-Bottom-{0..8}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-Top-{L,SL}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{7,8,13}-Top-All.in) \
                  $(shell echo Rivet/Star{,-Dipole}-{UE,Jets-{1..4}}.in ) \
	          $(shell echo Rivet/SppS{,-Dipole}-{200,500,900,546}-UE.in ) \
                  $(shell echo Rivet/LHC{,-Matchbox,-Matchbox-Powheg,-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{W-{e,mu},13-Z-{e,mu},8-Z-Mass{1..4}-{e,mu},Z-{e,mu,mu-SOPHTY},Z-LowMass-{e,mu},Z-MedMass-e,WZ,WW-{emu,ll},ZZ-{ll,lv},8-ZZ-lv,W-Z-{e,mu}}.in) \
                  $(shell echo Rivet/LHC{,-Dipole}-7-{W,Z}Gamma-{e,mu}.in) \
	          $(shell echo Rivet/LHC{,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{7-W-Jet-{1..3}-e,7-Z-Jet-{0..3}-e,7-Z-Jet-0-mu}.in) \
	          $(shell echo Rivet/LHC{-Matchbox,-Matchbox-Powheg,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{Z-b,Z-bb,W-b,8-Z-jj}.in) \
		  $(shell echo Rivet/LHC{,-Dipole}-{7,8}-PromptPhoton-{1..4}.in) Rivet/LHC-GammaGamma-7.in \
	          $(shell echo Rivet/LHC{,-Dipole,-Powheg}-7-{DiPhoton-GammaGamma,DiPhoton-GammaJet}.in) \
	          $(shell echo Rivet/LHC{,-Powheg,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{ggH,VBF,WH,ZH}.in) \
                  $(shell echo Rivet/LHC{,-Powheg,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-8-{ggH,VBF,WH,ZH}{,-GammaGamma}.in) \
                  $(shell echo Rivet/LHC{,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-ggHJet.in)
#                  $(shell echo Rivet/ISR-{30,44,53,62}-UE.in ) $(shell echo Rivet/SppS-{53,63}-UE.in )

Rivet-LEP: $(shell echo Rivet-LEP{,-Powheg,-Matchbox,-Dipole}-{9.4,12,13,17,27.6,29,30.2,30.7,30.75,30,31.3,34.8,43.6,50,52,55,56,57,60.8,60,61.4,10,12.8,22,26.8,35,44,48.0,91,93.0,130,133,136,161,172,177,183,189,192,196,197,200,202,206,91-nopi}) \
	   $(shell echo Rivet-LEP-{10.5,11.96,12.8,13.96,16.86,21.84,26.8,28.48,35.44,48.0,97.0}-gg) \
           $(shell echo Rivet-LEP{,-Powheg,-Matchbox-Powheg}-14.in)	
	rm -rf Rivet-LEP
	python/merge-LEP --with-gg LEP
	python/merge-LEP LEP-Powheg
	python/merge-LEP LEP-Matchbox
	python/merge-LEP LEP-Dipole
	rivet-mkhtml -o Rivet-LEP LEP.yoda:Hw LEP-Powheg.yoda:Hw-Powheg LEP-Matchbox.yoda:Hw-Matchbox LEP-Dipole.yoda:Hw-Dipole

Rivet-BFactory: $(shell echo Rivet-BFactory{,-Powheg,-Matchbox,-Dipole}-{10.52,10.52-sym,10.54,10.45}) \
                $(shell echo Rivet-BFactory-{Upsilon,Upsilon2,Upsilon4,Tau,10.58-res,10.58})
	rm -rf Rivet-BFactory
	python/merge-BFactory BFactory
	python/merge-BFactory BFactory-Powheg
	python/merge-BFactory BFactory-Matchbox
	python/merge-BFactory BFactory-Dipole
	rivet-mkhtml -o Rivet-BFactory BFactory.yoda:Hw BFactory-Powheg.yoda:Hw-Powheg BFactory-Matchbox.yoda:Hw-Matchbox BFactory-Dipole.yoda:Hw-Dipole

Rivet-DIS: $(shell echo Rivet-DIS{,-NoME,-Powheg,-Matchbox,-Dipole}-{e--LowQ2,e+-LowQ2,e+-HighQ2})
	rm -rf Rivet-DIS
	python/merge-DIS DIS 
	python/merge-DIS DIS-Powheg
	python/merge-DIS DIS-NoME
	python/merge-DIS DIS-Matchbox
	python/merge-DIS DIS-Dipole
	rivet-mkhtml -o Rivet-DIS DIS.yoda:Hw DIS-Powheg.yoda:Hw-Powheg DIS-NoME.yoda:Hw-NoME DIS-Matchbox.yoda:Hw-Matchbox DIS-Dipole.yoda:Hw-Dipole

Rivet-TVT-WZ:  $(shell echo Rivet-TVT{,-Powheg,-Matchbox,-Dipole}-{Run-I-Z,Run-I-W,Run-I-WZ,Run-II-Z-{e,{,LowMass-,HighMass-}mu},Run-II-W})
	rm -rf Rivet-TVT-WZ
	python/merge-TVT-EW TVT-Run-II-W.yoda TVT-Run-II-Z-{e,{,LowMass-,HighMass-}mu}.yoda\
                            TVT-Run-I-{W,Z,WZ}.yoda -o TVT-WZ.yoda
	python/merge-TVT-EW TVT-Powheg-Run-II-W.yoda   TVT-Powheg-Run-II-Z-{e,{,LowMass-,HighMass-}mu}.yoda\
                            TVT-Powheg-Run-I-{W,Z,WZ}.yoda -o TVT-Powheg-WZ.yoda
	python/merge-TVT-EW TVT-Matchbox-Run-II-W.yoda   TVT-Matchbox-Run-II-Z-{e,{,LowMass-,HighMass-}mu}.yoda\
                            TVT-Matchbox-Run-I-{W,Z,WZ}.yoda -o TVT-Matchbox-WZ.yoda
	python/merge-TVT-EW TVT-Dipole-Run-II-W.yoda   TVT-Dipole-Run-II-Z-{e,{,LowMass-,HighMass-}mu}.yoda\
                            TVT-Dipole-Run-I-{W,Z,WZ}.yoda -o TVT-Dipole-WZ.yoda
	rivet-mkhtml -o Rivet-TVT-WZ TVT-WZ.yoda:Hw TVT-Powheg-WZ.yoda:Hw-Powheg TVT-Matchbox-WZ.yoda:Hw-Matchbox TVT-Dipole-WZ.yoda:Hw-Dipole

Rivet-TVT-Photon: $(shell echo Rivet-TVT-Run-II-{DiPhoton-GammaGamma,DiPhoton-GammaJet,PromptPhoton}) \
                  $(shell echo Rivet-TVT-Powheg-Run-II-{DiPhoton-GammaGamma,DiPhoton-GammaJet})
	rm -rf Rivet-TVT-Photon
	python/merge-TVT-Photon TVT -o TVT-Photon.yoda
	python/merge-TVT-Photon TVT-Powheg -o TVT-Powheg-Photon.yoda
	rivet-mkhtml -o Rivet-TVT-Photon TVT-Photon.yoda:Hw TVT-Powheg-Photon.yoda:Hw-Powheg

Rivet-TVT-Jets:  $(shell echo Rivet-TVT-{Run-II-Jets-{0..11},Run-I-Jets-{1..8}} ) \
	         $(shell echo Rivet-TVT-{630-Jets-{1..3},300-Jets-1,900-Jets-1} ) \
                 $(shell echo Rivet-TVT-{Run-I,Run-II,300,630,900}-UE)
	rm -rf Rivet-TVT-Jets
	python/merge-TVT-Jets TVT
	rivet-mkhtml -o Rivet-TVT-Jets TVT-Jets.yoda:Hw

#python/merge-TVT-Energy TVT
#rivet-merge-CDF_2012_NOTE10874 TVT-300-Energy.yoda TVT-900-Energy.yoda TVT-1960-Energy.yoda
#flat2yoda RatioPlots.dat -o TVT-RatioPlots.yoda

Rivet-LHC-Jets: $(shell echo Rivet-LHC-7-DiJets-{1..7}-{A,B,C} ) \
	        $(shell echo Rivet-LHC-{7,8,13}-Jets-{0..10} ) \
	        $(shell echo Rivet-LHC-{900,2360,2760,7,8,13}-UE ) \
	        $(shell echo Rivet-LHC-{900,7,13}-UE-Long ) \
		$(shell echo Rivet-LHC-7-Charm-{1..5}) \
		$(shell echo Rivet-LHC-7-Bottom-{0..8}) \
		$(shell echo Rivet-LHC-7-Top-{L,SL})\
		$(shell echo Rivet-LHC-{7,8,13}-Top-All)
	rm -rf Rivet-LHC-Jets
	python/merge-LHC-Jets LHC
	rivet-mkhtml -o Rivet-LHC-Jets LHC-Jets.yoda:Hw

Rivet-Star: $(shell echo Rivet-Star-{UE,Jets-{1..4}} ) 
	rm -rf Rivet-Star
	python/merge-Star Star
	rivet-mkhtml -o Rivet-Star Star.yoda

Rivet-SppS: $(shell echo Rivet-ISR-{30,44,53,62}-UE ) \
	    $(shell echo Rivet-SppS-{53,63,200,500,900,546}-UE )
	rm -rf Rivet-SppS
	python/merge-SppS SppS
	rivet-mkhtml -o Rivet-SppS SppS.yoda

Rivet-LHC-EW: $(shell echo Rivet-LHC{,-Matchbox,-Powheg,-Dipole}-{13-Z-{e,mu},8-Z-Mass{1..4}-{e,mu},W-{e,mu},Z-{e,mu,mu-SOPHTY},Z-LowMass-{e,mu},Z-MedMass-e,WZ,WW-{emu,ll},ZZ-{ll,lv},8-ZZ-lv,W-Z-{e,mu}}) \
	      $(shell echo Rivet-LHC{,-Matchbox,-Dipole}-{7-W-Jet-{1..3}-e,7-Z-Jet-{0..3}-e,7-Z-Jet-0-mu}) \
	      $(shell echo Rivet-LHC{-Matchbox,-Dipole}-{Z-b,Z-bb,W-b,8-Z-jj}) \
              $(shell echo Rivet-LHC-7-{W,Z}Gamma-{e,mu}) \
	rm -rf Rivet-LHC-EW;
	python/merge-LHC-EW LHC-{13-Z-{e,mu},8-Z-Mass{1..4}-{e,mu},W-{e,mu},Z-{e,mu,mu-Short},Z-LowMass-{e,mu},Z-MedMass-e,W-Z-{e,mu},WW-{emu,ll},WZ,ZZ-{ll,lv},8-ZZ-lv}.yoda LHC-7-{W,Z}-Jet-{1,2,3}-e.yoda LHC-7-{W,Z}Gamma-{e,mu}.yoda -o LHC-EW.yoda;
	python/merge-LHC-EW LHC-Matchbox-{13-Z-{e,mu},8-Z-Mass{1..4}-{e,mu},W-{e,mu},Z-{e,mu,mu-Short},Z-LowMass-{e,mu},Z-MedMass-e,W-Z-{e,mu},WW-{emu,ll},WZ,ZZ-{ll,lv},8-ZZ-lv}.yoda LHC-Matchbox-7-{W,Z}-Jet-{1,2,3}-e.yoda -o LHC-Matchbox-EW.yoda;
	python/merge-LHC-EW LHC-Dipole-{13-Z-{e,mu},8-Z-Mass{1..4}-{e,mu},W-{e,mu},Z-{e,mu,mu-Short},Z-LowMass-{e,mu},Z-MedMass-e,W-Z-{e,mu},WW-{emu,ll},WZ,ZZ-{ll,lv},8-ZZ-lv}.yoda LHC-Dipole-7-{W,Z}-Jet-{1,2,3}-e.yoda -o LHC-Dipole-EW.yoda;
	python/merge-LHC-EW LHC-Powheg-{W-{e,mu},Z-{e,mu},Z-LowMass-{e,mu},Z-MedMass-e,W-Z-{e,mu},WW-{emu,ll},WZ,ZZ-{ll,lv},8-ZZ-lv}.yoda -o LHC-Powheg-EW.yoda;
	rivet-mkhtml -o Rivet-LHC-EW LHC-EW.yoda:Hw LHC-Powheg-EW.yoda:Hw-Powheg LHC-Matchbox-EW.yoda:Hw-Matchbox LHC-Matchbox-Z-b.yoda:Hw-Matchbox-Zb \
                                     LHC-Matchbox-Z-bb.yoda:Hw-Matchbox-Zbb LHC-Matchbox-W-b.yoda:Hw-Matchbox-W-bb LHC-Dipole-EW.yoda:Hw-Dipole \
                                     LHC-Dipole-Z-b.yoda:Hw-Dipole-Zb LHC-Dipole-Z-bb.yoda:Hw-Dipole-Zbb LHC-Dipole-W-b.yoda:Hw-Dipole-W-bb \
                                     LHC-Z-mu-SOPHTY.yoda:Hw LHC-Powheg-Z-mu-SOPHTY.yoda:Hw-Powheg LHC-Matchbox-Z-mu-SOPHTY.yoda:Hw-Matchbox
Rivet-LHC-Photon: $(shell echo Rivet-LHC-{7,8}-PromptPhoton-{1..4}) Rivet-LHC-GammaGamma-7 \
	          $(shell echo Rivet-LHC{,-Powheg}-7-{DiPhoton-GammaGamma,DiPhoton-GammaJet})
	rm -rf Rivet-LHC-Photon
	python/merge-LHC-Photon LHC -o LHC-Photon.yoda
	python/merge-LHC-Photon LHC-Powheg -o LHC-Powheg-Photon.yoda
	rivet-mkhtml -o Rivet-LHC-Photon LHC-Photon.yoda:Hw LHC-Powheg-Photon.yoda:Hw-Powheg

Rivet-LHC-Higgs:  $(shell echo Rivet-LHC{,-Powheg}-{ggH,VBF,WH,ZH})\
                  $(shell echo Rivet-LHC{,-Powheg}-8-{ggH,VBF,WH,ZH}{,-GammaGamma}) Rivet-LHC-ggHJet
	rm -rf Rivet-LHC-Higgs
	rivet-mkhtml -o Rivet-LHC-Higgs LHC-Powheg-ggH.yoda:gg-Powheg LHC-ggH.yoda:gg LHC-ggHJet.yoda:HJet \
                        LHC-Powheg-VBF.yoda:VBF-Powheg LHC-VBF.yoda:VBF LHC-WH.yoda:WH LHC-ZH.yoda:ZH \
                        LHC-Powheg-WH.yoda:WH-Powheg LHC-Powheg-ZH.yoda:ZH-Powheg

tests-Rivet : Rivet-LEP Rivet-BFactory Rivet-DIS Rivet-TVT-WZ Rivet-TVT-Photon Rivet-TVT-Jets Rivet-LHC-Jets Rivet-Star Rivet-SppS Rivet-LHC-EW Rivet-LHC-Photon Rivet-LHC-Higgs


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
