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
HWREAD = $(HERWIG) read -r $(REPO) -L $(builddir)/.libs
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

Rivet-LEP-% : Rivet/LEP-%.in
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

Rivet-LEP: Rivet-LEP-22  Rivet-LEP-35  Rivet-LEP-44  Rivet-LEP-91  \
           Rivet-LEP-133 Rivet-LEP-161 Rivet-LEP-172 Rivet-LEP-177 \
           Rivet-LEP-183 Rivet-LEP-189 Rivet-LEP-196 Rivet-LEP-197 \
           Rivet-LEP-200 Rivet-LEP-206 Rivet-LEP-14  Rivet-LEP-10\
           Rivet-LEP-Powheg-14  Rivet-LEP-Powheg-22 \
           Rivet-LEP-Powheg-35  Rivet-LEP-Powheg-44  \
           Rivet-LEP-Powheg-91  Rivet-LEP-Powheg-133 \
           Rivet-LEP-Powheg-161 Rivet-LEP-Powheg-172 \
           Rivet-LEP-Powheg-177 Rivet-LEP-Powheg-183 \
           Rivet-LEP-Powheg-189 Rivet-LEP-Powheg-196 \
           Rivet-LEP-Powheg-197 Rivet-LEP-Powheg-200 \
           Rivet-LEP-Powheg-206 Rivet-LEP-Powheg-10
	for i in LEP-*.aida; do rivet-rmgaps $$i; done;
	rivet-rmgaps LEP-91.aida;
	rivet-rmgaps LEP-Powheg-91.aida
	rm -rf Rivet-LEP
	python/merge-LEP LEP
	python/merge-LEP LEP-Powheg
	rivet-mkhtml -o Rivet-LEP LEP.aida:Hw++ LEP-Powheg.aida:Hw++-Powheg

Rivet-DIS: Rivet-DIS-e--LowQ2 \
	   Rivet-DIS-e+-LowQ2 Rivet-DIS-e+-HighQ2\
           Rivet-DIS-Powheg-e--LowQ2 \
	   Rivet-DIS-Powheg-e+-LowQ2 Rivet-DIS-Powheg-e+-HighQ2\
           Rivet-DIS-NoME-e--LowQ2 \
	   Rivet-DIS-NoME-e+-LowQ2 Rivet-DIS-NoME-e+-HighQ2
	rivet-rmgaps DIS-e+-LowQ2.aida
	rivet-rmgaps DIS-e--LowQ2.aida
	rivet-rmgaps DIS-e+-HighQ2.aida
	rivet-rmgaps DIS-Powheg-e+-LowQ2.aida
	rivet-rmgaps DIS-Powheg-e--LowQ2.aida
	rivet-rmgaps DIS-Powheg-e+-HighQ2.aida
	rivet-rmgaps DIS-NoME-e+-LowQ2.aida
	rivet-rmgaps DIS-NoME-e--LowQ2.aida
	rivet-rmgaps DIS-NoME-e+-HighQ2.aida
	python/merge-DIS DIS 
	python/merge-DIS DIS-Powheg
	python/merge-DIS DIS-NoME
	rivet-mkhtml -o Rivet-DIS DIS.aida:Hw++ DIS-Powheg.aida:Hw++-Powheg DIS-NoME.aida:Hw++-NoME

Rivet-TVT-WZ: Rivet-TVT-Run-I-Z     Rivet-TVT-Powheg-Run-I-Z \
              Rivet-TVT-Run-I-W     Rivet-TVT-Powheg-Run-I-W \
              Rivet-TVT-Run-I-WZ    Rivet-TVT-Powheg-Run-I-WZ\
              Rivet-TVT-Run-II-Z-e  Rivet-TVT-Powheg-Run-II-Z-e  \
              Rivet-TVT-Run-II-Z-mu Rivet-TVT-Powheg-Run-II-Z-mu \
              Rivet-TVT-Run-II-W    Rivet-TVT-Powheg-Run-II-W
	rivet-rmgaps TVT-Run-II-Z-e.aida;
	rivet-rmgaps TVT-Powheg-Run-II-Z-e.aida; 
	rm -rf Rivet-TVT-WZ
	python/merge-TVT-EW TVT-Run-II-W.aida TVT-Run-II-Z-{e,mu}.aida\
                            TVT-Run-I-{W,Z,WZ}.aida -o TVT-WZ.aida
	python/merge-TVT-EW TVT-Powheg-Run-II-W.aida   TVT-Powheg-Run-II-Z-{e,mu}.aida\
                            TVT-Powheg-Run-I-{W,Z,WZ}.aida -o TVT-Powheg-WZ.aida
	rivet-mkhtml -o Rivet-TVT-WZ TVT-WZ.aida:Hw++ TVT-Powheg-WZ.aida:Hw++-Powheg

Rivet-TVT-Photon: Rivet-TVT-Run-II-DiPhoton Rivet-TVT-Run-II-PromptPhoton
#                  Rivet-TVT-Run-I-PromptPhoton 
	rm -rf Rivet-TVT-Photon 
	python/merge-aida TVT-Run-II-DiPhoton.aida TVT-Run-II-PromptPhoton.aida\
                           -o TVT-Photon.aida
	rivet-mkhtml -o Rivet-TVT-Photon TVT-Photon.aida:Hw++


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
                Rivet-TVT-Run-I-UE\
                Rivet-TVT-630-UE        Rivet-TVT-630-Jets-1 \
                Rivet-TVT-630-Jets-2    Rivet-TVT-630-Jets-3
	rivet-rmgaps TVT-Run-I-Jets-4.aida	
	rm -rf Rivet-TVT-Jets
	python/merge-TVT-Jets TVT
	rivet-mkhtml -o Rivet-TVT-Jets TVT-Jets.aida:Hw++

Rivet-LHC-Jets: Rivet-LHC-7-Jets-0  Rivet-LHC-7-Jets-1  Rivet-LHC-7-Jets-2 \
                Rivet-LHC-7-Jets-3  Rivet-LHC-7-Jets-4 \
                Rivet-LHC-7-Jets-5  Rivet-LHC-7-Jets-6 \
                Rivet-LHC-7-Jets-7  Rivet-LHC-7-Jets-8 \
                Rivet-LHC-7-Jets-9  Rivet-LHC-7-Jets-10 \
                Rivet-LHC-7-Jets-11 Rivet-LHC-7-Jets-12 \
                Rivet-LHC-7-Jets-13 Rivet-LHC-7-UE \
	        Rivet-LHC-2360-UE Rivet-LHC-900-UE \
                Rivet-LHC-7-Heavy-1 Rivet-LHC-7-Heavy-2 \
                Rivet-LHC-7-Heavy-3 Rivet-LHC-7-Heavy-4 \
                Rivet-LHC-7-Heavy-5 Rivet-LHCX-7-Top \
	        Rivet-LHC-7-UE-Long Rivet-LHC-900-UE-Long
	rm -rf Rivet-LHC-Jets
	python/merge-LHC-Jets
	rivet-mkhtml -o Rivet-LHC-Jets LHC-Jets.aida:Hw++

Rivet-Star: Rivet-Star-UE     Rivet-Star-Jets-1 \
            Rivet-Star-Jets-2 Rivet-Star-Jets-3 \
            Rivet-Star-Jets-4
	rm -rf Rivet-Star
	rivet-rmgaps Star-UE.aida
	python/merge-Star Star
	rivet-mkhtml -o Rivet-Star Star.aida

Rivet-SppS: Rivet-ISR-44-UE Rivet-ISR-63-UE \
	    Rivet-SppS-63-UE  Rivet-ISR-30-UE \
	    Rivet-ISR-53-UE   Rivet-SppS-200-UE \
	    Rivet-SppS-500-UE Rivet-SppS-900-UE \
	    Rivet-SppS-546-UE Rivet-SppS-53-UE
	rm -rf Rivet-SppS
	python/merge-aida  ISR-44-UE.aida ISR-63-UE.aida SppS-63-UE.aida  ISR-30-UE.aida \
		           ISR-53-UE.aida SppS-200-UE.aida SppS-500-UE.aida SppS-900-UE.aida \
	                   SppS-546-UE.aida SppS-53-UE.aida -o SppS.aida
	rivet-mkhtml -o Rivet-SppS SppS.aida

Rivet-LHC-EW: Rivet-LHC-W-e  Rivet-LHC-Powheg-W-e  \
              Rivet-LHC-W-mu Rivet-LHC-Powheg-W-mu \
              Rivet-LHC-Z-e  Rivet-LHC-Powheg-Z-e  \
              Rivet-LHC-Z-mu Rivet-LHC-Powheg-Z-mu \
              Rivet-LHC-WW   Rivet-LHC-Powheg-WW\
              Rivet-LHC-ZZ   Rivet-LHC-Powheg-ZZ
	rm -rf Rivet-LHC-EW;
	python/merge-LHC-EW LHC-{W-e,W-mu,Z-e,Z-mu,WW,ZZ}.aida -o LHC-EW.aida;
	python/merge-LHC-EW LHC-Powheg-{W-e,W-mu,Z-e,Z-mu,WW,ZZ}.aida -o LHC-Powheg-EW.aida;
	rivet-mkhtml -o Rivet-LHC-EW LHC-EW.aida:Hw++ LHC-Powheg-EW.aida:Hw++-Powheg;

Rivet-LHC-Photon: Rivet-LHC-7-PromptPhoton Rivet-LHC-7-DiPhoton
	rm -rf Rivet-LHC-Photon 
	python/merge-aida LHC-7-PromptPhoton.aida LHC-7-DiPhoton.aida -o LHC-Photon.aida
	rivet-mkhtml -o Rivet-LHC-Photon LHC-Photon.aida:Hw++

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
	rm -f *.out *.log *.tex *.top *.run *.dump *.mult *.Bmult *.aida
