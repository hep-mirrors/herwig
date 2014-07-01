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

Rivet-LEP: Rivet-LEP-22  Rivet-LEP-35  Rivet-LEP-44  Rivet-LEP-91  \
           Rivet-LEP-130 Rivet-LEP-133 Rivet-LEP-136 Rivet-LEP-161 \
           Rivet-LEP-172 Rivet-LEP-177 Rivet-LEP-183 Rivet-LEP-189 \
           Rivet-LEP-192 Rivet-LEP-196 Rivet-LEP-197 Rivet-LEP-200 \
           Rivet-LEP-202 Rivet-LEP-206 Rivet-LEP-14  Rivet-LEP-10\
           Rivet-LEP-91-nopi \
           Rivet-LEP-Powheg-14  Rivet-LEP-Powheg-22 \
           Rivet-LEP-Powheg-35  Rivet-LEP-Powheg-44  \
           Rivet-LEP-Powheg-91  Rivet-LEP-Powheg-130 \
           Rivet-LEP-Powheg-133 Rivet-LEP-Powheg-136 \
           Rivet-LEP-Powheg-161 Rivet-LEP-Powheg-172 \
           Rivet-LEP-Powheg-177 Rivet-LEP-Powheg-183 \
           Rivet-LEP-Powheg-189 Rivet-LEP-Powheg-192 \
           Rivet-LEP-Powheg-196 Rivet-LEP-Powheg-197 \
           Rivet-LEP-Powheg-200 Rivet-LEP-Powheg-202 \
           Rivet-LEP-Powheg-206 Rivet-LEP-Powheg-10 \
           Rivet-LEP-Powheg-91-nopi
	rm -rf Rivet-LEP
	python/merge-LEP LEP
	python/merge-LEP LEP-Powheg
	rivet-mkhtml -o Rivet-LEP LEP.yoda:Hw++ LEP-Powheg.yoda:Hw++-Powheg

Rivet-BFactory: Rivet-BFactory-10.52 Rivet-BFactory-Powheg-10.52 \
                Rivet-BFactory-10.52-sym Rivet-BFactory-Powheg-10.52-sym \
                Rivet-BFactory-10.54 Rivet-BFactory-Powheg-10.54 \
                Rivet-BFactory-10.45 Rivet-BFactory-Powheg-10.45 \
                Rivet-BFactory-10.58 Rivet-BFactory-Powheg-10.58 \
                Rivet-BFactory-Upsilon Rivet-BFactory-Upsilon2 Rivet-BFactory-Upsilon4 \
                Rivet-BFactory-Tau
	rm -rf Rivet-BFactory
	python/merge-BFactory BFactory
	python/merge-BFactory BFactory-Powheg
	rivet-mkhtml -o Rivet-BFactory BFactory.yoda:Hw++ BFactory-Powheg.yoda:Hw++-Powheg

Rivet-DIS: Rivet-DIS-e--LowQ2 \
	   Rivet-DIS-e+-LowQ2 Rivet-DIS-e+-HighQ2\
           Rivet-DIS-Powheg-e--LowQ2 \
	   Rivet-DIS-Powheg-e+-LowQ2 Rivet-DIS-Powheg-e+-HighQ2\
           Rivet-DIS-NoME-e--LowQ2 \
	   Rivet-DIS-NoME-e+-LowQ2 Rivet-DIS-NoME-e+-HighQ2
	rm -rf Rivet-DIS
	python/merge-DIS DIS 
	python/merge-DIS DIS-Powheg
	python/merge-DIS DIS-NoME
	rivet-mkhtml -o Rivet-DIS DIS.yoda:Hw++ DIS-Powheg.yoda:Hw++-Powheg DIS-NoME.yoda:Hw++-NoME

Rivet-TVT-WZ: Rivet-TVT-Run-I-Z     Rivet-TVT-Powheg-Run-I-Z \
              Rivet-TVT-Run-I-W     Rivet-TVT-Powheg-Run-I-W \
              Rivet-TVT-Run-I-WZ    Rivet-TVT-Powheg-Run-I-WZ\
              Rivet-TVT-Run-II-Z-e  Rivet-TVT-Powheg-Run-II-Z-e  \
              Rivet-TVT-Run-II-Z-mu Rivet-TVT-Powheg-Run-II-Z-mu \
              Rivet-TVT-Run-II-W    Rivet-TVT-Powheg-Run-II-W
	rm -rf Rivet-TVT-WZ
	python/merge-TVT-EW TVT-Run-II-W.yoda TVT-Run-II-Z-{e,mu}.yoda\
                            TVT-Run-I-{W,Z,WZ}.yoda -o TVT-WZ.yoda
	python/merge-TVT-EW TVT-Powheg-Run-II-W.yoda   TVT-Powheg-Run-II-Z-{e,mu}.yoda\
                            TVT-Powheg-Run-I-{W,Z,WZ}.yoda -o TVT-Powheg-WZ.yoda
	rivet-mkhtml -o Rivet-TVT-WZ TVT-WZ.yoda:Hw++ TVT-Powheg-WZ.yoda:Hw++-Powheg

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
		Rivet-LHC-7-Jets-13 Rivet-LHC-7-UE \
		Rivet-LHC-2360-UE Rivet-LHC-900-UE Rivet-LHC-2760-UE \
		Rivet-LHC-7-Heavy-0 \
		Rivet-LHC-7-Heavy-1 Rivet-LHC-7-Heavy-2 \
		Rivet-LHC-7-Heavy-3 Rivet-LHC-7-Heavy-4 \
		Rivet-LHC-7-Heavy-5 Rivet-LHC-7-Heavy-6 \
		Rivet-LHC-7-Heavy-7 Rivet-LHC-7-Heavy-8
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

Rivet-SppS: Rivet-ISR-44-UE Rivet-ISR-63-UE \
	    Rivet-SppS-63-UE  Rivet-ISR-30-UE \
	    Rivet-ISR-53-UE   Rivet-SppS-200-UE \
	    Rivet-SppS-500-UE Rivet-SppS-900-UE \
	    Rivet-SppS-546-UE Rivet-SppS-53-UE
	rm -rf Rivet-SppS
	python/merge-SppS SppS
	rivet-mkhtml -o Rivet-SppS SppS.yoda

Rivet-LHC-EW: Rivet-LHC-W-e    Rivet-LHC-Powheg-W-e  \
              Rivet-LHC-W-mu   Rivet-LHC-Powheg-W-mu \
              Rivet-LHC-Z-e    Rivet-LHC-Powheg-Z-e  \
              Rivet-LHC-Z-mu   Rivet-LHC-Powheg-Z-mu \
              Rivet-LHC-WZ     Rivet-LHC-Powheg-WZ\
              Rivet-LHC-WW-emu Rivet-LHC-Powheg-WW-emu\
              Rivet-LHC-WW-ll  Rivet-LHC-Powheg-WW-ll\
              Rivet-LHC-ZZ-ll Rivet-LHC-Powheg-ZZ-ll\
              Rivet-LHC-ZZ-lv Rivet-LHC-Powheg-ZZ-lv\
              Rivet-LHC-7-W-Jet-1-e Rivet-LHC-7-W-Jet-2-e\
              Rivet-LHC-7-W-Jet-3-e Rivet-LHC-7-Z-Jet-1-e\
              Rivet-LHC-7-Z-Jet-2-e Rivet-LHC-7-Z-Jet-3-e
	rm -rf Rivet-LHC-EW;
	python/merge-LHC-EW LHC-{W-e,W-mu,Z-e,Z-mu,WW-emu,WW-ll,WZ,ZZ-ll,ZZ-lv}.yoda LHC-7-{W,Z}-Jet-{1,2,3}-e.yoda -o LHC-EW.yoda;
	python/merge-LHC-EW LHC-Powheg-{W-e,W-mu,Z-e,Z-mu,WW-emu,WW-ll,WZ,ZZ-ll,ZZ-lv}.yoda -o LHC-Powheg-EW.yoda;
	rivet-mkhtml -o Rivet-LHC-EW LHC-EW.yoda:Hw++ LHC-Powheg-EW.yoda:Hw++-Powheg;

Rivet-LHC-Photon: Rivet-LHC-7-PromptPhoton Rivet-LHC-7-DiPhoton Rivet-LHC-7-GammaGamma
	rm -rf Rivet-LHC-Photon 
	python/merge-LHC-Photon LHC-7-PromptPhoton.yoda LHC-7-DiPhoton.yoda LHC-7-GammaGamma.yoda -o LHC-Photon.yoda
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
