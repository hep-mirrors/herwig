AUTOMAKE_OPTIONS = -Wno-portability
AM_LDFLAGS += -module -avoid-version -rpath /dummy/path/not/used

EXTRA_DIST = Inputs

dist-hook:
	rm -rf $(distdir)/Inputs/.svn

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
tests-LEP : test-LEP-VV test-LEP-VH test-LEP-VBF test-LEP-BB test-LEP-Quarks test-LEP-Leptons test-LEP-default test-LEP-TopDecay
else
tests-LEP : test-LEP-VV test-LEP-VH test-LEP-VBF test-LEP-BB test-LEP-Quarks test-LEP-Leptons
endif

tests-DIS : test-DIS-Charged test-DIS-Neutral

if WANT_LIBFASTJET
tests-LHC : test-LHC-WW test-LHC-WZ test-LHC-ZZ test-LHC-ZGamma test-LHC-WGamma test-LHC-ZGamma \
	    test-LHC-ZH test-LHC-WH test-LHC-ZJet test-LHC-WJet test-LHC-Z test-LHC-W test-LHC-ZZVBF test-LHC-VBF \
	    test-LHC-WWVBF test-LHC-bbH test-LHC-ttH test-LHC-GammaGamma test-LHC-GammaJet test-LHC-Higgs \
	    test-LHC-HiggsJet test-LHC-QCDFast test-LHC-QCD test-LHC-Top \
            test-LHC-WHJet test-LHC-ZHJet test-LHC-HJet test-LHC-ZShower test-LHC-WShower\
            test-LHC-WHJet-Powheg test-LHC-ZHJet-Powheg test-LHC-HJet-Powheg \
            test-LHC-ZShower-Powheg test-LHC-WShower-Powheg
else
tests-LHC : test-LHC-WW test-LHC-WZ test-LHC-ZZ test-LHC-ZGamma test-LHC-WGamma test-LHC-ZGamma \
	    test-LHC-ZH test-LHC-WH test-LHC-ZJet test-LHC-WJet test-LHC-Z test-LHC-W test-LHC-ZZVBF test-LHC-VBF \
	    test-LHC-WWVBF test-LHC-bbH test-LHC-ttH test-LHC-GammaGamma test-LHC-GammaJet test-LHC-Higgs \
	    test-LHC-HiggsJet test-LHC-QCDFast test-LHC-QCD test-LHC-Top
endif

tests-Anomalous : test-LHC-ZZ-Anomalous

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

Rivet-LHC-% : Rivet/LHC-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

Rivet-Star-% : Rivet/Star-%.in
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

Rivet-LEP: Rivet-LEP-35 Rivet-LEP-44  Rivet-LEP-91  Rivet-LEP-133 \
           Rivet-LEP-161 Rivet-LEP-172 Rivet-LEP-177 Rivet-LEP-183 \
           Rivet-LEP-189 Rivet-LEP-196 Rivet-LEP-197 Rivet-LEP-200 \
           Rivet-LEP-206 
	export DPATH=`rivet-config --datadir`; \
	for i in LEP-*.aida; do echo $$i; \
	rivet-rmgaps $$DPATH/ALEPH_2004_S5765862.aida $$i; done; \
	rivet-rmgaps $$DPATH/ALEPH_1996_S3486095 LEP-91.aida
	rm -rf Rivet-LEP
	rivet-mkhtml --no-mc-errs -o Rivet-LEP LEP-*.aida

Rivet-TVT-WZ: Rivet-TVT-Run-I-Z     Rivet-TVT-Run-I-Powheg-Z \
              Rivet-TVT-Run-I-W     Rivet-TVT-Run-I-Powheg-W \
              Rivet-TVT-Run-I-WZ     Rivet-TVT-Run-I-Powheg-WZ\
              Rivet-TVT-Run-II-Z-e  Rivet-TVT-Run-II-Powheg-Z-e  \
              Rivet-TVT-Run-II-Z-mu Rivet-TVT-Run-II-Powheg-Z-mu \
              Rivet-TVT-Run-II-W    Rivet-TVT-Run-II-Powheg-W
	export DPATH=`rivet-config --datadir`; \
	rivet-rmgaps $$DPATH/CDF_2008_S8095620 TVT-Run-II-Z-e.aida; \
	rivet-rmgaps $$DPATH/CDF_2008_S8095620 TVT-Run-II-Powheg-Z-e.aida; 
	rm -rf Rivet-TVT-WZ
	rivet-mkhtml --no-mc-errs -o Rivet-TVT-WZ TVT-*-{W,Z,WZ}.aida \
              TVT-Run-II-Z{-e,-mu}.aida TVT-Run-II-Powheg-Z{-e,-mu}.aida

Rivet-TVT-Photon: Rivet-TVT-Run-II-DiPhoton Rivet-TVT-Run-II-PromptPhoton
#                  Rivet-TVT-Run-I-PromptPhoton 
	rm -rf Rivet-TVT-Photon
	rivet-mkhtml --no-mc-errs -o Rivet-Photon TVT-*Photon.aida

Rivet-TVT-Jets: Rivet-TVT-Run-II-Jets-1 Rivet-TVT-Run-II-Jets-2 \
                Rivet-TVT-Run-II-Jets-3 Rivet-TVT-Run-II-Jets-4 \
                Rivet-TVT-Run-II-Jets-5 Rivet-TVT-Run-II-Jets-6 \
                Rivet-TVT-Run-II-Jets-7 Rivet-TVT-Run-II-Jets-8 \
                Rivet-TVT-Run-II-UE \
                Rivet-TVT-Run-I-Jets-1  Rivet-TVT-Run-I-Jets-2 \
                Rivet-TVT-Run-I-Jets-3  Rivet-TVT-Run-I-Jets-4 \
                Rivet-TVT-Run-I-Jets-5  \
                Rivet-TVT-630-UE Rivet-TVT-Run-I-UE
	export DPATH=`rivet-config --datadir`; \
	rivet-rmgaps $$DPATH/CDF_1994_S2952106 TVT-Run-I-Jets-4.aida	
	rm -rf Rivet-TVT-Jets
	rivet-mkhtml --no-mc-errs -o Rivet-TVT-Jets TVT-*-Jets-?.aida TVT-*-UE*.aida

Rivet-Star: Rivet-Star-UE     Rivet-Star-Jets-1 \
            Rivet-Star-Jets-2 Rivet-Star-Jets-3 \
            Rivet-Star-Jets-4 Rivet-Star-Jets-5
	rm -rf Rivet-Star
	export DPATH=`rivet-config --datadir`; \
	rivet-rmgaps $$DPATH/STAR_2006_S6860818 Star-UE.aida
	rivet-mkhtml --no-mc-errs -o Rivet-Star Star-Jets-?.aida Star-UE.aida

Rivet-LHC-EW: Rivet-LHC-W  Rivet-LHC-Powheg-W \
              Rivet-LHC-Z  Rivet-LHC-Powheg-Z \
              Rivet-LHC-WW Rivet-LHC-Powheg-WW\
              Rivet-LHC-ZZ Rivet-LHC-Powheg-ZZ
	rm -rf Rivet-LHC-EW
	rivet-mkhtml --no-mc-errs -o Rivet-LHC-EW LHC{,-Powheg}-{W,Z,WW,ZZ}.aida
 
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
	rm -f *.out *.log *.tex *.top *.run *.dump *.mult *.Bmult
