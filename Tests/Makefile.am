AUTOMAKE_OPTIONS = -Wno-portability
AM_LDFLAGS += -module -avoid-version -rpath /dummy/path/not/used

EXTRA_DIST = Inputs

dist-hook:
	rm -rf $(distdir)/Inputs/.svn


EXTRA_LTLIBRARIES = LeptonTest.la GammaTest.la HadronTest.la DISTest.la

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

tests : tests-LEP test-Gamma

tests-LEP : test-LEP-VV test-LEP-VH test-LEP-VBF test-LEP-default test-LEP-BB test-LEP-Quarks test-LEP-Leptons

tests-DIS : test-DIS-Charged test-DIS-Neutral

tests-LHC : test-LHC-WW test-LHC-WZ test-LHC-ZZ test-LHC-ZGamma test-LHC-WGamma test-LHC-ZGamma \
	    test-LHC-ZH test-LHC-WH test-LHC-ZJet test-LHC-WJet test-LHC-Z test-LHC-W test-LHC-ZZVBF \
	    test-LHC-WWVBF test-LHC-bbH test-LHC-ttH test-LHC-GammaGamma test-LHC-GammaJet test-LHC-Higgs \
	    test-LHC-HiggsJet test-LHC-QCDFast test-LHC-QCD test-LHC-Top

tests-Gamma : test-Gamma-FF test-Gamma-WW test-Gamma-P

test-LEP-% : Inputs/LEP-%.in LeptonTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

test-Gamma-% : Inputs/Gamma-%.in GammaTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

test-DIS-% : Inputs/DIS-%.in DISTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

test-LHC-% : Inputs/LHC-%.in HadronTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<)) -N $${NUMEVENTS:-10000}

clean-local:
	rm -f *.out *.log *.tex *.top *.run *.dump *.mult *.Bmult
