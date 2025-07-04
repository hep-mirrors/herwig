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
HWBUILD = $(HERWIG) build --repo $(REPO) -L $(builddir)/.libs -i $(top_builddir)/src
HWINTEGRATE = $(HERWIG) integrate
HWRUN = $(HERWIG) run -N $${NUMEVENTS:-10000}


tests : tests-LEP tests-DIS tests-LHC tests-Gamma


LEPDEPS = \
test-LEP-VV \
test-LEP-VH \
test-LEP-VBF \
test-LEP-BB \
test-LEP-Quarks \
test-LEP-Leptons

if WANT_LIBFASTJET
LEPDEPS += test-LEP-TopDecay
endif

tests-LEP : $(LEPDEPS)

tests-DIS : test-DIS-Charged test-DIS-Neutral


LHCDEPS = \
test-LHC-WW test-LHC-WZ test-LHC-ZZ \
test-LHC-ZGamma test-LHC-WGamma \
test-LHC-ZH test-LHC-WH \
test-LHC-ZJet test-LHC-WJet \
test-LHC-Z test-LHC-W \
test-LHC-ZZVBF test-LHC-VBF \
test-LHC-WWVBF \
test-LHC-bbH test-LHC-ttH \
test-LHC-GammaGamma test-LHC-GammaJet \
test-LHC-Higgs test-LHC-HiggsJet \
test-LHC-QCDFast test-LHC-QCD \
test-LHC-Top


if WANT_LIBFASTJET
LHCDEPS += \
test-LHC-Bottom \
test-LHC-WHJet test-LHC-ZHJet test-LHC-HJet \
test-LHC-ZShower test-LHC-WShower \
test-LHC-WHJet-Powheg test-LHC-ZHJet-Powheg test-LHC-HJet-Powheg \
test-LHC-ZShower-Powheg test-LHC-WShower-Powheg
endif

tests-LHC : $(LHCDEPS)

tests-Gamma : test-Gamma-FF test-Gamma-WW test-Gamma-P



LEPLIBS = LeptonTest.la
HADLIBS = HadronTest.la

if WANT_LIBFASTJET
LEPLIBS += LeptonJetTest.la
HADLIBS += HadronJetTest.la
endif


test-LEP-% : Inputs/LEP-%.in $(LEPLIBS)
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<))

test-Gamma-% : Inputs/Gamma-%.in GammaTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<))

test-DIS-% : Inputs/DIS-%.in DISTest.la
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<))

test-LHC-% : Inputs/LHC-%.in GammaTest.la $(HADLIBS)
	$(HWREAD) $<
	$(HWRUN) $(notdir $(subst .in,.run,$<))



tests-Rivet : Rivet-EE Rivet-DIS Rivet-Fixed \
	      Rivet-TVT-EW Rivet-TVT-Photon Rivet-TVT-Jets \
	      Rivet-LHC-Jets Rivet-LHC-EW Rivet-LHC-Photon Rivet-LHC-Higgs

Rivet-%-UE-Cent.yoda : Rivet-%-UE-Cent.run Rivet-%-Cent.run
	rm -rf $(subst .yoda,,$@)
	mkdir $(subst .yoda,,$@)
	mv $(subst .yoda,.run,$@) $(subst UE-,,$(subst .yoda,.run,$@)) $(subst .yoda,,$@)
	cd $(subst .yoda,,$@); ../$(HWRUN) $(subst UE-,,$(subst .yoda,.run,$@)); mv $(subst UE-,,$@) ALICE_2015_CENT_PP.yoda
	cd $(subst .yoda,,$@); export RIVET_ANALYSIS_PATH=`pwd`; ../$(HWRUN)  $(subst .yoda,.run,$@)
	mv $(subst .yoda,,$@)/$@ .
	rm -rf $(subst .yoda,,$@)

Rivet-%.run : Rivet/%.in
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; $(HWBUILD) -c .cache/$(subst .run,,$@) $<

Rivet-Matchbox-%.yoda : Rivet-Matchbox-%.run
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; $(HWINTEGRATE) -c .cache/$(subst .run,,$<) $<
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; $(HWRUN)       -c .cache/$(subst .run,,$<) $<
	if [ -e $(subst .run,.hepmc,$<) ]; then gzip $(subst .run,.hepmc,$<); fi

Rivet-%.yoda : Rivet-%.run
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; $(HWRUN) $<
	if [ -e $(subst .run,.hepmc,$<) ]; then gzip $(subst .run,.hepmc,$<); fi

Rivet/%.in :
	python/make_input_files.py $(notdir $(subst .in,,$@))


Rivet-inputfiles: $(shell echo Rivet/EE{,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{12,12.8,13,14,14.8,17,21.5,22,22.5,25,26.8,27.6,27.7,29,30,30.2,30.3,30.5,30.7,30.8,31.2,31.3,31.6,34,34.5,34.8,35,36.2,41,41.5,42.1,42.6,43.5,43.6,44,45,48.0,50,52,53.3,55,55.3,56,57,58,59.5,60,60.8,61.4,65.4,66,75.7,76,82,85,91,91-nopi,93.0,130,130.1,133,136,136.1,161,161.3,172,172.3,177,182.8,183,188.6,189,192,194.4,196,197,200,200.2,202,205,206,206.2,207}.in) \
		  $(shell echo Rivet/EE-{183,189}-WW.in) \
	          $(shell echo Rivet/EE{,-Dipole}-{10.5,11.96,12.8,13.96,16.86,21.84,26.8,28.48,35.44,48.0,97.0}-gg.in) \
                  $(shell echo Rivet/EE{,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg}-{2.0,2.2,2.23,2.396,2.4,2.6,2.644,2.8,2.9,3.0,3.05,3.2,3.4,3.5,3.63,3.67,4.03,4.17,4.18,4.3,4.41,4.5,4.6,4.8,5.0,5.2,5.8,6.2,6.5,6.6,7.0,7.4,7.7,8.8,9.27,9.4,9.46,9.51,10,10.4,10.45,10.47,10.52,10.52-sym,10.54,10.55,10.58,10.58-sym,10.6}.in) \
                  $(shell echo Rivet/EE-{Bc,Ds,Eta,Etac,JPsi,JPsi-All,Lambdab,Lambdac,LL-10.58,LL-197,LL-29,LL-58,LL-91,Omega-Baryon,Omegac0,Omega-Meson,Phi,Psi2S,Psi2S-All,Psi3770,Tau,Upsilon,Upsilon2,Upsilon3,Upsilon4,Upsilon4-asym,Upsilon5,Xi0,Xic0,Xicp,Xim}.in) \
	          $(shell echo Rivet/DIS{,-NoME,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{225,251,300,318,318-CMS}-e+-{VeryLow,Low,Med,High}Q2.in) \
	          $(shell echo Rivet/DIS{,-NoME,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-318-{e+,e-}-CC-{VeryLow,Low,Med,High}Q2.in) \
		  $(shell echo Rivet/DIS{,-NoME,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{296,300,318}-e--{VeryLow,Low,Med,High}Q2.in) \
		  $(shell echo Rivet/Photo-{Direct,Resolved}-{300,318}-e+-Jets-{1..3}.in ) \
                  $(shell echo Rivet/TVT{,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{Run-I-Z,Run-I-W,Run-I-WZ,Run-II-Z-e,Run-II-Z-{,LowMass-,HighMass-}mu,Run-II-W}.in) \
	          $(shell echo Rivet/TVT{,-Dipole}-Run-II-{DiPhoton-GammaGamma,DiPhoton-GammaJet,PromptPhoton}.in) \
	          $(shell echo Rivet/TVT-Powheg-Run-II-{DiPhoton-GammaGamma,DiPhoton-GammaJet}.in) \
                  $(shell echo Rivet/TVT{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-Run-{I,II}-{Jets-{1..6},DiJets-{1..4}}.in ) \
	          $(shell echo Rivet/TVT{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{630-Jets-{1..3},300-Jets-1,900-Jets-1}.in ) \
                  $(shell echo Rivet/TVT{,-Dipole}-{Run-I,Run-II,300,630,900}-UE.in) \
                  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-DiJets-{1..7}-{A,B,C}.in ) \
                  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-13-DiJets-{{1..11}-A,{6..11}-B}.in ) \
                  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{7,8,13}-Jets-{0..10}.in ) \
	          $(shell echo Rivet/LHC{,-Dipole}-{900,2360,2760,7,8,13}-UE.in ) \
		  $(shell echo Rivet/LHC{,-Dipole}-2760-Jets-{1..3}.in ) \
	          $(shell echo Rivet/LHC{,-Dipole}-{900,7,13}-UE-Long.in ) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-Charm-{0..5}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{5,13}-Charm-0.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-Bottom-{0..9}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-13-Bottom-{0..6}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-Top-{L,SL}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{8,13}-Top-{All,L,SL}.in) \
                  $(shell echo Rivet/Star{,-Dipole}-{UE,Jets-{1..4}}.in ) \
	          $(shell echo Rivet/SppS{,-Dipole}-{53,63,200,500,546,900}-UE.in ) \
                  $(shell echo Rivet/LHC{,-Matchbox,-Matchbox-Powheg,-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{{,8-}W-{e,mu},13-W-mu,{2760-,,8-,13-}Z-{e,mu}-Mass{1..5},13-Z-nu,Z-mu-SOPHTY,WZ,WW-{emu,ll},13-WW-ll,13-ZZ-ll,ZZ-{ll,lv},{8,13}-WZ,8-ZZ-lv,8-WW-ll,Z-mu-Short}.in) \
                  $(shell echo Rivet/LHC{,-Dipole}-7-{W,Z}Gamma-{e,mu}.in) \
                  $(shell echo Rivet/LHC{,-Dipole}-8-ZGamma-{e,mu,nu}.in) \
	          $(shell echo Rivet/LHC{,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{7-W-Jet-{1..3}-e,7-Z-Jet-{0..3}-e,7-Z-Jet-0-mu}.in) \
	          $(shell echo Rivet/LHC{-Matchbox,-Matchbox-Powheg,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{Z-b,Z-bb,8-Z-b,8-Z-bb,13-Z-b,13-Z-bb,W-b,8-Z-jj}.in) \
		  $(shell echo Rivet/LHC{,-Dipole}-{7,8,13}-PromptPhoton-{1..5}.in) Rivet/LHC-GammaGamma-7.in \
	          $(shell echo Rivet/LHC{,-Powheg,-Dipole}-{7,8,13}-{DiPhoton-GammaGamma,DiPhoton-GammaJet}.in) \
	          $(shell echo Rivet/LHC{,-Powheg,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{ggH,VBF,WH,ZH}.in) \
                  $(shell echo Rivet/LHC{,-Powheg,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-8-{{ggH,VBF,WH,ZH}{,-GammaGamma},ggH-WW}.in) \
                  $(shell echo Rivet/LHC{,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-ggHJet.in) \
                  $(shell echo Rivet/ISR{,-Dipole}-{{30,44,53,62}-UE,{44,62}-Z-mu}.in Rivet/EHS{,-Dipole}{,-Pion}-UE.in) \
                  $(shell echo Rivet/SPS{,-Dipole}-{17.4-UE,200-Z-mu}.in ) \
                  $(shell echo Rivet/Fermilab{,-Dipole}-{27.4,38.8}-Z-mu.in ) \
	          $(shell echo Rivet/EE-Gamma-Direct-mumu-{161,172,183,189,196,206}.in ) \
                  $(shell echo Rivet/EE-Gamma-Direct-tautau-{189,196,206}.in ) \
                  $(shell echo Rivet/EE-Gamma-{pi0,Eta,EtaPrime,EtaC,ChiC2,ChiC0_2P,ChiC2_2P}-10.58.in ) \
                  $(shell echo Rivet/EE-Gamma-{Eta-{1,29},EtaPrime-{34.3,91.2}}.in ) \
	          $(shell echo Rivet/EE-Gamma-{Direct,Single-Resolved,Double-Resolved}-Jets-{198,206}.in )

#$(shell echo Rivet/LHC{,-Dipole}-{7,13}{,-UE}-Cent.in )

Analysis/RivetHerwig.so: Analysis/*.cc
	rivet-build Analysis/RivetHerwig.so Analysis/*.cc 

Rivet-GammaGamma: Rivet-GammaGamma/done
	touch $@

Rivet-GammaGamma/done: $(shell echo Rivet-GammaGamma-mumu-{3.5,4.5,5.5,6.5,7.5,9.0,12.5,17.5,30.0}.yoda )
	rm -rf Rivet-GammaGamma
	python/merge-GammaGamma GammaGamma
	rivet-mkhtml -o Rivet-GammaGamma GammaGamma.yoda:Hw
	touch $@

Rivet-EE-Gamma: Rivet-EE-Gamma/done
	touch $@

Rivet-EE-Gamma/done: $(shell echo Rivet-EE-Gamma-Direct-mumu-{161,172,183,189,196,206}.yoda ) \
	             $(shell echo Rivet-EE-Gamma-Direct-tautau-{189,196,206}.yoda ) \
	             $(shell echo Rivet-EE-Gamma-{pi0,Eta,EtaPrime,EtaC,ChiC2,ChiC0_2P,ChiC2_2P}-10.58.yoda ) \
                     $(shell echo Rivet-EE-Gamma-{Eta-{1,29},EtaPrime-{34.3,91.2}}.yoda ) \
	             $(shell echo Rivet-EE-Gamma-{Direct,Single-Resolved,Double-Resolved}-Jets-{196,198,206}.yoda )
	rm -rf Rivet-EE-Gamma
	python/merge-EE-Gamma EE-Gamma
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; rivet-mkhtml -o Rivet-EE-Gamma EE-Gamma.yoda:Hw
	touch $@

Rivet-EE : Rivet-EE/done
	touch $@

Rivet-EE/done : $(shell echo Rivet{,-Powheg}-EE-{12,12.8,13,14,14.8,17,21.5,22,22.5,25,26.8,27.6,27.7,29,30,30.2,30.3,30.5,30.7,30.8,31.2,31.3,31.6,34,34.5,34.8,35,36.2,41,41.5,42.1,42.6,43.5,43.6,44,45,48.0,50,52,53.3,55,55.3,56,57,58,59.5,60,60.8,61.4,65.4,66,75.7,76,82,85,91,91-nopi,93.0,130,130.1,133,136,136.1,161,161.3,172,172.3,177,182.8,183,188.6,189,192,194.4,196,197,200,200.2,202,205,206,206.2,207}.yoda) \
		$(shell echo Rivet-EE-{183,189}-WW.yoda) \
	  	$(shell echo Rivet-EE-{10.5,11.96,12.8,13.96,16.86,21.84,26.8,28.48,35.44,48.0,97.0}-gg.yoda) \
	  	$(shell echo Rivet-EE-{2.0,2.2,2.23,2.396,2.4,2.6,2.644,2.8,2.9,3.0,3.05,3.2,3.4,3.5,3.63,3.67,4.03,4.17,4.18,4.3,4.41,4.5,4.6,4.8,5.0,5.2,5.8,6.2,6.5,6.6,7.0,7.4,7.7,8.8,9.27,9.4,9.46,9.51,10,10.4,10.45,10.47,10.52,10.52-sym,10.54,10.55,10.58,10.58-sym,10.6}.yoda) \
	  	$(shell echo Rivet-EE-{Bc,Ds,Eta,Etac,JPsi,JPsi-All,Lambdab,Lambdac,LL-10.58,LL-197,LL-29,LL-58,LL-91,Omega-Baryon,Omegac0,Omega-Meson,Phi,Psi2S,Psi2S-All,Psi3770,Tau,Upsilon,Upsilon2,Upsilon3,Upsilon4,Upsilon4-asym,Upsilon5,Xi0,Xic0,Xicp,Xim}.yoda)
	rm -rf Rivet-EE
	python/merge-EE --with-gg --with-decay --with-WW EE
	python/merge-EE Powheg-EE
	rivet-mkhtml -o Rivet-EE EE.yoda:Hw Powheg-EE.yoda:Hw-Powheg
	python/plot-EE Rivet-EE
	touch $@

Onium-%.yoda:
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; $(HWBUILD) -c .cache/$(subst .yoda,,$@) Rivet/$(subst .yoda,.in,$@)
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; $(HWRUN)  $(subst .yoda,.run,$@)
	if [ -e $(subst .yoda,.hepmc,$@) ]; then gzip $(subst .yoda,.hepmc,$@); fi

Onium-Splitting: Analysis/RivetHerwig.so
	OUTPUT=`python/OniumSplitting.py --generate-input-files`; $(MAKE) $$OUTPUT NUMEVENTS=$${NUMEVENTS:-10000};
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; python/OniumSplitting.py --fast --analyse

Onium-Sigma: Analysis/RivetHerwig.so
	OUTPUT=`python/OniumSigma.py --generate-input-files | grep -v CT14 | grep -v LHAPDF | grep -v Eur.Phys`; $(MAKE) $$OUTPUT NUMEVENTS=$${NUMEVENTS:-10000};
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; python/OniumSigma.py --fast --analyse

Rivet-LowEnergy-%.yoda:
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; $(HWBUILD) -c .cache/$(subst .yoda,,$@) Rivet/$(subst .yoda,.in,$@)
	export RIVET_ANALYSIS_PATH=`pwd`/Analysis; $(HWRUN)  $(subst .yoda,.run,$@)
	if [ -e $(subst .yoda,.hepmc,$@) ]; then gzip $(subst .yoda,.hepmc,$@); fi

Rivet-LowEnergy-EE-%:
	args="--process "$(word 1,$(subst -, ,$(subst Rivet-LowEnergy-EE-,,$@))); if [ -n "$(strip $(word 2,$(subst -, ,$(subst Rivet-LowEnergy-EE-,,$@))))" ]; then args+=" --flavour "$(word 2,$(subst -, ,$(subst Rivet-LowEnergy-EE-,,$@))); fi; OUTPUT=`python/LowEnergy-EE.py $$args --non-perturbative --perturbative`; $(MAKE) $$OUTPUT NUMEVENTS=$${NUMEVENTS:-10000};
	args="--process "$(word 1,$(subst -, ,$(subst Rivet-LowEnergy-EE-,,$@))); plots=`python/LowEnergy-EE.py $$args --plots`;  plots2=`echo $$plots | sed s"A-eA-mA"g | sed s"A-fA-mA"g`; python/mergeLowEnergy.py $(subst Rivet-LowEnergy-,,$@) $$plots; if [ -e LowEnergy-NonPerturbative-EE-$(subst Rivet-LowEnergy-EE-,,$@).yoda ] && [ -e LowEnergy-Perturbative-EE-$(subst Rivet-LowEnergy-EE-,,$@).yoda ]; then rivet-mkhtml -o Rivet-LowEnergy-EE-$(subst Rivet-LowEnergy-EE-,,$@) LowEnergy-NonPerturbative-EE-$(subst Rivet-LowEnergy-EE-,,$@).yoda:"Non-Pert" LowEnergy-Perturbative-EE-$(subst Rivet-LowEnergy-EE-,,$@).yoda:"Pert" $$plots2; elif [ -e LowEnergy-NonPerturbative-EE-$(subst Rivet-LowEnergy-EE-,,$@).yoda ]; then rivet-mkhtml -o Rivet-LowEnergy-EE-$(subst Rivet-LowEnergy-EE-,,$@) LowEnergy-NonPerturbative-EE-$(subst Rivet-LowEnergy-EE-,,$@).yoda:"Non-Pert" $$plots2; elif [ -e  LowEnergy-Perturbative-EE-$(subst Rivet-LowEnergy-EE-,,$@).yoda ]; then  rivet-mkhtml -o Rivet-LowEnergy-EE-$(subst Rivet-LowEnergy-EE-,,$@)  LowEnergy-Perturbative-EE-$(subst Rivet-LowEnergy-EE-,,$@).yoda:"Pert" $$plots2; fi

Rivet-LowEnergy-Photon-%:
	args="--process "$(word 1,$(subst -, ,$(subst Rivet-LowEnergy-Photon-,,$@))); if [ -n "$(strip $(word 2,$(subst -, ,$(subst Rivet-LowEnergy-Photon-,,$@))))" ]; then args+=" --flavour "$(word 2,$(subst -, ,$(subst Rivet-LowEnergy-Photon-,,$@))); fi; OUTPUT=`python/LowEnergy-Photon.py $$args --non-perturbative --perturbative --resonance`; $(MAKE) $$OUTPUT NUMEVENTS=$${NUMEVENTS:-10000};
	args="--process "$(word 1,$(subst -, ,$(subst Rivet-LowEnergy-Photon-,,$@))); plots=`python/LowEnergy-Photon.py $$args --plots`; plots2=`echo $$plots | sed s"A-eA-mA"g | sed s"A-fA-mA"g`; yodaArgs=`python/mergeLowEnergy.py $(subst Rivet-LowEnergy-Photon-,GammaGamma-,$@) $$plots`; rivet-mkhtml -o Rivet-LowEnergy-Photon-$(subst Rivet-LowEnergy-Photon-,,$@) $$yodaArgs $$plots2;

Rivet-R:
	OUTPUT=`python/R.py --perturbative --non-perturbative`; $(MAKE) $$OUTPUT NUMEVENTS=$${NUMEVENTS:-10000};
	plots=`python/R.py --perturbative --non-perturbative --plots`; python/mergeLowEnergy.py R $$plots; rivet-mkhtml -o Rivet-R LowEnergy-Perturbative-R.yoda:"Pert" LowEnergy-NonPerturbative-R.yoda:"Non-Pert" $$plots

Rivet-DIS : Rivet-DIS/done
	touch $@

Rivet-DIS/done: $(shell echo Rivet{-DIS,-NoME-DIS,-Powheg-DIS}-{225,251,300,318,318-CMS}-e+-{VeryLow,Low,Med,High}Q2.yoda) \
		$(shell echo Rivet{-DIS,-NoME-DIS,-Powheg-DIS}-{296,300,318}-e--{VeryLow,Low,Med,High}Q2.yoda) \
                $(shell echo Rivet{-DIS,-NoME-DIS,-Powheg-DIS}-318-{e+,e-}-CC-{VeryLow,Low,Med,High}Q2.yoda)\
                $(shell echo Rivet-Photo-{Direct,Resolved}-{300,318}-e+-Jets-{1..3}.yoda)
	rm -rf Rivet-DIS
	python/merge-DIS DIS
	python/merge-DIS Powheg-DIS
	python/merge-DIS NoME-DIS
	rivet-mkhtml -o Rivet-DIS DIS.yoda:Hw Powheg-DIS.yoda:Hw-Powheg NoME-DIS.yoda:Hw-NoME
	touch $@

Rivet-TVT-EW : Rivet-TVT-EW/done
	touch $@

Rivet-TVT-EW/done:  $(shell echo Rivet{,-Powheg}-TVT-{Run-I-Z,Run-I-W,Run-I-WZ,Run-II-Z-{e,{,LowMass-,HighMass-}mu},Run-II-W}.yoda)
	rm -rf Rivet-TVT-EW
	python/merge-TVT-EW TVT
	python/merge-TVT-EW Powheg-TVT
	rivet-mkhtml -o Rivet-TVT-EW TVT-EW.yoda:Hw Powheg-TVT-EW.yoda:Hw-Powheg
	touch $@

Rivet-TVT-Photon : Rivet-TVT-Photon/done
	touch $@

Rivet-TVT-Photon/done: $(shell echo Rivet{,-Powheg}-TVT-Run-II-{DiPhoton-GammaGamma,DiPhoton-GammaJet}.yoda Rivet-TVT-Run-II-PromptPhoton.yoda)
	rm -rf Rivet-TVT-Photon
	python/merge-TVT-Photon TVT
	python/merge-TVT-Photon Powheg-TVT
	rivet-mkhtml -o Rivet-TVT-Photon TVT-Photon.yoda:Hw Powheg-TVT-Photon.yoda:Hw-Powheg
	touch $@



Rivet-TVT-Jets : Rivet-TVT-Jets/done
	touch $@

Rivet-TVT-Jets/done: $(shell echo Rivet-TVT-Run-{I,II}-{Jets-{1..6},DiJets-{1..4}}.yoda ) \
	             $(shell echo Rivet-TVT-{630-Jets-{1..3},300-Jets-1,900-Jets-1}.yoda ) \
                     $(shell echo Rivet-TVT-{Run-I,Run-II,300,630,900}-UE.yoda)
	rm -rf Rivet-TVT-Jets
	python/merge-TVT-Jets TVT
	rivet-mkhtml -o Rivet-TVT-Jets TVT-Jets.yoda:Hw
	touch $@

Rivet-Fixed : Rivet-Fixed/done
	touch $@


Rivet-Fixed/done : $(shell echo Rivet-SppS-{53,63,200,500,546,900}-UE.yoda ) \
	          $(shell echo Rivet-ISR-{{30,44,53,62}-UE,{44,62}-Z-mu}.yoda Rivet-EHS{,-Pion}-UE.yoda ) \
                  $(shell echo Rivet-Star-{UE,Jets-{1..4}}.yoda ) \
                  $(shell echo Rivet-SPS-{17.4-UE,200-Z-mu}.yoda ) \
                  $(shell echo Rivet-Fermilab-{27.4,38.8}-Z-mu.yoda )
	rm -rf Rivet-Fixed
	python/merge-Fixed Fixed
	rivet-mkhtml -o Rivet-Fixed Fixed.yoda:Hw
	touch $@

Rivet-LHC-Jets : Rivet-LHC-Jets/done
	touch $@

Rivet-LHC-Jets/done : \
	        $(shell echo Rivet-LHC-7-DiJets-{1..7}-{A,B,C}.yoda   ) \
	        $(shell echo Rivet-LHC-13-DiJets-{{1..11}-A,{6..11}-B}.yoda   ) \
	        $(shell echo Rivet-LHC-{7,8,13}-Jets-{0..10}.yoda     ) \
	        $(shell echo Rivet-LHC-2760-Jets-{1..3}.yoda     ) \
	        $(shell echo Rivet-LHC-{900,2360,2760,7,8,13}-UE.yoda ) \
	        $(shell echo Rivet-LHC-{900,7,13}-UE-Long.yoda        ) \
		$(shell echo Rivet-LHC-7-Charm-{0..5}.yoda            ) \
		$(shell echo Rivet-LHC-{5,13}-Charm-0.yoda            ) \
		$(shell echo Rivet-LHC-7-Bottom-{0..9}.yoda           ) \
		$(shell echo Rivet-LHC-13-Bottom-{0..6}.yoda           ) \
		$(shell echo Rivet-LHC-{7,8,13}-Top-{L,SL}.yoda ) \
		$(shell echo Rivet-LHC-{8,13}-Top-All.yoda )
#	        $(shell echo Rivet-LHC-{7,13}-UE-Cent.yoda )
	rm -rf Rivet-LHC-Jets
	python/merge-LHC-Jets LHC
	rivet-mkhtml -o Rivet-LHC-Jets LHC-Jets.yoda:Hw
	touch $@


Rivet-LHC-EW : Rivet-LHC-EW/done
	touch $@

Rivet-LHC-EW/done: \
		$(shell echo Rivet{,-Powheg}-LHC-{{,8-}W-{e,mu},13-W-mu,{2760-,,8-,13-}Z-{e,mu}-Mass{1..5},13-Z-nu,Z-mu-SOPHTY,WZ,WW-{emu,ll},13-WW-ll,ZZ-{ll,lv},{8,13}-WZ,13-ZZ-ll,8-ZZ-lv,8-WW-ll,Z-mu-Short}.yoda) \
		$(shell echo Rivet-LHC-{7-W-Jet-{1..3}-e,7-Z-Jet-{0..3}-e,7-Z-Jet-0-mu}.yoda) \
		$(shell echo Rivet-LHC-7-{W,Z}Gamma-{e,mu}.yoda) \
		$(shell echo Rivet-LHC-8-ZGamma-{e,mu,nu}.yoda)
	rm -rf Rivet-LHC-EW;
	python/merge-LHC-EW LHC
	python/merge-LHC-EW Powheg-LHC
	rivet-mkhtml -o Rivet-LHC-EW LHC-EW.yoda:Hw Powheg-LHC-EW.yoda:Hw-Powheg \
                                     Rivet-LHC-Z-mu-SOPHTY.yoda:Hw Rivet-Powheg-LHC-Z-mu-SOPHTY.yoda:Hw-Powheg
	touch $@




Rivet-LHC-Photon : Rivet-LHC-Photon/done
	touch $@

Rivet-LHC-Photon/done: \
		$(shell echo Rivet-LHC-{7,8,13}-PromptPhoton-{1..5}.yoda) \
		Rivet-LHC-GammaGamma-7.yoda \
	    $(shell echo Rivet{,-Powheg}-LHC-{7,8,13}-{DiPhoton-GammaGamma,DiPhoton-GammaJet}.yoda)
	rm -rf Rivet-LHC-Photon
	python/merge-LHC-Photon LHC
	python/merge-LHC-Photon Powheg-LHC
	rivet-mkhtml -o Rivet-LHC-Photon LHC-Photon.yoda:Hw Powheg-LHC-Photon.yoda:Hw-Powheg
	touch $@




Rivet-LHC-Higgs : Rivet-LHC-Higgs/done
	touch $@

Rivet-LHC-Higgs/done:  \
		$(shell echo Rivet{,-Powheg}-LHC-{ggH,VBF,WH,ZH}.yoda) \
        $(shell echo Rivet{,-Powheg}-LHC-8-{{ggH,VBF,WH,ZH}{,-GammaGamma},ggH-WW}.yoda) \
        Rivet-LHC-ggHJet.yoda
	rivet-merge Rivet-Powheg-LHC-8-{ggH{-GammaGamma,-WW,},{VBF,ZH,WH}{,-GammaGamma}}.yoda -o Powheg-LHC-Higgs.yoda
	rivet-merge Rivet-LHC-8-{ggH{-GammaGamma,-WW,},{VBF,ZH,WH}{,-GammaGamma}}.yoda -o LHC-Higgs.yoda
	rm -rf Rivet-LHC-Higgs
	rivet-mkhtml -o Rivet-LHC-Higgs Powheg-LHC-Higgs.yoda:Hw-Powheg LHC-Higgs.yoda:Hw\
	                Rivet-Powheg-LHC-ggH.yoda:gg-Powheg Rivet-LHC-ggH.yoda:gg Rivet-LHC-ggHJet.yoda:HJet \
                        Rivet-Powheg-LHC-VBF.yoda:VBF-Powheg Rivet-LHC-VBF.yoda:VBF Rivet-LHC-WH.yoda:WH Rivet-LHC-ZH.yoda:ZH \
                        Rivet-Powheg-LHC-WH.yoda:WH-Powheg Rivet-Powheg-LHC-ZH.yoda:ZH-Powheg
	touch $@





clean-local:
	rm -f *.out *.log *.tex *.top *.run *.dump *.mult *.Bmult *.yoda Rivet/*.in anatohepmc.txt hepmctoana.txt
	rm -rf Rivet-*
	rm -rf Onium-Splitting
	rm -rf Analysis/*.so

distclean-local:
	rm -rf .cache
