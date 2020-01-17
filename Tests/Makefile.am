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



tests-Rivet : Rivet-EE Rivet-DIS Rivet-Star Rivet-SppS \
	      Rivet-TVT-WZ Rivet-TVT-Photon Rivet-TVT-Jets \
	      Rivet-LHC-Jets Rivet-LHC-EW Rivet-LHC-Photon Rivet-LHC-Higgs

Rivet-%.run : Rivet/%.in
	$(HWBUILD) -c .cache/$(subst .run,,$@) $<

Rivet-Matchbox-%.yoda : Rivet-Matchbox-%.run
	$(HWINTEGRATE) -c .cache/$(subst .run,,$<) $<
	$(HWRUN)       -c .cache/$(subst .run,,$<) $<

Rivet-%.yoda : Rivet-%.run
	$(HWRUN) $<

Rivet/%.in :
	python/make_input_files.py $(notdir $(subst .in,,$@))

Rivet-inputfiles: $(shell echo Rivet/EE{,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{7.7,9.4,12,13,17,27.6,29,30.2,30.7,30,31.3,34,34.8,41,42.1,42.6,43.6,45,50,52,53.3,55,56,57,58,59.5,60.8,60,61.4,66,76,82,85,10,12.8,21.5,22,25,26.8,34.5,35,36.2,44,48.0,91,93.0,130,133,136,161,172,177,183,189,192,196,197,200,202,205,206,207,91-nopi}.in) \
                  $(shell echo Rivet/EE{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Powheg,-Matchbox-Powheg}-{14,14.8}.in) \
	          $(shell echo Rivet/EE{,-Dipole}-{10.5,11.96,12.8,13.96,16.86,21.84,26.8,28.48,35.44,48.0,97.0}-gg.in) \
                  $(shell echo Rivet/EE{,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg}-{2.2,2.6,3.0,3.2,4.6,4.8,3.63,4.03,4.5,9.46,10.00,10.52,10.52-sym,10.54,10.45,10.47}.in) \
                  $(shell echo Rivet/EE-{Upsilon,Upsilon2,Upsilon4,JPsi,Psi2S,Tau,Phi,10.58-res,10.58,Lambdac,Omega-Meson,Omega-Baryon,Eta,Xi0,Xic0,Xim}.in) \
                  $(shell echo Rivet/DIS{,-NoME,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{e--LowQ2,e+-LowQ2,e+-HighQ2}.in) \
                  $(shell echo Rivet/TVT{,-Powheg,-Matchbox,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox-Powheg,-Merging}-{Run-I-Z,Run-I-W,Run-I-WZ,Run-II-Z-e,Run-II-Z-{,LowMass-,HighMass-}mu,Run-II-W}.in) \
	          $(shell echo Rivet/TVT{,-Dipole}-Run-II-{DiPhoton-GammaGamma,DiPhoton-GammaJet,PromptPhoton}.in) \
	          $(shell echo Rivet/TVT-Powheg-Run-II-{DiPhoton-GammaGamma,DiPhoton-GammaJet}.in) \
                  $(shell echo Rivet/TVT{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{Run-II-Jets-{0..11},Run-I-Jets-{1..8}}.in ) \
	          $(shell echo Rivet/TVT{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{630-Jets-{1..3},300-Jets-1,900-Jets-1}.in ) \
                  $(shell echo Rivet/TVT{,-Dipole}-{Run-I,Run-II,300,630,900}-UE.in) \
                  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-DiJets-{1..7}-{A,B,C}.in ) \
                  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-13-DiJets-{{1..11}-A,{6..11}-B}.in ) \
                  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{7,8,13}-Jets-{0..10}.in ) \
	          $(shell echo Rivet/LHC{,-Dipole}-{900,2360,2760,7,8,13}-UE.in ) \
		  $(shell echo Rivet/LHC{,-Dipole}-2760-Jets-{1..3}.in ) \
	          $(shell echo Rivet/LHC{,-Dipole}-{900,7,13}-UE-Long.in ) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-Charm-{1..5}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-Bottom-{0..9}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-7-Top-{L,SL}.in) \
		  $(shell echo Rivet/LHC{,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Matchbox,-Matchbox-Powheg,-Merging}-{8,13}-Top-{All,L,SL}.in) \
                  $(shell echo Rivet/Star{,-Dipole}-{UE,Jets-{1..4}}.in ) \
	          $(shell echo Rivet/SppS{,-Dipole}-{53,63,200,500,546,900}-UE.in ) \
                  $(shell echo Rivet/LHC{,-Matchbox,-Matchbox-Powheg,-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{W-{e,mu},13-Z-{e,mu},Z-HighMass{1,2}-e,{8,13}-W-mu,8-Z-Mass{1..4}-{e,mu},Z-{e,mu,mu-SOPHTY},Z-LowMass-{e,mu},Z-MedMass-e,WZ,WW-{emu,ll},ZZ-{ll,lv},{8,13}-WZ,8-ZZ-lv,8-WW-ll,Z-mu-Short}.in) \
                  $(shell echo Rivet/LHC{,-Dipole}-7-{W,Z}Gamma-{e,mu}.in) \
                  $(shell echo Rivet/LHC{,-Dipole}-8-ZGamma-{e,mu}.in) \
	          $(shell echo Rivet/LHC{,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{7-W-Jet-{1..3}-e,7-Z-Jet-{0..3}-e,7-Z-Jet-0-mu}.in) \
	          $(shell echo Rivet/LHC{-Matchbox,-Matchbox-Powheg,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{Z-b,Z-bb,8-Z-b,8-Z-bb,W-b,8-Z-jj}.in) \
		  $(shell echo Rivet/LHC{,-Dipole}-{7,8,13}-PromptPhoton-{1..4}.in) Rivet/LHC-GammaGamma-7.in \
	          $(shell echo Rivet/LHC{,-Powheg,-Dipole}-{7,8}-{DiPhoton-GammaGamma,DiPhoton-GammaJet}.in) \
	          $(shell echo Rivet/LHC{,-Powheg,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-{ggH,VBF,WH,ZH}.in) \
                  $(shell echo Rivet/LHC{,-Powheg,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-8-{{ggH,VBF,WH,ZH}{,-GammaGamma},ggH-WW}.in) \
                  $(shell echo Rivet/LHC{,-Matchbox,-Matchbox-Powheg,-Dipole,-Dipole-MCatNLO,-Dipole-Matchbox-Powheg,-Merging}-ggHJet.in) \
                  $(shell echo Rivet/ISR{,-Dipole}-{30,44,53,62}-UE.in Rivet/EHS{,-Dipole}-UE.in)


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
	              $(shell echo Rivet-EE-Gamma-{Direct,Single-Resolved,Double-Resolved}-Jets-{198,206}.yoda )
	rm -rf Rivet-EE-Gamma
	python/merge-EE-Gamma EE-Gamma
	rivet-mkhtml -o Rivet-EE-Gamma EE-Gamma.yoda:Hw
	touch $@

Rivet-EE : Rivet-EE/done
	touch $@

Rivet-EE/done : $(shell echo Rivet{,-Powheg}-EE-{7.7,9.4,12,13,14,14.8,17,27.6,29,30.2,30.7,30,31.3,34,34.8,43.6,45,50,52,53.3,55,56,57,58,59.5,60.8,60,61.4,66,76,10,12.8,21.5,22,25,26.8,34.5,35,36.2,41,42.1,42.6,44,48.0,82,85,91,93.0,130,133,136,161,172,177,183,189,192,196,197,200,202,205,206,207,91-nopi}.yoda) \
	   $(shell echo Rivet-EE-{10.5,11.96,12.8,13.96,16.86,21.84,26.8,28.48,35.44,48.0,97.0}-gg.yoda) \
	   $(shell echo Rivet-EE-{10.52,10.52-sym,2.2,2.6,3.0,3.2,4.6,4.8,3.63,4.03,4.5,9.46,10.00,10.54,10.45,10.47,Upsilon,Upsilon2,Upsilon4,Tau,Phi,10.58-res,10.58,Lambdac,Omega-Meson,Omega-Baryon,Eta,Xi0,Xic0,Xim,JPsi,Psi2S}.yoda)
	rm -rf Rivet-EE
	python/merge-EE --with-gg EE
	python/merge-EE Powheg-EE
	rivet-mkhtml -o Rivet-EE EE.yoda:Hw Powheg-EE.yoda:Hw-Powheg
	python/plot-EE Rivet-EE
	touch $@

Rivet-DIS : Rivet-DIS/done
	touch $@

Rivet-DIS/done: $(shell echo Rivet{-DIS,-DIS-NoME,-Powheg-DIS}-{e--LowQ2,e+-LowQ2,e+-HighQ2}.yoda)
	rm -rf Rivet-DIS
	python/merge-DIS DIS 
	python/merge-DIS Powheg-DIS
	python/merge-DIS DIS-NoME
	rivet-mkhtml -o Rivet-DIS DIS.yoda:Hw Powheg-DIS.yoda:Hw-Powheg DIS-NoME.yoda:Hw-NoME
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

Rivet-TVT-Jets/done:  $(shell echo Rivet-TVT-{Run-II-Jets-{0..11},Run-I-Jets-{1..8}}.yoda ) \
	         $(shell echo Rivet-TVT-{630-Jets-{1..3},300-Jets-1,900-Jets-1}.yoda ) \
                 $(shell echo Rivet-TVT-{Run-I,Run-II,300,630,900}-UE.yoda)
	rm -rf Rivet-TVT-Jets
	python/merge-TVT-Jets TVT
	rivet-mkhtml -o Rivet-TVT-Jets TVT-Jets.yoda:Hw
	touch $@

Rivet-Star : Rivet-Star/done
	touch $@

Rivet-Star/done : $(shell echo Rivet-Star-{UE,Jets-{1..4}}.yoda ) 
	rm -rf Rivet-Star
	python/merge-Star Star
	rivet-mkhtml -o Rivet-Star Star.yoda:Hw
	touch $@



Rivet-SppS : Rivet-SppS/done
	touch $@


Rivet-SppS/done : $(shell echo Rivet-SppS-{53,63,200,500,546,900}-UE.yoda ) \
	          $(shell echo Rivet-ISR-{30,44,53,62}-UE.yoda ) Rivet-EHS-UE.yoda
	rm -rf Rivet-SppS
	python/merge-SppS SppS
	rivet-mkhtml -o Rivet-SppS SppS.yoda:Hw
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
		$(shell echo Rivet-LHC-7-Charm-{1..5}.yoda            ) \
		$(shell echo Rivet-LHC-7-Bottom-{0..9}.yoda           ) \
		$(shell echo Rivet-LHC-{7,8,13}-Top-{L,SL}.yoda ) \
		$(shell echo Rivet-LHC-{8,13}-Top-All.yoda )
	rm -rf Rivet-LHC-Jets
	python/merge-LHC-Jets LHC
	rivet-mkhtml -o Rivet-LHC-Jets LHC-Jets.yoda:Hw
	touch $@



Rivet-LHC-EW : Rivet-LHC-EW/done
	touch $@

Rivet-LHC-EW/done: \
		$(shell echo Rivet{,-Powheg}-LHC-{13-Z-{e,mu},{8,13}-W-mu,Z-HighMass{1,2}-e,8-Z-Mass{1..4}-{e,mu},W-{e,mu},Z-{e,mu,mu-SOPHTY},Z-LowMass-{e,mu},Z-MedMass-e,WZ,WW-{emu,ll},ZZ-{ll,lv},{8,13}-WZ,8-ZZ-lv,8-WW-ll,Z-mu-Short}.yoda) \
		$(shell echo Rivet-LHC-{7-W-Jet-{1..3}-e,7-Z-Jet-{0..3}-e,7-Z-Jet-0-mu}.yoda) \
		$(shell echo Rivet-LHC-7-{W,Z}Gamma-{e,mu}.yoda) \
		$(shell echo Rivet-LHC-8-ZGamma-{e,mu}.yoda) 
	rm -rf Rivet-LHC-EW;
	python/merge-LHC-EW LHC
	python/merge-LHC-EW Powheg-LHC
	rivet-mkhtml -o Rivet-LHC-EW LHC-EW.yoda:Hw Powheg-LHC-EW.yoda:Hw-Powheg \
                                     Rivet-LHC-Z-mu-SOPHTY.yoda:Hw Rivet-Powheg-LHC-Z-mu-SOPHTY.yoda:Hw-Powheg
	touch $@




Rivet-LHC-Photon : Rivet-LHC-Photon/done
	touch $@

Rivet-LHC-Photon/done: \
		$(shell echo Rivet-LHC-{7,8,13}-PromptPhoton-{1..4}.yoda) \
		Rivet-LHC-GammaGamma-7.yoda \
	    $(shell echo Rivet{,-Powheg}-LHC-{7,8}-{DiPhoton-GammaGamma,DiPhoton-GammaJet}.yoda)
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
	yodamerge --add Rivet-Powheg-LHC-8-{ggH{-GammaGamma,-WW,},{VBF,ZH,WH}{,-GammaGamma}}.yoda -o Powheg-LHC-Higgs.yoda
	yodamerge --add Rivet-LHC-8-{ggH{-GammaGamma,-WW,},{VBF,ZH,WH}{,-GammaGamma}}.yoda -o LHC-Higgs.yoda
	rm -rf Rivet-LHC-Higgs
	rivet-mkhtml -o Rivet-LHC-Higgs Powheg-LHC-Higgs.yoda:Hw-Powheg LHC-Higgs.yoda:Hw\
	                Rivet-Powheg-LHC-ggH.yoda:gg-Powheg Rivet-LHC-ggH.yoda:gg Rivet-LHC-ggHJet.yoda:HJet \
                        Rivet-Powheg-LHC-VBF.yoda:VBF-Powheg Rivet-LHC-VBF.yoda:VBF Rivet-LHC-WH.yoda:WH Rivet-LHC-ZH.yoda:ZH \
                        Rivet-Powheg-LHC-WH.yoda:WH-Powheg Rivet-Powheg-LHC-ZH.yoda:ZH-Powheg
	touch $@





clean-local:
	rm -f *.out *.log *.tex *.top *.run *.dump *.mult *.Bmult *.yoda Rivet/*.in anatohepmc.txt hepmctoana.txt
	rm -rf Rivet-*

distclean-local:
	rm -rf .cache
