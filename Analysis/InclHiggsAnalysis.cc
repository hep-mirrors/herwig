// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the InclHiggsAnalysis class.
//

#include "InclHiggsAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Herwig;

InclHiggsAnalysis::InclHiggsAnalysis() :
  _nJets_gt_20(-0.5,5.5,6),

  _ptHiggsHigh(0.,250.,25), _ptHiggsLow(0., 50.,25),

  _YHiggs(-5.25,5.25,21) ,

  _ptJet0(0.,250.,25)  , _ptJet1(0.,150.,15),
  _ptJet2(0.,100.,10)  , _ptJet3(0.,100.,10),

  _etaJet0(-5.25,5.25,21), _etaJet1(-5.25,5.25,20),
  _etaJet2(-5.25,5.25,21), _etaJet3(-5.25,5.25,20),

  _delta01(0.,10.,50)   , _delta12(0.,10.,50), _delta23(0.,10.,50),

  _htJets(0.,500.,50)  , _htJetsPlusHiggs(0.,500.,50),

  _y01(-6.,6.,50), _y12(-6.,6.,50), _y23(-6.,6.,50), _y34(-6.,6.,50)
{}

void InclHiggsAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {

  // Store all of the final state particles as leptons.
  tPVector leptons(event->getFinalState());

  // Store all non-higgses and non-leptonic particles in partons.
  tPVector partons;
  for(unsigned int ix=0; ix<leptons.size();++ix)
    if((abs(leptons[ix]->id())>16||abs(leptons[ix]->id())<11)
       &&leptons[ix]->id()!=25
       )
      partons.push_back(leptons[ix]);

  // Delete all non-higgses and non-leptonic particles from leptons.
  for(unsigned int ix=0; ix<leptons.size();++ix)
    if((abs(leptons[ix]->id())>16||abs(leptons[ix]->id())<11)
       &&leptons[ix]->id()!=25
       ) { 
      leptons.erase(leptons.begin()+ix);
      ix--;
    }

  // Store theHiggs
  tPPtr theHiggs;
  unsigned int position(0);
  for(unsigned int ix=0; ix<leptons.size();++ix)
    if(leptons[ix]->id()==25) {
      theHiggs=leptons[ix];
      position = ix;
    }
  leptons.erase(leptons.begin()+position);

  // Get the jets ordered by their pT (largest pT is first).
  vector<fastjet::PseudoJet> particlesToCluster;
  for(unsigned int jx=0; jx<partons.size();jx++) {
    fastjet::PseudoJet p(partons[jx]->momentum().x()/GeV,
			 partons[jx]->momentum().y()/GeV,
			 partons[jx]->momentum().z()/GeV,
			 partons[jx]->momentum().e()/GeV);
    p.set_user_index(jx);
    particlesToCluster.push_back(p);
  }
  fastjet::RecombinationScheme recombinationScheme = fastjet::E_scheme;
  fastjet::Strategy            strategy            = fastjet::Best;
  double R(0.7);
  fastjet::JetDefinition       jetDefinition(fastjet::kt_algorithm,
					     R,
					     recombinationScheme,
					     strategy);
  fastjet::ClusterSequence fastjetEvent(particlesToCluster,jetDefinition);
  vector<fastjet::PseudoJet> inclusiveJets = fastjetEvent.inclusive_jets();
  inclusiveJets = fastjet::sorted_by_pt(inclusiveJets);

  // How many jets were found?
  unsigned int nJets(inclusiveJets.size());

  // Jet multiplicity for jets with pT > 20 GeV:
  unsigned int nJets_gt_20(0);
  for(unsigned int ix=0;ix<nJets;ix++) 
    if(inclusiveJets[ix].perp()>20.&&fabs(inclusiveJets[ix].eta())<4.5) nJets_gt_20++;
  _nJets_gt_20.addWeighted(nJets_gt_20+0.5,1.);

  // Higgs pT spectrum in the high and low pT ranges:
  Energy ptHiggs(theHiggs->momentum().perp());
  _ptHiggsHigh.addWeighted(ptHiggs/GeV,1.);
  _ptHiggsLow.addWeighted(ptHiggs/GeV,1.);

  // Higgs pseudorapidity:
  double YHiggs(theHiggs->momentum().rapidity());
  _YHiggs.addWeighted(YHiggs,1.);

  // 1st hardest jet pT and pseudorapidity.
  if(nJets>0) {
    double ptJet0(inclusiveJets[0].perp());
    _ptJet0.addWeighted(ptJet0,1.);
    double etaJet0(inclusiveJets[0].eta());
    _etaJet0.addWeighted(etaJet0,1.);
  }

  // 2nd hardest jet pT and pseudorapidity, 
  // also DeltaR between 1st and 2nd hardest jets. 
  if(nJets>1) {
    double ptJet1(inclusiveJets[1].perp());
    _ptJet1.addWeighted(ptJet1,1.);
    double etaJet1(inclusiveJets[1].eta());
    _etaJet1.addWeighted(etaJet1,1.);
    double delta01 = sqrt( sqr(inclusiveJets[0].eta()-inclusiveJets[1].eta())
			  +sqr(inclusiveJets[0].phi()-inclusiveJets[1].phi()) );
    _delta01.addWeighted(delta01,1.);
  }

  // 3rd hardest jet pT and pseudorapidity, 
  // also DeltaR between 2nd and 3rd hardest jets. 
  if(nJets>2) {
    double ptJet2(inclusiveJets[2].perp());
    _ptJet2.addWeighted(ptJet2,1.);
    double etaJet2(inclusiveJets[2].eta());
    _etaJet2.addWeighted(etaJet2,1.);
    double delta12 = sqrt( sqr(inclusiveJets[1].eta()-inclusiveJets[2].eta())
			  +sqr(inclusiveJets[1].phi()-inclusiveJets[2].phi()) );
    _delta12.addWeighted(delta12,1.);
  }

  // 4th hardest jet pT and pseudorapidity, 
  // also DeltaR between 3rd and 4th hardest jets. 
  if(nJets>3) {
    double ptJet3(inclusiveJets[3].perp());
    _ptJet3.addWeighted(ptJet3,1.);
    double etaJet3(inclusiveJets[3].eta());
    _etaJet3.addWeighted(etaJet3,1.);
    double delta23 = sqrt( sqr(inclusiveJets[2].eta()-inclusiveJets[3].eta())
			  +sqr(inclusiveJets[2].phi()-inclusiveJets[3].phi()) );
    _delta23.addWeighted(delta23,1.);
  }

  // Scalar pT sum of all jets:
  double htJets(0.);
  for(unsigned int ix=0;ix<nJets;ix++) htJets += inclusiveJets[ix].perp();
  _htJets.addWeighted(htJets,1.);

  // Scalar pT sum of all jets plus the pT of the Higgs:
  _htJetsPlusHiggs.addWeighted(htJets+ptHiggs/GeV,1.);

//   _y01.addWeighted(,1.);
//   _y12.addWeighted(,1.);
//   _y23.addWeighted(,1.);
//   _y34.addWeighted(,1.);

}

void InclHiggsAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.

}

void InclHiggsAnalysis::analyze(tPPtr) {}

// ClassDescription<InclHiggsAnalysis> InclHiggsAnalysis::initInclHiggsAnalysis;
// // Definition of the static class description member.

NoPIOClassDescription<InclHiggsAnalysis> InclHiggsAnalysis::initInclHiggsAnalysis;
// Definition of the static class description member.

void InclHiggsAnalysis::Init() {

  static ClassDocumentation<InclHiggsAnalysis> documentation
    ("There is no documentation for the InclHiggsAnalysis class");

}

void InclHiggsAnalysis::dofinish() {
  AnalysisHandler::dofinish();

  using namespace HistogramOptions;

  ofstream topdrawFile;
  string topdrawFilename 
    = generator()->filename()+string("_")+name()+string("_mcatnlo.top");
  topdrawFile.open(topdrawFilename.c_str());

  _nJets_gt_20.topdrawOutput(topdrawFile,Frame,"RED","No.Jets pT>20 |eta|<4.5");
  _ptHiggsHigh.topdrawOutput(topdrawFile,Frame,"RED","pT Higgs");
  _ptHiggsLow.topdrawOutput(topdrawFile,Frame,"RED","pT Higgs");
  _YHiggs.topdrawOutput(topdrawFile,Frame,"RED","Y Higgs");
  _ptJet0.topdrawOutput(topdrawFile,Frame,"RED","pT 1st Jet");
  _ptJet1.topdrawOutput(topdrawFile,Frame,"RED","pT 2nd Jet");
  _ptJet2.topdrawOutput(topdrawFile,Frame,"RED","pT 3rd Jet");
  _ptJet3.topdrawOutput(topdrawFile,Frame,"RED","pT 4th Jet");
  _etaJet0.topdrawOutput(topdrawFile,Frame,"RED","eta 1st Jet");
  _etaJet1.topdrawOutput(topdrawFile,Frame,"RED","eta 2nd Jet");
  _etaJet2.topdrawOutput(topdrawFile,Frame,"RED","eta 3rd Jet");
  _etaJet3.topdrawOutput(topdrawFile,Frame,"RED","eta 4th Jet");
  _delta01.topdrawOutput(topdrawFile,Frame,"RED","Delta 01");
  _delta12.topdrawOutput(topdrawFile,Frame,"RED","Delta 12");
  _delta23.topdrawOutput(topdrawFile,Frame,"RED","Delta 23");
  _htJets.topdrawOutput(topdrawFile,Frame,"RED","HT Jets");
  _htJetsPlusHiggs.topdrawOutput(topdrawFile,Frame,"RED","HT Jets plus Higgs");
//   _y01.topdrawOutput(topdrawFile,Frame,"RED","y01");
//   _y12.topdrawOutput(topdrawFile,Frame,"RED","y12");
//   _y23.topdrawOutput(topdrawFile,Frame,"RED","y23");
//   _y34.topdrawOutput(topdrawFile,Frame,"RED","y34");

  _nJets_gt_20.normaliseToCrossSection();
  _nJets_gt_20.prefactor(_nJets_gt_20.prefactor()*1.e6);
  _nJets_gt_20.topdrawOutput(topdrawFile,Frame,"RED","No.Jets pT>20 |eta|<4.5 (fb)");
  _ptHiggsHigh.normaliseToCrossSection();
  _ptHiggsHigh.prefactor(_ptHiggsHigh.prefactor()*1.e6);
  _ptHiggsHigh.topdrawOutput(topdrawFile,Frame,"RED","pT Higgs (fb)");
  _ptHiggsLow.normaliseToCrossSection();
  _ptHiggsLow.prefactor(_ptHiggsLow.prefactor()*1.e6);
  _ptHiggsLow.topdrawOutput(topdrawFile,Frame,"RED","pT Higgs (fb)");
  _YHiggs.normaliseToCrossSection();
  _YHiggs.prefactor(_YHiggs.prefactor()*1.e6);
  _YHiggs.topdrawOutput(topdrawFile,Frame,"RED","Y Higgs (fb)");
  _ptJet0.normaliseToCrossSection();
  _ptJet0.prefactor(_ptJet0.prefactor()*1.e6);
  _ptJet0.topdrawOutput(topdrawFile,Frame,"RED","pT 1st Jet (fb)");
  _ptJet1.normaliseToCrossSection();
  _ptJet1.prefactor(_ptJet1.prefactor()*1.e6);
  _ptJet1.topdrawOutput(topdrawFile,Frame,"RED","pT 2nd Jet (fb)");
  _ptJet2.normaliseToCrossSection();
  _ptJet2.prefactor(_ptJet2.prefactor()*1.e6);
  _ptJet2.topdrawOutput(topdrawFile,Frame,"RED","pT 3rd Jet (fb)");
  _ptJet3.normaliseToCrossSection();
  _ptJet3.prefactor(_ptJet3.prefactor()*1.e6);
  _ptJet3.topdrawOutput(topdrawFile,Frame,"RED","pT 4th Jet (fb)");
  _etaJet0.normaliseToCrossSection();
  _etaJet0.prefactor(_etaJet0.prefactor()*1.e6);
  _etaJet0.topdrawOutput(topdrawFile,Frame,"RED","eta 1st Jet (fb)");
  _etaJet1.normaliseToCrossSection();
  _etaJet1.prefactor(_etaJet1.prefactor()*1.e6);
  _etaJet1.topdrawOutput(topdrawFile,Frame,"RED","eta 2nd Jet (fb)");
  _etaJet2.normaliseToCrossSection();
  _etaJet2.prefactor(_etaJet2.prefactor()*1.e6);
  _etaJet2.topdrawOutput(topdrawFile,Frame,"RED","eta 3rd Jet (fb)");
  _etaJet3.normaliseToCrossSection();
  _etaJet3.prefactor(_etaJet3.prefactor()*1.e6);
  _etaJet3.topdrawOutput(topdrawFile,Frame,"RED","eta 4th Jet (fb)");
  _delta01.normaliseToCrossSection();
  _delta01.prefactor(_delta01.prefactor()*1.e6);
  _delta01.topdrawOutput(topdrawFile,Frame,"RED","Delta 01 (fb)");
  _delta12.normaliseToCrossSection();
  _delta12.prefactor(_delta12.prefactor()*1.e6);
  _delta12.topdrawOutput(topdrawFile,Frame,"RED","Delta 12 (fb)");
  _delta23.normaliseToCrossSection();
  _delta23.prefactor(_delta23.prefactor()*1.e6);
  _delta23.topdrawOutput(topdrawFile,Frame,"RED","Delta 23 (fb)");
  _htJets.normaliseToCrossSection();
  _htJets.prefactor(_htJets.prefactor()*1.e6);
  _htJets.topdrawOutput(topdrawFile,Frame,"RED","HT Jets (fb)");
  _htJetsPlusHiggs.normaliseToCrossSection();
  _htJetsPlusHiggs.prefactor(_htJetsPlusHiggs.prefactor()*1.e6);
  _htJetsPlusHiggs.topdrawOutput(topdrawFile,Frame,"RED","HT Jets plus Higgs (fb)");
//   _y01.normaliseToCrossSection();
//   _y01.prefactor(_y01.prefactor()*1.*e6);
//   _y01.topdrawOutput(topdrawFile,Frame,"RED","y01 (fb)");
//   _y12.normaliseToCrossSection();
//   _y12.prefactor(_y12.prefactor()*1.*e6);
//   _y12.topdrawOutput(topdrawFile,Frame,"RED","y12 (fb)");
//   _y23.normaliseToCrossSection();
//   _y23.prefactor(_y23.prefactor()*1.*e6);
//   _y23.topdrawOutput(topdrawFile,Frame,"RED","y23 (fb)");
//   _y34.normaliseToCrossSection();
//   _y34.prefactor(_y34.prefactor()*1.*e6);
//   _y34.topdrawOutput(topdrawFile,Frame,"RED","y34 (fb)");

  topdrawFile.close();
}

