// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Higgspt class.
//

#include "Higgspt.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Histogram.h"


using namespace Herwig;
using namespace ThePEG;
using namespace std;

Histogram Higgspt_30(  0., 30.,50);
Histogram Higgspt_200( 0.,200.,50);
Histogram Higgspt_400( 0.,400.,50);
Histogram HiggsYjet_10(  -5.,  5.,50);
Histogram HiggsYjetYH_10(-5.,  5.,50);
Histogram HiggsYjet_40(  -5.,  5.,50);
Histogram HiggsYjetYH_40(-5.,  5.,50);
Histogram HiggsYjet_80(  -5.,  5.,50);
Histogram HiggsYjetYH_80(-5.,  5.,50);
Histogram HiggsYjet_10A(  -5.,  5.,50);
Histogram HiggsYjetYH_10A(-5.,  5.,50);
Histogram HiggsYjet_40A(  -5.,  5.,50);
Histogram HiggsYjetYH_40A(-5.,  5.,50);
Histogram HiggsYjet_80A(  -5.,  5.,50);
Histogram HiggsYjetYH_80A(-5.,  5.,50);
Histogram Njets_10(0.,10.,10);
Histogram Njets_40(0.,10.,10);
Histogram Njets_80(0.,10.,10);

void Higgspt::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void Higgspt::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<Higgspt> Higgspt::initHiggspt;
// Definition of the static class description member.

void Higgspt::Init() {

  static ClassDocumentation<Higgspt> documentation
    ("There is no documentation for the Higgspt class");

}

void Higgspt::analyze(tEventPtr event, long ieve, int loop, int state) {

  // Store all of the final state particles as particles.
  particles=event->getFinalState();
  // Extract the emitted parton
  tPVector::iterator pit;
  int Higgses(0);
  for(pit=particles.begin();pit!=particles.end();++pit) {
    if((*pit)->id()==25) { 
       Higgs = *pit;
       particles.erase(pit); // Delete the Higgs so it is not clustered
       Higgses++;
    }
  }
  if(Higgses>1||Higgses==0)
	throw Exception() << "Higgspt::analyze"  
			  << "\nMore / less than one Higgs in the final state"
			  << "\nNumberof Higgses is " << Higgses 
			  << Exception::warning;

  Energy pTh = Higgs->momentum().perp();
  Higgspt_30.addWeighted( pTh/GeV,1.);
  Higgspt_200.addWeighted(pTh/GeV,1.);
  Higgspt_400.addWeighted(pTh/GeV,1.);

  _kint.clearMap();
  KtJet::KtEvent ev = KtJet::KtEvent(_kint.convert(particles),4,2,1,0.7);
  // Get the two jets ordered by their Pt (largest Pt is first).
  vector<KtJet::KtLorentzVector> ktjets = ev.getJetsPt();
  if(ktjets.size()>0&&Higgs) {
  if(ktjets[0].perp()>=10000.) {
    HiggsYjet_10.addWeighted(ktjets[0].rapidity(),1.);
    HiggsYjetYH_10.addWeighted(
	ktjets[0].rapidity()-Higgs->momentum().rapidity()
	,1.);
  }
  if(ktjets[0].perp()>=40000.) {
    HiggsYjet_40.addWeighted(ktjets[0].rapidity(),1.);
    HiggsYjetYH_40.addWeighted(
	ktjets[0].rapidity()-Higgs->momentum().rapidity()
	,1.);
  }
  if(ktjets[0].perp()>=80000.) {
    HiggsYjet_80.addWeighted(ktjets[0].rapidity(),1.);
    HiggsYjetYH_80.addWeighted(
        ktjets[0].rapidity()-Higgs->momentum().rapidity()
	,1.);
  }
  unsigned int last(ktjets.size()-1);
  if(ktjets[last]>=10000.) {
    HiggsYjet_10A.addWeighted(ktjets[0].rapidity(),1.);
    HiggsYjetYH_10A.addWeighted(
	ktjets[0].rapidity()-Higgs->momentum().rapidity()
	,1.);
  }
  if(ktjets[last]>=40000.) {
    HiggsYjet_40A.addWeighted(ktjets[0].rapidity(),1.);
    HiggsYjetYH_40A.addWeighted(
	ktjets[0].rapidity()-Higgs->momentum().rapidity()
	,1.);
  }
  if(ktjets[last]>=80000.) {
    HiggsYjet_80A.addWeighted(ktjets[0].rapidity(),1.);
    HiggsYjetYH_80A.addWeighted(
	ktjets[0].rapidity()-Higgs->momentum().rapidity()
	,1.);
  }
  }
  unsigned int n10(0);
  unsigned int n40(0);
  unsigned int n80(0);
  for(unsigned int ix=0; ix<ktjets.size(); ++ix) {
    if(ktjets[ix].perp()>=10000.) n10++;
    if(ktjets[ix].perp()>=40000.) n40++;
    if(ktjets[ix].perp()>=80000.) n80++;
  }
  Njets_10.addWeighted(n10+0.001,1);
  Njets_40.addWeighted(n40+0.001,1);
  Njets_80.addWeighted(n80+0.001,1);
}

void Higgspt::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename()+string("-")+name()+string(".top");
  ofstream outfile(fname.c_str());

  using namespace HistogramOptions;
  Higgspt_30.normaliseToCrossSection();
  Higgspt_30.prefactor(Higgspt_30.prefactor()*1.e3);
  Higgspt_30.topdrawOutput(outfile,Frame,"RED","pT of Higgs (pb)");
  Higgspt_200.normaliseToCrossSection();
  Higgspt_200.prefactor(Higgspt_200.prefactor()*1.e3);
  Higgspt_200.topdrawOutput(outfile,Frame|Ylog,"RED","pT of Higgs (pb)");
  Higgspt_400.normaliseToCrossSection();
  Higgspt_400.prefactor(Higgspt_400.prefactor()*1.e3);
  Higgspt_400.topdrawOutput(outfile,Frame|Ylog,"RED","pT of Higgs (pb)");
  HiggsYjet_10.normaliseToCrossSection();
  HiggsYjet_10.prefactor(HiggsYjet_10.prefactor()*1.e3);
  HiggsYjet_10.topdrawOutput(outfile,Frame|Ylog,
			     "RED","Yjet (pb)");
  HiggsYjet_40.normaliseToCrossSection();
  HiggsYjet_40.prefactor(HiggsYjet_40.prefactor()*1.e3);
  HiggsYjet_40.topdrawOutput(outfile,Frame|Ylog,
			     "RED","Yjet (pb)");
  HiggsYjet_80.normaliseToCrossSection();
  HiggsYjet_80.prefactor(HiggsYjet_80.prefactor()*1.e3);
  HiggsYjet_80.topdrawOutput(outfile,Frame|Ylog,
			     "RED","Yjet (pb)");
  HiggsYjetYH_10.normaliseToCrossSection();
  HiggsYjetYH_10.prefactor(HiggsYjetYH_10.prefactor()*1.e3);
  HiggsYjetYH_10.topdrawOutput(outfile,Frame|Ylog,
			       "RED","Yjet-YH (pb)");
  HiggsYjetYH_40.normaliseToCrossSection();
  HiggsYjetYH_40.prefactor(HiggsYjetYH_40.prefactor()*1.e3);
  HiggsYjetYH_40.topdrawOutput(outfile,Frame|Ylog,
			       "RED","Yjet-YH (pb) pT>40 GeV");
  HiggsYjetYH_80.normaliseToCrossSection();
  HiggsYjetYH_80.prefactor(HiggsYjetYH_80.prefactor()*1.e3);
  HiggsYjetYH_80.topdrawOutput(outfile,Frame|Ylog,
			       "RED","Yjet-YH (pb) pT>80 GeV");
  HiggsYjet_10A.normaliseToCrossSection();
  HiggsYjet_10A.prefactor(HiggsYjet_10A.prefactor()*1.e3);
  HiggsYjet_10A.topdrawOutput(outfile,Frame|Ylog,
			      "RED","Yjet (pb) pT>10 GeV");
  HiggsYjet_40A.normaliseToCrossSection();
  HiggsYjet_40A.prefactor(HiggsYjet_40A.prefactor()*1.e3);
  HiggsYjet_40A.topdrawOutput(outfile,Frame|Ylog,
			      "RED","Yjet (pb) pT>40 GeV");
  HiggsYjet_80A.normaliseToCrossSection();
  HiggsYjet_80A.prefactor(HiggsYjet_80A.prefactor()*1.e3);
  HiggsYjet_80A.topdrawOutput(outfile,Frame|Ylog,
			      "RED","Yjet (pb) pT>80 GeV");
  HiggsYjetYH_10A.normaliseToCrossSection();
  HiggsYjetYH_10A.prefactor(HiggsYjetYH_10A.prefactor()*1.e3);
  HiggsYjetYH_10A.topdrawOutput(outfile,Frame|Ylog,
				"RED","Yjet-YH (pb) pT>10 GeV");
  HiggsYjetYH_40A.normaliseToCrossSection();
  HiggsYjetYH_40A.prefactor(HiggsYjetYH_40A.prefactor()*1.e3);
  HiggsYjetYH_40A.topdrawOutput(outfile,Frame|Ylog,
				"RED","Yjet-YH (pb) pT>40 GeV");
  HiggsYjetYH_80A.normaliseToCrossSection();
  HiggsYjetYH_80A.prefactor(HiggsYjetYH_80A.prefactor()*1.e3);
  HiggsYjetYH_80A.topdrawOutput(outfile,Frame|Ylog,
				"RED","Yjet-YH (pb) pT>80 GeV");
  Njets_10.normaliseToCrossSection();
  Njets_10.prefactor(Njets_10.prefactor()*1.e3);
  Njets_10.topdrawOutput(outfile,Frame,"RED","Njets (pb) pT>10 GeV");
  Njets_40.normaliseToCrossSection();
  Njets_40.prefactor(Njets_40.prefactor()*1.e3);
  Njets_40.topdrawOutput(outfile,Frame,"RED","Njets (pb) pT>40 GeV");
  Njets_80.normaliseToCrossSection();
  Njets_80.prefactor(Njets_80.prefactor()*1.e3);
  Njets_80.topdrawOutput(outfile,Frame,"RED","Njets (pb) pT>80 GeV");

  outfile.close();
}

