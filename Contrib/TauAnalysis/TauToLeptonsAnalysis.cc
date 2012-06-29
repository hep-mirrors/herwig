// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauToLeptonsAnalysis class.
//

#include "TauToLeptonsAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"


using namespace Herwig;

void TauToLeptonsAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  tPVector hadrons=event->getFinalState();
  map<tPPtr,ParticleVector> taus;
  for(unsigned int ix=0;ix<hadrons.size();++ix) {
    PPtr mother=hadrons[ix];
    do {
      if(!mother->parents().empty()) mother=mother->parents()[0];
      else                           mother=tPPtr();
    }
    while(mother&&abs(mother->id())!=ParticleID::tauminus);
    if(mother&&abs(mother->id())==ParticleID::tauminus) {
      if(taus.find(mother)==taus.end()) {
	taus.insert(make_pair(mother,ParticleVector()));
      }
      taus[mother].push_back(hadrons[ix]);
    }
  }
  map<tPPtr,ParticleVector>::const_iterator tit;
  for(tit=taus.begin();tit!=taus.end();++tit) {
    if(tit->second.size()!=3) continue;
    vector<Lorentz5Momentum> pdecay;
    int type(0);
    ParticleVector decay=tit->second;
    for(unsigned int ix=0;ix<decay.size();++ix) {
      long id = abs(decay[ix]->id());
      if(id>=11&&id<=14) {
	pdecay.push_back(decay[ix]->momentum());
	if(id==11) type=1;
	else if(id==13) type=2;
      }
    }
    if(pdecay.size()==2&&type!=0) {
      Lorentz5Momentum pw=pdecay[0]+pdecay[1];
      pw.rescaleMass();
      if(type==1) {
	*_emode+=pw.mass()/GeV;
      }
      else if(type==2) {
	*_mmode+=pw.mass()/GeV;
      }
    }
  }
}

NoPIOClassDescription<TauToLeptonsAnalysis> TauToLeptonsAnalysis::initTauToLeptonsAnalysis;
// Definition of the static class description member.

void TauToLeptonsAnalysis::Init() {

  static ClassDocumentation<TauToLeptonsAnalysis> documentation
    ("There is no documentation for the TauToLeptonsAnalysis class");

}

void TauToLeptonsAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _emode->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"e2-3N0O0e1 mass for T2-3RN0T1e2-3N0O0e1",
			" X XGUDX X          GX XWGXGX X XGUDX X",
			"1/SdS/dm0e2-3N0O0e11/GeV2-13",
			"  G G   X X XGUDX XX    X  X",
			"m0e2-3N0O0e11/GeV",
			" X X XGUDX XX    ");
  _mmode->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"M2-3N0O0M1 mass for T2-3RN0T1M2-3N0O0M1",
			"GX XGUDXGX          GX XWGXGXGX XGUDXGX",
			"1/SdS/dm0M2-3N0O0M11/GeV2-13",
			"  G G   XGX XGUDXGXX    X  X",
			"m0M2-3N0O0M11/GeV",
			" XGX XGUDXGXX    ");
}

void TauToLeptonsAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  _emode=new_ptr(Histogram(0.,1.8,200));
  _mmode=new_ptr(Histogram(0.,1.8,200));
}
