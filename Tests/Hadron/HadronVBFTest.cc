// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronVBFTest class.
//

#include "HadronVBFTest.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

void HadronVBFTest::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tPVector part = event->getFinalState();
  Lorentz5Momentum pjj;
  for(tPVector::const_iterator iter = part.begin(), end = part.end();
      iter!=end;++iter) {
    if((**iter).id()==ParticleID::h0) {
      *_mH     += (**iter).momentum().m()/GeV;
      *_yH     += (**iter).momentum().rapidity();
      *_phiH   += (**iter).momentum().phi()+Constants::pi;
      *_pTH[0] += (**iter).momentum().perp()/GeV;
      *_pTH[1] += (**iter).momentum().perp()/GeV;
    }
    else if((**iter).id()!=82) {
      *_yjet     += (**iter).momentum().rapidity();
      *_phijet   += (**iter).momentum().phi()+Constants::pi;
      *_pTjet[0] += (**iter).momentum().perp()/GeV;
      *_pTjet[1] += (**iter).momentum().perp()/GeV;
      pjj+=(**iter).momentum();
    }
  }
  *_mjj += pjj.m()/GeV;
}

IBPtr HadronVBFTest::clone() const {
  return new_ptr(*this);
}

IBPtr HadronVBFTest::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<HadronVBFTest> HadronVBFTest::initHadronVBFTest;
// Definition of the static class description member.

void HadronVBFTest::Init() {

  static ClassDocumentation<HadronVBFTest> documentation
    ("There is no documentation for the HadronVBFTest class");

}

void HadronVBFTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  string title,species;
  title = "mass of H";
  _mH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "rapidity of H";
  _yH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "pT of H";
  _pTH[0]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  _pTH[1]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  title = "azimuth of H";
  _phiH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "rapidity of jet";
  _yjet->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "pT of jet";
  _pTjet[0]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  _pTjet[1]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  title = "azimuth of jet";
  _phijet->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "mjj";
  _mjj->topdrawOutput(outfile,Frame,"BLACK",title);
}

void HadronVBFTest::doinitrun() {
  AnalysisHandler::doinitrun();
  if(getParticleData(ParticleID::h0)->mass()>200.*GeV) 
    _mH     = new_ptr(Histogram(200.,            400.,200));
  else
    _mH     = new_ptr(Histogram(114.,            116.0,200));
  _yH       = new_ptr(Histogram( -10.0,              10.0,200));
  _phiH     = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _pTH[0]      = new_ptr(Histogram(  0.0,1000.,1000));
  _pTH[1]      = new_ptr(Histogram(  0.0,1000.,100));
  _yjet     = new_ptr(Histogram( -10.0,              10.0,200));
  _phijet   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _pTjet[0] = new_ptr(Histogram(  0.0,1000.,1000));
  _pTjet[1] = new_ptr(Histogram(  0.0,1000.,100));
  _mjj = new_ptr(Histogram(0.0,2000.,100));
}
