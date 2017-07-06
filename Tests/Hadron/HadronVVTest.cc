// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronVVTest class.
//

#include "HadronVVTest.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

void HadronVVTest::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator stest =event->primaryCollision()->steps().end();
  StepVector::const_iterator send=sit;
  ++send;
  if(send==stest) --send;
  ++send;
  if(send==stest) --send;
  ++send;
  Lorentz5Momentum ptotal;
  for(;sit!=send;++sit) {
    ParticleSet part;
    (**sit).selectFinalState(inserter(part));
    ParticleSet::const_iterator iter = part.begin(), end = part.end();
    for( ;iter!=end;++iter) {
      if((**iter).id()==ParticleID::Wplus) {
	ptotal+=(**iter).momentum();
	_ptWp->addWeighted((**iter).momentum().perp()/GeV,event->weight());
	_rapWp->addWeighted((**iter).momentum().rapidity(),event->weight());
	_phiWp->addWeighted((**iter).momentum().phi()+Constants::pi,event->weight());
      }
      else if((**iter).id()==ParticleID::Wminus) {
	ptotal+=(**iter).momentum();
	_ptWm->addWeighted((**iter).momentum().perp()/GeV,event->weight());
	_rapWm->addWeighted((**iter).momentum().rapidity(),event->weight());
	_phiWm->addWeighted((**iter).momentum().phi()+Constants::pi,event->weight());
      }
      else if((**iter).id()==ParticleID::Z0) {
	ptotal+=(**iter).momentum();
	_ptZ->addWeighted((**iter).momentum().perp()/GeV,event->weight());
	_rapZ->addWeighted((**iter).momentum().rapidity(),event->weight());
	_phiZ->addWeighted((**iter).momentum().phi()+Constants::pi,event->weight());
      }
    }
  }
  _mass->addWeighted(ptotal.m()/GeV,event->weight());
}

IBPtr HadronVVTest::clone() const {
  return new_ptr(*this);
}

IBPtr HadronVVTest::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<HadronVVTest> HadronVVTest::initHadronVVTest;
// Definition of the static class description member.

void HadronVVTest::Init() {

  static ClassDocumentation<HadronVVTest> documentation
    ("There is no documentation for the HadronVVTest class");

}

void HadronVVTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _ptWp->topdrawOutput(outfile,Frame,"BLACK","Pt of W+");
  _rapWp->topdrawOutput(outfile,Frame,"BLACK","Rapidity of W+");
  _phiWp->topdrawOutput(outfile,Frame,"BLACK","Azimuthal angle for W+");
  _ptWm->topdrawOutput(outfile,Frame,"BLACK","Pt of W-");
  _rapWm->topdrawOutput(outfile,Frame,"BLACK","Rapidity of W-");
  _phiWm->topdrawOutput(outfile,Frame,"BLACK","Azimuthal angle for W-");
  _ptZ ->topdrawOutput(outfile,Frame,"BLACK","Pt of Z ");
  _rapZ ->topdrawOutput(outfile,Frame,"BLACK","Rapidity of Z ");
  _phiZ ->topdrawOutput(outfile,Frame,"BLACK","Azimuthal angle for Z ");
  _mass ->topdrawOutput(outfile,Frame,"BLACK","Pair mass");
}

void HadronVVTest::doinitrun() {
  AnalysisHandler::doinitrun();
  _ptWp  = new_ptr(Histogram(0.,1000.,200));
  _rapWp = new_ptr(Histogram(-10,10,200));
  _phiWp = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _ptWm  = new_ptr(Histogram(0.,1000.,200));
  _rapWm = new_ptr(Histogram(-10,10,200));
  _phiWm = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _ptZ   = new_ptr(Histogram(0.,1000.,200));
  _rapZ  = new_ptr(Histogram(-10,10,200));
  _phiZ  = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _mass  = new_ptr(Histogram(100.,1100.,200));
}
