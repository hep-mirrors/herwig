// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVTest class.
//

#include "VVTest.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"


using namespace Herwig;

void VVTest::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tPVector particles = event->getFinalState();
  tPVector::const_iterator iter = particles.begin(), end = particles.end();
  double sign=1;
  if(event->incoming().first->id()!=ParticleID::eminus) sign=-1;
  for( ;iter!=end;++iter) {
    if((**iter).id()==ParticleID::Wplus) {
      _cosWp->addWeighted(sign*(**iter).momentum().cosTheta(),event->weight());
      _phiWp->addWeighted((**iter).momentum().phi()+Constants::pi,event->weight());
    }
    else if((**iter).id()==ParticleID::Wminus) {
      _cosWm->addWeighted(sign*(**iter).momentum().cosTheta(),event->weight());
      _phiWm->addWeighted((**iter).momentum().phi()+Constants::pi,event->weight());
    }
    else if((**iter).id()==ParticleID::Z0) {
      _cosZ->addWeighted(sign*(**iter).momentum().cosTheta(),event->weight());
      _phiZ->addWeighted((**iter).momentum().phi()+Constants::pi,event->weight());
    }
  }
}

IBPtr VVTest::clone() const {
  return new_ptr(*this);
}

IBPtr VVTest::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<VVTest,AnalysisHandler>
describeHerwigVVTest("Herwig::VVTest", "LeptonTest.so");

void VVTest::Init() {

  static ClassDocumentation<VVTest> documentation
    ("There is no documentation for the VVTest class");

}

void VVTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _cosWp->topdrawOutput(outfile,Frame,"BLACK","Cos of polar angle for W+");
  _phiWp->topdrawOutput(outfile,Frame,"BLACK","Azimuthal angle for W+");
  _cosWm->topdrawOutput(outfile,Frame,"BLACK","Cos of polar angle for W-");
  _phiWm->topdrawOutput(outfile,Frame,"BLACK","Azimuthal angle for W-");
  _cosZ ->topdrawOutput(outfile,Frame,"BLACK","Cos of polar angle for Z ");
  _phiZ ->topdrawOutput(outfile,Frame,"BLACK","Azimuthal angle for Z ");
}

void VVTest::doinitrun() {
  AnalysisHandler::doinitrun();
  _cosWp    = new_ptr(Histogram( -1.0,              1.0,200));
  _phiWp    = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _cosWm    = new_ptr(Histogram( -1.0,              1.0,200));
  _phiWm    = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _cosZ     = new_ptr(Histogram( -1.0,              1.0,200));
  _phiZ     = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
}
