// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaPMETest class.
//

#include "GammaPMETest.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

void GammaPMETest::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tPVector particles = event->getFinalState();
  Lorentz5Momentum ptotal;
  for(tPVector::const_iterator it=particles.begin();it!=particles.end();++it) {
    int id = abs((**it).id());
    if(id<=5||id==ParticleID::g) {
      *_pt  += (**it).momentum().perp()/GeV;
      *_rap += (**it).momentum().rapidity();
      *_phi += (**it).momentum().phi()+Constants::pi;
      ptotal += (**it).momentum();
    }
  }
  *_mhat += ptotal.m()/GeV;
//   cerr << "testing at analysis " << ptotal.m()/GeV << "\n";
  *_yhat += ptotal.rapidity();
}

IBPtr GammaPMETest::clone() const {
  return new_ptr(*this);
}

IBPtr GammaPMETest::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<GammaPMETest> GammaPMETest::initGammaPMETest;
// Definition of the static class description member.

void GammaPMETest::Init() {

  static ClassDocumentation<GammaPMETest> documentation
    ("There is no documentation for the GammaPMETest class");

}

void GammaPMETest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _rap->topdrawOutput(outfile,Frame,"BLACK","Jet rapidity");
  _pt ->topdrawOutput(outfile,Frame,"BLACK","Jet pt");
  _phi->topdrawOutput(outfile,Frame,"BLACK","Jet azimuth");
  _mhat->topdrawOutput(outfile,Frame,"BLACK","CMS Mass");
  _yhat->topdrawOutput(outfile,Frame,"BLACK","CMS rapidity");
}

void GammaPMETest::doinitrun() {
  AnalysisHandler::doinitrun();
  _rap  = new_ptr(Histogram(-10.,10.,200));
  _pt   = new_ptr(Histogram(0.,1000.,1000));
  _phi  = new_ptr(Histogram(0,Constants::twopi,200));
  _mhat = new_ptr(Histogram(0.,1000.,1000));
  _yhat = new_ptr(Histogram(-10.,10.,200));
}
