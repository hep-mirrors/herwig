// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleVBFAnalysis class.
//

#include "SimpleVBFAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

void SimpleVBFAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tPVector part = event->getFinalState();
  tPVector jets;
  for(tPVector::const_iterator iter = part.begin(), end = part.end();
      iter!=end;++iter) {
    if((**iter).id()==ParticleID::h0) continue;
    if((**iter).momentum().perp()<20.*GeV) continue;
    if(abs((**iter).momentum().eta())>4.5) continue;
    jets.push_back(*iter);
  }
  if(jets.size()!=2) return;
  double eta[2]={jets[0]->momentum().eta(),
		 jets[1]->momentum().eta()};
  if(eta[0]*eta[1]>0.) return;
  if(abs(eta[0]-eta[1])<4.2) return;
  double deltaPhi = jets[0]->momentum().phi()-jets[1]->momentum().phi();
  if(deltaPhi<-Constants::pi) deltaPhi+=Constants::twopi;
  if(deltaPhi> Constants::pi) deltaPhi-=Constants::twopi;
  *_deltaPhi += abs(deltaPhi);
}

IBPtr SimpleVBFAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleVBFAnalysis::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<SimpleVBFAnalysis> SimpleVBFAnalysis::initSimpleVBFAnalysis;
// Definition of the static class description member.

void SimpleVBFAnalysis::Init() {

  static ClassDocumentation<SimpleVBFAnalysis> documentation
    ("There is no documentation for the SimpleVBFAnalysis class");

}

void SimpleVBFAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _deltaPhi->topdrawOutput(outfile,Frame,"BLACK","Delta Phi");
  outfile.close();
}

void SimpleVBFAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  _deltaPhi = new_ptr(Histogram(0,Constants::pi,200));
}
