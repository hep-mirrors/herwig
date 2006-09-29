// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiggsJetAnalysis class.
//

#include "HiggsJetAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HiggsJetAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HiggsJetAnalysis::~HiggsJetAnalysis() {}

namespace {
  inline Lorentz5Momentum getMomentum(tcPPtr particle) {
    return particle->momentum();
    //Lorentz5Momentum tmp = particle->children()[0]->next()->momentum();
    //tmp += particle->children()[1]->next()->momentum();
    //tmp.rescaleMass();
    //return tmp;

  }
}


void HiggsJetAnalysis::analyze(tEventPtr event, long, int, int) {
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z
  Lorentz5Momentum ph;  
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator send=event->primaryCollision()->steps().end();
  for(;sit!=send;++sit)
    {
      ParticleSet part=(**sit).all();
      ParticleSet::const_iterator iter=part.begin();
      ParticleSet::const_iterator end =part.end();
      for( ;iter!=end;++iter)
	{
	  if((**iter).id()==ParticleID::h0)
	    {
	      ph=getMomentum(*iter);
	      Energy pt = ph.perp()/GeV;
	      (_pth)+=(pt);
	      (_pthZoom)+=(pt);
	      double rap = 0.5*log((ph.e()+ph.z())/(ph.e()-ph.z()));
	      (_raph)+=(rap);
	      (_phih)+=ph.phi();
	    }
	}
    }
}

LorentzRotation HiggsJetAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void HiggsJetAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
}

void HiggsJetAnalysis::analyze(tPPtr) {}

void HiggsJetAnalysis::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void HiggsJetAnalysis::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<HiggsJetAnalysis> HiggsJetAnalysis::initHiggsJetAnalysis;
// Definition of the static class description member.

void HiggsJetAnalysis::Init() {

  static ClassDocumentation<HiggsJetAnalysis> documentation
    ("There is no documentation for the HiggsJetAnalysis class");

}

