// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaJetAnalysis class.
//

#include "GammaJetAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GammaJetAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

GammaJetAnalysis::~GammaJetAnalysis() {}

namespace {
  inline Lorentz5Momentum getMomentum(tcPPtr particle) {
    return particle->momentum();
    //Lorentz5Momentum tmp = particle->children()[0]->next()->momentum();
    //tmp += particle->children()[1]->next()->momentum();
    //tmp.rescaleMass();
    //return tmp;

  }
}


void GammaJetAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z
  Lorentz5Momentum pg;  
  bool foundphoton = false;
  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if((**it).id()==ParticleID::gamma) {
      // only book the hardest photon in the event
      if( (**it).momentum().perp() > pg.perp() ) {
	foundphoton = true;
	pg=getMomentum(*it);
      }
    }
  }

  if (foundphoton) {
    Energy pt = pg.perp();
    (_ptg)+=(pt)/GeV;
    (_Eg)+=pg.e()/GeV;
    (_ptgZoom)+=(pt)/GeV;
    double rap = 0.5*log((pg.e()+pg.z())/(pg.e()-pg.z()));
    (_rapg)+=(rap);
    (_phig)+=pg.phi();
  } else {
    cerr << "Analysis/GammaJetAnalysis: Found no hard photon in event " 
	 << event->number()  << ".\n";
    generator()->log() << "Analysis/GammaJetAnalysis: " 
		       << "Found no hard photon in event " 
		       << event->number()  << ".\n"
		       << *event;    
  }  
}

LorentzRotation GammaJetAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void GammaJetAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
}

void GammaJetAnalysis::analyze(tPPtr) {}

void GammaJetAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void GammaJetAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<GammaJetAnalysis> GammaJetAnalysis::initGammaJetAnalysis;
// Definition of the static class description member.

void GammaJetAnalysis::Init() {

  static ClassDocumentation<GammaJetAnalysis> documentation
    ("There is no documentation for the GammaJetAnalysis class");

}

