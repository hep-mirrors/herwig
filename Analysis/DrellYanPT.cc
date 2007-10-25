// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanPT class.
//

#include "DrellYanPT.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DrellYanPT.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DrellYanPT::~DrellYanPT() {}

void DrellYanPT::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator send=event->primaryCollision()->steps().end();
  for(;sit!=send;++sit) {
    ParticleSet part=(**sit).all();
    ParticleSet::const_iterator iter=part.begin();
    ParticleSet::const_iterator end =part.end();
    for( ;iter!=end;++iter) {
      if(((**iter).id()==ParticleID::Z0||(**iter).id()==ParticleID::gamma)
	 && (**iter).children().size()==2) {
	_Zpt.addWeighted((**iter).momentum().perp()/GeV,event->weight());
      } else if ((**iter).id()==ParticleID::Wplus && (**iter).children().size()==2) {
	_Wppt.addWeighted((**iter).momentum().perp()/GeV,event->weight());
	
      } else if ((**iter).id()==ParticleID::Wminus && (**iter).children().size()==2) {
	_Wmpt.addWeighted((**iter).momentum().perp()/GeV,event->weight());
      }
    }
  }
}

LorentzRotation DrellYanPT::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void DrellYanPT::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void DrellYanPT::analyze(tPPtr) {}

void DrellYanPT::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void DrellYanPT::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<DrellYanPT> DrellYanPT::initDrellYanPT;
// Definition of the static class description member.

void DrellYanPT::Init() {

  static ClassDocumentation<DrellYanPT> documentation
    ("Analyses the pt of weak bosons produces in Drell-Yan processes.");


}

