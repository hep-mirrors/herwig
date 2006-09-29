// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IdentifiedParticleAnalysis class.
//

#include "IdentifiedParticleAnalysis.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "IdentifiedParticleAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void IdentifiedParticleAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
}

LorentzRotation IdentifiedParticleAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void IdentifiedParticleAnalysis::analyze(const tPVector &) {
  // get the final-state
  tcEventPtr event=generator()->currentEvent();
  tPVector hadrons=event->getFinalState();
  // get the partons 
  tPVector partons=event->primaryCollision()->steps()[0]->getFinalState();
  int flav = getFlavour(partons);
  Energy Emax = 0.5*generator()->maximumCMEnergy();  
  for (tPVector::iterator it = hadrons.begin(); it != hadrons.end(); ++it ) 
    {
      // only looking at charged particles
      if(!(*it)->data().charged()) continue;
      // all particles
      double xp = _shapes->getX((*it)->momentum(), Emax); 
      *_xpa += xp;
      if(abs((*it)->id()) == ParticleID::piplus) 
	{
	  *_pipma += xp;
	  *_pipm += (*it)->momentum().vect().mag()/GeV; 
	} 
      else if(abs((*it)->id()) == ParticleID::Kplus) 
	{
	  *_kpma += xp; 
	  *_kpm += (*it)->momentum().vect().mag()/GeV;
	}
      else if(abs((*it)->id()) == ParticleID::pplus) 
	{
	  *_ppma += xp; 
	  *_ppm += (*it)->momentum().vect().mag()/GeV;
	}
      switch(flav) 
	{
	case 1:
	case 2:
	case 3:
	  *_xpl += xp; 
	  if(abs((*it)->id()) == ParticleID::piplus)     *_pipml += xp;
	  else if(abs((*it)->id()) == ParticleID::Kplus) *_kpml  += xp;
	  else if(abs((*it)->id()) == ParticleID::pplus) *_ppml  += xp;
	  *_udsxp += xp;
	  if (xp > 0) *_udsxip += -log(xp);
	  break;
	case 4: 
	  *_xpc += xp; 
	  if(abs((*it)->id()) == ParticleID::piplus) *_pipmc += xp; 
	  else if(abs((*it)->id()) == ParticleID::Kplus) *_kpmc += xp; 
	  else if(abs((*it)->id()) == ParticleID::pplus) *_ppmc += xp;
	  break;
	case 5:
	  *_xpb += xp; 
	  if(abs((*it)->id()) == ParticleID::piplus) *_pipmb += xp;
	  else if(abs((*it)->id()) == ParticleID::Kplus) *_kpmb += xp;
	  else if(abs((*it)->id()) == ParticleID::pplus) *_ppmb += xp;
	  break;
	default:
	  break;
	}
    }
  // finally for the lambda's
   set<tcPPtr> allparticles;
  StepVector steps = event->primaryCollision()->steps();
  for ( StepVector::const_iterator it = steps.begin()+2;
	it != steps.end(); ++it ) {
    (**it).select(inserter(allparticles), ThePEG::AllSelector());
  }

  for(set<tcPPtr>::const_iterator it = allparticles.begin(); 
      it != allparticles.end(); ++it) {
    if(abs( (*it)->id())== ParticleID::Lambda0)
      *_lpm+= (*it)->momentum().vect().mag()/Emax;
  }
}

void IdentifiedParticleAnalysis::analyze(tPPtr) {}

void IdentifiedParticleAnalysis::persistentOutput(PersistentOStream & os) const {
  os << _shapes;
}

void IdentifiedParticleAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> _shapes;
}

ClassDescription<IdentifiedParticleAnalysis> IdentifiedParticleAnalysis::initIdentifiedParticleAnalysis;
// Definition of the static class description member.

void IdentifiedParticleAnalysis::Init() {

  static ClassDocumentation<IdentifiedParticleAnalysis> documentation
    ("The IdentifiedParticleAnalysis class compares identified particle spectra with Z"
     " pole data");

  static Reference<IdentifiedParticleAnalysis,EventShapes> interfaceEventShapes
    ("EventShapes",
     "Pointer to the object which calculates the event shapes",
     &IdentifiedParticleAnalysis::_shapes, false, false, true, false, false);
}

