// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BranchingRatioReweighter class.
//

#include "BranchingRatioReweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "Herwig/Utilities/EnumParticles.h"

using namespace Herwig;

BranchingRatioReweighter::BranchingRatioReweighter() {}

BranchingRatioReweighter::~BranchingRatioReweighter() {}

void BranchingRatioReweighter::
handle(EventHandler & eh,const tPVector & ,const Hint & ) {
  tEventPtr event = eh.currentEvent();
  // weight
  double weight = 1.;
  // get all the particles
  set<tcPPtr> particles;
  event->select(inserter(particles),ThePEG::AllSelector());
  for(set<tcPPtr>::const_iterator it=particles.begin();it!=particles.end();++it) {
    // skip stable
    if((**it).dataPtr()->stable()) continue;
    // skip remnant and clusters
    if((**it).id()==ParticleID::Remnant || 
       (**it).id()==ParticleID::Cluster) continue;
    // if spacelike skip
    if((**it).mass()<ZERO) continue;
    if(*it == event->incoming().first || 
       *it == event->incoming().second ) continue;
    // find unique particles
    bool unique = true;
    for(unsigned int ix=0;ix<(**it).children().size();++ix) {
      if((**it).children()[ix]->id()==(**it).id()) {
	unique = false;
	break;
      }
    }
    if(!unique) continue;
    weight *= (**it).dataPtr()->decaySelector().sum();
  }
  // do the reweighting
  if ( dynamic_cast<StandardEventHandler*>(&eh) ) {
    StandardEventHandler& seh = 
      dynamic_cast<StandardEventHandler&>(eh);
    seh.reweight(weight);
  }
}

IBPtr BranchingRatioReweighter::clone() const {
  return new_ptr(*this);
}

IBPtr BranchingRatioReweighter::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type description system in ThePEG.
DescribeNoPIOClass<BranchingRatioReweighter,StepHandler>
describeHerwigBranchingRatioReweighter("Herwig::BranchingRatioReweighter", "Herwig.so");

void BranchingRatioReweighter::Init() {

  static ClassDocumentation<BranchingRatioReweighter> documentation
    ("The BranchingRatioReweighter class reweights events if some"
     " decay modes are switched off");

}

