// -*- C++ -*-
//
// HwDecayHandler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwDecayHandler class.
//

#include "HwDecayHandler.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "DecayIntegrator.h"
#include "DecayPhaseSpaceMode.h"
#include "Herwig++/Utilities/EnumParticles.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void HwDecayHandler::
handle(EventHandler &, const tPVector & tagged,
       const Hint &) throw(Veto, Stop, Exception) {
  // First go through the tagged particles for unstable ones
  tPVector parents;
  for(int i = 0, N = tagged.size(); i<N; ++i) {
    if(tagged[i]) {
      // add to parents if not stable
      if(!tagged[i]->data().stable() &&
	 tagged[i]->data().id() != ExtraParticleID::Remnant ) {
	parents.push_back(tagged[i]);
      }
      // if stable and has spinInfo set the developed flag
      else {
	develop(tagged[i]);
      }
    }
  }
  // if nothing to be decayed return
  if(parents.empty()) return;
  // Create a new step, decay all particles and add their children to the step
  StepPtr newstep = _newstep ? newStep() : currentStep();
  for(int i = 0, N = parents.size(); i<N; ++i) {
    performDecay(newstep->find(parents[i]), *newstep);
  }
}


// perform decay method including modifications for spin correlations
// and for the decayer to specify intermediate decay products
void HwDecayHandler::performDecay(tPPtr parent, Step & s) const
  throw(Veto, Exception) {
  long ntry = 0;
  tcSpinfoPtr hwspin;
  if ( maxLifeTime() >= 0.0*mm ) {
    if( ( lifeTimeOption() && parent->lifeLength().tau() > maxLifeTime())||
	(!lifeTimeOption() && parent->data().cTau()      > maxLifeTime()) ) {
      parent->setLifeLength(Distance());
      develop(parent);
      return;
    }
  }
  while ( true ) {
    // exit if fails
    if ( ++ntry >= maxLoop() ) 
      throw Exception() << "Too many tries " << maxLoop() << "to generate decay of "
			<< *parent << "in "
			<< "HwDecayHandler::performDecay" << Exception::eventerror;
    // select the decay mode
    tDMPtr   dm(parent->data().selectMode(*parent));
    // check we found a decay mode and it had a decayer
    if ( !dm ) {
      generator()->log() << *generator()->currentEvent() << "\n";
      generator()->log() << *parent << "\n";
      throw Exception() << "No DecayModes for " << parent->PDGName()
			<< " in HwDecayHandler::performDecay" 
			<< Exception::eventerror;
    }
    if ( !dm->decayer() ) throw Exception() << "No decayer for DecayMode of " 
					    << parent->PDGName()
					    << " in HwDecayHandler::performDecay" 
					    << Exception::eventerror;
    try {
      ParticleVector children = dm->decayer()->decay(*dm, *parent);
      if(children.empty()) continue;
      assert(parent->children().empty());
      // generate radiation in the decay
      tDecayIntegratorPtr hwdec=dynamic_ptr_cast<tDecayIntegratorPtr>(dm->decayer());
      if (hwdec && hwdec->canGeneratePhotons())
	children = hwdec->generatePhotons(*parent,children);
      // set up parent
      parent->decayMode(dm);
      // add children
      for ( int i = 0, N = children.size(); i < N; ++i ) {
	children[i]->setLabVertex(parent->labDecayVertex());
	if ( !s.addDecayProduct(parent, children[i]) ) 
	  throw Exception() << "Failed to add child " 
			    << children[i]->PDGName() 
			    << " in decay of " << parent->PDGName() 
			    << Exception::eventerror;
      }
      parent->scale(0.0*MeV2);
      // loop over the children
      for ( int i = 0, N = children.size(); i < N; ++i ) {
	// if the child has already been decayed add products to the record
	if(children[i]->decayed()) addDecayedParticle(children[i],s);
	// if not stable decay the child
	else if (!children[i]->data().stable()) {
	  performDecay(children[i], s);
	}
	// if stable and has spinInfo set up decay matrices etc.
	else {
	  develop(children[i]);
	}
      }
      // sort out the spinInfo for the parent after the decays
      if(parent->spinInfo()) {
	hwspin=dynamic_ptr_cast<tcSpinfoPtr>(parent->spinInfo());
	// if the parent has the right kind of spinInfo
	if(hwspin) {
	  // if the parent has been given a decay vertex
	  // calculate the decay matrix for the decay
	  if(hwspin->getDecayVertex()) hwspin->develop();
	  // if the particle was scalar then it doesn't matter that it
	  // doesn't have a decay vertex as there's no correlations
	  else if(hwspin->iSpin()==PDT::Spin0) hwspin->setDeveloped(true);
	}
      }
      return;
    }
    catch (Veto) 
      {}
  }
}

// method to add an intermediate which has already been decayed to the event record
void HwDecayHandler::addDecayedParticle(tPPtr parent, Step & s) const
  throw(Veto, Exception) 
{
  for ( int i = 0, N = parent->children().size(); i < N; ++i ) {
    parent->children()[i]->setLabVertex(parent->labDecayVertex());
    s.addDecayProduct(parent->children()[i]);
  }
  parent->scale(0.0*GeV2);
  for ( int i = 0, N = parent->children().size(); i < N; ++i ) {
    if((parent->children()[i])->decayed()) {
      for(unsigned int ix=0;ix<(parent->children()[i])->children().size();++ix)
	addDecayedParticle(parent->children()[i],s);
    }
    else if ( !(parent->children()[i])->data().stable() ) {
      performDecay(parent->children()[i], s);
    }
    else if(parent->children()[i]->data().stable()) {
      develop(parent->children()[i]);
    }
  }
  return;
}

void HwDecayHandler::persistentOutput(PersistentOStream & os) const {
  os << _newstep;
}

void HwDecayHandler::persistentInput(PersistentIStream & is, int)  {
  is >> _newstep;
}

ClassDescription<HwDecayHandler> HwDecayHandler::initHwDecayHandler;

void HwDecayHandler::Init() {

  static ClassDocumentation<HwDecayHandler> documentation
    ("This is the handler for decays in Herwig++.");

  static Switch<HwDecayHandler,bool> interfaceNewStep
    ("NewStep",
     "Add the particles in a new step",
     &HwDecayHandler::_newstep, true, false, false);
  static SwitchOption interfaceNewStepNew
    (interfaceNewStep,
     "Yes",
     "Add particles in a new step",
     true);
  static SwitchOption interfaceNewStepCurrent
    (interfaceNewStep,
     "No",
     "Add them in the current step",
     false);

}
