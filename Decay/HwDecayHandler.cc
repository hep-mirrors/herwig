// -*- C++ -*-
//
// HwDecayHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwDecayHandler class.
//

#include "HwDecayHandler.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "DecayIntegrator.h"
#include "DecayPhaseSpaceMode.h"
#include "ThePEG/PDT/MixedParticleData.h"
#include "Herwig/Utilities/EnumParticles.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void HwDecayHandler::
handle(EventHandler &, const tPVector & tagged,
       const Hint &) {
  // First go through the tagged particles for unstable ones
  tPVector parents;
  for(int i = 0, N = tagged.size(); i<N; ++i) {
    if(tagged[i]) {
      // add to parents if not stable
      if(!tagged[i]->data().stable() &&
	 tagged[i]->data().id() != ParticleID::Remnant &&
	 _excluded.find( tagged[i]->dataPtr() ) == _excluded.end() ) {
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
  useMe();
  // Create a new step, decay all particles and add their children to the step
  StepPtr newstep = _newstep ? newStep() : currentStep();
  for(int i = 0, N = parents.size(); i<N; ++i) {
    performDecay(newstep->find(parents[i]), *newstep);
  }
}

// perform decay method including modifications for spin correlations
// and for the decayer to specify intermediate decay products
void HwDecayHandler::performDecay(tPPtr parent, Step & s) const {
  long ntry = 0;
  tcSpinPtr hwspin;
  tcMixedParticleDataPtr 
    mixdata=dynamic_ptr_cast<tcMixedParticleDataPtr>(parent->dataPtr());
  if(mixdata) {
    pair<bool,Length> mixing = mixdata->generateLifeTime();
    develop(parent);
    parent->setLifeLength(Distance());
    PPtr newparent;
    if(mixing.first) {
      newparent = parent->dataPtr()->CC()->
	produceParticle(parent->momentum());
    }
    else {
      newparent = parent->dataPtr()      ->
	produceParticle(parent->momentum());
    }
    newparent->setLabVertex(parent->labDecayVertex());
    Lorentz5Distance lifeLength(mixing.second,
				parent->momentum().vect()*
				(mixing.second/parent->mass()));
    newparent->setLifeLength(lifeLength);
    s.addDecayProduct(parent, newparent);
    parent = newparent;
  }
  else if ( maxLifeTime() >= ZERO ) {
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
      parent->scale(ZERO);
      // loop over the children
      for ( int i = 0, N = children.size(); i < N; ++i ) {
	// if the child has already been decayed add products to the record
	if(children[i]->decayed()) addDecayedParticle(children[i],s);
	// if not stable decay the child
	else if (!children[i]->data().stable() &&
		 _excluded.find( children[i]->dataPtr() ) == _excluded.end() ) {
	  performDecay(children[i], s);
	}
	// if stable and has spinInfo set up decay matrices etc.
	else {
	  develop(children[i]);
	}
      }
      // sort out the spinInfo for the parent after the decays
      if(parent->spinInfo()) parent->spinInfo()->develop();
      return;
    }
    catch (Veto) 
      {}
  }
}

// method to add an intermediate which has already been decayed to the event record
void HwDecayHandler::addDecayedParticle(tPPtr parent, Step & s) const {
  for ( int i = 0, N = parent->children().size(); i < N; ++i ) {
    parent->children()[i]->setLabVertex(parent->labDecayVertex());
    s.addDecayProduct(parent->children()[i]);
  }
  parent->scale(ZERO);
  for ( int i = 0, N = parent->children().size(); i < N; ++i ) {
    if((parent->children()[i])->decayed()) {
      for(unsigned int ix=0;ix<(parent->children()[i])->children().size();++ix)
	addDecayedParticle(parent->children()[i],s);
    }
    else if ( ! parent->children()[i]->data().stable() &&
	      _excluded.find( parent->children()[i]->dataPtr() ) == _excluded.end() ) {
      performDecay(parent->children()[i], s);
    }
    else {
      develop(parent->children()[i]);
    }
  }
  return;
}

void HwDecayHandler::persistentOutput(PersistentOStream & os) const {
  os << _newstep << _excluded << _excludedVector;
}

void HwDecayHandler::persistentInput(PersistentIStream & is, int)  {
  is >> _newstep >> _excluded >> _excludedVector;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HwDecayHandler,DecayHandler>
describeHerwigHwDecayHandler("Herwig::HwDecayHandler", "Herwig.so");

void HwDecayHandler::Init() {

  static ClassDocumentation<HwDecayHandler> documentation
    ("This is the handler for decays in Herwig.",
     "Decays in Herwig include full spin correlations, based on \\cite{Richardson:2001df}.",
     "%\\cite{Richardson:2001df}\n"
     "\\bibitem{Richardson:2001df}\n"
     "  P.~Richardson,\n"
     "  ``Spin correlations in Monte Carlo simulations,''\n"
     "  JHEP {\\bf 0111}, 029 (2001)\n"
     "  [arXiv:hep-ph/0110108].\n"
     "  %%CITATION = JHEPA,0111,029;%%\n"
     );

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

  static RefVector<HwDecayHandler,ParticleData> interfaceExcluded
    ("Excluded",
     "Particles which should not be decayed",
     &HwDecayHandler::_excludedVector, -1, false, false, true, false, false);

}

void HwDecayHandler::doinit() {
  DecayHandler::doinit();
  _excluded = set<tcPDPtr>(_excludedVector.begin(),_excludedVector.end());
}
