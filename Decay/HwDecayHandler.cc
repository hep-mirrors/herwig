// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwDecayHandler class.
//

#include "HwDecayHandler.h"
#include "ThePEG/Handlers/CollisionHandler.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Timer.h"

using namespace ThePEG;
using namespace Herwig;

HwDecayHandler::~HwDecayHandler() {}

void HwDecayHandler::
handle(PartialCollisionHandler & ch, const tPVector & tagged,
       const Hint & hint) ThePEG_THROW_SPEC((Veto, Stop, Exception)) {
  // First go through the tagged particles for unstable ones
  Timer<46> timer("HwDecayHandler::handle");
  tPVector parents;
  for(int i = 0, N = tagged.size(); i<N; ++i)
    if(tagged[i] && !tagged[i]->data().stable())
      parents.push_back(tagged[i]);
  if(parents.empty()) return;

  // Create a new step, decay all particles and add their children
  // to the step
  tStepPtr newStep = ch.newStep();
  for(int i = 0, N = parents.size(); i<N; ++i)
    performDecay(newStep->find(parents[i]), *newStep);

  //cout << "Now checking for coloured particles\n";
  // Now lets get final state particles
  tPVector finalS = newStep->getFinalState();
  // and see if any are partonic (coloured)
  for(tPVector::iterator it = finalS.begin(); it!=finalS.end(); it++) {
    if((*it)->data().coloured()) {
      //cout <<  "Found some\n";
       // If coloured, add a new hadronization step to handle them
       ch.addStep(Group::main, Group::hadron, StepHdlPtr(), Hint::Default());
       break;
    }
  }
}
/*
void HwDecayHandler::
performDecay(tPPtr parent, Step & s) const
  ThePEG_THROW_SPEC((Veto, Exception)) {
  Timer<47> timer("HwDecayHandler::performDecay");
  long ntry = 0;
  while ( 1 ) {
    if ( ++ntry >= maxLoop() )
      throw DecHdlDecayFailed(parent->data(), maxLoop());
    tDMPtr dm = parent->data().selectMode(*parent);
    if ( !dm ) throw DecHdlNoDecayMode(parent->data());
    if ( !dm->decayer() ) throw DecHdlNoDecayer(parent->data(), *dm);
    try {
      ParticleVector children = dm->decayer()->decay(*dm, *parent);
      if(!children.empty() ) {
	parent->decayMode(dm);
	for(int i = 0, N = children.size(); i < N; ++i )
	  if(!s.addDecayProduct(parent, children[i]) )
	    throw DecHdlChildFail(parent->data(), children[i]->data());
	parent->scale(0.0*GeV2);
	for ( int i = 0, N = children.size(); i < N; ++i )
	  if(!children[i]->data().stable() ) performDecay(children[i], s);
	return;
      }
    }
    catch (DecHdlChildFail) {
      throw;
    }
    catch (Veto) {}
  }
}
*/

void HwDecayHandler::persistentOutput(PersistentOStream & os) const {}
void HwDecayHandler::persistentInput(PersistentIStream & is, int) {}
ClassDescription<HwDecayHandler> HwDecayHandler::initHwDecayHandler;
void HwDecayHandler::Init() {}

IBPtr HwDecayHandler::clone() const { return new_ptr(*this); }

IBPtr HwDecayHandler::fullclone() const { return new_ptr(*this); }

void HwDecayHandler::doupdate() throw(UpdateException) {
  StepHandler::doupdate();
}

void HwDecayHandler::doinit() throw(InitException) { StepHandler::doinit(); }

void HwDecayHandler::dofinish() { StepHandler::dofinish(); }

void HwDecayHandler::rebind(const TranslationMap & trans)
   throw(RebindException) {
  StepHandler::rebind(trans);
}

IVector HwDecayHandler::getReferences() {
  IVector ret = StepHandler::getReferences();
  return ret;
}

/*
DecHdlNoDecayMode::DecHdlNoDecayMode(const Interfaced & p) {
  theMessage << "An error occurred when trying to decay an unstable particle "
	     << "of type '" << p.name() << "'. No decay mode was found.";
  severity(runerror);
}

DecHdlNoDecayer::DecHdlNoDecayer(const Interfaced & p, const Interfaced & dm) {
  theMessage << "An error occurred when tryin to decay an unstable particle "
	     << "of type '" << p.name() << "'. The selected decay mode ("
	     << dm.name() << ") did not have a decayer associated with it.";
  severity(runerror);
}

DecHdlDecayFailed::DecHdlDecayFailed(const Interfaced & p, long n) {
  theMessage << "A possibly infinit loop was encountered while tryin to "
	     << "decay an unstable particle of type '" << p.name()
	     << "'. No acceptable decay mode found."
	     << "before reaching the limit of " << n << "iterations.";
  severity(eventerror);
}
  
DecHdlChildFail::DecHdlChildFail(const Interfaced & p, const Interfaced & c) {
  theMessage << "An error occurred when tryin to decay an unstable particle "
	     << "of type '" << p.name() << "'. One of the produced children "
	     << "(of type '" << c.name() << "') could not be added to the "
	     << "current step.";
  severity(abortnow);
}
*/
