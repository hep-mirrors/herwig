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
