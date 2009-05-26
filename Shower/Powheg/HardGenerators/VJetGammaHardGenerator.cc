// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VJetGammaHardGenerator class.
//

#include "VJetGammaHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

VJetGammaHardGenerator::VJetGammaHardGenerator() : pTmin_(2.*GeV)
{}

IBPtr VJetGammaHardGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr VJetGammaHardGenerator::fullclone() const {
  return new_ptr(*this);
}

void VJetGammaHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << alphaEM_ << ounit(pTmin_,GeV);
}

void VJetGammaHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> alphaEM_ >> iunit(pTmin_,GeV);
}

ClassDescription<VJetGammaHardGenerator> VJetGammaHardGenerator::initVJetGammaHardGenerator;
// Definition of the static class description member.

void VJetGammaHardGenerator::Init() {

  static ClassDocumentation<VJetGammaHardGenerator> documentation
    ("There is no documentation for the VJetGammaHardGenerator class");


  static Reference<VJetGammaHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the EM coupling constant",
     &VJetGammaHardGenerator::alphaEM_, false, false, true, false, false);

  static Parameter<VJetGammaHardGenerator, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &VJetGammaHardGenerator::pTmin_, GeV, 2.*GeV, ZERO, 100000.0*GeV,
     false, false, Interface::limited);
}

bool VJetGammaHardGenerator::canHandle(ShowerTreePtr tree) {
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // two outgoing particles
  if(tree->outgoingLines().size()!=2) return false;
  ShowerParticlePtr part[4];
  unsigned int ix=0;
  // extract the incoming particles
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    part[ix]=cit->first->progenitor();
    ++ix;
  }
  // extract the outgoing particles
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    part[ix]=cjt->first->progenitor();
    ++ix;
  }
  // check for W/Z
  if(abs(part[2]->id())!=ParticleID::Wplus&&part[2]->id()!=ParticleID::Z0)
    swap(part[2],part[3]);
  if(abs(part[2]->id())!=ParticleID::Wplus&&part[2]->id()!=ParticleID::Z0)
    return false;
  // check others all coloured
  if(!part[0]->dataPtr()->coloured()||!part[1]->dataPtr()->coloured()||
     !part[3]->dataPtr()->coloured()) return false;
  return true;
}


HardTreePtr VJetGammaHardGenerator::generateHardest(ShowerTreePtr) {
  return HardTreePtr();
}
