// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VGammaHardGenerator class.
//

#include "VGammaHardGenerator.h"
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
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

VGammaHardGenerator::VGammaHardGenerator() {}

IBPtr VGammaHardGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr VGammaHardGenerator::fullclone() const {
  return new_ptr(*this);
}

void VGammaHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS;
}

void VGammaHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS;
}

ClassDescription<VGammaHardGenerator> VGammaHardGenerator::initVGammaHardGenerator;
// Definition of the static class description member.

void VGammaHardGenerator::Init() {

  static ClassDocumentation<VGammaHardGenerator> documentation
    ("The VGammaHardGenerator class implements the generation of the "
     "hard QCD radiation in electroweak vector boson+photon processes");

  static Reference<VGammaHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &VGammaHardGenerator::_alphaS, false, false, true, false, false);
}

bool VGammaHardGenerator::canHandle(ShowerTreePtr tree) {
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // should be a quark and an antiquark
  unsigned int ix(0);
  ShowerParticlePtr part[2];
  // extract the incoming particles
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    part[ix]=cit->first->progenitor();
    ++ix;
  }
  // check incoming quark and antiquark
  if(!(QuarkMatcher::Check(part[0]->data())&&QuarkMatcher::Check(part[1]->data())&&
       double(part[0]->id())/double(part[1]->id())<0.)) return false;
  // two outgoing particles
  if(tree->outgoingLines().size()!=2) return false;
  // extract the outgoing particles
  ix=0;  
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    part[ix]=cjt->first->progenitor();
    ++ix;
  }
  // put photon first
  if(part[0]->id()!=ParticleID::gamma) swap(part[0],part[1]);
  // check photon
  if(part[0]->id()!=ParticleID::gamma) return false;
  // check gauge boson
  if(abs(part[1]->id())==ParticleID::Wplus ||
     part[1]->id() ==ParticleID::Z0) return true;
  else return false;
}

HardTreePtr VGammaHardGenerator::generateHardest(ShowerTreePtr tree) {
  // get the particles to be showered
  _beams.clear();
  _partons.clear();
  // find the incoming particles
  ShowerParticleVector incoming;
  _quarkplus = true;
  vector<ShowerProgenitorPtr> particlesToShower;
  //progenitor particles are produced in z direction.
  for( map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	 cit = tree->incomingLines().begin(); 
       cit != tree->incomingLines().end(); ++cit ) {
    incoming.push_back( cit->first->progenitor() );
    _beams.push_back( cit->first->beam() );
    _partons.push_back( cit->first->progenitor()->dataPtr() );
    // check that quark is along +ve z direction
    if(cit->first->progenitor()->id() > 0 && 
       cit->first->progenitor()->momentum().z() < ZERO ) 
      _quarkplus = false;
    particlesToShower.push_back( cit->first );
  }
  // find the outgoing particles
  tShowerParticlePtr boson,photon;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
 	cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(cjt->first->progenitor()->id()==ParticleID::gamma)
      photon = cjt->first->progenitor();
    else if(abs(cjt->first->progenitor()->id())==ParticleID::Wplus||
	    cjt->first->progenitor()->id()==ParticleID::Z0)
      boson  =  cjt->first->progenitor();
  }
  assert(boson&&photon);
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  if(_partons[0]->id()<_partons[1]->id()) {
    swap(_partons[0],_partons[1]);
    swap(_beams[0],_beams[1]);
  }
  cerr << "testing in generator hardest\n";
  cerr << *boson << "\n"<< *photon << "\n";
  return HardTreePtr();
}
