// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerHandler class.
//

#include "ShowerHandler.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Handlers/XComb.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include "Herwig++/Hadronization/Remnant.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include <cassert>

using namespace Herwig;

ShowerHandler::~ShowerHandler() {}

IBPtr ShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr ShowerHandler::fullclone() const {
  return new_ptr(*this);
}

ShowerHandler::ShowerHandler() : _maxtry(10) {
  _inputparticlesDecayInShower.push_back( 6 ); //  top
  _inputparticlesDecayInShower.push_back( 1000001 ); //  SUSY_d_L 
  _inputparticlesDecayInShower.push_back( 1000002 ); //  SUSY_u_L 
  _inputparticlesDecayInShower.push_back( 1000003 ); //  SUSY_s_L 
  _inputparticlesDecayInShower.push_back( 1000004 ); //  SUSY_c_L 
  _inputparticlesDecayInShower.push_back( 1000005 ); //  SUSY_b_1 
  _inputparticlesDecayInShower.push_back( 1000006 ); //  SUSY_t_1 
  _inputparticlesDecayInShower.push_back( 1000011 ); //  SUSY_e_Lminus 
  _inputparticlesDecayInShower.push_back( 1000012 ); //  SUSY_nu_eL 
  _inputparticlesDecayInShower.push_back( 1000013 ); //  SUSY_mu_Lminus 
  _inputparticlesDecayInShower.push_back( 1000014 ); //  SUSY_nu_muL 
  _inputparticlesDecayInShower.push_back( 1000015 ); //  SUSY_tau_1minus 
  _inputparticlesDecayInShower.push_back( 1000016 ); //  SUSY_nu_tauL 
  _inputparticlesDecayInShower.push_back( 1000021 ); //  SUSY_g 
  _inputparticlesDecayInShower.push_back( 1000022 ); //  SUSY_chi_10 
  _inputparticlesDecayInShower.push_back( 1000023 ); //  SUSY_chi_20 
  _inputparticlesDecayInShower.push_back( 1000024 ); //  SUSY_chi_1plus 
  _inputparticlesDecayInShower.push_back( 1000025 ); //  SUSY_chi_30 
  _inputparticlesDecayInShower.push_back( 1000035 ); //  SUSY_chi_40 
  _inputparticlesDecayInShower.push_back( 1000037 ); //  SUSY_chi_2plus 
  _inputparticlesDecayInShower.push_back( 1000039 ); //  SUSY_gravitino 
  _inputparticlesDecayInShower.push_back( 2000001 ); //  SUSY_d_R 
  _inputparticlesDecayInShower.push_back( 2000002 ); //  SUSY_u_R 
  _inputparticlesDecayInShower.push_back( 2000003 ); //  SUSY_s_R 
  _inputparticlesDecayInShower.push_back( 2000004 ); //  SUSY_c_R 
  _inputparticlesDecayInShower.push_back( 2000005 ); //  SUSY_b_2 
  _inputparticlesDecayInShower.push_back( 2000006 ); //  SUSY_t_2 
  _inputparticlesDecayInShower.push_back( 2000011 ); //  SUSY_e_Rminus 
  _inputparticlesDecayInShower.push_back( 2000012 ); //  SUSY_nu_eR 
  _inputparticlesDecayInShower.push_back( 2000013 ); //  SUSY_mu_Rminus 
  _inputparticlesDecayInShower.push_back( 2000014 ); //  SUSY_nu_muR 
  _inputparticlesDecayInShower.push_back( 2000015 ); //  SUSY_tau_2minus 
  _inputparticlesDecayInShower.push_back( 2000016 ); //  SUSY_nu_tauR 
  _inputparticlesDecayInShower.push_back( 25      ); //  h0
  _inputparticlesDecayInShower.push_back( 35      ); //  H0
  _inputparticlesDecayInShower.push_back( 36      ); //  A0
  _inputparticlesDecayInShower.push_back( 37      ); //  H+
  _inputparticlesDecayInShower.push_back( 23      ); // Z0
  _inputparticlesDecayInShower.push_back( 24      ); // W+/-
}

void ShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << _evolver << _maxtry << _inputparticlesDecayInShower
     << _particlesDecayInShower;
}

void ShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> _evolver >> _maxtry
     >> _inputparticlesDecayInShower
     >> _particlesDecayInShower;  
}

ClassDescription<ShowerHandler> ShowerHandler::initShowerHandler;
// Definition of the static class description member.

void ShowerHandler::Init() {

  static ClassDocumentation<ShowerHandler> documentation
    ("Main driver class for the showering.");

  static Reference<ShowerHandler,Evolver> 
    interfaceEvolver("Evolver", 
		     "A reference to the Evolver object", 
		     &Herwig::ShowerHandler::_evolver,
		     false, false, true, false);

  static Parameter<ShowerHandler,unsigned int> interfaceMaxTry
    ("MaxTry",
     "The maximum number of attempts for the main showering loop",
     &ShowerHandler::_maxtry, 10, 1, 100,
     false, false, Interface::limited);

  static ParVector<ShowerHandler,long> interfaceDecayInShower
    ("DecayInShower",
     "PDG codes of the particles to be decayed in the shower",
     &ShowerHandler::_inputparticlesDecayInShower, -1, 0l, -10000000l, 10000000l,
     false, false, Interface::limited);

}

void ShowerHandler::fillEventRecord() {
  // create a new step 
  StepPtr pstep = newStep();
  if(_done.empty()) throw Exception() << "Must have some showers to insert in "
				      << "ShowerHandler::fillEventRecord()" 
				      << Exception::runerror;
  if(!_done[0]->isHard()) throw Exception() << "Must start filling with hard process"
					    << " in ShowerHandler::fillEventRecord()" 
					    << Exception::runerror;
  // insert the steps
  for(unsigned int ix=0;ix<_done.size();++ix) {
    _done[ix]->fillEventRecord(pstep,
			       _evolver->isISRadiationON(),
			       _evolver->isFSRadiationON());
  }
} 

void ShowerHandler::findShoweringParticles() {
  // clear the storage
  _hard=ShowerTreePtr();
  _decay.clear();
  _done.clear();
  // temporary storage of the particles
  set<PPtr> hardParticles;
  // outgoing particles from the hard process
  ParticleVector outgoing=eventHandler()->currentCollision()->
    primarySubProcess()->outgoing();
  set<PPtr> outgoingset(outgoing.begin(),outgoing.end());
  // loop over the tagged particles
  tParticleVector::const_iterator taggedP = tagged().begin();
  bool isHard=false;
  for (;taggedP != tagged().end(); ++taggedP) {
    // if a remnant don't consider
    if(eventHandler()->currentCollision()->isRemnant(*taggedP))
      continue;
    // find the parent and if colourless s-channel resonance
    bool isDecayProd=false;
    tPPtr parent;
    if(!(*taggedP)->parents().empty()) {
      parent = (*taggedP)->parents()[0];
      // check if from s channel decaying colourless particle
      isDecayProd = decayProduct(parent);
    }
    // add to list of outgoing hard particles if needed
    isHard |=(outgoingset.find(*taggedP) != outgoingset.end());
    if(isDecayProd) hardParticles.insert(findParent(parent,isHard,outgoingset));
    else            hardParticles.insert(*taggedP);
  }
  // there must be something to shower
  if(hardParticles.empty()) 
    throw Exception() << "No particles to shower in "
		      << "ShowerHandler::fillShoweringParticles" 
		      << Exception::eventerror;
  if(!isHard)
    throw Exception() << "Starting on decay not yet implemented in "
		      << "ShowerHandler::findShoweringParticles()" 
		      << Exception::runerror;
  // create the hard process ShowerTree
  ParticleVector out(hardParticles.begin(),hardParticles.end());
  _hard=new_ptr(ShowerTree(eventHandler(),out,this,_decay));
  _hard->setParents();
}

void ShowerHandler::cascade() {
  // set the current step
  _current=currentStep();
  //  start of the try block for the whole showering process
  unsigned int countFailures=0;
  ShowerTreePtr hard;
  vector<ShowerTreePtr> decay;
  while (countFailures<_maxtry) {
    try {
      // find the particles in the hard process and the decayed particles to shower
      findShoweringParticles();
      // check if a hard process or decay
      bool isHard = _hard;
      // if a hard process perform the shower for the hard process
      if(isHard) {
	_evolver->showerHardProcess(_hard);
	_done.push_back(_hard);
	_hard->updateAfterShower(_decay,eventHandler());
      }
      // if no decaying particles to shower break out of the loop
      if(_decay.empty()) break;
      // if no hard process
      if(!isHard) 
	throw Exception() << "Shower starting with a decay is not yet implemented" 
			  << Exception::runerror;
      // shower the decay products
      while(!_decay.empty()) {
	multimap<Energy,ShowerTreePtr>::iterator dit=--_decay.end();
	while(!dit->second->parent()->hasShowered() && dit!=_decay.begin()) --dit;
	// get the particle and the width
	ShowerTreePtr decayingTree = dit->second;
	// 	    Energy largestWidthDecayingSystem=(*_decay.rbegin()).first;
	// remove it from the multimap
	_decay.erase(dit);
	// make sure the particle has been decayed
	decayingTree->decay(_decay,eventHandler());
	// now shower the decay
	_evolver->showerDecay(decayingTree);
	_done.push_back(decayingTree);
	decayingTree->updateAfterShower(_decay,eventHandler());
      }
      // suceeded break out of the loop
      break;
    }
    catch (Veto) {
      throw Exception() << "Problem with throwing Veto in ShowerHandler at the moment"
			<< Exception::eventerror;
      ++countFailures;
    }
  }
  // if loop exited because of too many tries, throw event away
  if (countFailures >= _maxtry) {
    throw Exception() << "Too many tries for main while loop "
		      << "in ShowerHandler::cascade()." 
		      << Exception::eventerror; 	
  }
  //enter the particles in the event record
  fillEventRecord();
  // remake the remnants (needs to be after the colours are sorted
  //                       out in the insertion into the event record)
  makeRemnants();
}

void ShowerHandler::makeRemnants() {
  // get the incoming particles
  PPair incoming=generator()->currentEvent()->incoming();
  ParticleVector in;
  in.push_back(incoming.first);
  in.push_back(incoming.second);
  // get the remnants
  tParticleSet remn=generator()->currentEvent()->primaryCollision()->getRemnants();
  // fix the momenta
  for(unsigned int ix=0;ix<in.size();++ix) {
    Lorentz5Momentum pnew;
    ParticleVector prem,pother;
    if(in[ix]->children().size()==1) continue;
    for(unsigned int iy=0;iy<in[ix]->children().size();++iy) {
      if(remn.find(in[ix]->children()[iy])==remn.end()) {
	pnew+=in[ix]->children()[iy]->momentum();
	pother.push_back(in[ix]->children()[iy]);
      }
      else
	prem.push_back(in[ix]->children()[iy]);
    }
    pnew=in[ix]->momentum()-pnew;
    pnew.rescaleMass();
    // throw exception if gone wrong
    if(prem.size()!=1||pother.size()!=1) 
      throw Exception() 
	<< "Must be one and only 1 remnant for beam in ShowerHandler::makeRemnants()"
	<< Exception::eventerror;
    // remake the remnant
    if(prem[0]->id()==ExtraParticleID::Remnant) {
      tRemnantPtr rem=dynamic_ptr_cast<tRemnantPtr>(prem[0]);
      if(rem) rem->regenerate(pother[0],pnew);
    }
  }
}

PPtr ShowerHandler::findParent(PPtr original, bool & isHard, 
			       set<PPtr> outgoingset) const {
  PPtr parent=original;
  isHard |=(outgoingset.find(original) != outgoingset.end());
  if(!original->parents().empty()) {
    PPtr orig=original->parents()[0];
    if(_current->find(orig)&&decayProduct(orig)) {
      parent=findParent(orig,isHard,outgoingset);
    }
  }
  return parent;
}

bool ShowerHandler::decayProduct(tPPtr particle) const{
  return 
    !(particle->dataPtr()->coloured()&&
      (particle->parents()[0]==eventHandler()->lastPartons().first||
       particle->parents()[0]==eventHandler()->lastPartons().second)) && 
    particle->momentum().m2()>0.0*GeV2&&
    particle != eventHandler()->lastPartons().first &&
    particle != eventHandler()->lastPartons().second;
}
