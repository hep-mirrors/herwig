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
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Utilities/Timer.h"
#include "Evolver.h"
#include "Herwig++/Shower/Kinematics/ShowerParticle.h"
#include <cassert>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerHandler.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ShowerTree.h"

using namespace Herwig;

void ShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << _evolver << _maxtry;
}

void ShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> _evolver >> _maxtry;  
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

}

void ShowerHandler::fillEventRecord() 
{
  // create a new step 
  StepPtr pstep = eventHandler()->newStep();
  if(_done.empty()) throw Exception() << "Must have some showers to insert in "
				      << "ShowerHandler::fillEventRecord()" 
				      << Exception::runerror;
  if(!_done[0]->isHard()) throw Exception() << "Must start filling with hard process"
					    << " in ShowerHandler::fillEventRecord()" 
					    << Exception::runerror;
  // insert the steps
  for(unsigned int ix=0;ix<_done.size();++ix)
    {
      _done[ix]->fillEventRecord(pstep,
				 _evolver->isISRadiationON(),
				 _evolver->isFSRadiationON());
    }
} 

void ShowerHandler::findShoweringParticles()
{
  Timer<1001> timer("ShowerHandler::findShoweringParticles");
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
  for (;taggedP != tagged().end(); ++taggedP) {
    // if a remnant don't consider
    if(eventHandler()->currentCollision()->isRemnant(*taggedP))
      continue;
    // find the parent and if colourless s-channel resonance
    bool isDecayProd=false;
    tPPtr parent;
    if(!(*taggedP)->parents().empty()) 
      {
	parent = (*taggedP)->parents()[0];
	// check if from s channel decaying colourless particle
	// (must be same as in findParent)
	isDecayProd = !parent->dataPtr()->coloured() && parent->momentum().m2()>0.;
      }
    // add to list of outgoing hard particles if needed
    if(outgoingset.find(*taggedP) != outgoingset.end())
      {
	if(isDecayProd) hardParticles.insert(findParent(parent));
	else            hardParticles.insert(*taggedP);
      }
    else
      {throw Exception() << "Starting on decay not yet implemented in "
			  << "ShowerHandler::findShoweringParticles()" 
			  << Exception::runerror;}
  }
  // there must be something to shower
  if(hardParticles.empty()) 
    throw Exception() << "No particles to shower in "
		      << "ShowerHandler::fillShoweringParticles" 
		      << Exception::eventerror;
  // create the hard process ShowerTree
  ParticleVector out(hardParticles.begin(),hardParticles.end());
  _hard=new_ptr(ShowerTree(eventHandler(),out,_evolver->showerVariables(),
			   _decay));
  _hard->setParents();
}

void ShowerHandler::cascade()
{
  Timer<1002> timer("ShowerHandler::cascade");
  // should we be doing anything
  if(!_evolver->showeringON()) return;
  //  start of the try block for the whole showering process
  unsigned int countFailures=0;
  ShowerTreePtr hard;
  vector<ShowerTreePtr> decay;
  while (countFailures<_maxtry) {
    try
      {
	// set the gluon mass to be used in the reconstruction
	_evolver->showerVariables()->setGluonMass(false);
	// find the particles in the hard process and the decayed particles to shower
	findShoweringParticles();
	// check if a hard process or decay
 	bool isHard = _hard;
 	// find the stopping scale for the shower if multi-scale shower is on
 	Energy largestWidth = Energy();
 	if(_evolver->showerVariables()->isMultiScaleShowerON()&&!_decay.empty())
 	  {
 	    largestWidth=(*_decay.rbegin()).first;
 	    if(largestWidth<
 	       _evolver->showerVariables()->globalParameters()->hadronizationScale())
 	      largestWidth = Energy();
 	  }
 	// if a hard process perform the shower for the hard process
 	if(isHard) 
	  {
	    _evolver->showerHardProcess(_hard);
	    _done.push_back(_hard);
	  }
 	// if no decaying particles to shower break out of the loop
 	if(_decay.empty()) break;
 	// if no hard process
 	if(!isHard) 
 	  throw Exception() << "Shower starting with a decay is not yet implemented" 
 			    << Exception::runerror;
 	// shower the decay products
 	while(!_decay.empty())
 	  {
	    multimap<Energy,ShowerTreePtr>::iterator dit=--_decay.end();
	    while(!dit->second->parent()->hasShowered() && dit!=_decay.begin())
	      --dit;
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
 	  }
	// suceeded break out of the loop
	break;
      }
    catch (Veto)
      {
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
  _evolver->makeRemnants(_hard);
}
