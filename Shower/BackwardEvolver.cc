// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BackwardEvolver class.
//

#include "BackwardEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
// #include "ThePEG/Interface/Parameter.h" 
#include "ThePEG/Interface/Reference.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ShowerParticle.h"
#include "ShowerKinematics.h"

using namespace Herwig;


BackwardEvolver::~BackwardEvolver() {}


void BackwardEvolver::persistentOutput(PersistentOStream & os) const {
  os << _splittingGenerator << _forwardEvolver;
}


void BackwardEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _splittingGenerator >> _forwardEvolver;
}


ClassDescription<BackwardEvolver> BackwardEvolver::initBackwardEvolver;
// Definition of the static class description member.

void BackwardEvolver::Init() {

  static ClassDocumentation<BackwardEvolver> documentation
    ("This class is responsible for the backward showering of space-like particles");

  static Reference<BackwardEvolver,SplittingGenerator> 
    interfaceSplitGen("SplittingGenerator", 
		      "A reference to the SplittingGenerator object", 
		      &Herwig::BackwardEvolver::_splittingGenerator,
		      false, false, true, false);
  static Reference<BackwardEvolver,ForwardEvolver> 
    interfaceForwardEvolver("ForwardEvolver", 
			    "A reference to the ForwardEvolver object", 
			    &Herwig::BackwardEvolver::_forwardEvolver,
			    false, false, true, false);

}

//------------------------------------------------------------------------------

bool BackwardEvolver::spaceLikeShower(tPartCollHdlPtr ch, 
				      const tShowerVarsPtr showerVars, 
				      //const tMECorrectionPtr meCorrectionPtr,
				      tShowerParticlePtr particle, 
				      ShowerParticleVector &allShowerParticles)
  throw (Veto, Stop, Exception) {
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "BackwardEvolver::spaceLikeShower "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl;
  }

  bool hasEmitted = false;
  tShowerParticlePtr part = particle;
  tShowerParticleVector particlesYetToShower;   // only time-like particles

  do {

    cout << "-- BackwardEvolver, part = " << part << "." << endl; 
    Branching bb = _splittingGenerator->
      chooseBackwardBranching(ch, *part);
    if(bb.first == ShoKinPtr() || bb.second == tSudakovPtr()) {      
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	generator()->log() << "-- no further backward branching."
			   << endl;
      }
      hasEmitted = false;
    } else {
      hasEmitted = true;
      
      // Assign the splitting function and the shower kinematics
      // to the emitting particle.
      part->setShowerKinematics(bb.first);
      part->setSplittingFn(bb.second->splittingFn()); 
      
      // For the time being we are considering only 1->2 branching
      tSplittingFnPtr splitF = bb.second->splittingFn();
      if(splitF) {	
	short id0, id2;
	id0 = bb.third[0];
	id2 = bb.third[2];
	if (bb.third[2] == ParticleID::g) {
	  id0 = part->id();
	} else if(bb.third[0] == ParticleID::g) {
	  id2 = -part->id();
	}	
	ShowerParticlePtr newParent = new_ptr(
	  ShowerParticle(getParticleData(id0)));
	ShowerParticlePtr otherChild = new_ptr(
          ShowerParticle(getParticleData(id2)));
	cout << "Branching " 
	     << id0 << " -> "
	     << part->id() << ", "
	     << id2 << endl;
	const ShowerIndex::InteractionType interaction = 
	  splitF->interactionType();
	const Energy scale = part->showerKinematics()->qtilde();
	// no z for angular ordering in backward branchings
	newParent->setEvolutionScale(interaction, scale);
	otherChild->setEvolutionScale(interaction, scale);
	
	// for the reconstruction of kinematics, parent/child
	// relationships are according to the branching process:
	// part -> (newParent, otherChild)
	ParticleVector theChildren; 
	theChildren.push_back(newParent); 
	theChildren.push_back(otherChild); 
	part->showerKinematics()->updateChildren(part, theChildren); 
	
	// *** set proper colour connections
	// *** set proper parent/child relationships

	newParent->addChild(part);
 	newParent->addChild(otherChild);
	newParent->x(part->x()/part->showerKinematics()->z());
	cout << "  new x = " 
	     << newParent->x() << endl;
	cout << "  | new parent   = " << newParent << endl;
	cout << "  |-- child 0    = " << newParent->children()[0] << endl;
	cout << "  |-- child 1    = " << newParent->children()[1] << endl;
	tPPtr hadron;
	if (part->parents().size() == 2)
	  hadron = part->parents()[0];
	else cerr << "not one parent!" << endl; 
	hadron->addChild(newParent);
	hadron->removeChild(part);
	part->removeParent(hadron);

	part = newParent;
      } // if ( splitF )
    }   // if (shoKin && sudakov) {...} else {

    cout << "BackwardEvolver." << endl;
  } while(hasEmitted);

  if (part->parents()[0]->id() == 2212) {
    cout << "Parent is proton!" << endl;
    if (part->id() != ParticleID::u &&
	part->id() != ParticleID::d &&
	part->id() != ParticleID::g) {
      cout << "particle is non-valence, forced splitting to gluon." << endl;

      // determine some new pair qtilde, z
      // set up new showerKinematics, store evolution variabls, 
      // ->updateChildren
      // create new particles
      // there is only one splitting possible
      // part = newParent
    }
    if (part->id() == ParticleID::g) {
      cout << "splitting gluon to valence." << endl;
      // basically same as above 
      // but select flavour u or d
    } 
  }

  // set remnant. 

  // do timelike evolution of new particles here? 
  // I think not before 1st reconstruction! 
//   while(!particlesYetToShower.empty()) {
//     tShowerParticlePtr part = particlesYetToShower.back();
//     particlesYetToShower.pop_back();
//     hasEmitted = hasEmitted || 
//      _forwardEvolver->timeLikeShower(ch, showerVars, part, allShowerParticles);
//   } 

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "BackwardEvolver::spaceLikeShower "
		       << " ===> END DEBUGGING <=== "
		       << endl;
  }
  return hasEmitted;
}


