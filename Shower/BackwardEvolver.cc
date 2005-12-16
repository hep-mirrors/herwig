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
#include "IS_QtildaShowerKinematics1to2.h"
#include "ForcedSplitting.h"

using namespace Herwig;

BackwardEvolver::~BackwardEvolver() {}


void BackwardEvolver::persistentOutput(PersistentOStream & os) const {
  os << _splittingGenerator << _forwardEvolver << _forcedSplitting;
}


void BackwardEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _splittingGenerator >> _forwardEvolver >> _forcedSplitting;
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

  static Reference<BackwardEvolver,ForcedSplitting> 
    interfaceForcedSplitting("ForcedSplitting", 
			    "A reference to the ForcedSplitting object", 
			    &Herwig::BackwardEvolver::_forcedSplitting,
			    false, false, true, false);
  
  
}

//------------------------------------------------------------------------------

int BackwardEvolver::spaceLikeShower(tEHPtr ch, 
				      const tShowerVarsPtr showerVars, 
				      tShowerParticlePtr particle, 
				      ShowerParticleVector &allShowerParticles)
  throw (Veto, Stop, Exception) {
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "BackwardEvolver::spaceLikeShower "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" 
		       << generator()->currentEventNumber() << endl;
  }

  int hasEmitted = 0;
  tShowerParticlePtr part = particle;
  tShowerParticleVector particlesYetToShower;   // only time-like particles
  Energy q0g = (showerVars->kinScale()-0.003*GeV)/2.3;

  do {
    bool cant = false;
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "  BackwardEvolver, part = " 
			 << part << "." << endl; 
    }
    //    cerr << "BackwardEvolver: calling SG->cBB... ";
    Branching bb = _splittingGenerator->chooseBackwardBranching(ch, *part);
    //    cerr << "done" << endl;
    if (bb.first != ShoKinPtr()) {
      double yy = 1.+sqr(q0g/bb.first->qtilde())/2.;
      double zm = yy - sqrt(sqr(yy)-1.); 
      double xp = part->x()/bb.first->z();
      // xp > zm^3 is super-safe enough for yet another emission!
      if (xp > zm*zm*zm) {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log() << "  BE: Can't split again! " << endl;
	}
	cant = true;
      } else {
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	    generator()->log() << "  check: xp = " << xp 
			       << ", zm = " << zm << endl
			       << "  pessimistic check: xp/zm = " 
			       << xp/zm << ", xp/zm^2 = " 
			       << xp/zm/zm <<  endl;
	  }
      }
    }
    if(bb.first == ShoKinPtr() || bb.second == tSudakovPtr() || cant) {      
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log()<< "  --- no further backward branching."
			  << endl;
      }
      hasEmitted = 0;
    } else {
      hasEmitted = 1;
      
      // Assign the splitting function and the shower kinematics
      // to the emitting particle.      
      part->setShowerKinematics(bb.first);
      part->setSplittingFn(bb.second->splittingFn()); 
      
      // check for additional angular ordering...
      static long calls=0, violated=0;
      if(part->children()[0]) {
	calls++;
	ShoKinPtr sk;
	if (dynamic_ptr_cast<ShowerParticlePtr>(part->children()[0])) {
	  sk = dynamic_ptr_cast<ShowerParticlePtr>(part->children()[0])->showerKinematics();
	  if ( bb.first->qtilde() > 
	       (1.-bb.first->z())/(1.-sk->z())*sk->qtilde()) {
	    //	    cout << "Angular Ordering Violated!" << endl;
	    violated++;
	    //	    cerr << double(violated)/double(calls) << " emissions violated ang ordering." << endl;
	  }
	}
      }

      // For the time being we are considering only 1->2 branching
      tSplittingFnPtr splitF = bb.second->splittingFn();
      if(splitF) {	
	short id0, id2;
	id0 = bb.third[0];
	id2 = bb.third[2];
	if(id2 == ParticleID::g) id0 = part->id();
	else if(id0 == ParticleID::g) id2 = -part->id();
	
	// Now create the actual particles, make the otherChild a final state
	// particle, while the newParent is not	
	ShowerParticlePtr newParent = new_ptr(
	  ShowerParticle(getParticleData(id0)));
	ShowerParticlePtr otherChild = new_ptr(
          ShowerParticle(getParticleData(id2)));
	otherChild->setFinalState(true);
	otherChild->setInitiatesTLS(true);
	newParent->setFinalState(false);
	//newParent->setFromHardSubprocess(true);

	// make sure, otherChild is included in TL shower.
	allShowerParticles.insert(allShowerParticles.end(), otherChild);
	allShowerParticles.insert(allShowerParticles.end(), newParent);

	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log() << "  Branching " 
			     << id0 << " -> "
			     << part->id() << ", "
			     << id2 << endl;
	}

	ShowerIndex::InteractionType inter = splitF->interactionType();
	Energy scale = part->showerKinematics()->qtilde();

	// Set up the colour connections and the parent/child relationships
	createBranching(part,newParent,otherChild,scale,inter);
	part = newParent;
      } // if ( splitF )
    }   // if (shoKin && sudakov) {...} else {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "  done one branching." << endl;
    }
  } while(hasEmitted > 0);

  

  /****
   * Now we need to force the final (up to two) splittings so that we are left
   * with only the valence quarks that make up the incoming hadron. If we have
   * terminated the shower on a sea quark, then we need two splittings, one
   * to a gluon and one to a valence quark. If we are on a gluon we just go
   * to a valence quark. We must also choose qtilda and z for each splitting so
   * that the kinematics can be reconstructed properly. This is no longer 
   * sampled according to the splitting functions, as they no longer have space
   * in the virtuality (since the shower has terminated). Instead we use a new
   * distribution.
   * NOTE: temporarily chosen linearly in z and logarithmically in qtilda, this
   * may be changed later.
   ****/
  hasEmitted = _forcedSplitting->split(part,allShowerParticles,ch);

  // Do we veto the whole shower after the final state showering or do we
  // seperately veto the initial state shower and final state shower?

  // do timelike evolution of new particles here? 
  // I think not before 1st reconstruction! 
//   while(!particlesYetToShower.empty()) {
//     tShowerParticlePtr part = particlesYetToShower.back();
//     particlesYetToShower.pop_back();
//     hasEmitted = hasEmitted || 
//       _forwardEvolver->timeLikeShower(ch, showerVars, part, allShowerParticles);
//   } 
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "BackwardEvolver::spaceLikeShower "
		       << " ===> END DEBUGGING <=== "
		       << endl;
  }
  return hasEmitted;
}



void BackwardEvolver::createBranching(ShowerParticlePtr part,
				      ShowerParticlePtr newParent,
				      ShowerParticlePtr otherChild,
				      Energy scale, 
				      ShowerIndex::InteractionType inter) {
  // no z for angular ordering in backward branchings
  newParent->setEvolutionScale(inter, scale);
  otherChild->setEvolutionScale(inter, scale);
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  Creating branching " 
		       << newParent->id() << "-->" << part->id() 
		       << " " << otherChild->id() << endl;
  }
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // part -> (newParent, otherChild)
  ParticleVector theChildren; 
  theChildren.push_back(newParent); 
  theChildren.push_back(otherChild); 
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "Updating children" << endl;
  }
  part->showerKinematics()->updateChildren(part, theChildren); 
  
  // *** set proper colour connections
  setColour(newParent,part,otherChild);
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  Colours set" << endl;
  }
  // *** set proper parent/child relationships
  newParent->addChild(part);
  newParent->addChild(otherChild);
  newParent->x(part->x()/part->showerKinematics()->z());
  
  // Print some messages
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  new x = " 
		       << newParent->x() << endl
		       << "  | new parent   = " << newParent->number() << " ("
		       << newParent << ")" << endl
		       << "  |-- child 0    = " 
		       << newParent->children()[0]->number() 
		       << " (" << newParent->children()[0] << ")" << endl
		       << "  |-- child 1    = " 
		       << newParent->children()[1]->number() 
		       << " (" << newParent->children()[1] << ")" << endl; 
  }
  // Now fix the hadrons connections
  tPPtr hadron;
  if(part->parents().size() == 2) hadron = part->parents()[0];
  else cerr << "Shower/BackwardEvolver::createBranching: not one parent!" 
	    << endl; 
  // If their is a ThePEGBase object, we must remove that from the 
  // hadron, and the extra call to remove the part is solely to remove
  // the hadron as a parent of part
  //if(part->getThePEGBase()) hadron->abandonChild(part->getThePEGBase());
  hadron->abandonChild(part);
  hadron->addChild(newParent);
}

void BackwardEvolver::setColour(ShowerParticlePtr &newParent,
				ShowerParticlePtr &oldParent,
				ShowerParticlePtr &otherChild) {

  /****
   * This method sets the colour connection for a backwards evolving
   * parton. The input is the created particles for the new parent (what
   * the parton came from), the original parton that backwards evolved,
   * and the other child which comes from the branching.
   ****/
  ShoColinePair parent = ShoColinePair();
  ShoColinePair child1 = ShoColinePair(oldParent->colourLine(),
				       oldParent->antiColourLine());
  ShoColinePair child2 = ShoColinePair();

  if(child1.first && child1.second) { // We had a colour octet
    if(newParent->id() == ParticleID::g) {
      if (UseRandom::rndbool()) {
	parent.first = child1.first;
	child2.first = child1.second;
	parent.second = new_ptr(ColourLine());
	child2.second = parent.second;
      } else {
	parent.second = child1.second;
	child2.second = child1.first;
	parent.first = new_ptr(ColourLine());
	child2.first = parent.first;
      }
    } else {
      if(newParent->id() < 0) { // a 3 bar state
	parent.second = child1.second;
	child2.second = child1.first;
      } else { // colour triplet
	parent.first = child1.first;
	child2.first = child1.second;
      }
    }
  } else if(child1.first) { // The child is a colour triplet
    if(newParent->hasColour() && newParent->hasAntiColour()) { // colour octet
      parent.first = child1.first;
      parent.second = child2.second = new_ptr(ColourLine());
    } else { // it must be a colour triplet, so child2 is a colour octet
      child2.second = child1.first;
      child2.first = parent.first = new_ptr(ColourLine());
    }
  } else if(child1.second) { // The child is a 3 bar state
    if(newParent->hasColour() && newParent->hasAntiColour()) { // colour octet
      parent.second = child1.second;
      parent.first = child2.first = new_ptr(ColourLine());
    } else { // it must be a colour triplet, so child2 is a colour octet
      child2.first = child1.second;
      child2.second = parent.second = new_ptr(ColourLine());
    }
  } else { // No colour info!
    cerr << "Shower/BackwardEvolver::setColour: "
	 << "No colour info for parton!\n";
    return;
  }
  if(parent.first) parent.first->addColoured(newParent);
  if(parent.second) parent.second->addAntiColoured(newParent);
  if(child2.first) child2.first->addColoured(otherChild);
  if(child2.second) child2.second->addAntiColoured(otherChild);
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  p  c = " << parent.first  << endl
		       << "  p  a = " << parent.second << endl
		       << "  c1 c = " << child1.first  << endl
		       << "  c1 a = " << child1.second << endl
		       << "  c2 c = " << child2.first  << endl
		       << "  c2 a = " << child2.second << endl;
  }
}

Ptr<SplittingGenerator>::pointer BackwardEvolver::splittingGenerator() {
  return _splittingGenerator;
}
