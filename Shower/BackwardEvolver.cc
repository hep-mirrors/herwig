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
				      tShowerParticlePtr particle, 
				      ShowerParticleVector &allShowerParticles)
  throw (Veto, Stop, Exception) {
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "BackwardEvolver::spaceLikeShower "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" 
		       << generator()->currentEventNumber() << endl;
  }
  cout << endl << "Process is " << *ch->currentEvent() << endl;

  bool hasEmitted = false;
  tShowerParticlePtr part = particle;
  tShowerParticleVector particlesYetToShower;   // only time-like particles

  do {

    cout << "-- BackwardEvolver, part = " << part << "." << endl; 
    Branching bb = _splittingGenerator->chooseBackwardBranching(ch, *part);
    if(bb.first == ShoKinPtr() || bb.second == tSudakovPtr()) {      
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	generator()->log() << "-- no further backward branching."
			   << endl;
      }
      hasEmitted = false;
      cout << "--- no further backward branching." << endl;
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
	if(id2 == ParticleID::g) id0 = part->id();
	else if(id0 == ParticleID::g) id2 = -part->id();
	
	// Now create the actual particles, make the otherChild a final state
	// particle, while the newParent is not	
	ShowerParticlePtr newParent = new_ptr(
	  ShowerParticle(getParticleData(id0)));
	ShowerParticlePtr otherChild = new_ptr(
          ShowerParticle(getParticleData(id2)));
	otherChild->setFinalState(true);
	newParent->setFinalState(false);
	//newParent->setFromHardSubprocess(true);

	cout << "Branching " 
	     << id0 << " -> "
	     << part->id() << ", "
	     << id2 << endl;

	ShowerIndex::InteractionType inter = splitF->interactionType();
	Energy scale = part->showerKinematics()->qtilde();

	// Set up the colour connections and the parent/child relationships
	createBranching(part,newParent,otherChild,scale,inter);
	part = newParent;
      } // if ( splitF )
    }   // if (shoKin && sudakov) {...} else {

    cout << "BackwardEvolver." << endl;
  } while(hasEmitted);

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
  long hadronId = part->parents()[0]->id();
  Energy oldQ;
  Energy minQ;
  long quarks[3];
  int maxIdx = 3;
  ShowerIndex::InteractionType inter = ShowerIndex::QCD;
  if(abs(hadronId) > 99) { // We have a hadron
    quarks[0] = hadronId % 10;
    quarks[1] = (hadronId/10)%10;
    quarks[2] = (hadronId/100)%10;
    if(quarks[2] == 0) maxIdx = 2; // we have a meson
    oldQ = part->evolutionScales()[ShowerIndex::QCD];

    // Look first at sea quarks, these must go to a gluon, we then handle
    // the gluons in the next step
    if(part->id() != quarks[0] || part->id() != quarks[1] || 
       part->id() != quarks[2] || part->id() != ParticleID::g) { 
      cout << "Particle is a sea quark\n"; 
      // determine some bounds for qtilde
      minQ = _splittingGenerator->showerVariables()->kinScale();
      cout << "CutOff = " << minQ << " and oldQ = " << oldQ << endl;
      // Create Shower Kinematics
      part->setSplittingFn(_splittingGenerator->
			   getSplittingFunction(part->id(),
						ParticleID::g));
      ShoKinPtr kinematics = forcedSplitting(*part,oldQ,minQ);
      part->setShowerKinematics(kinematics);

      // Create new particles, splitting is g->q qbar
      ShowerParticlePtr newParent = new_ptr(
		   ShowerParticle(getParticleData(ParticleID::g)));
      ShowerParticlePtr otherChild = new_ptr(
		   ShowerParticle(getParticleData(-part->id())));
      newParent->setFinalState(false);
      otherChild->setFinalState(true);
      //newParent->setFromHardSubprocess(true);

      createBranching(part,newParent,otherChild,kinematics->qtilde(),inter);

      // Store the old data so we can do the gluon splitting
      oldQ = kinematics->qtilde();
      part = newParent;

      // Put into list so it will be final showered
      //allShowerParticles.push_back(otherChild);
      //allShowerParticles.push_back(newParent);
      cout << "Created gluon splitting, gluon has scale " << oldQ << endl;
    }
    // We now handle the gluons, either it is where the shower terminated or
    // it has been created by splitting a sea quark
    int idx = 0;
    if(part->id() == ParticleID::g) { // gluon
      cout << "Particle is a gluon\n";
        // determine some bounds for qtilde
      minQ = _splittingGenerator->showerVariables()->kinScale();
      cout << "CutOff = " << minQ << " and oldQ = " << oldQ << endl;

      // Create new particles, splitting is q->g q
      // First choose which q
      idx = UseRandom::irnd(maxIdx);
      cout << "Chosen to split into " << quarks[idx] << endl;
      part->setSplittingFn(_splittingGenerator->
			   getSplittingFunction(ParticleID::g, quarks[idx]));

      // Create Shower Kinematics
      ShoKinPtr kinematics = forcedSplitting(*part,oldQ,minQ);
      part->setShowerKinematics(kinematics);

      ShowerParticlePtr newParent = new_ptr(
		   ShowerParticle(getParticleData(quarks[idx])));
      ShowerParticlePtr otherChild = new_ptr(
		   ShowerParticle(getParticleData(quarks[idx])));
      newParent->setFinalState(false);
      otherChild->setFinalState(true);
      //newParent->setFromHardSubprocess(true);

      // Set the colour and parent/child relationships
      createBranching(part,newParent,otherChild,kinematics->qtilde(),inter);

      // Add these so that they will be treated properly later
      //allShowerParticles.push_back(otherChild);
      allShowerParticles.push_back(newParent);
      //newParent->setFromHardSubprocess(true);
      part = newParent;
      cout << "Created final splitting " << kinematics->qtilde() 
	   << ", " << kinematics->z() << endl;    
    } else {
      // Otherwise figure out which particle we have ended on so we ignore it
      // in the remnant
      for(int i = 0; i<3; i++) if(part->id() == quarks[i]) idx = i;
      allShowerParticles.push_back(part);
      //part->setFromHardSubprocess(true);
    }
    // set remnant. 
    tPPtr hadron;
    if(part->parents().size() == 1) hadron = part->parents()[0];
    else cerr << "no remnant present!" << endl;

    // First decide what the remnant is
    long remId;
    int sign, spin;
    if(maxIdx == 2) { // Meson hadronic state
      remId = quarks[(idx+1)%2];
    } else { // Baryonic hadron
      // Get the other 2 elements of the array
      long id1 = quarks[(idx+1)%2];
      long id2 = quarks[(idx+2)%2];
      sign = (id1 < 0) ? -1 : 1; // Needed for the spin 0/1 part
      remId = id2*1000+id1*100;
      // Now decide if we have spin 0 diquark or spin 1 diquark
      if(id1 == id2 || UseRandom::rndbool()) spin = 3; // spin 1
      else spin = 1; // otherwise spin 0
      remId += sign*spin;

      // Create the remnant and set its momentum, also reset all of the decay 
      // products from the hadron
      PPtr newRemnant = new_ptr(Particle(getParticleData(remId)));
      newRemnant->setMomentum((1-part->x())*hadron->momentum());
      for(int i = hadron->children().size()-1; i!= -1; i--) {
	PPtr child = hadron->children()[i];
	hadron->removeChild(child);
	if(part != child) ch->currentStep()->removeParticle(child);
      }
      // Add the remnant to the step, this will be changed again if the
      // shower is vetoed. Set the colour connections as well
      hadron->addChild(newRemnant);
      //ch->currentStep()->addDecayProduct(hadron,newRemnant);
      hadron->addChild(part);
      if(part->id() > 0) part->antiColourNeighbour(newRemnant);
      else part->colourNeighbour(newRemnant);
    }
  }
  // Do we veto the whole shower after the final state showering or do we
  // seperately veto the initial state shower and final state shower?

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


ShoKinPtr BackwardEvolver::forcedSplitting(const ShowerParticle &particle,
					   Energy lastQ, Energy minQ) {
  // Now generate the new z and qtilde
  Energy newQ;
  double newZ,z0,z1;
  double randQ = UseRandom::rnd();
  double randZ = UseRandom::rnd();
  // Generate z with the same distributions as for regular splittings
  tSplittingFnPtr sf = particle.splitFun();
  // Bounds on z
  //if(particle.id() == ParticleID::g) z0 = (1.-sqrt(1.-4.*minQ/lastQ))/2.;
  //else z0 = minQ/lastQ/2.;
  z0 = minQ/lastQ/2.;
  z1 = 1-z0;
  if(!sf) {
    if(HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Shower) {
      generator()->log() << "The particle has no splitting function!"
			 << " Will use flat in z distribution." << endl;
    }
    newZ = (z1-z0)*randZ;
  } else {
    cout << "Trying to generate z we have a splitting fn of " << sf << endl;
    newZ = sf->invIntegOverP(sf->integOverP(z0) + randZ*(sf->integOverP(z1) - 
							 sf->integOverP(z0)));
  }
  // For the qtilde lets just start with a simple distribution weighted towards
  // the lower value: dP/dQ = 1/Q -> Q(R) = Q0^(1-R) Qmax^R
  newQ = pow(minQ,1-randQ)*pow(lastQ,randQ);

  Lorentz5Momentum p, n, ppartner, pcm;
  if(particle.isFromHardSubprocess()) {
    p = particle.momentum();
    ppartner = particle.partners()[ShowerIndex::QCD]->momentum();
    pcm = p; 
    pcm.boost((p + ppartner).findBoostToCM());	  
    n = Lorentz5Momentum( 0.0, -pcm.vect() ); 
    n.boost( -(p + ppartner).findBoostToCM() );
  } else {
    p = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
      ->showerKinematics()->getBasis()[0];
    n = dynamic_ptr_cast<ShowerParticlePtr>(particle.children()[0])
      ->showerKinematics()->getBasis()[1];
  } 
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  create ShowerKinematics with " 
		       << endl 
		       << "  p = " << p << endl 
		       << "  n = " << n << endl;
  }
  
  Ptr<IS_QtildaShowerKinematics1to2>::pointer showerKin = 
    new_ptr(IS_QtildaShowerKinematics1to2(p, n));

  // Phi is uniform
  ShowerVarsPtr vars = _splittingGenerator->showerVariables();;
  showerKin->qtilde(newQ);
  showerKin->setResScale(vars->cutoffQScale(ShowerIndex::QCD));
  showerKin->setKinScale(vars->kinScale()); 
  showerKin->z(newZ);
  showerKin->phi(2.*pi*UseRandom::rnd());

  return showerKin;
}

void BackwardEvolver::createBranching(ShowerParticlePtr part,
				      ShowerParticlePtr newParent,
				      ShowerParticlePtr otherChild,
				      Energy scale, 
				      ShowerIndex::InteractionType inter) {
  // no z for angular ordering in backward branchings
  newParent->setEvolutionScale(inter, scale);
  otherChild->setEvolutionScale(inter, scale);
  
  cout << "Creating branching " << newParent->id() << "-->" << part->id() 
       << " " << otherChild->id() << endl;
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // part -> (newParent, otherChild)
  ParticleVector theChildren; 
  theChildren.push_back(newParent); 
  theChildren.push_back(otherChild); 
  cout << "Updating children" << endl;
  part->showerKinematics()->updateChildren(part, theChildren); 
  
  cout << "Calling set Colour" << endl;
  // *** set proper colour connections
  setColour(newParent,part,otherChild);
  cout << "Colours set" << endl;

  // *** set proper parent/child relationships
  newParent->addChild(part);
  newParent->addChild(otherChild);
  newParent->x(part->x()/part->showerKinematics()->z());
  
  // Print some messages
  cout << "  new x = " 
       << newParent->x() << endl;
  cout << "  | new parent   = " << newParent << endl;
  cout << "  |-- child 0    = " << newParent->children()[0] << endl;
  cout << "  |-- child 1    = " << newParent->children()[1] << endl; 

  // Now fix the hadrons connections
  tPPtr hadron;
  if(part->parents().size() == 2) hadron = part->parents()[0];
  else cerr << "not one parent!" << endl; 
  hadron->addChild(newParent);
  hadron->removeChild(part);
  part->removeParent(hadron);
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
    if(newParent->id() < 0) { // a 3 bar state
      parent.second = child1.second;
      child2.second = child1.first;
    } else { // colour triplet
      parent.first = child1.first;
      child2.first = child1.second;
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
    cout << "No colour info for parton!\n";
    return;
  }
  if(parent.first) parent.first->addColoured(newParent);
  if(parent.second) parent.second->addAntiColoured(newParent);
  if(child2.first) child2.first->addColoured(otherChild);
  if(child2.second) child2.second->addAntiColoured(otherChild);
}
