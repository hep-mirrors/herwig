// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QEDEvolver class.
//

#include "QEDEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/ShowerHandler.h"

using namespace Herwig;

QEDEvolver::QEDEvolver() {}

IBPtr QEDEvolver::clone() const {
  return new_ptr(*this);
}

IBPtr QEDEvolver::fullclone() const {
  return new_ptr(*this);
}

void QEDEvolver::persistentOutput(PersistentOStream & os) const {
}

void QEDEvolver::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<QEDEvolver> QEDEvolver::initQEDEvolver;
// Definition of the static class description member.

void QEDEvolver::Init() {

  static ClassDocumentation<QEDEvolver> documentation
    ("There is no documentation for the QEDEvolver class");

}

void QEDEvolver::showerHardProcess(ShowerTreePtr hard) {
//   PowhegEvolver::showerHardProcess(hard);
  // set the current tree
  currentTree(hard);
  // extract particles to shower
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(true);
  // setup the maximum scales for the shower, given by the hard process
  if (hardVetoOn()) 
    setupMaximumScales(currentTree(), particlesToShower);
  // generate the intrinsic p_T once and for all
  generateIntrinsicpT(particlesToShower);
  // main showering loop
  for(ShowerInteraction::Type inter=ShowerInteraction::QCD;
      inter<=ShowerInteraction::QED;inter=ShowerInteraction::Type(int(inter)+1)) {
    unsigned int ntry(0);
    do {
      // clear results of last attempt
      if(ntry!=0) {
	currentTree()->clear();
	setEvolutionPartners(true,ShowerInteraction::QCD);
      }
      // generate the shower
      // pick random starting point 
      unsigned int istart=UseRandom::irnd(particlesToShower.size());
      unsigned int istop = particlesToShower.size();
      // loop over particles with random starting point
      for(unsigned int ix=istart;ix<=istop;++ix) {
	if(ix==particlesToShower.size()) {
	  if(istart!=0) {
	    istop = istart-1;
	    ix=0;
	  }
	  else break;
	}
	// set the progenitor
	progenitor(particlesToShower[ix]);
	// initial-state
	if(!progenitor()->progenitor()->isFinalState()) {
	  if(!isISRadiationON()) continue;
	  // get the PDF
	  setBeamParticle(progenitor()->beam());
	  assert(beamParticle());
	  // perform the shower
	  // set the beam particle
	  tPPtr beamparticle=progenitor()->original();
	  if(!beamparticle->parents().empty()) 
	    beamparticle=beamparticle->parents()[0];
	  // generate the shower
	  progenitor()->hasEmitted(startSpaceLikeShower(beamparticle,inter));
	}
	// final-state
	else {
	  if(!isFSRadiationON()) continue;
	  // perform shower
	  progenitor()->hasEmitted(startTimeLikeShower(inter));
	}
      }
    }
    while(!showerModel()->kinematicsReconstructor()->
	  reconstructHardJets(hard,intrinsicpT())&&
	  maximumTries()>++ntry);
    if(maximumTries()==ntry) throw ShowerHandler::ShowerTriesVeto(ntry);
    if(inter==ShowerInteraction::QCD) {
      vector<HardBranchingPtr> spaceBranchings,allBranchings;
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
 	if(particlesToShower[ix]->progenitor()->isFinalState()) {
	  HardBranchingPtr newBranch;
	  if(particlesToShower[ix]->hasEmitted()) {
	    newBranch = 
	      new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				    particlesToShower[ix]->progenitor()->
				    showerKinematics()->SudakovFormFactor(),
				    HardBranchingPtr(),false));
	    constructTimeLikeLine(newBranch,particlesToShower[ix]->progenitor());
	  }
	  else {
	    newBranch = 
	      new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				    SudakovPtr(),HardBranchingPtr(),false));
	  }
	  allBranchings.push_back(newBranch);
	}
	else {
	  HardBranchingPtr first,last;
	  constructSpaceLikeLine(particlesToShower[ix]->progenitor(),first,last,
				 SudakovPtr(),
				 particlesToShower[ix]->original()->parents()[0]);
	  spaceBranchings.push_back(first);
	  allBranchings.push_back(last);
 	}
      }
      HardTreePtr QCDTree = new_ptr(HardTree(allBranchings,spaceBranchings));
      // do the inverse recon
      if(!showerModel()->kinematicsReconstructor()->
	 deconstructHardJets(QCDTree,this,ShowerInteraction::QED))
	throw Exception() << "Can't to shower deconstruction for QED shower in"
			  << "QEDEvolver::showerHard" << Exception::eventerror;
      // clear the old shower
      currentTree()->clear();
      // set the hard tree
      hardTree(QCDTree);
      // set the charge partners
      setEvolutionPartners(false,ShowerInteraction::QED);
      // get the particles to be showered
      map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
      map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
      particlesToShower.clear();
      // incoming particles
      for(cit=currentTree()->incomingLines().begin();
	  cit!=currentTree()->incomingLines().end();++cit)
	particlesToShower.push_back(((*cit).first));
      assert(particlesToShower.size()==2);
      // outgoing particles
      for(cjt=currentTree()->outgoingLines().begin();
	  cjt!=currentTree()->outgoingLines().end();++cjt)
	particlesToShower.push_back(((*cjt).first));
    }
  }
  // tree has now showered
  currentTree()->hasShowered(true);
}

void QEDEvolver::showerDecay(ShowerTreePtr decay) {
  // set the ShowerTree to be showered
  currentTree(decay);
  // extract particles to be shower, set scales and perform hard matrix element 
  // correction
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(false);
  setupMaximumScales(currentTree(), particlesToShower);
  // compute the minimum mass of the final-state
  Energy minmass(ZERO);
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if(particlesToShower[ix]->progenitor()->isFinalState())
      minmass+=particlesToShower[ix]->progenitor()->mass();
  }
  // main showering loop
  for(ShowerInteraction::Type inter=ShowerInteraction::QCD;
      inter<=ShowerInteraction::QED;inter=ShowerInteraction::Type(int(inter)+1)) {
    unsigned int ntry(0);
    do {
      // clear results of last attempt
      if(ntry!=0) {
	currentTree()->clear();
	setEvolutionPartners(false,inter);
      }
      unsigned int istart=UseRandom::irnd(particlesToShower.size());
      unsigned int istop = particlesToShower.size();
      // loop over particles with random starting point
      for(unsigned int ix=istart;ix<=istop;++ix) {
	if(ix==particlesToShower.size()) {
	  if(istart!=0) {
	    istop = istart-1;
	    ix=0;
	  }
	  else break;
	}
	// extract the progenitor
	progenitor(particlesToShower[ix]);
	// final-state radiation
	if(progenitor()->progenitor()->isFinalState()) {
	  if(!isFSRadiationON()) continue;
	  // perform shower
	  progenitor()->hasEmitted(startTimeLikeShower(inter));
	}
	// initial-state radiation
	else {
	  if(!isISRadiationON()) continue;
	  // perform shower
	  // set the scales correctly. The current scale is the maximum scale for
	  // emission not the starting scale
	  Energy maxscale=progenitor()->progenitor()->evolutionScale();
	  Energy startScale=progenitor()->progenitor()->mass();
	  progenitor()->progenitor()->setEvolutionScale(startScale);
	  // perform the shower
	  progenitor()->hasEmitted(spaceLikeDecayShower(progenitor()->progenitor(),
							maxscale,minmass,inter)); 
	}
      }
    }
    while(!showerModel()->kinematicsReconstructor()->reconstructDecayJets(decay)&&
	  maximumTries()>++ntry);
    if(maximumTries()==ntry) 
      throw Exception() << "Failed to generate the shower after "
			<< ntry << " attempts in Evolver::showerDecay()"
			<< Exception::eventerror;
    if(inter==ShowerInteraction::QCD) {
      vector<HardBranchingPtr> spaceBranchings,allBranchings;
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	if(particlesToShower[ix]->progenitor()->isFinalState()) {
	  HardBranchingPtr newBranch;
	  if(particlesToShower[ix]->hasEmitted()) {
	    newBranch = 
	      new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				    particlesToShower[ix]->progenitor()->
				    showerKinematics()->SudakovFormFactor(),
				    HardBranchingPtr(),false));
	    constructTimeLikeLine(newBranch,particlesToShower[ix]->progenitor());
	  }
	  else {
	    newBranch = 
	      new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				    SudakovPtr(),HardBranchingPtr(),false));
	  }
	  allBranchings.push_back(newBranch);
	}
	else {
	  spaceBranchings.
	    push_back(new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
					    SudakovPtr(),HardBranchingPtr(),true)));
	  allBranchings.push_back(spaceBranchings.back());
	}
      }
      HardTreePtr QCDTree = new_ptr(HardTree(allBranchings,spaceBranchings));
      // do the inverse recon
      if(!showerModel()->kinematicsReconstructor()->
	 deconstructDecayJets(QCDTree,this,ShowerInteraction::QED))
	throw Exception() << "Can't to shower deconstruction for QED shower in"
			  << "QEDEvolver::showerDecay" << Exception::eventerror;
      // clear the old shower
      currentTree()->clear();
      // set the hard tree
      hardTree(QCDTree);
      // set the charge partners
      setEvolutionPartners(false,ShowerInteraction::QED);
      // get the particles to be showered
      map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
      map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
      particlesToShower.clear();
      // incoming particles
      for(cit=currentTree()->incomingLines().begin();
	  cit!=currentTree()->incomingLines().end();++cit)
	particlesToShower.push_back(((*cit).first));
      assert(particlesToShower.size()==1);
      // outgoing particles
      for(cjt=currentTree()->outgoingLines().begin();
	  cjt!=currentTree()->outgoingLines().end();++cjt)
	particlesToShower.push_back(((*cjt).first));
    }
  }
  // tree has now showered
  currentTree()->hasShowered(true);
  checkShowerMomentum( particlesToShower );
}

void QEDEvolver::constructTimeLikeLine(tHardBranchingPtr branch,
				       tShowerParticlePtr particle) {
  for(unsigned int ix=0;ix<particle->children().size();++ix) {
    tShowerParticlePtr child = dynamic_ptr_cast<ShowerParticlePtr>(particle->children()[ix]);
    if(child->children().empty()) {
      HardBranchingPtr newBranch = 
	new_ptr(HardBranching(child,SudakovPtr(),branch,false));
      branch->addChild(newBranch);
    }
    else {
      HardBranchingPtr newBranch = 
	new_ptr(HardBranching(child,child->showerKinematics()->SudakovFormFactor(),
			      branch,false));
      constructTimeLikeLine(newBranch,child);
      branch->addChild(newBranch);
    }
  }
}

void QEDEvolver::constructSpaceLikeLine(tShowerParticlePtr particle,
					HardBranchingPtr & first,
					HardBranchingPtr & last,
					SudakovPtr sud,PPtr beam) {
  if(!particle) return;
  if(!particle->parents().empty()) {
    tShowerParticlePtr parent = 
      dynamic_ptr_cast<ShowerParticlePtr>(particle->parents()[0]);
    SudakovPtr newSud=particle->showerKinematics()->SudakovFormFactor();
    constructSpaceLikeLine(parent,first,last,newSud,beam);
  }
  HardBranchingPtr newBranch = 
    new_ptr(HardBranching(particle,sud,last,true));
  newBranch->beam(beam);
  if(!first) {
    first=newBranch;
    last =newBranch;
    return;
  }
  last->addChild(newBranch);
  tShowerParticlePtr timeChild = 
    dynamic_ptr_cast<ShowerParticlePtr>(particle->parents()[0]->children()[1]);
  HardBranchingPtr timeBranch;
  if(!timeChild->children().empty()) {
    timeBranch = 
      new_ptr(HardBranching(timeChild,
			    timeChild->showerKinematics()->SudakovFormFactor(),
			    last,false));
    constructTimeLikeLine(timeBranch,timeChild);
  }
  else {
    timeBranch = 
      new_ptr(HardBranching(timeChild,SudakovPtr(),last,false));
  }
  last->addChild(timeBranch);
  last=newBranch;
}
