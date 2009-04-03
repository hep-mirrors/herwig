// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QEDEvolver class.
//

#include "QEDEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Shower/ShowerHandler.h"

using namespace Herwig;

QEDEvolver::QEDEvolver() : QCDFirst_(true) {
  interactions_.push_back(ShowerInteraction::QCD);
  interactions_.push_back(ShowerInteraction::QED);
}

void QEDEvolver::doinit() {
  PowhegEvolver::doinit();
  if(!QCDFirst_) swap(interactions_[0],interactions_[1]);
}
IBPtr QEDEvolver::clone() const {
  return new_ptr(*this);
}

IBPtr QEDEvolver::fullclone() const {
  return new_ptr(*this);
}

void QEDEvolver::persistentOutput(PersistentOStream & os) const {
  os << QCDFirst_ << interactions_.size();
  for(unsigned int ix=0;ix<interactions_.size();++ix) os << oenum(interactions_[ix]); 
}

void QEDEvolver::persistentInput(PersistentIStream & is, int) {
  unsigned int isize;
  is >> QCDFirst_ >> isize;
  interactions_.resize(isize);
  for(unsigned int ix=0;ix<interactions_.size();++ix) is >> ienum(interactions_[ix]);
}

ClassDescription<QEDEvolver> QEDEvolver::initQEDEvolver;
// Definition of the static class description member.

void QEDEvolver::Init() {

  static ClassDocumentation<QEDEvolver> documentation
    ("The QEDEvolver class implements the shower including QED radiation.");

  static Switch<QEDEvolver,bool> interfaceFirstInteraction
    ("FirstInteraction",
     "The interaction to use first in the shower ",
     &QEDEvolver::QCDFirst_, true, false, false);
  static SwitchOption interfaceFirstInteractionQCD
    (interfaceFirstInteraction,
     "QCD",
     "Do QCD first",
     true);
  static SwitchOption interfaceFirstInteractionQED
    (interfaceFirstInteraction,
     "QED",
     "Do QED first",
     false);

}

void QEDEvolver::showerHardProcess(ShowerTreePtr hard) {
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
  for(unsigned int inter=0;inter<interactions_.size();++inter) {
    // zero pt so only added first time round
    if(inter!=0) intrinsicpT().clear();
    // set up for second pass if required
    if(inter!=0) constructHardTree(particlesToShower,interactions_[inter]);
    // main shower loop
    unsigned int ntry(0);
    do {
      // clear results of last attempt if needed
      if(ntry!=0) {
	currentTree()->clear();
	setEvolutionPartners(true,interactions_[inter]);
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
	  progenitor()->hasEmitted(startSpaceLikeShower(beamparticle,interactions_[inter]));
	}
	// final-state
	else {
	  if(!isFSRadiationON()) continue;
	  // perform shower
	  progenitor()->hasEmitted(startTimeLikeShower(interactions_[inter]));
	}
      }
    }
    while(!showerModel()->kinematicsReconstructor()->
	  reconstructHardJets(hard,intrinsicpT())&&
	  maximumTries()>++ntry);
    if(maximumTries()==ntry) throw ShowerHandler::ShowerTriesVeto(ntry);
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
  // loop over possible interactions
  for(unsigned int inter=0;inter<interactions_.size();++inter) {
    unsigned int ntry(0);
    // set up for second pass if required
    if(inter!=0) constructDecayTree(particlesToShower,interactions_[inter]);
    // main showering loop
    do {
      // clear results of last attempt
      if(ntry!=0) {
	currentTree()->clear();
	setEvolutionPartners(false,interactions_[inter]);
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
	  progenitor()->hasEmitted(startTimeLikeShower(interactions_[inter]));
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
	  progenitor()->hasEmitted(startSpaceLikeDecayShower(maxscale,minmass,
							     interactions_[inter])); 
	}
      }
    }
    while(!showerModel()->kinematicsReconstructor()->reconstructDecayJets(decay)&&
	  maximumTries()>++ntry);
    if(maximumTries()==ntry) 
      throw Exception() << "Failed to generate the shower after "
			<< ntry << " attempts in Evolver::showerDecay()"
			<< Exception::eventerror;
  }
  // tree has now showered
  currentTree()->hasShowered(true);
  checkShowerMomentum( particlesToShower );
}

void QEDEvolver::constructTimeLikeLine(tHardBranchingPtr branch,
				       tShowerParticlePtr particle) {
  for(unsigned int ix=0;ix<particle->children().size();++ix) {
//     HardBranching::Status status = 
//       particle->children()[ix]->id()==particle->id() ?
//       branch->status() : HardBranching::Outgoing;
    HardBranching::Status status = branch->status();
    tShowerParticlePtr child = 
      dynamic_ptr_cast<ShowerParticlePtr>(particle->children()[ix]);
    if(child->children().empty()) {
      HardBranchingPtr newBranch = 
	new_ptr(HardBranching(child,SudakovPtr(),branch,status));
      branch->addChild(newBranch);
    }
    else {
      HardBranchingPtr newBranch = 
	new_ptr(HardBranching(child,child->showerKinematics()->SudakovFormFactor(),
			      branch,status));
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
    new_ptr(HardBranching(particle,sud,last,HardBranching::Incoming));
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
			    last,HardBranching::Outgoing));
    constructTimeLikeLine(timeBranch,timeChild);
  }
  else {
    timeBranch = 
      new_ptr(HardBranching(timeChild,SudakovPtr(),last,HardBranching::Outgoing));
  }
  last->addChild(timeBranch);
  last=newBranch;
}

void QEDEvolver::constructHardTree(vector<ShowerProgenitorPtr> & particlesToShower,
				   ShowerInteraction::Type inter) {
  bool noEmission = true;
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if(particlesToShower[ix]->progenitor()->isFinalState()) {
      HardBranchingPtr newBranch;
      if(particlesToShower[ix]->hasEmitted()) {
	noEmission = false;
	newBranch = 
	  new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				particlesToShower[ix]->progenitor()->
				showerKinematics()->SudakovFormFactor(),
				HardBranchingPtr(),HardBranching::Outgoing));
	constructTimeLikeLine(newBranch,particlesToShower[ix]->progenitor());
      }
      else {
	newBranch = 
	  new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				SudakovPtr(),HardBranchingPtr(),
				HardBranching::Outgoing));
      }
      allBranchings.push_back(newBranch);
    }
    else {
      HardBranchingPtr first,last;
      if(!particlesToShower[ix]->progenitor()->parents().empty()) {
	noEmission = false;
	constructSpaceLikeLine(particlesToShower[ix]->progenitor(),
			       first,last,SudakovPtr(),
			       particlesToShower[ix]->original()->parents()[0]);
      }
      else {
	first = new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				      SudakovPtr(),HardBranchingPtr(),
				      HardBranching::Incoming));
	last = first;
      }
      spaceBranchings.push_back(first);
      allBranchings.push_back(last);
    }
  }
  if(!noEmission) {
    HardTreePtr QCDTree = new_ptr(HardTree(allBranchings,spaceBranchings));
    // do the inverse recon
    if(!showerModel()->kinematicsReconstructor()->
       deconstructHardJets(QCDTree,this,inter))
      throw Exception() << "Can't to shower deconstruction for QED shower in"
			<< "QEDEvolver::showerHard" << Exception::eventerror;
    // set the hard tree
    hardTree(QCDTree);
  }
  // clear the old shower
  currentTree()->clear();
  // set the charge partners
  setEvolutionPartners(true,inter);
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

void QEDEvolver::constructDecayTree(vector<ShowerProgenitorPtr> & particlesToShower,
				    ShowerInteraction::Type inter) {
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if(particlesToShower[ix]->progenitor()->isFinalState()) {
      HardBranchingPtr newBranch;
      if(particlesToShower[ix]->hasEmitted()) {
	newBranch = 
	  new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				particlesToShower[ix]->progenitor()->
				showerKinematics()->SudakovFormFactor(),
				HardBranchingPtr(),HardBranching::Outgoing));
	constructTimeLikeLine(newBranch,particlesToShower[ix]->progenitor());
      }
      else {
	newBranch = 
	  new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				SudakovPtr(),HardBranchingPtr(),
				HardBranching::Outgoing));
      }
      allBranchings.push_back(newBranch);
    }
    else {
      HardBranchingPtr newBranch;
      if(particlesToShower[ix]->hasEmitted()) {
	newBranch = 
	  new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				particlesToShower[ix]->progenitor()->
				showerKinematics()->SudakovFormFactor(),
				HardBranchingPtr(),HardBranching::Decay));
	constructTimeLikeLine(newBranch,particlesToShower[ix]->progenitor());
	HardBranchingPtr last=newBranch;
	do {
	  for(unsigned int ix=0;ix<last->children().size();++ix) {
	    if(last->children()[ix]->branchingParticle()->id()==
	       particlesToShower[ix]->id()) {
	      last = last->children()[ix];
	      continue;
	    }
	  }
	}
	while(!last->children().empty());
	last->status(HardBranching::Incoming);
	spaceBranchings.push_back(newBranch);
	allBranchings  .push_back(last);
      }
      else {
	newBranch = 
	  new_ptr(HardBranching(particlesToShower[ix]->progenitor(),
				SudakovPtr(),HardBranchingPtr(),
				HardBranching::Incoming));
	spaceBranchings.push_back(newBranch);
	allBranchings  .push_back(newBranch);
      }
    }
  }
  HardTreePtr QCDTree = new_ptr(HardTree(allBranchings,spaceBranchings));
  // set the charge partners
  ShowerParticleVector particles;
  particles.push_back(spaceBranchings.back()->branchingParticle());
  for(set<HardBranchingPtr>::iterator cit=QCDTree->branchings().begin();
      cit!=QCDTree->branchings().end();++cit) {
    if((*cit)->status()==HardBranching::Outgoing)
      particles.push_back((*cit)->branchingParticle());
  }
//   showerModel()->partnerFinder()->setInitialEvolutionScales(particles,true,inter,true);
  // test fix to get same partners again
  showerModel()->partnerFinder()->setInitialEvolutionScales(particles,true,inter,true);
  // do the inverse recon
  if(!showerModel()->kinematicsReconstructor()->
     deconstructDecayJets(QCDTree,this,inter))
    throw Exception() << "Can't to shower deconstruction for QED shower in"
		      << "QEDEvolver::showerDecay" << Exception::eventerror;
  // clear the old shower
  currentTree()->clear();
  // set the hard tree
  hardTree(QCDTree);
  // set the charge partners
  setEvolutionPartners(false,inter);
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
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit=hardTree()->particles().end(),
      mit = hardTree()->particles().find(particlesToShower[ix]->progenitor());
    if( mit != eit) {
      if(mit->second->status()==HardBranching::Outgoing)
	particlesToShower[ix]->progenitor()->set5Momentum(mit->second->pVector());
    }
  }
}
