// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NasonEvolver class.
//

#include "NasonEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/ShowerHandler.h"
#include "DefaultEmissionGenerator.h"
#include "Herwig++/Shower/Default/FS_QtildaShowerKinematics1to2.h"
#include "NasonTree.h"

using namespace Herwig;

void NasonEvolver::persistentOutput(PersistentOStream & os) const {
  os << _hardgenerator;
}

void NasonEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _hardgenerator;
}

ClassDescription<NasonEvolver> NasonEvolver::initNasonEvolver;
// Definition of the static class description member.

void NasonEvolver::Init() {

  static ClassDocumentation<NasonEvolver> documentation
    ("The NasonEvolver implements the Nason approach to MC@NLO");

  static RefVector<NasonEvolver,HardestEmissionGenerator> interfaceHardGenerator
    ("HardGenerator",
     "The objects responsible for generating the hardestr emission",
     &NasonEvolver::_hardgenerator, -1, false, false, true, false, false);

}

void NasonEvolver::showerDecay(ShowerTreePtr tree) {
  // set the tree
  currentTree(tree);
  // set up the shower
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(false);
  // main showering loop
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      currentTree()->clear();
      setColourPartners(false);
    }
    // initial-state radiation
//     if(_splittingGenerator->isISRadiationON()) {
//       // compute the minimum mass of the final-state
//       Energy minmass(0.);
//       for(unsigned int ix=0;ix<particlesToShower.size();++ix)
// 	{if(particlesToShower[ix]->progenitor()->isFinalState())
// 	    minmass+=particlesToShower[ix]->progenitor()->mass();}
//       for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
// 	// only consider initial-state particles
// 	if(particlesToShower[ix]->progenitor()->isFinalState()) continue;
// 	// perform shower
// 	_progenitor=particlesToShower[ix];
// 	// set the scales correctly. The current scale is the maximum scale for
// 	// emission not the starting scale
// 	vector<Energy> maxscale=progenitor()->progenitor()->evolutionScales();
// 	Energy startScale=progenitor()->progenitor()->mass();
// 	progenitor()->progenitor()->setEvolutionScale(ShowerIndex::QCD,startScale);
// 	progenitor()->progenitor()->setEvolutionScale(ShowerIndex::QED,startScale);
// 	progenitor()->progenitor()->setEvolutionScale(ShowerIndex::EWK,startScale);
// 	// perform the shower
// 	progenitor()->hasEmitted(spaceLikeDecayShower(progenitor()->progenitor(),
// 						     maxscale,minmass)); 
//       }
//     }
    // final-state radiation
    map<ShowerParticlePtr,tNasonBranchingPtr>::const_iterator mit,eit;
    if(_nasontree) eit=_nasontree->particles().end();
    if(isFSRadiationON()) {
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
 	// only consider final-state particles
 	if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
 	// perform shower
 	progenitor(particlesToShower[ix]);
	if(_nasontree) mit=_nasontree->particles().find(progenitor()->progenitor());
	if(_nasontree&&mit!=eit&&!mit->second->_children.empty()) {
	  progenitor()->
	    hasEmitted(truncatedTimeLikeShower(particlesToShower[ix]->progenitor(),
					       mit->second));
	}
	else {
	  progenitor()->hasEmitted(timeLikeShower(particlesToShower[ix]->progenitor()));
	}
      }
    }
  }
  while(!showerModel()->kinematicsReconstructor()->reconstructDecayJets(tree)&&
	maximumTries()>++ntry);
  if(maximumTries()==ntry) 
    throw Exception() << "Failed to generate the shower after "
		      << ntry << " attempts in NasonEvolver::showerDecay()"
		      << Exception::eventerror;
  currentTree()->hasShowered(true);
}

void NasonEvolver::showerHardProcess(ShowerTreePtr tree) {
  // set the tree
  currentTree(tree);
  // set up the shower
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(true);
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      currentTree()->clear();
      setColourPartners(true);
    }
    // initial-state radiation
    map<ShowerParticlePtr,tNasonBranchingPtr>::const_iterator mit,eit;
    if(_nasontree) eit=_nasontree->particles().end();
    if(isISRadiationON()) {
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	// only consider initial-state particles
	if(particlesToShower[ix]->progenitor()->isFinalState()) continue;
	// get the PDF
	setBeamParticle(particlesToShower[ix]->beam());
	if(!beamParticle()) throw Exception() << "The Beam particle does not have"
					      << " BeamParticleData in Evolver::" 
					      << "showerhardProcess()" 
					      << Exception::runerror;
	// perform the shower
	progenitor(particlesToShower[ix]);
	if(_nasontree) mit=_nasontree->particles().find(progenitor()->progenitor());
	if(_nasontree&&mit!=eit&&mit->second->_parent) {
	  progenitor()->
	    hasEmitted(truncatedSpaceLikeShower(particlesToShower[ix]->progenitor(),
						particlesToShower[ix]->
						original()->parents()[0],
						mit->second->_parent));
	}
	else {
	  progenitor()->
	    hasEmitted(spaceLikeShower(particlesToShower[ix]->progenitor(),
				       particlesToShower[ix]->
				       original()->parents()[0]));
	}
      }
    }
    // final-state radiation
    if(isFSRadiationON()) {
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	// only consider final-state particles
	if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
	// perform shower
	progenitor(particlesToShower[ix]);
	progenitor()->hasEmitted(timeLikeShower(particlesToShower[ix]->progenitor()));
      }
    }
  }
  while(!showerModel()->kinematicsReconstructor()->reconstructHardJets(tree)&&
 	maximumTries()>++ntry);
  if(maximumTries()==ntry) throw Exception() << "Failed to generate the shower after "
 				      << ntry 
 				      << " attempts in NasonEvolver::showerHardProcess()"
 				      << Exception::eventerror;
  currentTree()->hasShowered(true);
}

vector<ShowerProgenitorPtr> NasonEvolver::setupShower(bool hard) {
  // set the colour partners
  setColourPartners(hard);
  // generate the hardest emission
  hardestEmission();
  // get the particles to be showered
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=currentTree()->incomingLines().begin();
      cit!=currentTree()->incomingLines().end();++cit)
    particlesToShower.push_back(((*cit).first));
  assert((particlesToShower.size()==1&&!hard)||(particlesToShower.size()==2&&hard));
  // outgoing particles
  for(cjt=currentTree()->outgoingLines().begin();
      cjt!=currentTree()->outgoingLines().end();++cjt)
    particlesToShower.push_back(((*cjt).first));
  // remake the colour partners if needed
  if(currentTree()->hardMatrixElementCorrection()) {
    setColourPartners(hard);
    currentTree()->resetShowerProducts();
  }
  return particlesToShower;
}

void NasonEvolver::hardestEmission() {
  // see if there is an appropriate hard emission generator
  HardestEmissionGeneratorPtr currenthard=HardestEmissionGeneratorPtr();
  for(unsigned int ix=0;ix<_hardgenerator.size();++ix) {
    if(!_hardgenerator[ix]->canHandle(currentTree())) continue;
    if(currenthard) {
      if(dynamic_ptr_cast<DefaultEmissionGeneratorPtr>(currenthard)) {
	currenthard=_hardgenerator[ix];
      }
      else if(!dynamic_ptr_cast<DefaultEmissionGeneratorPtr>(_hardgenerator[ix])) {
	ostringstream output;
	output << "There is more than one possible hard emission generator"
	       << " which could be used for ";
	map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
	map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
	for(cit=currentTree()->incomingLines().begin();
	    cit!=currentTree()->incomingLines().end();++cit)
	  {output << cit->first->progenitor()->PDGName() << " ";}
	output << " -> ";
	for(cjt=currentTree()->outgoingLines().begin();
	    cjt!=currentTree()->outgoingLines().end();++cjt)
	  {output << cjt->first->progenitor()->PDGName() << " ";}
	output << "in NasonEvolver::hardestEmission()\n";
	throw Exception() << output << Exception::runerror;
      }
    }
    else {
      currenthard=_hardgenerator[ix];
    }
  }
  // if no suitable generator return
  _nasontree=NasonTreePtr();
  if(!currenthard) return;
  // generate the hardest emission
  _nasontree = currenthard->generateHardest(currentTree());
}

bool NasonEvolver::truncatedTimeLikeShower(tShowerParticlePtr particle,
					   NasonBranchingPtr branch) {
  bool vetoed = true;
  Branching fb;
  unsigned int iout;
  tcPDPtr pdata[2];
  while (vetoed) {
    iout=0;
    vetoed = false;
    fb=splittingGenerator()->chooseForwardBranching(*particle,1.);
    // check haven't evolved too far
    if(!fb.kinematics||fb.kinematics->scale()<branch->_scale) {
      fb=Branching();
      break;
    }
    // get the particle data objects
    for(unsigned int ix=0;ix<2;++ix) pdata[ix]=getParticleData(fb.ids[ix+1]);
    if(particle->id()!=fb.ids[0]) {
      for(unsigned int ix=0;ix<2;++ix) {
	tPDPtr cc(pdata[ix]->CC());
	if(cc) pdata[ix]=cc;
      }
    }
    // find the truncated line
    if(pdata[0]->id()!=pdata[1]->id()) {
      if(pdata[0]->id()==particle->id())       iout=1;
      else if (pdata[1]->id()==particle->id()) iout=2;
    }
    else if(pdata[0]->id()==particle->id()) {
      if(fb.kinematics->z()>0.5) iout=1;
      else                       iout=2;
    }
    // apply the vetos for the truncated shower
    if(iout==0) vetoed=true;
    else {
      if(iout==1&&fb.kinematics->z()<0.5)      vetoed=true;
      else if(iout==2&&fb.kinematics->z()>0.5) vetoed=true;
    }
    // pt veto
    if(fb.kinematics->pT()>progenitor()->maximumpT()) vetoed = true;
    if(vetoed) particle->setEvolutionScale(ShowerIndex::QCD,fb.kinematics->scale());
  }
  // if no branching set decay matrix and return
  if(!fb.kinematics) {
    // construct the kinematics for the hard emission
    ShoKinPtr showerKin=
      branch->_sudakov->createFinalStateBranching(branch->_scale,
						  branch->_children[0]->_z,
						  branch->_phi,
						  branch->_children[0]->_pt);
    particle->setEvolutionScale(ShowerIndex::QCD,branch->_scale);
    showerKin->initialize(*particle,PPtr());
    IdList idlist(3);
    idlist[0]=particle->id();
    idlist[1]=branch->_children[0]->_particle->id();
    idlist[2]=branch->_children[1]->_particle->id();
    fb=Branching(showerKin,idlist,branch->_sudakov);
    // Assign the shower kinematics to the emitting particle.
    particle->setShowerKinematics(fb.kinematics);
    // Assign the splitting function to the emitting particle. 
    // For the time being we are considering only 1->2 branching
    // Create the ShowerParticle objects for the two children of
    // the emitting particle; set the parent/child relationship
    // if same as definition create particles, otherwise create cc
    ShowerParticleVector theChildren;
    theChildren.push_back(new_ptr(ShowerParticle(branch->_children[0]->
						 _particle->dataPtr(),true)));
    theChildren.push_back(new_ptr(ShowerParticle(branch->_children[1]->
						 _particle->dataPtr(),true)));
    particle->showerKinematics()->updateChildren(particle, theChildren);
    // update the history if needed
    if(particle==currentTree()->getFinalStateShowerProduct(progenitor()))
      currentTree()->updateFinalStateShowerProduct(progenitor(),
						   particle,theChildren);
    currentTree()->addFinalStateBranching(particle,theChildren);
    // shower the first  particle
    //
    //  need to set rho here
    //
    if(branch->_children[0]->_children.empty()) {
      timeLikeShower(theChildren[0]);
    }
    else {
      truncatedTimeLikeShower(theChildren[0],branch->_children[0]);
    } 
    // shower the second particle
    //
    //   need to set rho here
    //
    if(branch->_children[1]->_children.empty()) {
      timeLikeShower(theChildren[1]);
    }
    else {
      truncatedTimeLikeShower(theChildren[1],branch->_children[1]);
    }
    return true;
  }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->setShowerKinematics(fb.kinematics);
  // Assign the splitting function to the emitting particle. 
  // For the time being we are considering only 1->2 branching
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  // if same as definition create particles, otherwise create cc
  ShowerParticleVector theChildren; 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[0],true))); 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[1],true)));
  particle->showerKinematics()->updateChildren(particle, theChildren);
  // update the history if needed
  if(particle==currentTree()->getFinalStateShowerProduct(progenitor()))
    currentTree()->updateFinalStateShowerProduct(progenitor(),
 						particle,theChildren);
  currentTree()->addFinalStateBranching(particle,theChildren);
  // shower the first  particle
  //
  //  need to set rho here
  //
  if(iout==1) truncatedTimeLikeShower(theChildren[0],branch);
  else        timeLikeShower(theChildren[0]);
  // shower the second particle
  //
  //   need to set rho here
  //
  if(iout==2) truncatedTimeLikeShower(theChildren[1],branch);
  else        timeLikeShower(theChildren[1]);
  // branching has happened
  return true;
}

bool NasonEvolver::truncatedSpaceLikeShower(tShowerParticlePtr particle, PPtr beam,
					    NasonBranchingPtr branch) {
  bool vetoed(true);
  Branching bb;
  // generate branching
  tcPDPtr part[2];
  while (vetoed) {
    vetoed=false;
    bb=splittingGenerator()->chooseBackwardBranching(*particle,beam,1.,beamParticle());
    if(!bb.kinematics||bb.kinematics->scale()<branch->_scale) {
      bb=Branching();
      break;
    }
    // particles as in Sudakov form factor
    part[0]=getParticleData(bb.ids[0]);
    part[1]=getParticleData(bb.ids[2]);
    if(particle->id()!=bb.ids[1]) {
      if(part[0]->CC()) part[0]=part[0]->CC();
      if(part[1]->CC()) part[1]=part[1]->CC();
    }
    // apply the vetos for the truncated shower
    // if particle changes type
    if(part[0]->id()!=particle->id()) vetoed=true;
    // if doesn't carry most of momentum
    if(bb.kinematics->z()<0.5) vetoed=true;
    // pt veto
    if(bb.kinematics->pT()>progenitor()->maximumpT()) vetoed = true;
    if(vetoed) particle->setEvolutionScale(ShowerIndex::QCD,bb.kinematics->scale());
  }
  if(!bb.kinematics) {
    double z(0.);
    NasonBranchingPtr timelike;
    for(unsigned int ix=0;ix<branch->_children.size();++ix) {
      if(!branch->_children[ix]->_incoming) {
	timelike=branch->_children[ix];
      }
      if( branch->_children[ix]->_incoming) z   = branch->_children[ix]->_z;
    }
    ShoKinPtr kinematics=
      branch->_sudakov->createInitialStateBranching(branch->_scale,z,branch->_phi,
						    branch->_children[0]->_pt);
    kinematics->initialize(*particle,beam);
   // assign the splitting function and shower kinematics
    particle->setShowerKinematics(kinematics);
    // For the time being we are considering only 1->2 branching
    // Now create the actual particles, make the otherChild a final state
    // particle, while the newParent is not
    ShowerParticlePtr newParent=new_ptr(ShowerParticle(branch->_particle->dataPtr(),
						       false));
    ShowerParticlePtr otherChild = new_ptr(ShowerParticle(timelike->_particle->dataPtr(),
							  true,true));
    ShowerParticleVector theChildren;
    theChildren.push_back(particle); 
    theChildren.push_back(otherChild);
    particle->showerKinematics()->updateParent(newParent, theChildren);
    // update the history if needed
    currentTree()->updateInitialStateShowerProduct(progenitor(),newParent);
    currentTree()->addInitialStateBranching(particle,newParent,otherChild);
    // for the reconstruction of kinematics, parent/child
    // relationships are according to the branching process:
    // now continue the shower
    bool emitted;
    if(branch->_parent) {
      emitted=truncatedSpaceLikeShower(newParent,beam,branch->_parent);
    }
    else {
      emitted=spaceLikeShower(newParent,beam);
    }
    if(!emitted) {
      kinematics->updateLast(newParent);
    }
    particle->showerKinematics()->updateChildren(newParent, theChildren);
    // perform the shower of the final-state particle
    if(timelike->_children.empty()) {
      timeLikeShower(otherChild);
    }
    else {
      truncatedTimeLikeShower(otherChild,timelike);
    }
    // return the emitted
    return true;
  }
  // assign the splitting function and shower kinematics
  particle->setShowerKinematics(bb.kinematics);
  // For the time being we are considering only 1->2 branching
  // Now create the actual particles, make the otherChild a final state
  // particle, while the newParent is not
  ShowerParticlePtr newParent=new_ptr(ShowerParticle(part[0],false));
  ShowerParticlePtr otherChild = new_ptr(ShowerParticle(part[1],true,true));
  ShowerParticleVector theChildren; 
  theChildren.push_back(particle); 
  theChildren.push_back(otherChild);
  particle->showerKinematics()->updateParent(newParent, theChildren);
  // update the history if needed
  currentTree()->updateInitialStateShowerProduct(progenitor(),newParent);
  currentTree()->addInitialStateBranching(particle,newParent,otherChild);
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // now continue the shower
  bool emitted=truncatedSpaceLikeShower(newParent,beam,branch);
  // now reconstruct the momentum
  if(!emitted) bb.kinematics->updateLast(newParent);
  particle->showerKinematics()->updateChildren(newParent, theChildren);
  // perform the shower of the final-state particle
  timeLikeShower(otherChild);
  // return the emitted
  return true;
}
