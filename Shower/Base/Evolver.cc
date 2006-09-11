// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Evolver class.
//
#include "Evolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "Herwig++/Shower/Default/QtildaShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Utilities/Timer.h"
#include "ShowerTree.h"
#include "ShowerProgenitor.h"
#include "KinematicsReconstructor.h"
#include "PartnerFinder.h"
#include "MECorrectionBase.h"

using namespace Herwig;

Evolver::~Evolver() {}

void Evolver::persistentOutput(PersistentOStream & os) const {
  os << _model << _splittingGenerator << _maxtry << _meCorrMode;
}

void Evolver::persistentInput(PersistentIStream & is, int) {
  is >> _model >> _splittingGenerator >> _maxtry >> _meCorrMode;
}

void Evolver::doinitrun() {
  Interfaced::doinitrun();
  for(unsigned int ix=0;ix<_model->meCorrections().size();++ix) {
    _model->meCorrections()[ix]->evolver(this);
  }
}

ClassDescription<Evolver> Evolver::initEvolver;
// Definition of the static class description member.

void Evolver::Init() {

  static ClassDocumentation<Evolver> documentation
    ("This class is responsible for carrying out the showering,",
     "including the kinematics reconstruction, in a given scale range.");

  static Reference<Evolver,SplittingGenerator> 
    interfaceSplitGen("SplittingGenerator", 
		      "A reference to the SplittingGenerator object", 
		      &Herwig::Evolver::_splittingGenerator,
		      false, false, true, false);

  static Reference<Evolver,ShowerModel> interfaceShowerModel
    ("ShowerModel",
     "The pointer to the object which defines the shower evolution model.",
     &Evolver::_model, false, false, true, false, false);

  static Parameter<Evolver,unsigned int> interfaceMaxTry
    ("MaxTry",
     "The maximum number of attempts to generate the shower from a"
     " particular ShowerTree",
     &Evolver::_maxtry, 100, 1, 1000,
     false, false, Interface::limited);

  static Switch<Evolver, unsigned int> ifaceMECorrMode
    ("MECorrMode",
     "Choice of the ME Correction Mode",
     &Evolver::_meCorrMode, 1, false, false);
  static SwitchOption off
    (ifaceMECorrMode,"MEC-off","MECorrections off", 0);
  static SwitchOption on
    (ifaceMECorrMode,"MEC-on","hard+soft on", 1);
  static SwitchOption hard
    (ifaceMECorrMode,"MEC-hard","only hard on", 2);
  static SwitchOption soft
    (ifaceMECorrMode,"MEC-soft","only soft on", 3);

}

void Evolver::showerHardProcess(ShowerTreePtr hard) {
  // set the current tree
  _currenttree=hard;
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(true);
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      _currenttree->clear();
      setColourPartners(true);
    }
    // initial-state radiation
    if(_splittingGenerator->isISRadiationON()) {
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	// only consider initial-state particles
	if(particlesToShower[ix]->progenitor()->isFinalState()) continue;
	// get the PDF
	_beam=dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>
	  (particlesToShower[ix]->original()->parents()[0]->dataPtr());
	if(!_beam) throw Exception() << "The Beam particle does not have"
				     << " BeamParticleData in Evolver::" 
				     << "showerhardProcess()" 
				     << Exception::runerror;
	// perform the shower
	_progenitor=particlesToShower[ix];
	_progenitor->hasEmitted(spaceLikeShower(particlesToShower[ix]->progenitor()));
      }
    }
    // final-state radiation
    if(_splittingGenerator->isFSRadiationON()) {
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	// only consider final-state particles
	if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
	// perform shower
	_progenitor=particlesToShower[ix];
	_progenitor->hasEmitted(timeLikeShower(particlesToShower[ix]->progenitor()));
      }
    }
  }
  while(!_model->kinematicsReconstructor()->reconstructHardJets(hard)&&_maxtry>++ntry);
  if(_maxtry==ntry) throw Exception() << "Failed to generate the shower after "
				      << ntry 
				      << " attempts in Evolver::showerHardProcess()"
				      << Exception::eventerror;
  _currenttree->hasShowered(true);
}

void Evolver::hardMatrixElementCorrection() {
  // set me correction to null pointer
  _currentme=MECorrectionPtr();
  // set the initial enhancement factors for the soft correction
  _initialenhance=1.;
  _finalenhance  =1.;
  // if hard matrix element switched off return
  if(!MECOn()) return;
  // see if there is an appropraite matrix element correction
  for(unsigned int ix=0;ix<_model->meCorrections().size();++ix) {
    double initial,final;
    if(!_model->meCorrections()[ix]->canHandle(_currenttree,
					       initial,final)) continue;
    if(_currentme) {
      ostringstream output;
      output << "There is more than one possible matrix"
	     << "element which could be applied for ";
      map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
      for(cit=_currenttree->incomingLines().begin();
	  cit!=_currenttree->incomingLines().end();++cit)
	{output << cit->first->progenitor()->PDGName() << " ";}
      output << " -> ";
      for(cit=_currenttree->outgoingLines().begin();
	  cit!=_currenttree->outgoingLines().end();++cit)
	{output << cit->first->progenitor()->PDGName() << " ";}
      output << "in Evolver::hardMatrixElementCorrection()\n";
      throw Exception() << output << Exception::runerror;
    }
    else {
      _currentme=_model->meCorrections()[ix];
      _initialenhance = initial;
      _finalenhance   = final;
    }
  }
  // if no suitable me correction
  if(!_currentme) return; 
  // now apply the hard correction
  if(hardMEC()) _currentme->applyHardMatrixElementCorrection(_currenttree);
}

bool Evolver::timeLikeShower(tShowerParticlePtr particle) {
  Timer<1005> timer("Evolver::timeLikeShower");
  bool vetoed = true;
  Branching fb;
  while (vetoed) {
    vetoed = false; 
    fb=_splittingGenerator->chooseForwardBranching(*particle,_finalenhance);
      // apply vetos if needed
    if(fb.kinematics && _currentme && softMEC())
      vetoed=_currentme->softMatrixElementVeto(_progenitor,particle,fb);
  }
  // if no branching set decay matrix and return
  if(!fb.kinematics) {
    //
    // add decay matrix stuff here
    //
    return false;
  }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->setShowerKinematics(fb.kinematics);
  // Assign the splitting function to the emitting particle. 
  // For the time being we are considering only 1->2 branching
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  // if same as definition create particles, otherwise create cc
  tcPDPtr pdata[2];
  for(unsigned int ix=0;ix<2;++ix) pdata[ix]=getParticleData(fb.ids[ix+1]);
  if(particle->id()!=fb.ids[0]) {
    for(unsigned int ix=0;ix<2;++ix) {
      tPDPtr cc(pdata[ix]->CC());
      if(cc) pdata[ix]=cc;
    }
  }
  ShowerParticleVector theChildren; 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[0]))); 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[1]))); 
  // some code moved to updateChildren
  particle->showerKinematics()->updateChildren(particle, theChildren);
  // update the history if needed
  if(particle==_currenttree->getFinalStateShowerProduct(_progenitor))
    _currenttree->updateFinalStateShowerProduct(_progenitor,
						particle,theChildren);
  _currenttree->addFinalStateBranching(particle,theChildren);
  // shower the first  particle
  //
  //  need to set rho here
  //
  timeLikeShower(theChildren[0]);
  // shower the second particle
  //
  //   need to set rho here
  //
  timeLikeShower(theChildren[1]);
  // branching has happened
  return true;
}

bool Evolver::spaceLikeShower(tShowerParticlePtr particle) {
  Timer<1006> timer("Evolver::spaceLikeShower");
  bool vetoed(true);
  Branching bb;
  // generate branching
  while (vetoed) {
    vetoed=false;
    bb=_splittingGenerator->chooseBackwardBranching(*particle,_initialenhance,_beam);
    // apply the soft correction
    if(bb.kinematics && _currentme && softMEC())
      vetoed=_currentme->softMatrixElementVeto(_progenitor,particle,bb);
  }
  if(!bb.kinematics) return false;
  // assign the splitting function and shower kinematics
  particle->setShowerKinematics(bb.kinematics);
  // For the time being we are considering only 1->2 branching
  // particles as in Sudakov form factor
  tcPDPtr part[2]={getParticleData(bb.ids[0]),
		   getParticleData(bb.ids[2])};
  if(particle->id()!=bb.ids[1]) {
    if(part[0]->CC()) part[0]=part[0]->CC();
    if(part[1]->CC()) part[1]=part[1]->CC();
  }
  // Now create the actual particles, make the otherChild a final state
  // particle, while the newParent is not
  ShowerParticlePtr newParent=new_ptr(ShowerParticle(part[0]));
  newParent->setFinalState(false);
  ShowerParticlePtr otherChild = new_ptr(ShowerParticle(part[1]));
  otherChild->setFinalState(true);
  otherChild->setInitiatesTLS(true);
  ShowerParticleVector theChildren; 
  theChildren.push_back(particle); 
  theChildren.push_back(otherChild); 
  particle->showerKinematics()->updateChildren(newParent, theChildren);
  // update the history if needed
  _currenttree->updateInitialStateShowerProduct(_progenitor,newParent);
  _currenttree->addInitialStateBranching(particle,newParent,otherChild);
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // now continue the shower
  bool emitted=spaceLikeShower(newParent);
  //bool emitted=false;
  // now reconstruct the momentum
  if(!emitted) bb.kinematics->updateLast(newParent,0);
  particle->showerKinematics()->updateChildren(newParent, theChildren);
  // perform the shower of the final-state particle
  timeLikeShower(otherChild);
  // return the emitted
  return true;
}

void Evolver::showerDecay(ShowerTreePtr decay) {
  // set the ShowerTree to be showered
  _currenttree=decay;
  // extract particles to be shower, set scales and perform hard matrix element 
  // correction
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(false);
  // main showering loop
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      _currenttree->clear();
      setColourPartners(false);
    }
    // initial-state radiation
    if(_splittingGenerator->isISRadiationON()) {
      // compute the minimum mass of the final-state
      Energy minmass(0.);
      for(unsigned int ix=0;ix<particlesToShower.size();++ix)
	{if(particlesToShower[ix]->progenitor()->isFinalState())
	    minmass+=particlesToShower[ix]->progenitor()->mass();}
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	// only consider initial-state particles
	if(particlesToShower[ix]->progenitor()->isFinalState()) continue;
	// perform shower
	_progenitor=particlesToShower[ix];
	// set the scales correctly. The current scale is the maximum scale for
	// emission not the starting scale
	vector<Energy> maxscale=_progenitor->progenitor()->evolutionScales();
	Energy startScale=_progenitor->progenitor()->mass();
	_progenitor->progenitor()->setEvolutionScale(ShowerIndex::QCD,startScale);
	_progenitor->progenitor()->setEvolutionScale(ShowerIndex::QED,startScale);
	_progenitor->progenitor()->setEvolutionScale(ShowerIndex::EWK,startScale);
	// perform the shower
	_progenitor->hasEmitted(spaceLikeDecayShower(_progenitor->progenitor(),
						     maxscale,minmass)); 
      }
    }
    // final-state radiation
    if(_splittingGenerator->isFSRadiationON()) {
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	// only consider final-state particles
	if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
	// perform shower
	_progenitor=particlesToShower[ix];
	_progenitor->hasEmitted(timeLikeShower(particlesToShower[ix]->progenitor()));
      }
    }
  }
  while(!_model->kinematicsReconstructor()->reconstructDecayJets(decay)&&_maxtry>++ntry);
  if(_maxtry==ntry) throw Exception() << "Failed to generate the shower after "
				      << ntry << " attempts in Evolver::showerDecay()"
				      << Exception::eventerror;
  _currenttree->hasShowered(true);
}

bool Evolver::spaceLikeDecayShower(tShowerParticlePtr particle,vector<Energy> maxscale,
				   Energy minmass) {
  bool vetoed = true;
  Branching fb;
  while (vetoed) 
    {
      vetoed = false;
      fb=_splittingGenerator->chooseDecayBranching(*particle,maxscale,minmass,
						   _initialenhance);
      // apply the soft correction
      if(fb.kinematics && _currentme && softMEC())
	vetoed=_currentme->softMatrixElementVeto(_progenitor,particle,fb);
    }
  // if no branching set decay matrix and return
  if(!fb.kinematics) {
    //
    // add decay matrix stuff here
    //
    return false;
  }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->setShowerKinematics(fb.kinematics);
  // For the time being we are considering only 1->2 branching
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  // if same as definition create particles, otherwise create cc
  tcPDPtr pdata[2];
  for(unsigned int ix=0;ix<2;++ix) pdata[ix]=getParticleData(fb.ids[ix+1]);
  if(particle->id()!=fb.ids[0]) {
    for(unsigned int ix=0;ix<2;++ix) {
	tPDPtr cc(pdata[ix]->CC());
	if(cc) pdata[ix]=cc;
      }
  }
  ShowerParticleVector theChildren; 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[0]))); 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[1]))); 
  // some code moved to updateChildren
  particle->showerKinematics()->updateChildren(particle, theChildren);
  // In the case of splittings which involves coloured particles,
  // set properly the colour flow of the branching.
  // update the history if needed
  _currenttree->updateInitialStateShowerProduct(_progenitor,theChildren[0]);
  _currenttree->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
  // shower the first  particle
  //
  //  need to set rho here
  //
  spaceLikeDecayShower(theChildren[0],maxscale,minmass);
  // shower the second particle
  //
  //   need to set rho here
  //
  timeLikeShower(theChildren[1]);
  // branching has happened
  return true;
}

vector<ShowerProgenitorPtr> Evolver::setupShower(bool hard) {
  // put all the particles into a data structure which has the particles
  // and the maximum pt for emission from the particle
  // set the initial colour partners
  setColourPartners(hard);
  // generate the hard matrix element correction if needed
  hardMatrixElementCorrection();
  // get the particles to be showered
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=_currenttree->incomingLines().begin();
      cit!=_currenttree->incomingLines().end();++cit)
    particlesToShower.push_back(((*cit).first));
  assert((particlesToShower.size()==1&&!hard)||(particlesToShower.size()==2&&hard));
  // outgoing particles
  for(cit=_currenttree->outgoingLines().begin();
      cit!=_currenttree->outgoingLines().end();++cit)
    particlesToShower.push_back(((*cit).first));
  // remake the colour partners if needed
  if(_currenttree->hardMatrixElementCorrection()) {
    setColourPartners(hard);
    _currenttree->resetShowerProducts();
  }
  return particlesToShower;
}

void Evolver::setColourPartners(bool hard) {
  vector<ShowerParticlePtr> particles;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  for(cit=_currenttree->incomingLines().begin();
      cit!=_currenttree->incomingLines().end();++cit)
    particles.push_back(cit->first->progenitor());
  assert((particles.size()==1&&!hard)||(particles.size()==2&&hard));
  // outgoing particles
  for(cit=_currenttree->outgoingLines().begin();
      cit!=_currenttree->outgoingLines().end();++cit)
    particles.push_back(cit->first->progenitor());
  // Set the initial evolution scales
  if(_splittingGenerator->isInteractionON(ShowerIndex::QCD))
    _model->partnerFinder()->setQCDInitialEvolutionScales(particles,!hard);
  if(_splittingGenerator->isInteractionON(ShowerIndex::QED))
    _model->partnerFinder()->setQEDInitialEvolutionScales(particles,!hard);
  if(_splittingGenerator->isInteractionON(ShowerIndex::EWK))
    _model->partnerFinder()->setEWKInitialEvolutionScales(particles,!hard);
}
