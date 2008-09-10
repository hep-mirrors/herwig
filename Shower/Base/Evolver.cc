// -*- C++ -*-
//
// Evolver.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Evolver class.
//
#include "Evolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ShowerKinematics.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Utilities/Throw.h"
#include "ShowerTree.h"
#include "ShowerProgenitor.h"
#include "KinematicsReconstructor.h"
#include "PartnerFinder.h"
#include "MECorrectionBase.h"
#include "ThePEG/Handlers/XComb.h"

using namespace Herwig;

IBPtr Evolver::clone() const {
  return new_ptr(*this);
}

IBPtr Evolver::fullclone() const {
  return new_ptr(*this);
}

void Evolver::persistentOutput(PersistentOStream & os) const {
  os << _model << _splittingGenerator << _maxtry 
     << _meCorrMode << _hardVetoMode << _hardVetoRead << _limitEmissions << _ptVetoDefinition
     << ounit(_iptrms,GeV) << _beta << ounit(_gamma,GeV) << ounit(_iptmax,GeV) << _vetoes;
}

void Evolver::persistentInput(PersistentIStream & is, int) {
  is >> _model >> _splittingGenerator >> _maxtry 
     >> _meCorrMode >> _hardVetoMode >> _hardVetoRead >> _limitEmissions >> _ptVetoDefinition
     >> iunit(_iptrms,GeV) >> _beta >> iunit(_gamma,GeV) >> iunit(_iptmax,GeV) >> _vetoes;
}

void Evolver::doinitrun() {
  Interfaced::doinitrun();
  for(unsigned int ix=0;ix<showerModel()->meCorrections().size();++ix) {
    showerModel()->meCorrections()[ix]->evolver(this);
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
    (ifaceMECorrMode,"No","MECorrections off", 0);
  static SwitchOption on
    (ifaceMECorrMode,"Yes","hard+soft on", 1);
  static SwitchOption hard
    (ifaceMECorrMode,"Hard","only hard on", 2);
  static SwitchOption soft
    (ifaceMECorrMode,"Soft","only soft on", 3);

  static Switch<Evolver, unsigned int> ifaceHardVetoMode
    ("HardVetoMode",
     "Choice of the Hard Veto Mode",
     &Evolver::_hardVetoMode, 1, false, false);
  static SwitchOption HVoff
    (ifaceHardVetoMode,"No","hard vetos off", 0);
  static SwitchOption HVon
    (ifaceHardVetoMode,"Yes","hard vetos on", 1);
  static SwitchOption HVIS
    (ifaceHardVetoMode,"Initial", "only IS emissions vetoed", 2);
  static SwitchOption HVFS
    (ifaceHardVetoMode,"Final","only FS emissions vetoed", 3);

  static Switch<Evolver, unsigned int> ifaceHardVetoRead
    ("HardVetoScaleSource",
     "If hard veto scale is to be read",
     &Evolver::_hardVetoRead, 0, false, false);
  static SwitchOption HVRcalc
    (ifaceHardVetoRead,"Calculate","Calculate from hard process", 0);
  static SwitchOption HVRread
    (ifaceHardVetoRead,"Read","Read from XComb->lastScale", 1);

  static Parameter<Evolver, Energy> ifaceiptrms
    ("IntrinsicPtGaussian",
     "RMS of intrinsic pT of Gaussian distribution:\n"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &Evolver::_iptrms, GeV, 0*GeV, 0*GeV, 1000000.0*GeV,
     false, false, Interface::limited);

  static Parameter<Evolver, double> ifacebeta
    ("IntrinsicPtBeta",
     "Proportion of inverse quadratic distribution in generating intrinsic pT.\n"
     "(1-Beta) is the proportion of Gaussian distribution",
     &Evolver::_beta, 0, 0, 1,
     false, false, Interface::limited);

  static Parameter<Evolver, Energy> ifacegamma
    ("IntrinsicPtGamma",
     "Parameter for inverse quadratic:\n"
     "2*Beta*Gamma/(sqr(Gamma)+sqr(intrinsicpT))",
     &Evolver::_gamma,GeV, 0*GeV, 0*GeV, 100000.0*GeV,
     false, false, Interface::limited);

  static Parameter<Evolver, Energy> ifaceiptmax
    ("IntrinsicPtIptmax",
     "Upper bound on intrinsic pT for inverse quadratic",
     &Evolver::_iptmax,GeV, 0*GeV, 0*GeV, 100000.0*GeV,
     false, false, Interface::limited);

  static RefVector<Evolver,ShowerVeto> ifaceVetoes
    ("Vetoes",
     "The vetoes to be checked during showering",
     &Evolver::_vetoes, -1,
     false,false,true,true,false);

  static Switch<Evolver,unsigned int> interfaceLimitEmissions
    ("LimitEmissions",
     "Limit the number and type of emissions for testing",
     &Evolver::_limitEmissions, 0, false, false);
  static SwitchOption interfaceLimitEmissionsNoLimit
    (interfaceLimitEmissions,
     "NoLimit",
     "Allow an arbitrary number of emissions",
     0);
  static SwitchOption interfaceLimitEmissionsOneInitialStateEmission
    (interfaceLimitEmissions,
     "OneInitialStateEmission",
     "Allow one emission in the initial state and none in the final state",
     1);
  static SwitchOption interfaceLimitEmissionsOneFinalStateEmission
    (interfaceLimitEmissions,
     "OneFinalStateEmission",
     "Allow one emission in the final state and none in the initial state",
     2);

  static Switch<Evolver, unsigned int> ifaceJetMeasureMode
    ("JetMeasure",
     "Choice of the jet measure algorithm",
     &Evolver::_ptVetoDefinition, 1, false, false);
  
  static SwitchOption Durham
    (ifaceJetMeasureMode,"Durham","Durham jet algorithm", 0);
  
  static SwitchOption LUCLUS
    (ifaceJetMeasureMode,"LUCLUS","LUCLUS jet algorithm", 1);

  static Parameter<Evolver,double> interfacePtCut
    ("JetCut",
     "The jet cut (in durham or shower definition)",
     &Evolver::_y_cut, 1.1, 0., 1.1,
     false, false, Interface::limited );

  /*
  static Switch<Evolver, unsigned int> ifacePtVetoDefinition
    ("PtVetoDefinition",
     "Choice of the pt deinition used in pt veto",
     &Evolver::_ptVetoDefinition, 0, false, false);
  static SwitchOption PtDefShower
    (ifacePtVetoDefinition,"Shower","shower pt definition", 0);
  static SwitchOption PtDefDurham
    (ifacePtVetoDefinition,"Durham","durham pt", 1);
  */
}

void Evolver::generateIntrinsicpT(vector<ShowerProgenitorPtr> particlesToShower) {
  _intrinsic.clear();
  if ( !ipTon() || !isISRadiationON() ) return;
  // don't do anything for the moment for secondary scatters
  if( !ShowerHandler::currentHandler()->firstInteraction() ) return;
  // generate intrinsic pT
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    // only consider initial-state particles
    if(particlesToShower[ix]->progenitor()->isFinalState()) continue;
    if(!particlesToShower[ix]->progenitor()->dataPtr()->coloured()) continue;
    Energy ipt;
    if(UseRandom::rnd() > _beta) {
      ipt=_iptrms*sqrt(-log(UseRandom::rnd()));
    }
    else {
      ipt=_gamma*sqrt(pow(1.+sqr(_iptmax/_gamma), UseRandom::rnd())-1.);
    }
    pair<Energy,double> pt = make_pair(ipt,UseRandom::rnd(Constants::twopi));
    _intrinsic[particlesToShower[ix]] = pt;
  }
}

void Evolver::setupMaximumScales(ShowerTreePtr hard, 
				 vector<ShowerProgenitorPtr> p) {
  // find out if hard partonic subprocess.
  bool isPartonic(false); 
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
    cit = _currenttree->incomingLines().begin();
  Lorentz5Momentum pcm;
  for(; cit!=hard->incomingLines().end(); ++cit) {
    pcm += cit->first->progenitor()->momentum();
    isPartonic |= cit->first->progenitor()->coloured();
  }
  // find maximum pt from hard process, the maximum pt from all outgoing
  // coloured lines (this is simpler and more general than
  // 2stu/(s^2+t^2+u^2)).  Maximum scale for scattering processes will
  // be transverse mass.
  Energy ptmax = -1.0*GeV;
  // general case calculate the scale  
  if (!hardVetoXComb()) {
    // scattering process
    if(hard->isHard()) {
      // coloured incoming particles
      if (isPartonic) {
	map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cjt = hard->outgoingLines().begin();
	for(; cjt!=hard->outgoingLines().end(); ++cjt) {
	  if (cjt->first->progenitor()->coloured())
	    ptmax = max(ptmax,cjt->first->progenitor()->momentum().mt());
	}
      }
      if (ptmax < 0.0*GeV) ptmax = pcm.m();
    } 
    // decay, incoming() is the decaying particle.
    else { 
      ptmax = hard->incomingLines().begin()->first
	->progenitor()->momentum().mass(); 
    }
  }
  // hepeup.SCALUP is written into the lastXComb by the
  // LesHouchesReader itself - use this by user's choice. 
  // Can be more general than this. 
  else {
    ptmax = sqrt( ShowerHandler::currentHandler()
		  ->lastXCombPtr()->lastScale() );
  }
  // set maxHardPt for all progenitors.  For partonic processes this
  // is now the max pt in the FS, for non-partonic processes or
  // processes with no coloured FS the invariant mass of the IS
  vector<ShowerProgenitorPtr>::const_iterator ckt = p.begin();
  for (; ckt != p.end(); ckt++) (*ckt)->maxHardPt(ptmax);
}

void Evolver::showerHardProcess(ShowerTreePtr hard) {
  // set the current tree
  currentTree(hard);
  // zero number of emissions
  _nis = _nfs = 0;
  // extract particles to shower
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(true);
  // setup the maximum scales for the shower, given by the hard process
  if (hardVetoOn()) 
    setupMaximumScales(_currenttree, particlesToShower);
  // generate the intrinsic p_T once and for all
  generateIntrinsicpT(particlesToShower);
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      currentTree()->clear();
      setColourPartners(true);
      _nis = 0;
      _nfs = 0;
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
      _progenitor=particlesToShower[ix];
      // initial-state
      if(!_progenitor->progenitor()->isFinalState()) {
	if(!isISRadiationON()) continue;
	// get the PDF
	setBeamParticle(_progenitor->beam());
	assert(beamParticle());
	// perform the shower
	// set the beam particle
	tPPtr beamparticle=_progenitor->original();
	if(!beamparticle->parents().empty()) 
	  beamparticle=beamparticle->parents()[0];
	// generate the shower
	_progenitor->hasEmitted(startSpaceLikeShower(beamparticle));
      }
      // final-state
      else {
	if(!isFSRadiationON()) continue;
	// perform shower
	_progenitor->hasEmitted(startTimeLikeShower());
      }
    }
  }
  while(!_model->kinematicsReconstructor()->reconstructHardJets(hard,_intrinsic)&&
	_maxtry>++ntry);
  if(_maxtry==ntry) throw ShowerHandler::ShowerTriesVeto(ntry);
  // the tree has now showered
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
  // see if there is an appropriate matrix element correction
  for(unsigned int ix=0;ix<showerModel()->meCorrections().size();++ix) {
    double initial,final;
    if(!showerModel()->meCorrections()[ix]->canHandle(_currenttree,
						      initial,final,this)) continue;
    if(_currentme) {
      ostringstream output;
      output << "There is more than one possible matrix"
	     << "element which could be applied for ";
      map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
      for(cit=_currenttree->incomingLines().begin();
	  cit!=_currenttree->incomingLines().end();++cit)
	{output << cit->first->progenitor()->PDGName() << " ";}
      output << " -> ";
      map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
      for(cjt=_currenttree->outgoingLines().begin();
	  cjt!=_currenttree->outgoingLines().end();++cjt)
	{output << cjt->first->progenitor()->PDGName() << " ";}
      output << "in Evolver::hardMatrixElementCorrection()\n";
      throw Exception() << output << Exception::runerror;
    }
    else {
      _currentme=showerModel()->meCorrections()[ix];
      _initialenhance = initial;
      _finalenhance   = final;
    }
  }
  // now apply the hard correction
  if(_currentme && hardMEC()) 
    _currentme->applyHardMatrixElementCorrection(_currenttree);
}

bool Evolver::timeLikeShower(tShowerParticlePtr particle) {
  // don't do anything if not needed
  if(_limitEmissions == 1 || ( _limitEmissions == 2 && _nfs != 0) ) return false;
  // generate the emission
  Branching fb;
  while (true) {
    fb=_splittingGenerator->chooseForwardBranching(*particle,_finalenhance);
    // no emission return
    if(!fb.kinematics) return false;
    // if emission OK break
    if(!timeLikeVetoed(fb,particle)) break;
    // otherwise reset scale and continue
    particle->setEvolutionScale(fb.kinematics->scale());
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
  theChildren.push_back(new_ptr(ShowerParticle(pdata[0],true))); 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[1],true))); 
  // update the children
  particle->showerKinematics()->updateChildren(particle, theChildren);
  // update the history if needed
  if(particle==_currenttree->getFinalStateShowerProduct(_progenitor))
    _currenttree->updateFinalStateShowerProduct(_progenitor,
						particle,theChildren);
  _currenttree->addFinalStateBranching(particle,theChildren);
  // update number of emissions
  ++_nfs;
  if(_limitEmissions!=0) return true;
  // shower the first  particle
  timeLikeShower(theChildren[0]);
  // shower the second particle
  timeLikeShower(theChildren[1]);
  // branching has happened
  return true;
}

bool 
Evolver::spaceLikeShower(tShowerParticlePtr particle, PPtr beam) {
  // don't do anything if not needed
  if(_limitEmissions == 2  || 
     ( _limitEmissions == 1 && _nis != 0 ) ) return false;
  Branching bb;
  // generate branching
  while (true) {
    bb=_splittingGenerator->chooseBackwardBranching(*particle,beam,
						    _initialenhance,_beam);
    // return if no emission
    if(!bb.kinematics) return false;
    // if not vetoed break
    if(!spaceLikeVetoed(bb,particle)) break;
    // otherwise reset scale and continue
    particle->setEvolutionScale(bb.kinematics->scale());
  }
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
  ShowerParticlePtr newParent=new_ptr(ShowerParticle(part[0],false));
  ShowerParticlePtr otherChild = new_ptr(ShowerParticle(part[1],true,true));
  ShowerParticleVector theChildren;
  theChildren.push_back(particle); 
  theChildren.push_back(otherChild); 
  particle->showerKinematics()->updateParent(newParent, theChildren);
  // update the history if needed
  _currenttree->updateInitialStateShowerProduct(_progenitor,newParent);
  _currenttree->addInitialStateBranching(particle,newParent,otherChild);
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // now continue the shower
  ++_nis;
  bool emitted = _limitEmissions==0 ? spaceLikeShower(newParent,beam) : false;
  // now reconstruct the momentum
  if(!emitted) {
    if(_intrinsic.find(_progenitor)==_intrinsic.end()) {
      bb.kinematics->updateLast(newParent,0.*MeV,0.*MeV);
    }
    else {
      pair<Energy,double> kt=_intrinsic[_progenitor];
      bb.kinematics->updateLast(newParent,
				kt.first*cos(kt.second),
				kt.first*sin(kt.second));
    }
  }
  particle->showerKinematics()->updateChildren(newParent, theChildren);
  if(_limitEmissions!=0) return true;
  // perform the shower of the final-state particle
  timeLikeShower(otherChild);
  // return the emitted
  return true;
}

void Evolver::showerDecay(ShowerTreePtr decay) {
  // set the ShowerTree to be showered
  currentTree(decay);
  // extract particles to be shower, set scales and perform hard matrix element 
  // correction
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(false);
  setupMaximumScales(_currenttree, particlesToShower);
  // compute the minimum mass of the final-state
  Energy minmass(0.*MeV);
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if(particlesToShower[ix]->progenitor()->isFinalState())
      minmass+=particlesToShower[ix]->progenitor()->mass();
  }
  // main showering loop
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      _currenttree->clear();
      setColourPartners(false);
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
      _progenitor=particlesToShower[ix];
      // final-state radiation
      if(_progenitor->progenitor()->isFinalState()) {
	if(!isFSRadiationON()) continue;
	// perform shower
	_progenitor->hasEmitted(startTimeLikeShower());
      }
      // initial-state radiation
      else {
	if(!isISRadiationON()) continue;
	// perform shower
	// set the scales correctly. The current scale is the maximum scale for
	// emission not the starting scale
	Energy maxscale=_progenitor->progenitor()->evolutionScale();
	Energy startScale=_progenitor->progenitor()->mass();
	_progenitor->progenitor()->setEvolutionScale(startScale);
	// perform the shower
	_progenitor->hasEmitted(spaceLikeDecayShower(_progenitor->progenitor(),
						     maxscale,minmass)); 
      }
    }
  }
  while(!showerModel()->kinematicsReconstructor()->reconstructDecayJets(decay)&&
	maximumTries()>++ntry);
  if(maximumTries()==ntry) throw Exception() << "Failed to generate the shower after "
				      << ntry << " attempts in Evolver::showerDecay()"
				      << Exception::eventerror;
  // tree has now showered
  _currenttree->hasShowered(true);
}

bool Evolver::spaceLikeDecayShower(tShowerParticlePtr particle,
				   Energy maxscale,
				   Energy minmass) {
  Branching fb;
  while (true) {
    fb=_splittingGenerator->chooseDecayBranching(*particle,maxscale,minmass,
						 _initialenhance);
    // return if no radiation
    if(!fb.kinematics) return false;
    // if not vetoed break
    if(!spaceLikeDecayVetoed(fb,particle)) break;
    // otherwise reset scale and continue
    particle->setEvolutionScale(fb.kinematics->scale());
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
  theChildren.push_back(new_ptr(ShowerParticle(pdata[0],true))); 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[1],true))); 
  // some code moved to updateChildren
  particle->showerKinematics()->updateChildren(particle, theChildren);
  // In the case of splittings which involves coloured particles,
  // set properly the colour flow of the branching.
  // update the history if needed
  _currenttree->updateInitialStateShowerProduct(_progenitor,theChildren[0]);
  _currenttree->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
  // shower the first  particle
  spaceLikeDecayShower(theChildren[0],maxscale,minmass);
  // shower the second particle
  timeLikeShower(theChildren[1]);
  // branching has happened
  return true;
}

vector<ShowerProgenitorPtr> Evolver::setupShower(bool hard) {
  // put all the particles into a data structure which has the particles
  // and the maximum pt for emission from the particle
  // set the initial colour partners
  setColourPartners(hard);
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  // generate the hard matrix element correction if needed
  hardMatrixElementCorrection();
  // get the particles to be showered
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=_currenttree->incomingLines().begin();
      cit!=_currenttree->incomingLines().end();++cit)
    particlesToShower.push_back((*cit).first);
  // outgoing particles
  for(cjt=_currenttree->outgoingLines().begin();
      cjt!=_currenttree->outgoingLines().end();++cjt)
    particlesToShower.push_back(((*cjt).first));
  // remake the colour partners if needed
  if(_currenttree->hardMatrixElementCorrection()) {
    setColourPartners(hard);
    _currenttree->resetShowerProducts();
  }
  return particlesToShower;
}

void Evolver::setColourPartners(bool hard) {
  vector<ShowerParticlePtr> particles;
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cit=_currenttree->incomingLines().begin();
      cit!=_currenttree->incomingLines().end();++cit)
    particles.push_back(cit->first->progenitor());
  assert((particles.size()==1&&!hard)||(particles.size()==2&&hard));
  // outgoing particles
  for(cjt=_currenttree->outgoingLines().begin();
      cjt!=_currenttree->outgoingLines().end();++cjt)
    particles.push_back(cjt->first->progenitor());
  // Set the initial evolution scales
  showerModel()->partnerFinder()->setInitialEvolutionScales(particles,!hard);
}

bool Evolver::startTimeLikeShower() {
  return timeLikeShower(_progenitor->progenitor());
}

bool Evolver::startSpaceLikeShower(PPtr parent) {
  return spaceLikeShower(_progenitor->progenitor(),parent);
}

bool Evolver::timeLikeVetoed(const Branching & fb,
			     ShowerParticlePtr particle) {
  // check whether emission was harder than largest pt of hard subprocess
  if (hardVetoFS() && fb.kinematics->pT() > _progenitor->maxHardPt()) 
    return true;
  // soft matrix elementr correction veto
  if(_currentme && softMEC() &&
     _currentme->softMatrixElementVeto(_progenitor,particle,fb))
    return true;
  // veto on maximum pt
  Energy ptVeto;
  // if set use the y_cut pt veto with the appropriate pt_veto definition
  if( _y_cut > 1. ) ptVeto = _progenitor->maximumpT();
  else              ptVeto = sqrt(ShowerHandler::currentHandler()
				  ->lastXCombPtr()->lastS() * _y_cut );
  if( fb.kinematics && _ptVetoDefinition == 0 ) 
    ptVeto *= max( fb.kinematics->z(), 1. - fb.kinematics->z() );
  if(fb.kinematics->pT() > ptVeto)
    return true;
  // general vetos
  if (fb.kinematics && !_vetoes.empty()) {
    bool vetoed=false;
    for (vector<ShowerVetoPtr>::iterator v = _vetoes.begin();
	 v != _vetoes.end(); ++v) {
      bool test = (**v).vetoTimeLike(_progenitor,particle,fb);
      switch((**v).vetoType()) {
      case ShowerVeto::Emission:
	vetoed |= test;
	break;
      case ShowerVeto::Shower:
	if(test) throw VetoShower();
	break;
      case ShowerVeto::Event:
	if(test) throw Veto();
	break;
      }
    }
    if(vetoed) return true;
  }
  return false;
}

bool Evolver::spaceLikeVetoed(const Branching & bb,ShowerParticlePtr particle) {
  // check whether emission was harder than largest pt of hard subprocess
  if (hardVetoIS() && bb.kinematics->pT() > _progenitor->maxHardPt())
    return true;
  // apply the soft correction
  if(_currentme && softMEC() &&
     _currentme->softMatrixElementVeto(_progenitor,particle,bb))
    return true;
  // check vs max pt for the shower
  Energy ptVeto;
  //if set use the y_cut pt veto with the appropriate pt_veto definition
  if( _y_cut > 1. ) ptVeto = _progenitor->maximumpT();
  else              ptVeto = sqrt( ShowerHandler::currentHandler()
				   ->lastXCombPtr()->lastS() * _y_cut );
  if( bb.kinematics && _ptVetoDefinition == 0 ) 
    ptVeto *= max( bb.kinematics->z(), 1. - bb.kinematics->z() );
  if(bb.kinematics->pT()> ptVeto )
    return true;
  // the more general vetos
  if (!_vetoes.empty()) {
    bool vetoed=false;
    for (vector<ShowerVetoPtr>::iterator v = _vetoes.begin();
	 v != _vetoes.end(); ++v) {
      bool test = (**v).vetoSpaceLike(_progenitor,particle,bb);
      switch ((**v).vetoType()) {
      case ShowerVeto::Emission:
	vetoed |= test;
	break;
      case ShowerVeto::Shower:
	if(test) throw VetoShower();
	break;
      case ShowerVeto::Event:
	if(test) throw Veto();
	break;
      }
    }
    if (vetoed) return true;
  }
  return false;
}

bool Evolver::spaceLikeDecayVetoed(const Branching & fb,
				     ShowerParticlePtr particle) {
  // apply the soft correction
  if(_currentme && softMEC() &&
     _currentme->softMatrixElementVeto(_progenitor,particle,fb))
    return true;
  // veto on hardest pt in the shower
  Energy ptVeto = _progenitor->maximumpT();
  if( fb.kinematics && _ptVetoDefinition == 0 ) 
    ptVeto *= max( fb.kinematics->z(), 1. - fb.kinematics->z() );
  if(fb.kinematics->pT()> ptVeto) return true;
  // general vetos
  if (!_vetoes.empty()) {
    bool vetoed=false;
    for (vector<ShowerVetoPtr>::iterator v = _vetoes.begin();
	 v != _vetoes.end(); ++v) {
      bool test = (**v).vetoSpaceLike(_progenitor,particle,fb);
      switch((**v).vetoType()) {
      case ShowerVeto::Emission:
	vetoed |= test;
	break;
      case ShowerVeto::Shower:
	if(test) VetoShower();
	break;
      case ShowerVeto::Event:
	if(test) throw Veto();
	break;
      }
      if (vetoed) return true;
    }
  }
  return false;
}
