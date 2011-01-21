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
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Shower/ShowerHandler.h" 

using namespace Herwig;

IBPtr Evolver::clone() const {
  return new_ptr(*this);
}

IBPtr Evolver::fullclone() const {
  return new_ptr(*this);
}

void Evolver::persistentOutput(PersistentOStream & os) const {
  os << _model << _splittingGenerator << _maxtry 
     << _meCorrMode << _hardVetoMode << _hardVetoRead << _limitEmissions
     << ounit(_iptrms,GeV) << _beta << ounit(_gamma,GeV) << ounit(_iptmax,GeV) 
     << _vetoes << _hardonly << _trunc_Mode << _hardEmissionMode;
}

void Evolver::persistentInput(PersistentIStream & is, int) {
  is >> _model >> _splittingGenerator >> _maxtry 
     >> _meCorrMode >> _hardVetoMode >> _hardVetoRead >> _limitEmissions
     >> iunit(_iptrms,GeV) >> _beta >> iunit(_gamma,GeV) >> iunit(_iptmax,GeV) 
     >> _vetoes >> _hardonly >> _trunc_Mode >> _hardEmissionMode;
}

ClassDescription<Evolver> Evolver::initEvolver;
// Definition of the static class description member.

void Evolver::Init() {
  
  static ClassDocumentation<Evolver> documentation
    ("This class is responsible for carrying out the showering,",
     "including the kinematics reconstruction, in a given scale range,"
     "including the option of the POWHEG approach to simulated next-to-leading order"
     " radiation\\cite{Nason:2004rx}.",
     "%\\cite{Nason:2004rx}\n"
     "\\bibitem{Nason:2004rx}\n"
     "  P.~Nason,\n"
     "  ``A new method for combining NLO QCD with shower Monte Carlo algorithms,''\n"
     "  JHEP {\\bf 0411} (2004) 040\n"
     "  [arXiv:hep-ph/0409146].\n"
     "  %%CITATION = JHEPA,0411,040;%%\n");

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
     &Evolver::_iptrms, GeV, ZERO, ZERO, 1000000.0*GeV,
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
     &Evolver::_gamma,GeV, ZERO, ZERO, 100000.0*GeV,
     false, false, Interface::limited);

  static Parameter<Evolver, Energy> ifaceiptmax
    ("IntrinsicPtIptmax",
     "Upper bound on intrinsic pT for inverse quadratic",
     &Evolver::_iptmax,GeV, ZERO, ZERO, 100000.0*GeV,
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
  static SwitchOption interfaceLimitEmissionsHardOnly
    (interfaceLimitEmissions,
     "HardOnly",
     "Only allow radiation from the hard ME correction",
     3);

  static Switch<Evolver,bool> interfaceHardOnly
    ("HardOnly",
     "Only generate the emission supplied by the hardest emission"
     " generator, for testing only.",
     &Evolver::_hardonly, false, false, false);
  static SwitchOption interfaceHardOnlyNo
    (interfaceHardOnly,
     "No",
     "Generate full shower",
     false);
  static SwitchOption interfaceHardOnlyYes
    (interfaceHardOnly,
     "Yes",
     "Only the hardest emission",
     true);

  static Switch<Evolver,bool> interfaceTruncMode
    ("TruncatedShower", "Include the truncated shower?", 
     &Evolver::_trunc_Mode, 1, false, false);
  static SwitchOption interfaceTruncMode0
    (interfaceTruncMode,"No","Truncated Shower is OFF", 0);
  static SwitchOption interfaceTruncMode1
    (interfaceTruncMode,"Yes","Truncated Shower is ON", 1);

  static Switch<Evolver,unsigned int> interfaceHardEmissionMode
    ("HardEmissionMode",
     "Whether to use ME corrections or POWHEG for the hardest emission",
     &Evolver::_hardEmissionMode, 0, false, false);
  static SwitchOption interfaceHardEmissionModeMECorrection
    (interfaceHardEmissionMode,
     "MECorrection",
     "Old fashioned ME correction",
     0);
  static SwitchOption interfaceHardEmissionModePOWHEG
    (interfaceHardEmissionMode,
     "POWHEG",
     "Powheg style hard emission",
     1);

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
      if (ptmax < ZERO) ptmax = pcm.m();
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

void Evolver::showerHardProcess(ShowerTreePtr hard, XCPtr xcomb) {
  _hardme = HwMEBasePtr();
  // extract the matrix element
  tStdXCombPtr lastXC = dynamic_ptr_cast<tStdXCombPtr>(xcomb);
  if(lastXC) {
    _hardme = dynamic_ptr_cast<HwMEBasePtr>(lastXC->matrixElement());
  }
  _decayme = HwDecayerBasePtr();
  // set the current tree
  currentTree(hard);
  // zero number of emissions
  _nis = _nfs = 0;
  // extract particles to shower
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(true);
  // setup the maximum scales for the shower, given by the hard process
  if (hardVetoOn()) setupMaximumScales(currentTree(), particlesToShower);
  // generate the intrinsic p_T once and for all
  generateIntrinsicpT(particlesToShower);
  // main shower loop
  unsigned int ntry(0);
  do {
    // clear results of last attempt if needed
    if(ntry!=0) {
      currentTree()->clear();
      setEvolutionPartners(true,ShowerInteraction::QCD);
      _nis = _nfs = 0;
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
	tPPtr beamparticle=progenitor()->original();
	if(!beamparticle->parents().empty()) 
	  beamparticle=beamparticle->parents()[0];
	// generate the shower
	progenitor()->hasEmitted(startSpaceLikeShower(beamparticle,
						      ShowerInteraction::QCD));
      }
      // final-state
      else {
	if(!isFSRadiationON()) continue;
	// perform shower
	progenitor()->hasEmitted(startTimeLikeShower(ShowerInteraction::QCD));
      }
    }
  }
  while(!showerModel()->kinematicsReconstructor()->
	reconstructHardJets(hard,intrinsicpT())&&
	maximumTries()>++ntry);
  _hardme=HwMEBasePtr();
  _hardtree=HardTreePtr();
  if(_maxtry==ntry) throw ShowerHandler::ShowerTriesVeto(ntry);
  
  // the tree has now showered
  _currenttree->hasShowered(true);
}

void Evolver::hardMatrixElementCorrection(bool hard) {
  // set the initial enhancement factors for the soft correction
  _initialenhance = 1.;
  _finalenhance   = 1.;
  // if hard matrix element switched off return
  if(!MECOn()) return;
  // see if we can get the correction from the matrix element
  // or decayer
  if(hard) {
    if(_hardme&&_hardme->hasMECorrection()) {
      _hardme->initializeMECorrection(_currenttree,
				      _initialenhance,_finalenhance);
      if(hardMEC())
	_hardme->applyHardMatrixElementCorrection(_currenttree);
    }
  }
  else {
    if(_decayme&&_decayme->hasMECorrection()) {
      _decayme->initializeMECorrection(_currenttree,
				       _initialenhance,_finalenhance);
      if(hardMEC())
	_decayme->applyHardMatrixElementCorrection(_currenttree);
    }
  }
}

bool Evolver::timeLikeShower(tShowerParticlePtr particle, 
			     ShowerInteraction::Type type) {
  // don't do anything if not needed
  if(_limitEmissions == 1 || _limitEmissions == 3 || 
     ( _limitEmissions == 2 && _nfs != 0) ) return false;
  // generate the emission
  Branching fb;
  while (true) {
    fb=_splittingGenerator->chooseForwardBranching(*particle,_finalenhance,type);
    // no emission return
    if(!fb.kinematics) return false;
    // if emission OK break
    if(!timeLikeVetoed(fb,particle)) break;
    // otherwise reset scale and continue - SO IS involved in veto algorithm
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
  particle->showerKinematics()->updateChildren(particle, theChildren,true);
  // update the history if needed
  if(particle==_currenttree->getFinalStateShowerProduct(_progenitor))
    _currenttree->updateFinalStateShowerProduct(_progenitor,
						particle,theChildren);
  _currenttree->addFinalStateBranching(particle,theChildren);
  // update number of emissions
  ++_nfs;
  if(_limitEmissions!=0) return true;
  // shower the first  particle
  timeLikeShower(theChildren[0],type);
  // shower the second particle
  timeLikeShower(theChildren[1],type);
  // branching has happened
  return true;
}

bool 
Evolver::spaceLikeShower(tShowerParticlePtr particle, PPtr beam,
			 ShowerInteraction::Type type) {
  //using the pdf's associated with the ShowerHandler assures, that
  //modified pdf's are used for the secondary interactions via 
  //CascadeHandler::resetPDFs(...)
  tcPDFPtr pdf;
  if(ShowerHandler::currentHandler()->firstPDF().particle() == _beam)
    pdf = ShowerHandler::currentHandler()->firstPDF().pdf();
  if(ShowerHandler::currentHandler()->secondPDF().particle() == _beam)
    pdf = ShowerHandler::currentHandler()->secondPDF().pdf();
  Energy freeze = ShowerHandler::currentHandler()->pdfFreezingScale();
  // don't do anything if not needed
  if(_limitEmissions == 2  || _limitEmissions == 3  ||
     ( _limitEmissions == 1 && _nis != 0 ) ) return false;
  Branching bb;
  // generate branching
  while (true) {
    bb=_splittingGenerator->chooseBackwardBranching(*particle,beam,
						    _initialenhance,
						    _beam,type,
						    pdf,freeze);
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
  //this updates the evolution scale
  particle->showerKinematics()->updateParent(newParent, theChildren,true);
  // update the history if needed
  _currenttree->updateInitialStateShowerProduct(_progenitor,newParent);
  _currenttree->addInitialStateBranching(particle,newParent,otherChild);
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // now continue the shower
  ++_nis;
  bool emitted = _limitEmissions==0 ? 
    spaceLikeShower(newParent,beam,type) : false;
  // now reconstruct the momentum
  if(!emitted) {
    if(_intrinsic.find(_progenitor)==_intrinsic.end()) {
      bb.kinematics->updateLast(newParent,ZERO,ZERO);
    }
    else {
      pair<Energy,double> kt=_intrinsic[_progenitor];
      bb.kinematics->updateLast(newParent,
				kt.first*cos(kt.second),
				kt.first*sin(kt.second));
    }
  }
  particle->showerKinematics()->updateChildren(newParent, theChildren,true);
  if(_limitEmissions!=0) return true;
  // perform the shower of the final-state particle
  timeLikeShower(otherChild,type);
  // return the emitted
  return true;
}

void Evolver::showerDecay(ShowerTreePtr decay) {
  _decayme = HwDecayerBasePtr();
  _hardme  = HwMEBasePtr();
  // find the decayer
  // try the normal way if possible
  tDMPtr dm = decay->incomingLines().begin()->first->original()   ->decayMode();
  if(!dm) dm = decay->incomingLines().begin()->first->copy()      ->decayMode();
  if(!dm) dm = decay->incomingLines().begin()->first->progenitor()->decayMode();
  // otherwise make a string and look it up
  if(!dm) {
    string tag = decay->incomingLines().begin()->first->original()->dataPtr()->name() 
      + "->";
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  it=decay->outgoingLines().begin();it!=decay->outgoingLines().end();++it) {
      if(it!=decay->outgoingLines().begin()) tag += ",";
      tag += it->first->original()->dataPtr()->name();
    }
    tag += ";";
    dm = generator()->findDecayMode(tag);
  }
  if(dm) _decayme = dynamic_ptr_cast<HwDecayerBasePtr>(dm->decayer());
  // set the ShowerTree to be showered
  currentTree(decay);
  decay->applyTransforms();
  // extract particles to be shower, set scales and 
  // perform hard matrix element correction
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(false);
  setupMaximumScales(currentTree(), particlesToShower);
  // compute the minimum mass of the final-state
  Energy minmass(ZERO);
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if(particlesToShower[ix]->progenitor()->isFinalState())
      minmass+=particlesToShower[ix]->progenitor()->mass();
  }
  // main showering loop
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      currentTree()->clear();
      setEvolutionPartners(false,ShowerInteraction::QCD);
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
	progenitor()->hasEmitted(startTimeLikeShower(ShowerInteraction::QCD));
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
							   ShowerInteraction::QCD)); 
      }
    }
  }
  while(!showerModel()->kinematicsReconstructor()->reconstructDecayJets(decay)&&
	maximumTries()>++ntry);
  _decayme = HwDecayerBasePtr();
  if(maximumTries()==ntry) 
    throw Exception() << "Failed to generate the shower after "
		      << ntry << " attempts in Evolver::showerDecay()"
		      << Exception::eventerror;
  
  // tree has now showered
  _currenttree->hasShowered(true);
}

bool Evolver::spaceLikeDecayShower(tShowerParticlePtr particle,
				   Energy maxscale,
				   Energy minmass,ShowerInteraction::Type type) {
  Branching fb;
  while (true) {
    fb=_splittingGenerator->chooseDecayBranching(*particle,maxscale,minmass,
						 _initialenhance,type);
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
  particle->showerKinematics()->updateChildren(particle, theChildren,true);
  // In the case of splittings which involves coloured particles,
  // set properly the colour flow of the branching.
  // update the history if needed
  _currenttree->updateInitialStateShowerProduct(_progenitor,theChildren[0]);
  _currenttree->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
  // shower the first  particle
  spaceLikeDecayShower(theChildren[0],maxscale,minmass,type);
  // shower the second particle
  timeLikeShower(theChildren[1],type);
  // branching has happened
  return true;
}

vector<ShowerProgenitorPtr> Evolver::setupShower(bool hard) {
  // generate POWHEG hard emission if needed
  if(_hardEmissionMode==1) hardestEmission(hard);
  // set the initial colour partners
  setEvolutionPartners(hard,ShowerInteraction::QCD);
  // get the particles to be showered
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  // generate hard me if needed
  if(_hardEmissionMode==0) hardMatrixElementCorrection(hard);
  // get the particles to be showered
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=currentTree()->incomingLines().begin();
      cit!=currentTree()->incomingLines().end();++cit)
    particlesToShower.push_back((*cit).first);
  assert((particlesToShower.size()==1&&!hard)||
	 (particlesToShower.size()==2&&hard));
  // outgoing particles
  for(cjt=currentTree()->outgoingLines().begin();
      cjt!=currentTree()->outgoingLines().end();++cjt)
    particlesToShower.push_back(((*cjt).first));
  // remake the colour partners if needed
  if(_hardEmissionMode==0 && _currenttree->hardMatrixElementCorrection()) {
    setEvolutionPartners(hard,ShowerInteraction::QCD);
    _currenttree->resetShowerProducts();
  }
  return particlesToShower;
}

void Evolver::setEvolutionPartners(bool hard,ShowerInteraction::Type ) {
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  vector<ShowerParticlePtr> particles;
  // match the particles in the ShowerTree and hardTree
  if(hardTree() && !hardTree()->connect(currentTree()))
    throw Exception() << "Can't match trees in "
		      << "Evolver::setEvolutionPartners()"
		      << Exception::eventerror;
  // sort out the colour partners
  for(cit=currentTree()->incomingLines().begin();
      cit!=currentTree()->incomingLines().end();++cit)
    particles.push_back(cit->first->progenitor());
  assert((particles.size()==1&&!hard)||(particles.size()==2&&hard));
  // outgoing particles
  for(cjt=currentTree()->outgoingLines().begin();
      cjt!=currentTree()->outgoingLines().end();++cjt)
    particles.push_back(cjt->first->progenitor());
  if(hardTree()) {
    // find the partner
    for(unsigned int ix=0;ix<particles.size();++ix) {
      tHardBranchingPtr partner = 
	hardTree()->particles()[particles[ix]]->colourPartner();
      if(!partner) continue;
      for(map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator
	    it=hardTree()->particles().begin();
	  it!=hardTree()->particles().end();++it) {
	if(it->second==partner) particles[ix]->setPartner(it->first);
      }
      if(!particles[ix]->partner()) 
	throw Exception() << "Can't match partners in "
			  << "Evolver::setEvolutionPartners()"
			  << Exception::eventerror;
    }
  }
  // Set the initial evolution scales
  showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,!hard,ShowerInteraction::QCD,!_hardtree);
}

bool Evolver::startTimeLikeShower(ShowerInteraction::Type type) {
  if(hardTree()) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit=hardTree()->particles().end(),
      mit = hardTree()->particles().find(progenitor()->progenitor());
    if( mit != eit && !mit->second->children().empty() ) {
      return truncatedTimeLikeShower(progenitor()->progenitor(), mit->second ,type);
    }
  }
  return  hardOnly() ? false :
    timeLikeShower(progenitor()->progenitor() ,type) ;
}

bool Evolver::startSpaceLikeShower(PPtr parent, ShowerInteraction::Type type) {
  if(hardTree()) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit =hardTree()->particles().end(),
      mit = hardTree()->particles().find(progenitor()->progenitor());
    if( mit != eit && mit->second->parent() ) {
      return truncatedSpaceLikeShower( progenitor()->progenitor(),
				       parent, mit->second->parent(), type );
    } 
  }
  return  hardOnly() ? false :
    spaceLikeShower(progenitor()->progenitor(),parent,type);
}

bool Evolver::startSpaceLikeDecayShower(Energy maxscale,Energy minimumMass,
					ShowerInteraction::Type type) {
  if(hardTree()) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit =hardTree()->particles().end(),
      mit = hardTree()->particles().find(progenitor()->progenitor());
    if( mit != eit && mit->second->parent() ) {
      HardBranchingPtr branch=mit->second;
      while(branch->parent()) branch=branch->parent();
      return truncatedSpaceLikeDecayShower(progenitor()->progenitor(),maxscale,
					   minimumMass, branch ,type);
    }
  }
  return  hardOnly() ? false :
    spaceLikeDecayShower(progenitor()->progenitor(),maxscale,minimumMass,type);
}

bool Evolver::timeLikeVetoed(const Branching & fb,
			     ShowerParticlePtr particle) {
  // check whether emission was harder than largest pt of hard subprocess
  if ( hardVetoFS() && fb.kinematics->pT() > _progenitor->maxHardPt() ) 
    return true;
  // soft matrix element correction veto
  if( softMEC()) {
    if(_hardme && _hardme->hasMECorrection()) {
      if(_hardme->softMatrixElementVeto(_progenitor,particle,fb))
	return true;
    }
    else if(_decayme && _decayme->hasMECorrection()) {
      if(_decayme->softMatrixElementVeto(_progenitor,particle,fb))
	return true;
    }
  }
  // veto on maximum pt
  if(fb.kinematics->pT()>_progenitor->maximumpT()) return true;
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
  if( softMEC() && _hardme && _hardme->hasMECorrection() ) {
    if(_hardme->softMatrixElementVeto(_progenitor,particle,bb))
      return true;
  }
  // the more general vetos

  // check vs max pt for the shower
  if(bb.kinematics->pT()>_progenitor->maximumpT()) return true;

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

bool Evolver::spaceLikeDecayVetoed( const Branching & fb,
				    ShowerParticlePtr particle ) {
  // apply the soft correction
  if( softMEC() && _decayme && _decayme->hasMECorrection() ) {
    if(_decayme->softMatrixElementVeto(_progenitor,particle,fb))
      return true;
  }
  // veto on hardest pt in the shower
  if(fb.kinematics->pT()> _progenitor->maximumpT()) return true;
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
	if(test) throw VetoShower();
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

void Evolver::hardestEmission(bool hard) {
  if( ( _hardme &&  _hardme->hasPOWHEGCorrection()) ||
      (_decayme && _decayme->hasPOWHEGCorrection())) {
    if(_hardme)
      _hardtree =  _hardme->generateHardest( currentTree() );
    else
      _hardtree = _decayme->generateHardest( currentTree() );
    if(!_hardtree) return;
    // join up the two tree
    connectTrees(currentTree(),_hardtree,hard);
  }
  else {
    _hardtree = ShowerHandler::currentHandler()->generateCKKW(currentTree());
  }
}

bool Evolver::truncatedTimeLikeShower(tShowerParticlePtr particle,
				      HardBranchingPtr branch,
				      ShowerInteraction::Type type) {
  Branching fb;
  unsigned int iout=0;
  tcPDPtr pdata[2];
  while (true) {
    // no truncated shower break
    if(!isTruncatedShowerON()||hardOnly()) break;
    // generate emission
    fb=splittingGenerator()->chooseForwardBranching(*particle,1.,type);
    // no emission break
    if(!fb.kinematics) break;
    // check haven't evolved too far
    if(fb.kinematics->scale() < branch->scale()) {
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
    iout=0;
    if(pdata[0]->id()!=pdata[1]->id()) {
      if(pdata[0]->id()==particle->id())       iout=1;
      else if (pdata[1]->id()==particle->id()) iout=2;
    }
    else if(pdata[0]->id()==particle->id()) {
      if(fb.kinematics->z()>0.5) iout=1;
      else                       iout=2;
    }
    // apply the vetos for the truncated shower
    // no flavour changing branchings
    if(iout==0) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    double zsplit = iout==1 ? fb.kinematics->z() : 1-fb.kinematics->z();
    // only if same interaction for forced branching 
    // and evolution
    if(type==branch->sudakov()->interactionType()) {
      if(zsplit < 0.5 || // hardest line veto
	 fb.kinematics->scale()*zsplit < branch->scale() ) { // angular ordering veto
	particle->setEvolutionScale(fb.kinematics->scale());
	continue;
      }
    }
    // pt veto
    if(fb.kinematics->pT() > progenitor()->maximumpT()) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    // should do base class vetos as well
    if(timeLikeVetoed(fb,particle)) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    break;
  }
  // if no branching set decay matrix and return
  if(!fb.kinematics) {
    // construct the kinematics for the hard emission
    ShoKinPtr showerKin=
      branch->sudakov()->createFinalStateBranching(branch->scale(),
						   branch->children()[0]->z(),
						   branch->phi(),
						   branch->children()[0]->pT());
    particle->setEvolutionScale(branch->scale() );
    showerKin->initialize( *particle,PPtr() );
    IdList idlist(3);
    idlist[0] = particle->id();
    idlist[1] = branch->children()[0]->branchingParticle()->id();
    idlist[2] = branch->children()[1]->branchingParticle()->id();
    fb = Branching( showerKin, idlist, branch->sudakov() );
    // Assign the shower kinematics to the emitting particle.
    particle->setShowerKinematics( fb.kinematics );
    // Assign the splitting function to the emitting particle. 
    // For the time being we are considering only 1->2 branching
    // Create the ShowerParticle objects for the two children of
    // the emitting particle; set the parent/child relationship
    // if same as definition create particles, otherwise create cc
    ShowerParticleVector theChildren;
    theChildren.push_back(new_ptr(ShowerParticle(branch->children()[0]->
						 branchingParticle()->dataPtr(),true)));
    theChildren.push_back(new_ptr(ShowerParticle(branch->children()[1]->
						 branchingParticle()->dataPtr(),true)));
    particle->showerKinematics()->
      updateChildren(particle, theChildren,type==branch->sudakov()->interactionType());
    // update the history if needed
    if(particle==currentTree()->getFinalStateShowerProduct(progenitor()))
      currentTree()->updateFinalStateShowerProduct(progenitor(),
						   particle,theChildren);
    currentTree()->addFinalStateBranching(particle,theChildren);
    // shower the first  particle
    if( branch->children()[0]->children().empty() ) {
      if( ! hardOnly() )
	timeLikeShower(theChildren[0],type);
    }
    else {
      truncatedTimeLikeShower( theChildren[0],branch->children()[0],type);
    } 
    // shower the second particle
    if( branch->children()[1]->children().empty() ) {
      if( ! hardOnly() )
	timeLikeShower( theChildren[1] , type);
    }
    else {
      truncatedTimeLikeShower( theChildren[1],branch->children()[1] ,type);
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
  theChildren.push_back( new_ptr( ShowerParticle( pdata[0], true ) ) ); 
  theChildren.push_back( new_ptr( ShowerParticle( pdata[1], true ) ) );
  particle->showerKinematics()->
    updateChildren( particle, theChildren , true);
  // update the history if needed
  if( particle == currentTree()->getFinalStateShowerProduct( progenitor() ) )
    currentTree()->updateFinalStateShowerProduct( progenitor(),
						  particle, theChildren );
  currentTree()->addFinalStateBranching( particle, theChildren );
  // shower the first  particle
  if( iout == 1 ) truncatedTimeLikeShower( theChildren[0], branch , type );
  else            timeLikeShower( theChildren[0]  , type);
  // shower the second particle
  if( iout == 2 ) truncatedTimeLikeShower( theChildren[1], branch , type );
  else            timeLikeShower( theChildren[1]  , type);
  // branching has happened
  return true;
}

bool Evolver::truncatedSpaceLikeShower(tShowerParticlePtr particle, PPtr beam,
				       HardBranchingPtr branch,
				       ShowerInteraction::Type type) {
  tcPDFPtr pdf;
  if(ShowerHandler::currentHandler()->firstPDF().particle()  == beamParticle())
    pdf = ShowerHandler::currentHandler()->firstPDF().pdf();
  if(ShowerHandler::currentHandler()->secondPDF().particle() == beamParticle())
    pdf = ShowerHandler::currentHandler()->secondPDF().pdf();
  Energy freeze = ShowerHandler::currentHandler()->pdfFreezingScale();
  Branching bb;
  // generate branching
  tcPDPtr part[2];
  while (true) {
    if( !isTruncatedShowerON() || hardOnly() ) break;
    bb = splittingGenerator()->chooseBackwardBranching( *particle, 
							beam, 1., beamParticle(), 
							type , pdf,freeze);
    if( !bb.kinematics || bb.kinematics->scale() < branch->scale() ) {
      bb = Branching();
      break;
    }
    // particles as in Sudakov form factor
    part[0] = getParticleData( bb.ids[0] );
    part[1] = getParticleData( bb.ids[2] );
    
    //is emitter anti-particle
    if( particle->id() != bb.ids[1]) {
      if( part[0]->CC() ) part[0] = part[0]->CC();
      if( part[1]->CC() ) part[1] = part[1]->CC();
    }
    double zsplit = bb.kinematics->z();
    // apply the vetos for the truncated shower
    // if doesn't carry most of momentum
    if(type==branch->sudakov()->interactionType() &&
       zsplit < 0.5) {
      particle->setEvolutionScale(bb.kinematics->scale() );
      continue;
    }
    // others
    if( part[0]->id() != particle->id() || // if particle changes type
	bb.kinematics->pT() > progenitor()->maximumpT() ||   // pt veto
	bb.kinematics->scale() < branch->scale()) { // angular ordering veto
      particle->setEvolutionScale(bb.kinematics->scale() );
      continue;
    }
    // and those from the base class
    if(spaceLikeVetoed(bb,particle)) {
      particle->setEvolutionScale(bb.kinematics->scale() );
      continue;
    }
    break;
  }
  if( !bb.kinematics ) {
    //do the hard emission
    double z(0.);
    HardBranchingPtr timelike;
    for( unsigned int ix = 0; ix < branch->children().size(); ++ix ) {
      if( branch->children()[ix]->status() ==HardBranching::Outgoing) {
	timelike = branch->children()[ix];
      }
      if( branch->children()[ix]->status() ==HardBranching::Incoming )
	z = branch->children()[ix]->z();
    }
    ShoKinPtr kinematics =
      branch->sudakov()->createInitialStateBranching( branch->scale(), z, branch->phi(),
						      branch->children()[0]->pT() );
    kinematics->initialize( *particle, beam );
    // assign the splitting function and shower kinematics
    particle->setShowerKinematics( kinematics );
    // For the time being we are considering only 1->2 branching
    // Now create the actual particles, make the otherChild a final state
    // particle, while the newParent is not
    ShowerParticlePtr newParent = 
      new_ptr( ShowerParticle( branch->branchingParticle()->dataPtr(), false ) );
    ShowerParticlePtr otherChild = 
      new_ptr( ShowerParticle( timelike->branchingParticle()->dataPtr(),
			       true, true ) );
    ShowerParticleVector theChildren;
    theChildren.push_back( particle ); 
    theChildren.push_back( otherChild );
    particle->showerKinematics()->
      updateParent( newParent, theChildren, type==branch->sudakov()->interactionType() );
    // update the history if needed
    currentTree()->updateInitialStateShowerProduct( progenitor(), newParent );
    currentTree()->addInitialStateBranching( particle, newParent, otherChild );
    // for the reconstruction of kinematics, parent/child
    // relationships are according to the branching process:
    // now continue the shower
    bool emitted=false;
    if(!hardOnly()) {
      if( branch->parent() ) {
	emitted = truncatedSpaceLikeShower( newParent, beam, branch->parent() , type);
      }
      else {
	emitted = spaceLikeShower( newParent, beam , type);
      }
    }
    if( !emitted ) {
      if( intrinsicpT().find( progenitor() ) == intrinsicpT().end() ) {
	kinematics->updateLast( newParent, ZERO, ZERO );
      }
      else {
	pair<Energy,double> kt = intrinsicpT()[progenitor()];
	kinematics->updateLast( newParent,
				kt.first*cos( kt.second ),
				kt.first*sin( kt.second ) );
      }
    }
    particle->showerKinematics()->
      updateChildren( newParent, theChildren,
		      type==branch->sudakov()->interactionType() );
    if(hardOnly()) return true;
    // perform the shower of the final-state particle
    if( timelike->children().empty() ) {
      timeLikeShower( otherChild , type);
    }
    else {
      truncatedTimeLikeShower( otherChild, timelike , type);
    }
    // return the emitted
    return true;
  }
  // assign the splitting function and shower kinematics
  particle->setShowerKinematics( bb.kinematics );
  // For the time being we are considering only 1->2 branching
  // Now create the actual particles, make the otherChild a final state
  // particle, while the newParent is not
  ShowerParticlePtr newParent = new_ptr( ShowerParticle( part[0], false ) );
  ShowerParticlePtr otherChild = new_ptr( ShowerParticle( part[1], true, true ) );
  ShowerParticleVector theChildren; 
  theChildren.push_back( particle ); 
  theChildren.push_back( otherChild );
  particle->showerKinematics()->updateParent( newParent, theChildren , true);
  // update the history if needed
  currentTree()->updateInitialStateShowerProduct( progenitor(), newParent );
  currentTree()->addInitialStateBranching( particle, newParent, otherChild );
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // now continue the shower
  bool emitted = truncatedSpaceLikeShower( newParent, beam, branch,type);
  // now reconstruct the momentum
  if( !emitted ) {
    if( intrinsicpT().find( progenitor() ) == intrinsicpT().end() ) {
      bb.kinematics->updateLast( newParent, ZERO, ZERO );
    }
    else {
      pair<Energy,double> kt = intrinsicpT()[ progenitor() ];
      bb.kinematics->updateLast( newParent,
				 kt.first*cos( kt.second ),
				 kt.first*sin( kt.second ) );
    }
  }
  particle->showerKinematics()->updateChildren( newParent, theChildren , true);
  // perform the shower of the final-state particle
  timeLikeShower( otherChild , type);
  // return the emitted
  return true;
}

bool Evolver::
truncatedSpaceLikeDecayShower(tShowerParticlePtr particle, Energy maxscale,
			      Energy minmass, HardBranchingPtr branch,
			      ShowerInteraction::Type type) {
  Branching fb;
  unsigned int iout=0;
  tcPDPtr pdata[2];
  while (true) {
    // no truncated shower break
    if(!isTruncatedShowerON()||hardOnly()) break;
    fb=splittingGenerator()->chooseDecayBranching(*particle,maxscale,minmass,1.,type);
    // return if no radiation
    if(!fb.kinematics) break;
    // check haven't evolved too far
    if(fb.kinematics->scale() < branch->scale()) {
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
    iout=0;
    if(pdata[0]->id()!=pdata[1]->id()) {
      if(pdata[0]->id()==particle->id())       iout=1;
      else if (pdata[1]->id()==particle->id()) iout=2;
    }
    else if(pdata[0]->id()==particle->id()) {
      if(fb.kinematics->z()>0.5) iout=1;
      else                       iout=2;
    }
    // apply the vetos for the truncated shower
    // no flavour changing branchings
    if(iout==0) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    double zsplit = iout==1 ? fb.kinematics->z() : 1-fb.kinematics->z();
    if(type==branch->sudakov()->interactionType()) {
      if(zsplit < 0.5 || // hardest line veto
	 fb.kinematics->scale()*zsplit < branch->scale() ) { // angular ordering veto
	particle->setEvolutionScale(fb.kinematics->scale());
	continue;
      }
    }
    // pt veto
    if(fb.kinematics->pT() > progenitor()->maximumpT()) {
      particle->setEvolutionScale(fb.kinematics->scale());
      continue;
    }
    // should do base class vetos as well
    // if not vetoed break
    if(!spaceLikeDecayVetoed(fb,particle)) break;
    // otherwise reset scale and continue
    particle->setEvolutionScale(fb.kinematics->scale());
  }
  // if no branching set decay matrix and return
  if(!fb.kinematics) {
    // construct the kinematics for the hard emission
    ShoKinPtr showerKin=
      branch->sudakov()->createDecayBranching(branch->scale(),
					      branch->children()[0]->z(),
					      branch->phi(),
					      branch->children()[0]->pT());
    particle->setEvolutionScale(branch->scale() );
    showerKin->initialize( *particle,PPtr() );
    IdList idlist(3);
    idlist[0] = particle->id();
    idlist[1] = branch->children()[0]->branchingParticle()->id();
    idlist[2] = branch->children()[1]->branchingParticle()->id();
    fb = Branching( showerKin, idlist, branch->sudakov() );
    // Assign the shower kinematics to the emitting particle.
    particle->setShowerKinematics( fb.kinematics );
    // Assign the splitting function to the emitting particle. 
    // For the time being we are considering only 1->2 branching
    // Create the ShowerParticle objects for the two children of
    // the emitting particle; set the parent/child relationship
    // if same as definition create particles, otherwise create cc
    ShowerParticleVector theChildren;
    theChildren.push_back(new_ptr(ShowerParticle(branch->children()[0]->
						 branchingParticle()->dataPtr(),true)));
    theChildren.push_back(new_ptr(ShowerParticle(branch->children()[1]->
						 branchingParticle()->dataPtr(),true)));
    particle->showerKinematics()->
      updateChildren(particle, theChildren,
		     type==branch->sudakov()->interactionType());
    if(theChildren[0]->id()==particle->id()) {
      // update the history if needed
      currentTree()->updateInitialStateShowerProduct(progenitor(),theChildren[0]);
      currentTree()->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
      // shower the space-like particle
      if( branch->children()[0]->children().empty() ) {
	if( ! hardOnly() ) spaceLikeDecayShower(theChildren[0],maxscale,minmass,type);
      }
      else {
	truncatedSpaceLikeDecayShower( theChildren[0],maxscale,minmass,
				       branch->children()[0],type);
      }
      // shower the second particle
      if( branch->children()[1]->children().empty() ) {
	if( ! hardOnly() ) timeLikeShower( theChildren[1] , type);
      }
      else {
	truncatedTimeLikeShower( theChildren[1],branch->children()[1] ,type);
      }
    }
    else {
      // update the history if needed
      currentTree()->updateInitialStateShowerProduct(progenitor(),theChildren[1]);
      currentTree()->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
      // shower the space-like particle
      if( branch->children()[1]->children().empty() ) {
	if( ! hardOnly() ) spaceLikeDecayShower(theChildren[1],maxscale,minmass,type);
      }
      else {
	truncatedSpaceLikeDecayShower( theChildren[1],maxscale,minmass,
				       branch->children()[1],type);
      }
      // shower the second particle
      if( branch->children()[0]->children().empty() ) {
	if( ! hardOnly() ) timeLikeShower( theChildren[0] , type);
      }
      else {
	truncatedTimeLikeShower( theChildren[0],branch->children()[0] ,type);
      }
    }
    return true;
  }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->setShowerKinematics(fb.kinematics);
  // For the time being we are considering only 1->2 branching
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  // if same as definition create particles, otherwise create cc
  ShowerParticleVector theChildren; 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[0],true))); 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[1],true)));
  particle->showerKinematics()->updateChildren(particle, theChildren,true);
  // In the case of splittings which involves coloured particles,
  // set properly the colour flow of the branching.
  // update the history if needed
  currentTree()->updateInitialStateShowerProduct(progenitor(),theChildren[0]);
  currentTree()->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
  // shower the first  particle
  truncatedSpaceLikeDecayShower(theChildren[0],maxscale,minmass,branch,type);
  // shower the second particle
  timeLikeShower(theChildren[1],type);
  // branching has happened
  return true;
}

void Evolver::connectTrees(ShowerTreePtr showerTree, HardTreePtr hardTree, bool hard )const {
  ShowerParticleVector particles;
  // find the Sudakovs
  for(set<HardBranchingPtr>::iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    // Sudakovs for ISR
    if((**cit).parent()&&(**cit).status()==HardBranching::Incoming) {
      IdList br(3);
      br[0] = (**cit).parent()->branchingParticle()->id();
      br[1] = (**cit).          branchingParticle()->id();
      br[2] = (**cit).parent()->children()[0]==*cit ?
	(**cit).parent()->children()[1]->branchingParticle()->id() :
	(**cit).parent()->children()[0]->branchingParticle()->id();
      BranchingList branchings = splittingGenerator()->initialStateBranchings();
      if(br[1]<0&&br[0]==br[1]) {
	br[0] = abs(br[0]);
	br[1] = abs(br[1]);
      }
      else if(br[1]<0) {
	br[1] = -br[1];
	br[2] = -br[2];
      }
      long index = abs(br[1]);
      SudakovPtr sudakov;
      for(BranchingList::const_iterator cjt = branchings.lower_bound(index); 
	  cjt != branchings.upper_bound(index); ++cjt ) {
	IdList ids = cjt->second.second;
	if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
	  sudakov=cjt->second.first;
	  break;
	}
      }
      if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				     << "Evolver::connectTrees() for ISR" 
				     << Exception::runerror;
      (**cit).parent()->sudakov(sudakov);
    }
    // Sudakovs for FSR
    else if(!(**cit).children().empty()) {
      IdList br(3);
      br[0] = (**cit)               .branchingParticle()->id();
      br[1] = (**cit).children()[0]->branchingParticle()->id();
      br[2] = (**cit).children()[1]->branchingParticle()->id();
      BranchingList branchings = splittingGenerator()->finalStateBranchings();
      if(br[0]<0) {
	br[0] = abs(br[0]);
	br[1] = abs(br[1]);
	br[2] = abs(br[2]);
      }
      long index = br[0];
      SudakovPtr sudakov;
      for(BranchingList::const_iterator cjt = branchings.lower_bound(index); 
	  cjt != branchings.upper_bound(index); ++cjt ) {
	IdList ids = cjt->second.second;
	if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
	  sudakov=cjt->second.first;
	  break;
	}
      }
      if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				     << "Evolver::generateHardest()" 
				     << Exception::runerror;
      (**cit).sudakov(sudakov);
    }
  }
  // calculate the evolution scale
  for(set<HardBranchingPtr>::iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    particles.push_back((*cit)->branchingParticle());
  }
  showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,!hard,hardTree->interaction(),true);
  // inverse reconstruction
  if(hard)
    showerModel()->kinematicsReconstructor()->
      deconstructHardJets(hardTree,ShowerHandler::currentHandler()->evolver(),
			  hardTree->interaction());
  else
    showerModel()->kinematicsReconstructor()->
      deconstructDecayJets(hardTree,ShowerHandler::currentHandler()->evolver(),
			   hardTree->interaction());
  // now reset the momenta of the showering particles
  vector<ShowerProgenitorPtr> particlesToShower;
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator
	cit=showerTree->incomingLines().begin();
      cit!=showerTree->incomingLines().end();++cit )
    particlesToShower.push_back(cit->first);
  // extract the showering particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator
	  cit=showerTree->outgoingLines().begin();
	cit!=showerTree->outgoingLines().end();++cit )
    particlesToShower.push_back(cit->first);
  // match them
  vector<bool> matched(particlesToShower.size(),false);
  for(set<HardBranchingPtr>::const_iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    Energy2 dmin( 1e30*GeV2 );
    int iloc(-1);
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      if(matched[ix]) continue;
      if( (**cit).branchingParticle()->id() !=  particlesToShower[ix]->progenitor()->id() ) continue;
      if( (**cit).branchingParticle()->isFinalState() !=
	  particlesToShower[ix]->progenitor()->isFinalState() ) continue;
      Energy2 dtest = 
	sqr( particlesToShower[ix]->progenitor()->momentum().x() - (**cit).showerMomentum().x() ) +
	sqr( particlesToShower[ix]->progenitor()->momentum().y() - (**cit).showerMomentum().y() ) +
	sqr( particlesToShower[ix]->progenitor()->momentum().z() - (**cit).showerMomentum().z() ) +
	sqr( particlesToShower[ix]->progenitor()->momentum().t() - (**cit).showerMomentum().t() );
      // add mass difference for identical particles (e.g. Z0 Z0 production)
      dtest += 1e10*sqr(particlesToShower[ix]->progenitor()->momentum().m()-
			(**cit).showerMomentum().m());
      if( dtest < dmin ) {
	iloc = ix;
	dmin = dtest;
      }
    }
    if(iloc<0) throw Exception() << "Failed to match shower and hard trees in Evolver::hardestEmission"
				 << Exception::eventerror;
    particlesToShower[iloc]->progenitor()->set5Momentum((**cit).showerMomentum());
    matched[iloc] = true;
  }
  // correction boosts for daughter trees
  for(map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator 
	tit  = showerTree->treelinks().begin();
      tit != showerTree->treelinks().end();++tit) {
    ShowerTreePtr decayTree = tit->first;
    map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
      cit = decayTree->incomingLines().begin();
    // reset the momentum of the decay particle
    Lorentz5Momentum oldMomentum = cit->first->progenitor()->momentum();
    Lorentz5Momentum newMomentum = tit->second.second->momentum();
    LorentzRotation boost( oldMomentum.findBoostToCM(),oldMomentum.e()/oldMomentum.mass());
    boost.boost          (-newMomentum.findBoostToCM(),newMomentum.e()/newMomentum.mass());
    decayTree->transform(boost,true);
  }
}
