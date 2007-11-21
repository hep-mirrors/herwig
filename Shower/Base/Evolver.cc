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

#include "Herwig++/Shower/CKKW/Clustering/CascadeReconstructor.h"
#include "Herwig++/Shower/CKKW/Reweighting/DefaultReweighter.h"
#include "Herwig++/Shower/CKKW/Reweighting/DefaultCKKWVeto.h"

using namespace Herwig;

void Evolver::persistentOutput(PersistentOStream & os) const {
  os << _model << _splittingGenerator << _maxtry 
     << _meCorrMode << _hardVetoMode
     << ounit(_iptrms,GeV) << _beta << ounit(_gamma,GeV) << _vetoes
     << _reconstructor << _reweighter << _ckkwVeto << _useCKKW;
}

void Evolver::persistentInput(PersistentIStream & is, int) {
  is >> _model >> _splittingGenerator >> _maxtry 
     >> _meCorrMode >> _hardVetoMode
     >> iunit(_iptrms,GeV) >> _beta >> iunit(_gamma,GeV) >> _vetoes
     >> _reconstructor >> _reweighter >> _ckkwVeto >> _useCKKW;
}

void Evolver::doinitrun() {
  Interfaced::doinitrun();
  for(unsigned int ix=0;ix<showerModel()->meCorrections().size();++ix) {
    showerModel()->meCorrections()[ix]->evolver(this);
  }
#ifdef HERWIG_CHECK_VETOES
  _vetoed_points.timelike.open("vetoed_timelike.dat");
  _vetoed_points.spacelike.open("vetoed_spacelike.dat");
  _vetoed_points.spacelike_decay.open("vetoed_spacelike_decay.dat");
#endif
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

  static RefVector<Evolver,ShowerVeto> ifaceVetoes
    ("Vetoes",
     "The vetoes to be checked during showering",
     &Evolver::_vetoes, -1,
     false,false,true,true,false);


 }


void Evolver::useCKKW (CascadeReconstructorPtr cr, ReweighterPtr rew) {

  DefaultReweighterPtr drew = dynamic_ptr_cast<DefaultReweighterPtr>(rew);
  if (!drew) Throw<InitException>()
    << "Shower : Evolver::useCKKW : DefaultReweighter needed by Evolver, found Reweighter.";
  _reconstructor = cr;
  _reweighter = drew;

  // alpha_s

  // alpha_s is set as a reference in Reweighter

  // now register the splitting functions with
  // the reweighter

  if (!_splittingGenerator->isISRadiationON(ShowerIndex::QCD) && 
      !_splittingGenerator->isFSRadiationON(ShowerIndex::QCD))
    Throw<InitException>() << "Shower : Evolver::useCKKW : Attempt to use CKKW with QCD radiation switched off.";

  if(_splittingGenerator->isISRadiationON(ShowerIndex::QCD)) {
    BranchingList bb = _splittingGenerator->initialStateBranchings();
    for(BranchingList::iterator b = bb.begin(); b != bb.end(); ++b) {
      if (b->second.first->interactionType() == ShowerIndex::QCD)
	_reweighter->insertSplitting(b->second.second,
				    b->second.first->splittingFn(),
				    true);
    }
  }

  if(_splittingGenerator->isFSRadiationON(ShowerIndex::QCD)) {
    BranchingList fb = _splittingGenerator->finalStateBranchings();
    for(BranchingList::iterator b = fb.begin(); b != fb.end(); ++b) {
      if (b->second.first->interactionType() == ShowerIndex::QCD)
	_reweighter->insertSplitting(b->second.second,
				    b->second.first->splittingFn(),
				    false);
    }
  }

  // setup the reweighter

  _reweighter->setup(_reconstructor);

  // and insert an appropriate veto

  if (!dynamic_ptr_cast<DefaultJetMeasurePtr>(_reweighter->resolution()))
    Throw<InitException>() << "Shower : Evolver::useCKKW : DefaultJetMeasure needed by Evolver, found JetMeasure.";

  _ckkwVeto =
    new_ptr(DefaultCKKWVeto(dynamic_ptr_cast<DefaultJetMeasurePtr>(_reweighter->resolution())));

  _ckkwVeto->enable();

  addVeto(_ckkwVeto);

  // indicate that we are doing CKKW

  _useCKKW = true;

}

void Evolver::initCKKWShower (unsigned int currentMult, unsigned int maxMult) {

  _ckkwVeto->eventGenerator(generator());

  // disable the veto for maximum multiplicity,
  // if wanted.

#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "== Evolver::initCKKWShower (" << currentMult << ", " << maxMult << ")" << endl;
#endif

  if(!_reweighter->vetoHighest() && currentMult == maxMult) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "max multiplicity and no vetoes on max multiplicity" << endl;
#endif
    _ckkwVeto->disable();
  }
  else {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "veto will be applied." << endl;
#endif
    _ckkwVeto->enable();
  }

}


void Evolver::generateIntrinsicpT(vector<ShowerProgenitorPtr> particlesToShower) {
  _intrinsic.clear();
  if ( !ipTon() || !isISRadiationON() ) return;
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    // only consider initial-state particles
    if(particlesToShower[ix]->progenitor()->isFinalState()) continue;
    if(!particlesToShower[ix]->progenitor()->dataPtr()->coloured()) continue;
    Energy ipt;
    if(UseRandom::rnd() > _beta) {
      ipt=_iptrms*sqrt(-log(UseRandom::rnd()));
    }
    else {
      ipt=_gamma*tan(Constants::pi*UseRandom::rnd()/2.);
    }
    pair<Energy,double> pt = make_pair(ipt,UseRandom::rnd()*Constants::twopi);
    _intrinsic[particlesToShower[ix]] = pt;
  }
}

void Evolver::setupMaximumScales(ShowerTreePtr hard, 
				 vector<ShowerProgenitorPtr> p) {

  // find out if hard partonic subprocess.
  bool isPartonic(false); 
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  cit = _currenttree->incomingLines().begin();
  Lorentz5Momentum pcm;
  for(; cit!=hard->incomingLines().end(); ++cit) {
    pcm += cit->first->progenitor()->momentum(); 
    if (isPartonic || cit->first->progenitor()->coloured()) {
      isPartonic = true; 
    }
  }

  // find maximum pt from hard process, the maximum pt from all
  // outgoing coloured lines (this is simpler and more general than
  // 2stu/(s^2+t^2+u^2)).
  Energy ptmax = -1.0*GeV, pt = 0.0*GeV;

  if (isPartonic) {
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
    cjt = hard->outgoingLines().begin();
    for(; cjt!=hard->outgoingLines().end(); ++cjt) {
      if (cjt->first->progenitor()->coloured()) {
	pt = cjt->first->progenitor()->momentum().perp();
	if (pt > ptmax) {
	  ptmax = pt;      
	}
      }
    }
    // if there are no coloured FS particles, use shat as maximum
    if (ptmax < 0.0*GeV) {
      ptmax = pcm.m(); 
    }
  } else {
    // if no coloured IS use shat as well
    ptmax = pcm.m();
  }
  
  // set maxHardPt for all progenitors.  For partonic processes this
  // is now the max pt in the FS, for non-partonic processes or
  // processes with no coloured FS the invariant mass of the IS
  vector<ShowerProgenitorPtr>::const_iterator ckt = p.begin();
  for (; ckt != p.end(); ckt++) {
    (*ckt)->maxHardPt(ptmax);
  }
}

void Evolver::showerHardProcess(ShowerTreePtr hard) {
  // set the current tree
  currentTree(hard);
  vector<ShowerProgenitorPtr> particlesToShower=setupShower(true);
  // setup the maximum scales for the shower, given by the hard process
  if (hardVetoOn()) {
    setupMaximumScales(_currenttree, particlesToShower);
  }

  // generate the intrinsic p_T once and for all
  generateIntrinsicpT(particlesToShower);
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      currentTree()->clear();
      setColourPartners(true);
    }
    // initial-state radiation
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
	_progenitor=particlesToShower[ix];

	tPPtr beamparticle; 
	if ( particlesToShower[ix]->original()->parents().empty() ) {
	  beamparticle=particlesToShower[ix]->original();
	} 
	else {
	  beamparticle=particlesToShower[ix]->original()->parents()[0];
	}
	_progenitor->hasEmitted(spaceLikeShower(particlesToShower[ix]->progenitor(),
						beamparticle));
      }
    }
    // final-state radiation
    if(isFSRadiationON()) {
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	// only consider final-state particles
	if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
	// perform shower
	_progenitor=particlesToShower[ix];
	_progenitor->hasEmitted(timeLikeShower(particlesToShower[ix]->progenitor()));
      }
    }
  }
  while(!_model->kinematicsReconstructor()->reconstructHardJets(hard,_intrinsic)&&
	_maxtry>++ntry);
  if(_maxtry==ntry) throw ShowerHandler::ShowerTriesVeto(ntry);

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
  // if no suitable me correction
  if(!_currentme) return; 
  // now apply the hard correction
  if(hardMEC()) _currentme->applyHardMatrixElementCorrection(_currenttree);
}

bool Evolver::timeLikeShower(tShowerParticlePtr particle) {
  bool vetoed = true;
  Branching fb;
  while (vetoed) {
    vetoed = false; 
    fb=_splittingGenerator->chooseForwardBranching(*particle,_finalenhance);
    if(!fb.kinematics) break;
    // check whether emission was harder than largest pt of hard subprocess
    if (hardVetoFS() && fb.kinematics->pT() > _progenitor->maxHardPt()) {
      vetoed = true;
      particle->setEvolutionScale(ShowerIndex::QCD, fb.kinematics->scale());
      continue;
    }

    // apply vetos if needed
    if(_currentme && softMEC() && !_theUseCKKW) {
      vetoed=_currentme->softMatrixElementVeto(_progenitor,particle,fb);
      if (vetoed) continue;
    }
    if(fb.kinematics->pT()>_progenitor->maximumpT()) vetoed=true;

    if (fb.kinematics && !_vetoes.empty()) {
      for (vector<ShowerVetoPtr>::iterator v = _vetoes.begin();
	   v != _vetoes.end(); ++v) {
	if ((**v).vetoType() == ShowerVetoType::Emission) {
	  vetoed |= (**v).vetoTimeLike(_progenitor,particle,fb);
	}
	if ((**v).vetoType() == ShowerVetoType::Shower && (**v).vetoTimeLike(_progenitor,particle,fb)) {
	  throw VetoShower ();
	}
	if ((**v).vetoType() == ShowerVetoType::Event && (**v).vetoTimeLike(_progenitor,particle,fb))
	  throw Veto ();
      }
      if (vetoed) {
#ifdef HERWIG_CHECK_VETOES
	_vetoed_points.timelike << fb.kinematics->scale()/GeV << "\t"
				<< fb.kinematics->z() << endl;
#endif
	particle->setEvolutionScale(ShowerIndex::QCD, fb.kinematics->scale());
	continue;
      }
    }

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
  theChildren.push_back(new_ptr(ShowerParticle(pdata[0],true))); 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[1],true))); 
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

bool Evolver::spaceLikeShower(tShowerParticlePtr particle, PPtr beam) {
  bool vetoed(true);
  Branching bb;
  // generate branching
  while (vetoed) {
    vetoed=false;
    bb=_splittingGenerator->chooseBackwardBranching(*particle,beam,
						    _initialenhance,_beam);

    // check whether emission was harder than largest pt of hard subprocess
    if (hardVetoIS() && bb.kinematics 
	&& bb.kinematics->pT() > _progenitor->maxHardPt()) {
      vetoed = true;
      particle->setEvolutionScale(ShowerIndex::QCD, bb.kinematics->scale());
      continue;
    }

    // apply the soft correction
    if(bb.kinematics && _currentme && softMEC() && !_theUseCKKW) {
      vetoed=_currentme->softMatrixElementVeto(_progenitor,particle,bb);
      if (vetoed) continue;
    }
    if(bb.kinematics && bb.kinematics->pT()>_progenitor->maximumpT()) vetoed=true;

    if (bb.kinematics && _vetoes.size()) {
      for (vector<ShowerVetoPtr>::iterator v = _vetoes.begin();
	   v != _vetoes.end(); ++v) {
	if ((**v).vetoType() == ShowerVetoType::Emission) {
	  vetoed |= (**v).vetoSpaceLike(_progenitor,particle,bb);
	}
	if ((**v).vetoType() == ShowerVetoType::Shower && (**v).vetoSpaceLike(_progenitor,particle,bb))
	  throw VetoShower ();
	if ((**v).vetoType() == ShowerVetoType::Event && (**v).vetoSpaceLike(_progenitor,particle,bb))
	  throw Veto ();
      }
      if (vetoed) {
#ifdef HERWIG_CHECK_VETOES
	_vetoed_points.spacelike << bb.kinematics->scale()/GeV << "\t"
				 << bb.kinematics->z() << endl;
#endif
	particle->setEvolutionScale(ShowerIndex::QCD, bb.kinematics->scale());
	continue;
      }
    }
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
  bool emitted=spaceLikeShower(newParent,beam);
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
  // main showering loop
  unsigned int ntry(0);
  do {
    // clear results of last attempt
    if(ntry!=0) {
      _currenttree->clear();
      setColourPartners(false);
    }
    // initial-state radiation
    if(isISRadiationON()) {
      // compute the minimum mass of the final-state
      Energy minmass(0.*MeV);
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
    if(isFSRadiationON()) {
      for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
	// only consider final-state particles
	if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
	// perform shower
	_progenitor=particlesToShower[ix];
	_progenitor->hasEmitted(timeLikeShower(particlesToShower[ix]->progenitor()));
      }
    }
  }
  while(!showerModel()->kinematicsReconstructor()->reconstructDecayJets(decay)&&
	maximumTries()>++ntry);
  if(maximumTries()==ntry) throw Exception() << "Failed to generate the shower after "
				      << ntry << " attempts in Evolver::showerDecay()"
				      << Exception::eventerror;
  _currenttree->hasShowered(true);
}

bool Evolver::spaceLikeDecayShower(tShowerParticlePtr particle,vector<Energy> maxscale,
				   Energy minmass) {
  bool vetoed = true;
  Branching fb;
  while (vetoed) {
    vetoed = false;
    fb=_splittingGenerator->chooseDecayBranching(*particle,maxscale,minmass,
						 _initialenhance);
    
    // SG: here could be the veto of too hard radiation.  I think the
    // initial conditions don't allow this anyway.

    // apply the soft correction
    if(fb.kinematics && _currentme && softMEC() && !_theUseCKKW) {
      vetoed=_currentme->softMatrixElementVeto(_progenitor,particle,fb);
      if (vetoed) continue;
    }
    if(fb.kinematics && fb.kinematics->pT()>_progenitor->maximumpT()) vetoed=true;

    if (fb.kinematics && _vetoes.size()) {
      for (vector<ShowerVetoPtr>::iterator v = _vetoes.begin();
	   v != _vetoes.end(); ++v) {
	if ((**v).vetoType() == ShowerVetoType::Emission) {
	  vetoed |= (**v).vetoSpaceLike(_progenitor,particle,fb);
	}
	if ((**v).vetoType() == ShowerVetoType::Shower && (**v).vetoSpaceLike(_progenitor,particle,fb))
	  throw VetoShower ();
	if ((**v).vetoType() == ShowerVetoType::Event && (**v).vetoSpaceLike(_progenitor,particle,fb))
	  throw Veto ();
      }
      if (vetoed) {
#ifdef HERWIG_CHECK_VETOES
	_vetoed_points.spacelike_decay << fb.kinematics->scale()/GeV << "\t"
				       << fb.kinematics->z() << endl;
#endif
	particle->setEvolutionScale(ShowerIndex::QCD, fb.kinematics->scale());
	continue;
      }
    }
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

  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;

  // If we use CKKW and there is information available from
  // the last parton shower history for each coloured line
  // in the current tree, we reset scales and apply vetoes
  _theUseCKKW = false;
  if (_useCKKW) {
    _theUseCKKW = true;
    for(cit=_currenttree->incomingLines().begin();
	cit!=_currenttree->incomingLines().end();++cit) {
      if ((*cit).first->original()->coloured())
	if (!_reconstructor->clusteringParticle((*cit).first->original())) {
	  _theUseCKKW = false;
	  break;
	}
    }
    if(_theUseCKKW)
      for(cjt=_currenttree->outgoingLines().begin();
	  cjt!=_currenttree->outgoingLines().end();++cjt) {
	if ((*cjt).first->original()->coloured())
	  if (!_reconstructor->clusteringParticle((*cjt).first->original())) {
	    _theUseCKKW = false;
	    break;
	  }
      }
    // if we failed here, disable the veto
    if(!_theUseCKKW) _ckkwVeto->disable();
  }

  // generate the hard matrix element correction if needed
  // don't do this if we are doing CKKW
  if (!_theUseCKKW)
    hardMatrixElementCorrection();
  // get the particles to be showered
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=_currenttree->incomingLines().begin();
      cit!=_currenttree->incomingLines().end();++cit)
    particlesToShower.push_back((*cit).first);
  assert((particlesToShower.size()==1&&!hard)
	 ||(particlesToShower.size()==2&&hard));
  // outgoing particles
  for(cjt=_currenttree->outgoingLines().begin();
      cjt!=_currenttree->outgoingLines().end();++cjt)
    particlesToShower.push_back(((*cjt).first));
  // remake the colour partners if needed
  if(_currenttree->hardMatrixElementCorrection() && !_theUseCKKW) {
    setColourPartners(hard);
    _currenttree->resetShowerProducts();
  }
  // if CKKW, reset the initial scales
  if (_theUseCKKW) {

    for (vector<ShowerProgenitorPtr>::iterator p = particlesToShower.begin();
	 p != particlesToShower.end(); ++p)
      // if not coloured, don't consider
      if ((**p).original()->coloured()) {
	tClusteringParticlePtr historyP = _reconstructor->clusteringParticle((**p).original());
	if (!historyP)
	  throw Exception() << "Shower : Evolver::setupShower (CKKW): No external leg for particle"
			    << " found in cascade history" << Exception::eventerror;
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "resetting hard scale for " << (**p).original()
			   << " PDGId = " << (**p).original()->id()
			   << " to " << sqrt(historyP->showerScale())/GeV << " GeV "
			   << " using clustering partilce " << endl;
	historyP->debugDump(generator()->log());
#endif
	(**p).progenitor()->setEvolutionScale(ShowerIndex::QCD,sqrt(historyP->showerScale()));
      }

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
  if(_splittingGenerator->isInteractionON(ShowerIndex::QCD))
    showerModel()->partnerFinder()->setQCDInitialEvolutionScales(particles,!hard);
  if(_splittingGenerator->isInteractionON(ShowerIndex::QED))
    showerModel()->partnerFinder()->setQEDInitialEvolutionScales(particles,!hard);
  if(_splittingGenerator->isInteractionON(ShowerIndex::EWK))
    showerModel()->partnerFinder()->setEWKInitialEvolutionScales(particles,!hard);
}
