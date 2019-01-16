// -*- C++ -*-
//
// QTildeShowerHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeShowerHandler class.
//

#include "QTildeShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include "Herwig/PDF/MPIPDF.h"
#include "Herwig/PDF/MinBiasPDF.h"
#include "Herwig/Shower/Core/Base/ShowerTree.h"
#include "Herwig/Shower/Core/Base/HardTree.h"
#include "Herwig/Shower/QTilde/Base/KinematicsReconstructor.h"
#include "Herwig/Shower/QTilde/Base/PartnerFinder.h"
#include "Herwig/PDF/HwRemDecayer.h"
#include "Herwig/Shower/Core/Base/ShowerVertex.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Base/SubtractedME.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "Herwig/Shower/RealEmissionProcess.h"

using namespace Herwig;

bool QTildeShowerHandler::_hardEmissionWarn = true;
bool QTildeShowerHandler::_missingTruncWarn = true;

QTildeShowerHandler::QTildeShowerHandler() :
  _maxtry(100), _meCorrMode(1), _reconOpt(0),
  _hardVetoReadOption(false),
  _iptrms(ZERO), _beta(0.), _gamma(ZERO), _iptmax(),
  _limitEmissions(0), _initialenhance(1.), _finalenhance(1.),
  _nReWeight(100), _reWeight(false),
  interaction_(ShowerInteraction::Both),
  _trunc_Mode(true), _hardEmission(1),
  _spinOpt(1), _softOpt(2), _hardPOWHEG(false), muPt(ZERO),
  _maxTryFSR(100000), _maxFailFSR(100), _fracFSR(0.001),
  _nFSR(0), _nFailedFSR(0)
{}

QTildeShowerHandler::~QTildeShowerHandler() {}

IBPtr QTildeShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr QTildeShowerHandler::fullclone() const {
  return new_ptr(*this);
}

void QTildeShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << _model << _splittingGenerator << _maxtry 
     << _meCorrMode << _hardVetoReadOption
     << _limitEmissions << _spinOpt << _softOpt << _hardPOWHEG
     << ounit(_iptrms,GeV) << _beta << ounit(_gamma,GeV) << ounit(_iptmax,GeV) 
     << _vetoes << _fullShowerVetoes << _nReWeight << _reWeight
     << _trunc_Mode << _hardEmission << _reconOpt 
     << ounit(muPt,GeV)
     << oenum(interaction_) << _maxTryFSR << _maxFailFSR << _fracFSR;
}

void QTildeShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> _model >> _splittingGenerator >> _maxtry 
     >> _meCorrMode >> _hardVetoReadOption
     >> _limitEmissions >> _spinOpt >> _softOpt >> _hardPOWHEG
     >> iunit(_iptrms,GeV) >> _beta >> iunit(_gamma,GeV) >> iunit(_iptmax,GeV)
     >> _vetoes >> _fullShowerVetoes >> _nReWeight >> _reWeight
     >> _trunc_Mode >> _hardEmission >> _reconOpt
     >> iunit(muPt,GeV)
     >> ienum(interaction_) >> _maxTryFSR >> _maxFailFSR >> _fracFSR;
}


// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<QTildeShowerHandler,ShowerHandler>
describeHerwigQTildeShowerHandler("Herwig::QTildeShowerHandler", "HwShower.so");

void QTildeShowerHandler::Init() {

  static ClassDocumentation<QTildeShowerHandler> documentation
    ("TheQTildeShowerHandler class is the main class"
     " for the angular-ordered parton shower",
     "The Shower evolution was performed using an algorithm described in "
     "\\cite{Marchesini:1983bm,Marchesini:1987cf,Gieseke:2003rz,Bahr:2008pv}.",
     "%\\cite{Marchesini:1983bm}\n"
     "\\bibitem{Marchesini:1983bm}\n"
     "  G.~Marchesini and B.~R.~Webber,\n"
     "  ``Simulation Of QCD Jets Including Soft Gluon Interference,''\n"
     "  Nucl.\\ Phys.\\  B {\\bf 238}, 1 (1984).\n"
     "  %%CITATION = NUPHA,B238,1;%%\n"
     "%\\cite{Marchesini:1987cf}\n"
     "\\bibitem{Marchesini:1987cf}\n"
     "  G.~Marchesini and B.~R.~Webber,\n"
     "   ``Monte Carlo Simulation of General Hard Processes with Coherent QCD\n"
     "  Radiation,''\n"
     "  Nucl.\\ Phys.\\  B {\\bf 310}, 461 (1988).\n"
     "  %%CITATION = NUPHA,B310,461;%%\n"
     "%\\cite{Gieseke:2003rz}\n"
     "\\bibitem{Gieseke:2003rz}\n"
     "  S.~Gieseke, P.~Stephens and B.~Webber,\n"
     "  ``New formalism for QCD parton showers,''\n"
     "  JHEP {\\bf 0312}, 045 (2003)\n"
     "  [arXiv:hep-ph/0310083].\n"
     "  %%CITATION = JHEPA,0312,045;%%\n"
     );

  static Reference<QTildeShowerHandler,SplittingGenerator> 
    interfaceSplitGen("SplittingGenerator", 
		      "A reference to the SplittingGenerator object", 
		      &Herwig::QTildeShowerHandler::_splittingGenerator,
		      false, false, true, false);

  static Reference<QTildeShowerHandler,ShowerModel> interfaceShowerModel
    ("ShowerModel",
     "The pointer to the object which defines the shower evolution model.",
     &QTildeShowerHandler::_model, false, false, true, false, false);

  static Parameter<QTildeShowerHandler,unsigned int> interfaceMaxTry
    ("MaxTry",
     "The maximum number of attempts to generate the shower from a"
     " particular ShowerTree",
     &QTildeShowerHandler::_maxtry, 100, 1, 100000,
     false, false, Interface::limited);

  static Parameter<QTildeShowerHandler,unsigned int> interfaceNReWeight
    ("NReWeight",
     "The number of attempts for the shower when reweighting",
     &QTildeShowerHandler::_nReWeight, 100, 10, 10000,
     false, false, Interface::limited);

  static Switch<QTildeShowerHandler, unsigned int> ifaceMECorrMode
    ("MECorrMode",
     "Choice of the ME Correction Mode",
     &QTildeShowerHandler::_meCorrMode, 1, false, false);
  static SwitchOption on
    (ifaceMECorrMode,"HardPlusSoft","hard+soft on", 1);
  static SwitchOption hard
    (ifaceMECorrMode,"Hard","only hard on", 2);
  static SwitchOption soft
    (ifaceMECorrMode,"Soft","only soft on", 3);

  static Switch<QTildeShowerHandler, bool> ifaceHardVetoReadOption
    ("HardVetoReadOption",
     "Apply read-in scale veto to all collisions or just the primary one?",
     &QTildeShowerHandler::_hardVetoReadOption, false, false, false);
  static SwitchOption AllCollisions
    (ifaceHardVetoReadOption,
     "AllCollisions",
     "Read-in pT veto applied to primary and secondary collisions.",
     false);
  static SwitchOption PrimaryCollision
    (ifaceHardVetoReadOption,
     "PrimaryCollision",
     "Read-in pT veto applied to primary but not secondary collisions.",
     true);

  static Parameter<QTildeShowerHandler, Energy> ifaceiptrms
    ("IntrinsicPtGaussian",
     "RMS of intrinsic pT of Gaussian distribution:\n"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &QTildeShowerHandler::_iptrms, GeV, ZERO, ZERO, 1000000.0*GeV,
     false, false, Interface::limited);

  static Parameter<QTildeShowerHandler, double> ifacebeta
    ("IntrinsicPtBeta",
     "Proportion of inverse quadratic distribution in generating intrinsic pT.\n"
     "(1-Beta) is the proportion of Gaussian distribution",
     &QTildeShowerHandler::_beta, 0, 0, 1,
     false, false, Interface::limited);

  static Parameter<QTildeShowerHandler, Energy> ifacegamma
    ("IntrinsicPtGamma",
     "Parameter for inverse quadratic:\n"
     "2*Beta*Gamma/(sqr(Gamma)+sqr(intrinsicpT))",
     &QTildeShowerHandler::_gamma,GeV, ZERO, ZERO, 100000.0*GeV,
     false, false, Interface::limited);

  static Parameter<QTildeShowerHandler, Energy> ifaceiptmax
    ("IntrinsicPtIptmax",
     "Upper bound on intrinsic pT for inverse quadratic",
     &QTildeShowerHandler::_iptmax,GeV, ZERO, ZERO, 100000.0*GeV,
     false, false, Interface::limited);

  static RefVector<QTildeShowerHandler,ShowerVeto> ifaceVetoes
    ("Vetoes",
     "The vetoes to be checked during showering",
     &QTildeShowerHandler::_vetoes, -1,
     false,false,true,true,false);
 
  static RefVector<QTildeShowerHandler,FullShowerVeto> interfaceFullShowerVetoes
    ("FullShowerVetoes",
     "The vetos to be appliede on the full final state of the shower",
     &QTildeShowerHandler::_fullShowerVetoes, -1, false, false, true, false, false);

  static Switch<QTildeShowerHandler,unsigned int> interfaceLimitEmissions
    ("LimitEmissions",
     "Limit the number and type of emissions for testing",
     &QTildeShowerHandler::_limitEmissions, 0, false, false);
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
  static SwitchOption interfaceLimitEmissionsOneEmission
    (interfaceLimitEmissions,
     "OneEmission",
     "Allow one emission in either the final state or initial state, but not both",
     4);

  static Switch<QTildeShowerHandler,bool> interfaceTruncMode
    ("TruncatedShower", "Include the truncated shower?", 
     &QTildeShowerHandler::_trunc_Mode, 1, false, false);
  static SwitchOption interfaceTruncMode0
    (interfaceTruncMode,"No","Truncated Shower is OFF", 0);
  static SwitchOption interfaceTruncMode1
    (interfaceTruncMode,"Yes","Truncated Shower is ON", 1);

  static Switch<QTildeShowerHandler,int> interfaceHardEmission
    ("HardEmission",
     "Whether to use ME corrections or POWHEG for the hardest emission",
     &QTildeShowerHandler::_hardEmission, 0, false, false);
  static SwitchOption interfaceHardEmissionNone
    (interfaceHardEmission,
     "None",
     "No Corrections",
     0);
  static SwitchOption interfaceHardEmissionMECorrection
    (interfaceHardEmission,
     "MECorrection",
     "Old fashioned ME correction",
     1);
  static SwitchOption interfaceHardEmissionPOWHEG
    (interfaceHardEmission,
     "POWHEG",
     "Powheg style hard emission",
     2);

  static Switch<QTildeShowerHandler,ShowerInteraction> interfaceInteractions
    ("Interactions",
     "The interactions to be used in the shower",
     &QTildeShowerHandler::interaction_, ShowerInteraction::Both, false, false);
  static SwitchOption interfaceInteractionsQCD
    (interfaceInteractions,
     "QCD",
     "Only QCD radiation",
     ShowerInteraction::QCD);
  static SwitchOption interfaceInteractionsQED
    (interfaceInteractions,
     "QED",
     "Only QEd radiation",
     ShowerInteraction::QED);
  static SwitchOption interfaceInteractionsQCDandQED
    (interfaceInteractions,
     "QCDandQED",
     "Both QED and QCD radiation",
     ShowerInteraction::Both);

  static Switch<QTildeShowerHandler,unsigned int> interfaceReconstructionOption
    ("ReconstructionOption",
     "Treatment of the reconstruction of the transverse momentum of "
     "a branching from the evolution scale.",
     &QTildeShowerHandler::_reconOpt, 0, false, false);
  static SwitchOption interfaceReconstructionOptionCutOff
    (interfaceReconstructionOption,
     "CutOff",
     "Use the cut-off masses in the calculation",
     0);
  static SwitchOption interfaceReconstructionOptionOffShell
    (interfaceReconstructionOption,
     "OffShell",
     "Use the off-shell masses in the calculation veto the emission of the parent,"
     " no veto in generation of emissions from children",
     1);
  static SwitchOption interfaceReconstructionOptionOffShell2
    (interfaceReconstructionOption,
     "OffShell2",
     "Use the off-shell masses in the calculation veto the emissions from the children."
     " no veto in generation of emissions from children",
     2);
  static SwitchOption interfaceReconstructionOptionOffShell3
    (interfaceReconstructionOption,
     "OffShell3",
     "Use the off-shell masses in the calculation veto the emissions from the children."
     " veto in generation of emissions from children using cut-off for second parton",
     3);
  static SwitchOption interfaceReconstructionOptionOffShell4
    (interfaceReconstructionOption,
     "OffShell4",
     "As OffShell3 but with a restriction on the mass of final-state"
     " jets produced via backward evolution.",
     4);
  static SwitchOption interfaceReconstructionOptionOffShell5
    (interfaceReconstructionOption,
     "OffShell5",
     "Try and preserve q2 but if pt negative just zero it",
     5);

  static Switch<QTildeShowerHandler,unsigned int> interfaceSpinCorrelations
    ("SpinCorrelations",
     "Treatment of spin correlations in the parton shower",
     &QTildeShowerHandler::_spinOpt, 1, false, false);
  static SwitchOption interfaceSpinCorrelationsNo
    (interfaceSpinCorrelations,
     "No",
     "No spin correlations",
     0);
  static SwitchOption interfaceSpinCorrelationsSpin
    (interfaceSpinCorrelations,
     "Yes",
     "Include the azimuthal spin correlations only",
     1);

  static Switch<QTildeShowerHandler,unsigned int> interfaceSoftCorrelations
    ("SoftCorrelations",
     "Option for the treatment of soft correlations in the parton shower",
     &QTildeShowerHandler::_softOpt, 2, false, false);
  static SwitchOption interfaceSoftCorrelationsNone
    (interfaceSoftCorrelations,
     "No",
     "No soft correlations",
     0);
  static SwitchOption interfaceSoftCorrelationsFull
    (interfaceSoftCorrelations,
     "Full",
     "Use the full eikonal",
     1);
  static SwitchOption interfaceSoftCorrelationsSingular
    (interfaceSoftCorrelations,
     "Singular",
     "Use original Webber-Marchisini form",
     2);

  static Switch<QTildeShowerHandler,bool> interfaceHardPOWHEG
    ("HardPOWHEG",
     "Treatment of powheg emissions which are too hard to have a shower interpretation",
     &QTildeShowerHandler::_hardPOWHEG, false, false, false);
  static SwitchOption interfaceHardPOWHEGAsShower
    (interfaceHardPOWHEG,
     "AsShower",
     "Still interpret as shower emissions",
     false);
  static SwitchOption interfaceHardPOWHEGRealEmission
    (interfaceHardPOWHEG,
     "RealEmission",
     "Generate shower from the real emmission configuration",
     true);

  static Parameter<QTildeShowerHandler,unsigned int> interfaceMaxTryFSR
    ("MaxTryFSR",
     "The maximum number of attempted FSR emissions in"
     " the generation of the FSR",
     &QTildeShowerHandler::_maxTryFSR, 100000, 10, 100000000,
     false, false, Interface::limited);

  static Parameter<QTildeShowerHandler,unsigned int> interfaceMaxFailFSR
    ("MaxFailFSR",
     "Maximum number of failures generating the FSR",
     &QTildeShowerHandler::_maxFailFSR, 100, 1, 100000000,
     false, false, Interface::limited);

  static Parameter<QTildeShowerHandler,double> interfaceFSRFailureFraction
    ("FSRFailureFraction",
     "Maximum fraction of events allowed to fail due to too many FSR emissions",
     &QTildeShowerHandler::_fracFSR, 0.001, 1e-10, 1,
     false, false, Interface::limited);

}

tPPair QTildeShowerHandler::cascade(tSubProPtr sub,
				    XCPtr xcomb) {
  // use me for reference in tex file etc
  useMe();
  prepareCascade(sub);
  // set things up in the base class
  resetWeights();
  hard_=ShowerTreePtr();
  decay_.clear();
  done_.clear();
  // start of the try block for the whole showering process
  unsigned int countFailures=0;
  while (countFailures<maxtry()) {
    try {
      decay_.clear();
      done_.clear();
      PerturbativeProcessPtr hard;
      DecayProcessMap decay;
      splitHardProcess(firstInteraction() ? tagged() :
      		       tPVector(currentSubProcess()->outgoing().begin(),
      				currentSubProcess()->outgoing().end()),
      		       hard,decay);
      ShowerTree::constructTrees(hard_,decay_,hard,decay);
      // if no hard process
      if(!hard_)  throw Exception() << "Shower starting with a decay"
				    << "is not implemented" 
				    << Exception::runerror;
      // perform the shower for the hard process
      showerHardProcess(hard_,xcomb);
      done_.push_back(hard_);
      hard_->updateAfterShower(decay_);
      // if no decaying particles to shower break out of the loop
      if(decay_.empty()) break;
      // shower the decay products
      while(!decay_.empty()) {
	// find particle whose production process has been showered
	ShowerDecayMap::iterator dit = decay_.begin();
	while(!dit->second->parent()->hasShowered() && dit!=decay_.end()) ++dit;
	assert(dit!=decay_.end());
	// get the particle
	ShowerTreePtr decayingTree = dit->second;
	// remove it from the multimap
	decay_.erase(dit);
	// make sure the particle has been decayed
	QTildeShowerHandler::decay(decayingTree,decay_);
	// now shower the decay
	showerDecay(decayingTree);
	done_.push_back(decayingTree);
	decayingTree->updateAfterShower(decay_);
      }
      // suceeded break out of the loop
      break;
    }
    catch (KinematicsReconstructionVeto) {
      resetWeights();
      ++countFailures;
    }
    catch ( ... ) {
      hard_=ShowerTreePtr();
      decay_.clear();
      done_.clear();
      throw;
    }
  }
  // if loop exited because of too many tries, throw event away
  if (countFailures >= maxtry()) {
    resetWeights();
    hard_=ShowerTreePtr();
    decay_.clear();
    done_.clear();
    throw Exception() << "Too many tries for main while loop "
		      << "in QTildeShowerHandler::cascade()." 
		      << Exception::eventerror; 	
  }
  //enter the particles in the event record
  fillEventRecord();
  // clear storage
  hard_=ShowerTreePtr();
  decay_.clear();
  done_.clear();
  // non hadronic case return
  if (!isResolvedHadron(incomingBeams().first ) && 
      !isResolvedHadron(incomingBeams().second) )
    return incomingBeams();
  // remake the remnants (needs to be after the colours are sorted
  //                       out in the insertion into the event record)
  if ( firstInteraction() ) return remakeRemnant(sub->incoming());
  //Return the new pair of incoming partons. remakeRemnant is not
  //necessary here, because the secondary interactions are not yet
  //connected to the remnants.
  return make_pair(findFirstParton(sub->incoming().first ),
		   findFirstParton(sub->incoming().second));
}

void QTildeShowerHandler::fillEventRecord() {
  // create a new step 
  StepPtr pstep = newStep();
  assert(!done_.empty());
  assert(done_[0]->isHard());
  // insert the steps
  for(unsigned int ix=0;ix<done_.size();++ix) {
    done_[ix]->fillEventRecord(pstep,doISR(),doFSR());
  }
}

HardTreePtr QTildeShowerHandler::generateCKKW(ShowerTreePtr ) const {
  return HardTreePtr();
}

void QTildeShowerHandler::doinit() {
  ShowerHandler::doinit();
  // interactions may have been changed through a setup file so we
  // clear it up here
  // calculate max no of FSR vetos
  _maxFailFSR = max(int(_maxFailFSR), int(_fracFSR*double(generator()->N())));
  // check on the reweighting
  for(unsigned int ix=0;ix<_fullShowerVetoes.size();++ix) {
    if(_fullShowerVetoes[ix]->behaviour()==1) {
      _reWeight = true;
      break;
    }
  }
  if(_reWeight && maximumTries()<_nReWeight) {
    throw Exception() << "Reweight being performed in the shower but the number of attempts for the"
		      << "shower is less than that for the reweighting.\n"
		      << "Maximum number of attempt for the shower "
		      << fullName() << ":MaxTry is " << maximumTries() << "\nand for reweighting is "
		      << fullName() << ":NReWeight is " << _nReWeight << "\n"
		      << "we recommend the number of attempts is 10 times the number for reweighting\n"
		      << Exception::runerror;
  }
}

void QTildeShowerHandler::generateIntrinsicpT(vector<ShowerProgenitorPtr> particlesToShower) {
  _intrinsic.clear();
  if ( !ipTon() || !doISR() ) return;
  // don't do anything for the moment for secondary scatters
  if( !firstInteraction() ) return;
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

void QTildeShowerHandler::setupMaximumScales(const vector<ShowerProgenitorPtr> & p,
				 XCPtr xcomb) {
  // let POWHEG events radiate freely
  if(_hardEmission==2&&hardTree()) {
    vector<ShowerProgenitorPtr>::const_iterator ckt = p.begin();
    for (; ckt != p.end(); ckt++) (*ckt)->maxHardPt(Constants::MaxEnergy);
    return;
  }
  // return if no vetos
  if (!restrictPhasespace()) return; 
  // find out if hard partonic subprocess.
  bool isPartonic(false); 
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
    cit = _currenttree->incomingLines().begin();
  Lorentz5Momentum pcm;
  for(; cit!=currentTree()->incomingLines().end(); ++cit) {
    pcm += cit->first->progenitor()->momentum();
    isPartonic |= cit->first->progenitor()->coloured();
  }
  // find minimum pt from hard process, the maximum pt from all outgoing
  // coloured lines (this is simpler and more general than
  // 2stu/(s^2+t^2+u^2)).  Maximum scale for scattering processes will
  // be transverse mass.
  Energy ptmax = generator()->maximumCMEnergy();
  // general case calculate the scale  
  if ( !hardScaleIsMuF() || (hardVetoReadOption()&&!firstInteraction()) ) {
    // scattering process
    if(currentTree()->isHard()) {
      assert(xcomb);
      // coloured incoming particles
      if (isPartonic) {
	map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cjt = currentTree()->outgoingLines().begin();
	for(; cjt!=currentTree()->outgoingLines().end(); ++cjt) {
	  if (cjt->first->progenitor()->coloured())
	    ptmax = min(ptmax,cjt->first->progenitor()->momentum().mt());
	}
      }
      if (ptmax == generator()->maximumCMEnergy() ) ptmax = pcm.m();
      if(hardScaleIsMuF()&&hardVetoReadOption()&&
	 !firstInteraction()) {
	ptmax=min(ptmax,sqrt(xcomb->lastShowerScale()));
      }
    } 
    // decay, incoming() is the decaying particle.
    else { 
      ptmax = currentTree()->incomingLines().begin()->first
	->progenitor()->momentum().mass(); 
    }
  }
  // hepeup.SCALUP is written into the lastXComb by the
  // LesHouchesReader itself - use this by user's choice. 
  // Can be more general than this. 
  else {
    if(currentTree()->isHard()) {
      assert(xcomb);
      ptmax = sqrt( xcomb->lastShowerScale() );
    }
    else {
      ptmax = currentTree()->incomingLines().begin()->first
	->progenitor()->momentum().mass(); 
    }
  }
  ptmax *= hardScaleFactor();
  // set maxHardPt for all progenitors.  For partonic processes this
  // is now the max pt in the FS, for non-partonic processes or
  // processes with no coloured FS the invariant mass of the IS
  vector<ShowerProgenitorPtr>::const_iterator ckt = p.begin();
  for (; ckt != p.end(); ckt++) (*ckt)->maxHardPt(ptmax);
}

void QTildeShowerHandler::setupHardScales(const vector<ShowerProgenitorPtr> & p,
					  XCPtr xcomb) {
  if ( hardScaleIsMuF() &&
       (!hardVetoReadOption() || firstInteraction()) ) {
    Energy hardScale = ZERO;
    if(currentTree()->isHard()) {
      assert(xcomb);
      hardScale = sqrt( xcomb->lastShowerScale() );
    }
    else {
      hardScale = currentTree()->incomingLines().begin()->first
	->progenitor()->momentum().mass(); 
    }
    hardScale *= hardScaleFactor();
    vector<ShowerProgenitorPtr>::const_iterator ckt = p.begin();
    for (; ckt != p.end(); ckt++) (*ckt)->hardScale(hardScale);
    muPt = hardScale;
  }
}

void QTildeShowerHandler::showerHardProcess(ShowerTreePtr hard, XCPtr xcomb) {
  _hardme = HwMEBasePtr();
  // extract the matrix element
  tStdXCombPtr lastXC = dynamic_ptr_cast<tStdXCombPtr>(xcomb);
  if(lastXC) {
    _hardme = dynamic_ptr_cast<HwMEBasePtr>(lastXC->matrixElement());
  }
  _decayme = HwDecayerBasePtr();
  // set the current tree
  currentTree(hard);
  hardTree(HardTreePtr());
  // work out the type of event
  currentTree()->xcombPtr(dynamic_ptr_cast<StdXCombPtr>(xcomb));
  currentTree()->identifyEventType();
  checkFlags();
  // generate the showering
  doShowering(true,xcomb);
}

RealEmissionProcessPtr QTildeShowerHandler::hardMatrixElementCorrection(bool hard) {
  // set the initial enhancement factors for the soft correction
  _initialenhance = 1.;
  _finalenhance   = 1.;
  // see if we can get the correction from the matrix element
  // or decayer
  RealEmissionProcessPtr real;
  if(hard) {
    if(_hardme&&_hardme->hasMECorrection()) {
      _hardme->initializeMECorrection(_currenttree->perturbativeProcess(),
     				      _initialenhance,_finalenhance);
      if(hardMEC())
     	real = 
	  _hardme->applyHardMatrixElementCorrection(_currenttree->perturbativeProcess());
    }
  }
  else {
    if(_decayme&&_decayme->hasMECorrection()) {
      _decayme->initializeMECorrection(_currenttree->perturbativeProcess(),
   				       _initialenhance,_finalenhance);
      if(hardMEC())
   	real = _decayme->applyHardMatrixElementCorrection(_currenttree->perturbativeProcess());
    }
  }
  return real;
}

ShowerParticleVector QTildeShowerHandler::createTimeLikeChildren(tShowerParticlePtr, IdList ids) {
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  // if same as definition create particles, otherwise create cc
  ShowerParticleVector children;
  for(unsigned int ix=0;ix<2;++ix) {
    children.push_back(new_ptr(ShowerParticle(ids[ix+1],true)));
    if(children[ix]->id()==_progenitor->id()&&!ids[ix+1]->stable()&&abs(ids[ix+1]->id())!=ParticleID::tauminus)
      children[ix]->set5Momentum(Lorentz5Momentum(_progenitor->progenitor()->mass()));
    else
      children[ix]->set5Momentum(Lorentz5Momentum(ids[ix+1]->mass()));
  }
  return children;
}

bool QTildeShowerHandler::timeLikeShower(tShowerParticlePtr particle, 
					 ShowerInteraction type,
					 Branching fb, bool first) {
  // don't do anything if not needed
  if(_limitEmissions == 1 || hardOnly() || 
     ( _limitEmissions == 2 && _nfs != 0) ||
     ( _limitEmissions == 4 && _nfs + _nis != 0) ) {
    if(particle->spinInfo()) particle->spinInfo()->develop();
    return false;
  }
  // too many tries
  if(_nFSR>=_maxTryFSR) {
    ++_nFailedFSR;
    // too many failed events
    if(_nFailedFSR>=_maxFailFSR)
      throw Exception() << "Too many events have failed due to too many shower emissions, in\n"
			<< "QTildeShowerHandler::timeLikeShower(). Terminating run\n"
			<< Exception::runerror;
    throw Exception() << "Too many attempted emissions in QTildeShowerHandler::timeLikeShower()\n"
		      << Exception::eventerror;
  }
  // generate the emission
  ShowerParticleVector children;
  int ntry=0;
  // generate the emission
  if(!fb.kinematics) 
    fb = selectTimeLikeBranching(particle,type,HardBranchingPtr());
  // no emission, return
  if(!fb.kinematics) {
    if(particle->spinInfo()) particle->spinInfo()->develop();
    return false;
  }
  Branching fc[2];
  bool setupChildren = true;
  while (ntry<50) {
    fc[0] = Branching();
    fc[1] = Branching();
    ++ntry;
    assert(fb.kinematics);
    // has emitted
    // Assign the shower kinematics to the emitting particle.
    if(setupChildren) {
      ++_nFSR;
      particle->showerKinematics(fb.kinematics);
      // check highest pT
      if(fb.kinematics->pT()>progenitor()->highestpT())
	progenitor()->highestpT(fb.kinematics->pT());
      // create the children
      children = createTimeLikeChildren(particle,fb.ids);
      // update the children
      particle->showerKinematics()->
	updateChildren(particle, children,fb.type,_reconOpt==3 || _reconOpt==4);
      // update number of emissions
      ++_nfs;
      if(_limitEmissions!=0) {
	if(children[0]->spinInfo()) children[0]->spinInfo()->develop();
	if(children[1]->spinInfo()) children[1]->spinInfo()->develop();
	if(particle->spinInfo()) particle->spinInfo()->develop();
	return true;
      }
      setupChildren = false;
    }
    // select branchings for children
    fc[0] = selectTimeLikeBranching(children[0],type,HardBranchingPtr());
    fc[1] = selectTimeLikeBranching(children[1],type,HardBranchingPtr());
    // old default
    if(_reconOpt==0||_reconOpt==5) {
      break;
    }
    // all other options
    else {
      // cut-off masses for the branching
      const vector<Energy> & virtualMasses = fb.sudakov->virtualMasses(fb.ids);
      // compute the masses of the children
      Energy masses[3];
      for(unsigned int ix=0;ix<2;++ix) {
   	if(fc[ix].kinematics) {
   	  const vector<Energy> & vm = fc[ix].sudakov->virtualMasses(fc[ix].ids);
   	  Energy2 q2 = 
   	    fc[ix].kinematics->z()*(1.-fc[ix].kinematics->z())*sqr(fc[ix].kinematics->scale());
   	  if(fc[ix].ids[0]->id()!=ParticleID::g) q2 += sqr(vm[0]);
   	  masses[ix+1] = sqrt(q2);
   	}
   	else {
   	  masses[ix+1] = virtualMasses[ix+1];
   	}
      }
      masses[0] = fb.ids[0]->id()!=ParticleID::g ? virtualMasses[0] : ZERO;
      double z = fb.kinematics->z();
      Energy2 pt2 = z*(1.-z)*(z*(1.-z)*sqr(fb.kinematics->scale()) + sqr(masses[0]))
   	- sqr(masses[1])*(1.-z) - sqr(masses[2])*z;
      if(pt2>=ZERO) break;
      // clean up the vetoed emission
      if(_reconOpt==1) {
	particle->showerKinematics(ShoKinPtr());
	for(unsigned int ix=0;ix<children.size();++ix)
	  particle->abandonChild(children[ix]);
	children.clear();
	if(particle->spinInfo()) particle->spinInfo()->decayVertex(VertexPtr());
	particle->vetoEmission(fb.type,fb.kinematics->scale());
	// generate the new emission
	fb = selectTimeLikeBranching(particle,type,HardBranchingPtr());
	// no emission, return
	if(!fb.kinematics) {
	  if(particle->spinInfo()) particle->spinInfo()->develop();
	  return false;
	}
	setupChildren = true;
	continue;
      }
      // clean up vetoed children
      else if(_reconOpt>=2) {
	// reset the scales for the children
	for(unsigned int ix=0;ix<2;++ix) {
	  if(fc[ix].kinematics)
	    children[ix]->vetoEmission(fc[ix].type,fc[ix].kinematics->scale());
	  else
	    children[ix]->vetoEmission(ShowerPartnerType::QCDColourLine,ZERO);
	  children[ix]->virtualMass(ZERO);
	} 
      }
    }
  };
  // shower the first  particle
  if(fc[0].kinematics) timeLikeShower(children[0],type,fc[0],false);
  if(children[0]->spinInfo()) children[0]->spinInfo()->develop();
  // shower the second particle
  if(fc[1].kinematics) timeLikeShower(children[1],type,fc[1],false);
  if(children[1]->spinInfo()) children[1]->spinInfo()->develop();
  if(_reconOpt>=1)
    particle->showerKinematics()->updateParent(particle, children,fb.type);
  // branching has happened
  if(first&&!children.empty())
    particle->showerKinematics()->resetChildren(particle,children);
  if(particle->spinInfo()) particle->spinInfo()->develop();
  return true;
}

bool 
QTildeShowerHandler::spaceLikeShower(tShowerParticlePtr particle, PPtr beam,
				     ShowerInteraction type) {
  //using the pdf's associated with the ShowerHandler assures, that
  //modified pdf's are used for the secondary interactions via 
  //CascadeHandler::resetPDFs(...)
  tcPDFPtr pdf;
  if(firstPDF().particle() == _beam)
    pdf = firstPDF().pdf();
  if(secondPDF().particle() == _beam)
    pdf = secondPDF().pdf();
  Energy freeze = pdfFreezingScale();
  // don't do anything if not needed
  if(_limitEmissions == 2  || hardOnly() ||
     ( _limitEmissions == 1 && _nis != 0 ) ||
     ( _limitEmissions == 4 && _nis + _nfs != 0 ) ) {
    if(particle->spinInfo()) particle->spinInfo()->develop();
    return false;
  }
  Branching bb;
  // generate branching
  while (true) {
    bb=_splittingGenerator->chooseBackwardBranching(*particle,beam,
						    _initialenhance,
						    _beam,type,
						    pdf,freeze);
    // return if no emission
    if(!bb.kinematics) {
      if(particle->spinInfo()) particle->spinInfo()->develop();
      return false;
    }
    // if not vetoed break
    if(!spaceLikeVetoed(bb,particle)) break;
    // otherwise reset scale and continue
    particle->vetoEmission(bb.type,bb.kinematics->scale());
    if(particle->spinInfo()) particle->spinInfo()->decayVertex(VertexPtr());
  }
  // assign the splitting function and shower kinematics
  particle->showerKinematics(bb.kinematics);
  if(bb.kinematics->pT()>progenitor()->highestpT())
    progenitor()->highestpT(bb.kinematics->pT());
  // For the time being we are considering only 1->2 branching
  // particles as in Sudakov form factor
  tcPDPtr part[2]={bb.ids[0],bb.ids[2]};
  // Now create the actual particles, make the otherChild a final state
  // particle, while the newParent is not
  ShowerParticlePtr newParent  = new_ptr(ShowerParticle(part[0],false));
  ShowerParticlePtr otherChild = new_ptr(ShowerParticle(part[1],true,true));
  ShowerParticleVector theChildren;
  theChildren.push_back(particle); 
  theChildren.push_back(otherChild);
  //this updates the evolution scale
  particle->showerKinematics()->
    updateParent(newParent, theChildren,bb.type);
  // update the history if needed
  _currenttree->updateInitialStateShowerProduct(_progenitor,newParent);
  _currenttree->addInitialStateBranching(particle,newParent,otherChild);
  // for the reconstruction of kinematics, parent/child
  // relationships are according to the branching process:
  // now continue the shower
  ++_nis;
  bool emitted = _limitEmissions==0 ? 
    spaceLikeShower(newParent,beam,type) : false;
  if(newParent->spinInfo()) newParent->spinInfo()->develop();
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
  particle->showerKinematics()->
    updateChildren(newParent, theChildren,bb.type,_reconOpt>=4);
  if(_limitEmissions!=0) {
    if(particle->spinInfo()) particle->spinInfo()->develop();
    return true;
  }
  // perform the shower of the final-state particle
  timeLikeShower(otherChild,type,Branching(),true);
  updateHistory(otherChild);
  if(theChildren[1]->spinInfo()) theChildren[1]->spinInfo()->develop();
  // return the emitted
  if(particle->spinInfo()) particle->spinInfo()->develop();
  return true;
}

void QTildeShowerHandler::showerDecay(ShowerTreePtr decay) {
  // work out the type of event
  currentTree()->xcombPtr(StdXCombPtr());
  currentTree()->identifyEventType();
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
    OrderedParticles outgoing;
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  it=decay->outgoingLines().begin();it!=decay->outgoingLines().end();++it) {
      if(abs(decay->incomingLines().begin()->first->original()->id()) == ParticleID::t &&
	 abs(it->first->original()->id())==ParticleID::Wplus &&
	 decay->treelinks().size() == 1) {
	ShowerTreePtr Wtree = decay->treelinks().begin()->first;
	for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator
	      it2=Wtree->outgoingLines().begin();it2!=Wtree->outgoingLines().end();++it2) {
	  outgoing.insert(it2->first->original()->dataPtr());
	}
      }
      else {
	outgoing.insert(it->first->original()->dataPtr());
      }
    }
    for(OrderedParticles::const_iterator it=outgoing.begin(); it!=outgoing.end();++it) {
      if(it!=outgoing.begin()) tag += ",";
      tag +=(**it).name();
    }
    tag += ";";
    dm = findDecayMode(tag);
  }
  if(dm) _decayme = dynamic_ptr_cast<HwDecayerBasePtr>(dm->decayer());
  // set the ShowerTree to be showered
  currentTree(decay);
  decay->applyTransforms();
  hardTree(HardTreePtr());
  // generate the showering
  doShowering(false,XCPtr());
  // if no vetos 
  // force calculation of spin correlations
  SpinPtr spInfo = decay->incomingLines().begin()->first->progenitor()->spinInfo();
  if(spInfo) {
    if(!spInfo->developed()) spInfo->needsUpdate();
    spInfo->develop();
  }
}

bool QTildeShowerHandler::spaceLikeDecayShower(tShowerParticlePtr particle,
				   const ShowerParticle::EvolutionScales & maxScales,
				   Energy minmass,ShowerInteraction type,
				   Branching fb) {
  // too many tries
  if(_nFSR>=_maxTryFSR) {
    ++_nFailedFSR;
    // too many failed events
    if(_nFailedFSR>=_maxFailFSR)
      throw Exception() << "Too many events have failed due to too many shower emissions, in\n"
			<< "QTildeShowerHandler::timeLikeShower(). Terminating run\n"
			<< Exception::runerror;
    throw Exception() << "Too many attempted emissions in QTildeShowerHandler::timeLikeShower()\n"
		      << Exception::eventerror;
  }
  // generate the emission
  ShowerParticleVector children;
  int ntry=0;
  // generate the emission
  if(!fb.kinematics) 
    fb = selectSpaceLikeDecayBranching(particle,maxScales,minmass,type,
				       HardBranchingPtr());
  // no emission, return
  if(!fb.kinematics) return false;
  Branching fc[2];
  bool setupChildren = true;
  while (ntry<50) {
    if(particle->virtualMass()==ZERO) 
      particle->virtualMass(_progenitor->progenitor()->mass());
    fc[0] = Branching();
    fc[1] = Branching();
    ++ntry;
    assert(fb.kinematics);
    // has emitted
    // Assign the shower kinematics to the emitting particle.
    if(setupChildren) {
      ++_nFSR;
      // Assign the shower kinematics to the emitting particle.
      particle->showerKinematics(fb.kinematics);
      if(fb.kinematics->pT()>progenitor()->highestpT())
	progenitor()->highestpT(fb.kinematics->pT());
      // create the ShowerParticle objects for the two children
      children = createTimeLikeChildren(particle,fb.ids);
      // updateChildren the children
      particle->showerKinematics()->
	updateChildren(particle, children, fb.type,_reconOpt>=3);
      setupChildren = false;
    }
    // select branchings for children
    fc[0] = selectSpaceLikeDecayBranching(children[0],maxScales,minmass,
					  type,HardBranchingPtr());
    fc[1] = selectTimeLikeBranching      (children[1],type,HardBranchingPtr());
    // old default
    if(_reconOpt==0) {
      // shower the first  particle
      _currenttree->updateInitialStateShowerProduct(_progenitor,children[0]);
      _currenttree->addInitialStateBranching(particle,children[0],children[1]);
      if(fc[0].kinematics) spaceLikeDecayShower(children[0],maxScales,minmass,type,Branching());
      // shower the second particle
      if(fc[1].kinematics) timeLikeShower(children[1],type,fc[1],true);
      updateHistory(children[1]);
      // branching has happened
      break;
    }
    // Herwig default
    else if(_reconOpt==1) {
      // shower the first  particle
      _currenttree->updateInitialStateShowerProduct(_progenitor,children[0]);
      _currenttree->addInitialStateBranching(particle,children[0],children[1]);
      if(fc[0].kinematics) spaceLikeDecayShower(children[0],maxScales,minmass,type,Branching());
      // shower the second particle
      if(fc[1].kinematics) timeLikeShower(children[1],type,fc[1],true);
      updateHistory(children[1]);
      // branching has happened
      particle->showerKinematics()->updateParent(particle, children,fb.type);
      // clean up the vetoed emission
      if(particle->virtualMass()==ZERO) {
    	particle->showerKinematics(ShoKinPtr());
   	for(unsigned int ix=0;ix<children.size();++ix)
   	  particle->abandonChild(children[ix]);
    	children.clear();
 	particle->vetoEmission(fb.type,fb.kinematics->scale());
 	// generate the new emission
  	fb = selectSpaceLikeDecayBranching(particle,maxScales,minmass,type,
				       HardBranchingPtr());
 	// no emission, return
  	if(!fb.kinematics) {
  	  return false;
   	}
  	setupChildren = true;
   	continue;
      }
      else
   	break;
    }
    else if(_reconOpt>=2) {
      // cut-off masses for the branching
      const vector<Energy> & virtualMasses = fb.sudakov->virtualMasses(fb.ids);
      // compute the masses of the children
      Energy masses[3];
      // space-like children
      masses[1] = children[0]->virtualMass();
      // time-like child
      if(fc[1].kinematics) {
	const vector<Energy> & vm = fc[1].sudakov->virtualMasses(fc[1].ids);
	Energy2 q2 = 
	  fc[1].kinematics->z()*(1.-fc[1].kinematics->z())*sqr(fc[1].kinematics->scale());
	if(fc[1].ids[0]->id()!=ParticleID::g) q2 += sqr(vm[0]);
	masses[2] = sqrt(q2);
      }
      else {
	masses[2] = virtualMasses[2];
      } 
      masses[0]=particle->virtualMass();
      double z = fb.kinematics->z();
      Energy2 pt2 = (1.-z)*(z*sqr(masses[0])-sqr(masses[1])-z/(1.-z)*sqr(masses[2]));
      if(pt2>=ZERO) {
  	break;
      }
      else {
  	// reset the scales for the children
  	for(unsigned int ix=0;ix<2;++ix) {
  	  if(fc[ix].kinematics)
  	    children[ix]->vetoEmission(fc[ix].type,fc[ix].kinematics->scale());
  	  else {
	    if(ix==0) 
	      children[ix]->vetoEmission(ShowerPartnerType::QCDColourLine,Constants::MaxEnergy);
	    else
	      children[ix]->vetoEmission(ShowerPartnerType::QCDColourLine,ZERO);
	  }
   	} 
	children[0]->virtualMass(_progenitor->progenitor()->mass());
	children[1]->virtualMass(ZERO);
      }
    }
  };
  if(_reconOpt>=2) {
    // In the case of splittings which involves coloured particles,
    // set properly the colour flow of the branching.
    // update the history if needed
    _currenttree->updateInitialStateShowerProduct(_progenitor,children[0]);
    _currenttree->addInitialStateBranching(particle,children[0],children[1]);
    // shower the first  particle
    if(fc[0].kinematics) spaceLikeDecayShower(children[0],maxScales,minmass,type,Branching());
    // shower the second particle
    if(fc[1].kinematics) timeLikeShower(children[1],type,fc[1],true);
    updateHistory(children[1]);
    // branching has happened
    particle->showerKinematics()->updateParent(particle, children,fb.type);
  }
  // branching has happened
  return true;
}

vector<ShowerProgenitorPtr> QTildeShowerHandler::setupShower(bool hard) {
  RealEmissionProcessPtr real;
  // generate hard me if needed
  if(_hardEmission==1) {
    real = hardMatrixElementCorrection(hard);
    if(real&&!real->outgoing().empty()) setupMECorrection(real);
  }
  // generate POWHEG hard emission if needed
  else if(_hardEmission==2) 
    hardestEmission(hard);
  // set the initial colour partners
  setEvolutionPartners(hard,interaction_,false);
  // get the particles to be showered
  vector<ShowerProgenitorPtr> particlesToShower = 
    currentTree()->extractProgenitors();
  // return the answer
  return particlesToShower;
}

void QTildeShowerHandler::setEvolutionPartners(bool hard,ShowerInteraction type,
					       bool clear) {
  // match the particles in the ShowerTree and hardTree
  if(hardTree() && !hardTree()->connect(currentTree()))
    throw Exception() << "Can't match trees in "
		      << "QTildeShowerHandler::setEvolutionPartners()"
		      << Exception::eventerror;
  // extract the progenitors
  vector<ShowerParticlePtr> particles = 
    currentTree()->extractProgenitorParticles();
  // clear the partners if needed
  if(clear) {
    for(unsigned int ix=0;ix<particles.size();++ix) {
      particles[ix]->partner(ShowerParticlePtr());
      particles[ix]->clearPartners();
    }
  }
  // sort out the colour partners
  if(hardTree()) {
    // find the partner
    for(unsigned int ix=0;ix<particles.size();++ix) {
      tShowerParticlePtr partner = hardTree()->particles()[particles[ix]]->branchingParticle()->partner();
      if(!partner) continue;
      for(map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator
	    it=hardTree()->particles().begin();
	  it!=hardTree()->particles().end();++it) {
	if(it->second->branchingParticle()==partner) {
	  particles[ix]->partner(it->first);
	  break;
	}
      }
      if(!particles[ix]->partner()) 
	throw Exception() << "Can't match partners in "
			  << "QTildeShowerHandler::setEvolutionPartners()"
			  << Exception::eventerror;
    }
  }

  // Set the initial evolution scales
  showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,!hard,interaction_,!_hardtree);
  if(hardTree() && _hardPOWHEG) {
    bool tooHard=false;
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit=hardTree()->particles().end();
    for(unsigned int ix=0;ix<particles.size();++ix) {
      map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
	mit = hardTree()->particles().find(particles[ix]);
      Energy hardScale(ZERO);
      ShowerPartnerType type(ShowerPartnerType::Undefined);
      // final-state
      if(particles[ix]->isFinalState()) {
	if(mit!= eit && !mit->second->children().empty()) {
	  hardScale = mit->second->scale();
	  type = mit->second->type();
	}
      }
      // initial-state
      else {
	if(mit!= eit && mit->second->parent()) {
	  hardScale = mit->second->parent()->scale();
	  type = mit->second->parent()->type();
	}
      }
      if(type!=ShowerPartnerType::Undefined) {
	if(type==ShowerPartnerType::QED) {
	  tooHard |= particles[ix]->scales().QED_noAO<hardScale;
	}
	else if(type==ShowerPartnerType::QCDColourLine) {
	  tooHard |= particles[ix]->scales().QCD_c_noAO<hardScale;
	}
	else if(type==ShowerPartnerType::QCDAntiColourLine) {
	  tooHard |= particles[ix]->scales().QCD_ac_noAO<hardScale;
	}
      }
    }
    if(tooHard) convertHardTree(hard,type);
  }
}

void QTildeShowerHandler::updateHistory(tShowerParticlePtr particle) {
  if(!particle->children().empty()) {
    ShowerParticleVector theChildren;
    for(unsigned int ix=0;ix<particle->children().size();++ix) {
      ShowerParticlePtr part = dynamic_ptr_cast<ShowerParticlePtr>
	(particle->children()[ix]);
      theChildren.push_back(part);
    }
    // update the history if needed
    if(particle==_currenttree->getFinalStateShowerProduct(_progenitor))
      _currenttree->updateFinalStateShowerProduct(_progenitor,
						  particle,theChildren);
    _currenttree->addFinalStateBranching(particle,theChildren);
    for(unsigned int ix=0;ix<theChildren.size();++ix)
      updateHistory(theChildren[ix]);
  }
}

bool QTildeShowerHandler::startTimeLikeShower(ShowerInteraction type) {
  _nFSR = 0;
  // initialize basis vectors etc
  if(!progenitor()->progenitor()->partner()) return false;
  progenitor()->progenitor()->initializeFinalState();
  if(hardTree()) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit=hardTree()->particles().end(),
      mit = hardTree()->particles().find(progenitor()->progenitor());
    if( mit != eit && !mit->second->children().empty() ) {
      bool output=truncatedTimeLikeShower(progenitor()->progenitor(),
					  mit->second ,type,Branching(),true);
      if(output) updateHistory(progenitor()->progenitor());
      return output;
    }
  }
  // do the shower
  bool output = hardOnly() ? false :
    timeLikeShower(progenitor()->progenitor() ,type,Branching(),true) ;
  if(output) updateHistory(progenitor()->progenitor());
  return output;
}

bool QTildeShowerHandler::startSpaceLikeShower(PPtr parent, ShowerInteraction type) {
  // initialise the basis vectors
  if(!progenitor()->progenitor()->partner()) return false;
  progenitor()->progenitor()->initializeInitialState(parent);
  if(hardTree()) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit =hardTree()->particles().end(),
      mit = hardTree()->particles().find(progenitor()->progenitor());
    if( mit != eit && mit->second->parent() ) {
      return truncatedSpaceLikeShower( progenitor()->progenitor(),
				       parent, mit->second->parent(), type );
    } 
  }
  // perform the shower
  return  hardOnly() ? false :
    spaceLikeShower(progenitor()->progenitor(),parent,type);
}

bool QTildeShowerHandler::
startSpaceLikeDecayShower(const ShowerParticle::EvolutionScales & maxScales,
			  Energy minimumMass,ShowerInteraction type) {
  _nFSR = 0;
  // set up the particle basis vectors
  if(!progenitor()->progenitor()->partner()) return false;
  progenitor()->progenitor()->initializeDecay();
  if(hardTree()) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit =hardTree()->particles().end(),
      mit = hardTree()->particles().find(progenitor()->progenitor());
    if( mit != eit && mit->second->parent() ) {
      HardBranchingPtr branch=mit->second;
      while(branch->parent()) branch=branch->parent();
      return truncatedSpaceLikeDecayShower(progenitor()->progenitor(),maxScales,
					   minimumMass, branch ,type, Branching());
    }
  }
  // perform the shower
  return  hardOnly() ? false :
    spaceLikeDecayShower(progenitor()->progenitor(),maxScales,minimumMass,type,Branching());
}

bool QTildeShowerHandler::timeLikeVetoed(const Branching & fb,
			     ShowerParticlePtr particle) {
  // work out type of interaction
  ShowerInteraction type = convertInteraction(fb.type);
  // check whether emission was harder than largest pt of hard subprocess
  if ( restrictPhasespace() && fb.kinematics->pT() > _progenitor->maxHardPt() ) 
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
  if(fb.kinematics->pT()>_progenitor->maximumpT(type)) return true;
  // general vetos
  if (fb.kinematics && !_vetoes.empty()) {
    bool vetoed=false;
    for (vector<ShowerVetoPtr>::iterator v = _vetoes.begin();
	 v != _vetoes.end(); ++v) {
      bool test = (**v).vetoTimeLike(_progenitor,particle,fb,currentTree());
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
  if ( firstInteraction() &&
       profileScales() ) {
    double weight = 
      profileScales()->
      hardScaleProfile(_progenitor->hardScale(),fb.kinematics->pT());
    if ( UseRandom::rnd() > weight )
      return true;
  }
  return false;
}

bool QTildeShowerHandler::spaceLikeVetoed(const Branching & bb,
			      ShowerParticlePtr particle) {
  // work out type of interaction
  ShowerInteraction type = convertInteraction(bb.type);
  // check whether emission was harder than largest pt of hard subprocess
  if (restrictPhasespace() && bb.kinematics->pT() > _progenitor->maxHardPt())
    return true;
  // apply the soft correction
  if( softMEC() && _hardme && _hardme->hasMECorrection() ) {
    if(_hardme->softMatrixElementVeto(_progenitor,particle,bb))
      return true;
  }
  // the more general vetos

  // check vs max pt for the shower
  if(bb.kinematics->pT()>_progenitor->maximumpT(type)) return true;

  if (!_vetoes.empty()) {
    bool vetoed=false;
    for (vector<ShowerVetoPtr>::iterator v = _vetoes.begin();
	 v != _vetoes.end(); ++v) {
      bool test = (**v).vetoSpaceLike(_progenitor,particle,bb,currentTree());
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
  if ( firstInteraction() &&
       profileScales() ) {
    double weight = 
      profileScales()->
      hardScaleProfile(_progenitor->hardScale(),bb.kinematics->pT());
    if ( UseRandom::rnd() > weight )
      return true;
  }
  return false;
}

bool QTildeShowerHandler::spaceLikeDecayVetoed( const Branching & fb,
				    ShowerParticlePtr particle) {
  // work out type of interaction
  ShowerInteraction type = convertInteraction(fb.type);
  // apply the soft correction
  if( softMEC() && _decayme && _decayme->hasMECorrection() ) {
    if(_decayme->softMatrixElementVeto(_progenitor,particle,fb))
      return true;
  }
  // veto on hardest pt in the shower
  if(fb.kinematics->pT()> _progenitor->maximumpT(type)) return true;
  // general vetos
  if (!_vetoes.empty()) {
    bool vetoed=false;
    for (vector<ShowerVetoPtr>::iterator v = _vetoes.begin();
	 v != _vetoes.end(); ++v) {
      bool test = (**v).vetoSpaceLike(_progenitor,particle,fb,currentTree());
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

void QTildeShowerHandler::hardestEmission(bool hard) {
  HardTreePtr ISRTree;
  // internal POWHEG in production or decay
  if( (( _hardme  &&  _hardme->hasPOWHEGCorrection()!=0 ) ||
       ( _decayme && _decayme->hasPOWHEGCorrection()!=0 ) ) ) {
    RealEmissionProcessPtr real;
    unsigned int type(0);
    // production
    if(_hardme) {
      assert(hard);
      real = _hardme->generateHardest( currentTree()->perturbativeProcess(),
				       interaction_);
      type = _hardme->hasPOWHEGCorrection();
    }
    // decay
    else {
      assert(!hard);
      real = _decayme->generateHardest( currentTree()->perturbativeProcess() );
      type = _decayme->hasPOWHEGCorrection();
    }
    if(real) {
      // set up ther hard tree
      if(!real->outgoing().empty()) _hardtree = new_ptr(HardTree(real));
      // set up the vetos
      currentTree()->setVetoes(real->pT(),type);
    }
    // store initial state POWHEG radiation
    if(_hardtree && _hardme && _hardme->hasPOWHEGCorrection()==1) 
      ISRTree = _hardtree;
  }
  else if (hard) {
    // Get minimum pT cutoff used in shower approximation
    Energy maxpt = 1.*GeV;

    if ( currentTree()->showerApproximation() ) {
      int colouredIn  = 0;
      int colouredOut = 0;
      for( map< ShowerProgenitorPtr, tShowerParticlePtr >::iterator it
	     = currentTree()->outgoingLines().begin();
	   it != currentTree()->outgoingLines().end(); ++it ) {
	if( it->second->coloured() ) ++colouredOut;
      }
      for( map< ShowerProgenitorPtr, ShowerParticlePtr >::iterator it
	     = currentTree()->incomingLines().begin();
	   it != currentTree()->incomingLines().end(); ++it ) {
	if( it->second->coloured() ) ++colouredIn;
      }
      if ( currentTree()->showerApproximation()->ffPtCut() == currentTree()->showerApproximation()->fiPtCut() &&
	   currentTree()->showerApproximation()->ffPtCut() == currentTree()->showerApproximation()->iiPtCut() )
	maxpt = currentTree()->showerApproximation()->ffPtCut();
      else if ( colouredIn == 2 && colouredOut == 0 )
	maxpt = currentTree()->showerApproximation()->iiPtCut();
      else if ( colouredIn == 0 && colouredOut > 1 )
	maxpt = currentTree()->showerApproximation()->ffPtCut();
      else if ( colouredIn == 2 && colouredOut == 1 )
	maxpt = min(currentTree()->showerApproximation()->iiPtCut(), currentTree()->showerApproximation()->fiPtCut());
      else if ( colouredIn == 1 && colouredOut > 1 )
	maxpt = min(currentTree()->showerApproximation()->ffPtCut(), currentTree()->showerApproximation()->fiPtCut());
      else
	maxpt = min(min(currentTree()->showerApproximation()->iiPtCut(), currentTree()->showerApproximation()->fiPtCut()), 
		    currentTree()->showerApproximation()->ffPtCut());
    }

    // Generate hardtree from born and real emission subprocesses
    _hardtree = generateCKKW(currentTree());

    // Find transverse momentum of hardest emission
    if (_hardtree){
      for(set<HardBranchingPtr>::iterator it=_hardtree->branchings().begin();
     	  it!=_hardtree->branchings().end();++it) {
      	if ((*it)->parent() && (*it)->status()==HardBranching::Incoming)
      	  maxpt=(*it)->branchingParticle()->momentum().perp();
      	if ((*it)->children().size()==2 && (*it)->status()==HardBranching::Outgoing){
	  if ((*it)->branchingParticle()->id()!=21 &&
	      abs((*it)->branchingParticle()->id())>5 ){
	    if ((*it)->children()[0]->branchingParticle()->id()==21 ||
		abs((*it)->children()[0]->branchingParticle()->id())<6)
	      maxpt=(*it)->children()[0]->branchingParticle()->momentum().perp();
	    else if ((*it)->children()[1]->branchingParticle()->id()==21 ||
		     abs((*it)->children()[1]->branchingParticle()->id())<6)
	      maxpt=(*it)->children()[1]->branchingParticle()->momentum().perp();
	  }
	  else {
	    if ( abs((*it)->branchingParticle()->id())<6){
	      if (abs((*it)->children()[0]->branchingParticle()->id())<6)
		maxpt = (*it)->children()[1]->branchingParticle()->momentum().perp();
	      else 
		maxpt = (*it)->children()[0]->branchingParticle()->momentum().perp();
	    }
	    else maxpt = (*it)->children()[1]->branchingParticle()->momentum().perp();
	  }
      	}
      } 
    }
     
    
    // Hardest (pt) emission should be the first powheg emission.
    maxpt=min(sqrt(lastXCombPtr()->lastShowerScale()),maxpt);

    // set maximum pT for subsequent emissions from S events
    if ( currentTree()->isPowhegSEvent() ) {
      for( map< ShowerProgenitorPtr, tShowerParticlePtr >::iterator it
	     = currentTree()->outgoingLines().begin(); 
	   it != currentTree()->outgoingLines().end(); ++it ) {
	if( ! it->second->coloured() ) continue;
	it->first->maximumpT(maxpt, ShowerInteraction::QCD  );
      }  
      for( map< ShowerProgenitorPtr, ShowerParticlePtr >::iterator it
	     = currentTree()->incomingLines().begin(); 
	   it != currentTree()->incomingLines().end(); ++it ) {
	if( ! it->second->coloured() ) continue;
	it->first->maximumpT(maxpt, ShowerInteraction::QCD );
      }
    }
  }
  else 
    _hardtree = generateCKKW(currentTree());

  // if hard me doesn't have a FSR powheg 
  // correction use decay powheg correction
  if (_hardme && _hardme->hasPOWHEGCorrection()<2) {
    addFSRUsingDecayPOWHEG(ISRTree);
  }
  // connect the trees
  if(_hardtree) {
    connectTrees(currentTree(),_hardtree,hard); 
  }
}

void QTildeShowerHandler::addFSRUsingDecayPOWHEG(HardTreePtr ISRTree) {
  // check for intermediate colour singlet resonance
  const ParticleVector inter =  _hardme->subProcess()->intermediates();
  if (inter.size()!=1 || inter[0]->momentum().m2()/GeV2 < 0 || 
      inter[0]->dataPtr()->iColour()!=PDT::Colour0) {
    return;
  }
   
  // ignore cases where outgoing particles are not coloured
  map<ShowerProgenitorPtr, tShowerParticlePtr > out = currentTree()->outgoingLines();
  if (out.size() != 2 ||
      out. begin()->second->dataPtr()->iColour()==PDT::Colour0 ||
      out.rbegin()->second->dataPtr()->iColour()==PDT::Colour0) {
    return;
  }

  // look up decay mode
  tDMPtr dm;
  string tag;
  string inParticle = inter[0]->dataPtr()->name() + "->";
  vector<string> outParticles;
  outParticles.push_back(out.begin ()->first->progenitor()->dataPtr()->name());
  outParticles.push_back(out.rbegin()->first->progenitor()->dataPtr()->name());
  for (int it=0; it<2; ++it){
    tag = inParticle + outParticles[it] + "," + outParticles[(it+1)%2] + ";";
    dm = generator()->findDecayMode(tag);
    if(dm) break;
  }

  // get the decayer
  HwDecayerBasePtr decayer;   
  if(dm) decayer = dynamic_ptr_cast<HwDecayerBasePtr>(dm->decayer());
  // check if decayer has a FSR POWHEG correction 
  if (!decayer || decayer->hasPOWHEGCorrection()<2) {
    return;
  }
  // generate the hardest emission
  // create RealEmissionProcess
  PPtr in = new_ptr(*inter[0]);
  RealEmissionProcessPtr newProcess(new_ptr(RealEmissionProcess()));
  newProcess->bornIncoming().push_back(in);
  newProcess->bornOutgoing().push_back(out.begin ()->first->progenitor());
  newProcess->bornOutgoing().push_back(out.rbegin()->first->progenitor());
  // generate the FSR
  newProcess = decayer->generateHardest(newProcess);
  HardTreePtr FSRTree;
  if(newProcess) {
    // set up ther hard tree
    if(!newProcess->outgoing().empty()) FSRTree = new_ptr(HardTree(newProcess));
    // set up the vetos
    currentTree()->setVetoes(newProcess->pT(),2);
  }
  if(!FSRTree) return;
  
  // if there is no ISRTree make _hardtree from FSRTree
  if (!ISRTree){
    vector<HardBranchingPtr> inBranch,hardBranch;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator
	  cit =currentTree()->incomingLines().begin();
	cit!=currentTree()->incomingLines().end();++cit ) {
      inBranch.push_back(new_ptr(HardBranching(cit->second,SudakovPtr(),
					       HardBranchingPtr(),
					       HardBranching::Incoming)));
      inBranch.back()->beam(cit->first->original()->parents()[0]);
      hardBranch.push_back(inBranch.back());
    }
    if(inBranch[0]->branchingParticle()->dataPtr()->coloured()) {
      inBranch[0]->colourPartner(inBranch[1]);
      inBranch[1]->colourPartner(inBranch[0]);
    }
    for(set<HardBranchingPtr>::iterator it=FSRTree->branchings().begin();
	it!=FSRTree->branchings().end();++it) {
      if((**it).branchingParticle()->id()!=in->id()) 
	hardBranch.push_back(*it);
    } 
    hardBranch[2]->colourPartner(hardBranch[3]);
    hardBranch[3]->colourPartner(hardBranch[2]);
    HardTreePtr newTree = new_ptr(HardTree(hardBranch,inBranch,
					   ShowerInteraction::QCD));            
    _hardtree = newTree;    
  }
  
  // Otherwise modify the ISRTree to include the emission in FSRTree
  else {
    vector<tShowerParticlePtr> FSROut, ISROut;   
    set<HardBranchingPtr>::iterator itFSR, itISR;
    // get outgoing particles 
    for(itFSR =FSRTree->branchings().begin();
	itFSR!=FSRTree->branchings().end();++itFSR){
      if ((**itFSR).status()==HardBranching::Outgoing) 
	FSROut.push_back((*itFSR)->branchingParticle());
    }     
    for(itISR =ISRTree->branchings().begin();
	itISR!=ISRTree->branchings().end();++itISR){
      if ((**itISR).status()==HardBranching::Outgoing) 
	ISROut.push_back((*itISR)->branchingParticle());
    }
    
    // find COM frame formed by outgoing particles
    LorentzRotation eventFrameFSR, eventFrameISR;
    eventFrameFSR = ((FSROut[0]->momentum()+FSROut[1]->momentum()).findBoostToCM());  
    eventFrameISR = ((ISROut[0]->momentum()+ISROut[1]->momentum()).findBoostToCM());

    // find rotation between ISR and FSR frames
    int j=0;
    if (ISROut[0]->id()!=FSROut[0]->id()) j=1;
    eventFrameISR.rotateZ( (eventFrameFSR*FSROut[0]->momentum()).phi()-
			   (eventFrameISR*ISROut[j]->momentum()).phi() );
    eventFrameISR.rotateY( (eventFrameFSR*FSROut[0]->momentum()).theta()-
			   (eventFrameISR*ISROut[j]->momentum()).theta() );
    eventFrameISR.invert();
    
    for (itFSR=FSRTree->branchings().begin();
	 itFSR!=FSRTree->branchings().end();++itFSR){
      if ((**itFSR).branchingParticle()->id()==in->id()) continue;
      for (itISR =ISRTree->branchings().begin();
	   itISR!=ISRTree->branchings().end();++itISR){
	if ((**itISR).status()==HardBranching::Incoming) continue;
	if ((**itFSR).branchingParticle()->id()==
	    (**itISR).branchingParticle()->id()){
	  // rotate FSRTree particle to ISRTree event frame
	  (**itISR).branchingParticle()->setMomentum(eventFrameISR*
						     eventFrameFSR*
						     (**itFSR).branchingParticle()->momentum());
	  (**itISR).branchingParticle()->rescaleMass();
	  // add the children of the FSRTree particles to the ISRTree
	  if(!(**itFSR).children().empty()){
	    (**itISR).addChild((**itFSR).children()[0]);
	    (**itISR).addChild((**itFSR).children()[1]);
	    // rotate momenta to ISRTree event frame
	    (**itISR).children()[0]->branchingParticle()->setMomentum(eventFrameISR*
								      eventFrameFSR*
								      (**itFSR).children()[0]->branchingParticle()->momentum());
	    (**itISR).children()[1]->branchingParticle()->setMomentum(eventFrameISR*
								      eventFrameFSR*
								      (**itFSR).children()[1]->branchingParticle()->momentum());
	  }
	}
      }	
    }
    _hardtree = ISRTree;
  }
}

bool QTildeShowerHandler::truncatedTimeLikeShower(tShowerParticlePtr particle,
				      HardBranchingPtr branch,
				      ShowerInteraction type,
				      Branching fb, bool first) {
  // select a branching if we don't have one
  if(!fb.kinematics)
    fb = selectTimeLikeBranching(particle,type,branch);
  // must be an emission, the forced one it not a truncated one
  assert(fb.kinematics);
  ShowerParticleVector children;
  int ntry=0;
  Branching fc[2];
  bool setupChildren = true;
  while (ntry<50) {
    if(!fc[0].hard) fc[0] = Branching();
    if(!fc[1].hard) fc[1] = Branching();
    ++ntry;
    // Assign the shower kinematics to the emitting particle.
    if(setupChildren) {
      ++_nFSR;
      // Assign the shower kinematics to the emitting particle.
      particle->showerKinematics(fb.kinematics);
      if(fb.kinematics->pT()>progenitor()->highestpT())
	progenitor()->highestpT(fb.kinematics->pT());
      // create the children
      children = createTimeLikeChildren(particle,fb.ids);
      // update the children
      particle->showerKinematics()->
   	updateChildren(particle, children,fb.type,_reconOpt>=3);
      setupChildren = false;
    }
    // select branchings for children
    if(!fc[0].kinematics) {
      // select branching for first particle
      if(!fb.hard && fb.iout ==1 )
	fc[0] = selectTimeLikeBranching(children[0],type,branch);
      else if(fb.hard && !branch->children()[0]->children().empty() )
	fc[0] = selectTimeLikeBranching(children[0],type,branch->children()[0]);
      else
	fc[0] = selectTimeLikeBranching(children[0],type,HardBranchingPtr());
    }
    // select branching for the second particle
    if(!fc[1].kinematics) {
      // select branching for first particle
      if(!fb.hard && fb.iout ==2 )
	fc[1] = selectTimeLikeBranching(children[1],type,branch);
      else if(fb.hard && !branch->children()[1]->children().empty() )
	fc[1] = selectTimeLikeBranching(children[1],type,branch->children()[1]);
      else
	fc[1] = selectTimeLikeBranching(children[1],type,HardBranchingPtr());
    }
    // old default
    if(_reconOpt==0 || (_reconOpt==1 && fb.hard) ) {
      // shower the first  particle
      if(fc[0].kinematics) {
	// the parent has truncated emission and following line
	if(!fb.hard && fb.iout == 1)
	  truncatedTimeLikeShower(children[0],branch,type,fc[0],false);
	// hard emission and subsquent hard emissions
	else if(fb.hard && !branch->children()[0]->children().empty() )
	  truncatedTimeLikeShower(children[0],branch->children()[0],type,fc[0],false);
	// normal shower
	else
	  timeLikeShower(children[0],type,fc[0],false);
      }
      if(children[0]->spinInfo()) children[0]->spinInfo()->develop();
      // shower the second particle
      if(fc[1].kinematics) {
	// the parent has truncated emission and following line
	if(!fb.hard && fb.iout == 2)
	  truncatedTimeLikeShower(children[1],branch,type,fc[1],false);
	// hard emission and subsquent hard emissions
	else if(fb.hard && !branch->children()[1]->children().empty() )
	  truncatedTimeLikeShower(children[1],branch->children()[1],type,fc[1],false);
	else
	  timeLikeShower(children[1],type,fc[1],false);
      }
      if(children[1]->spinInfo()) children[1]->spinInfo()->develop();
      // branching has happened
      particle->showerKinematics()->updateParent(particle, children,fb.type);
      break;
    }
    // H7 default
    else if(_reconOpt==1) {
      // shower the first  particle
      if(fc[0].kinematics) {
	// the parent has truncated emission and following line
	if(!fb.hard && fb.iout == 1)
	  truncatedTimeLikeShower(children[0],branch,type,fc[0],false);
	else
	  timeLikeShower(children[0],type,fc[0],false);
      }
      if(children[0]->spinInfo()) children[0]->spinInfo()->develop();
      // shower the second particle
      if(fc[1].kinematics) {
	// the parent has truncated emission and following line
	if(!fb.hard && fb.iout == 2)
	  truncatedTimeLikeShower(children[1],branch,type,fc[1],false);
	else
	  timeLikeShower(children[1],type,fc[1],false);
      }
      if(children[1]->spinInfo()) children[1]->spinInfo()->develop();
      // branching has happened
      particle->showerKinematics()->updateParent(particle, children,fb.type);
      // clean up the vetoed emission
      if(particle->virtualMass()==ZERO) {
   	particle->showerKinematics(ShoKinPtr());
    	for(unsigned int ix=0;ix<children.size();++ix)
   	  particle->abandonChild(children[ix]);
    	children.clear();
   	if(particle->spinInfo()) particle->spinInfo()->decayVertex(VertexPtr());
  	particle->vetoEmission(fb.type,fb.kinematics->scale());
   	// generate the new emission
   	fb = selectTimeLikeBranching(particle,type,branch);
	// must be at least hard emission
	assert(fb.kinematics);
   	setupChildren = true;
  	continue;
      }
      else
   	break;
    }
    else if(_reconOpt>=2) {
      // cut-off masses for the branching
      const vector<Energy> & virtualMasses = fb.sudakov->virtualMasses(fb.ids);
      // compute the masses of the children
      Energy masses[3];
      for(unsigned int ix=0;ix<2;++ix) {
   	if(fc[ix].kinematics) {
   	  const vector<Energy> & vm = fc[ix].sudakov->virtualMasses(fc[ix].ids);
   	  Energy2 q2 = 
   	    fc[ix].kinematics->z()*(1.-fc[ix].kinematics->z())*sqr(fc[ix].kinematics->scale());
   	  if(fc[ix].ids[0]->id()!=ParticleID::g) q2 += sqr(vm[0]);
   	  masses[ix+1] = sqrt(q2);
   	}
   	else {
   	  masses[ix+1] = virtualMasses[ix+1];
   	}
      }
      masses[0] = fb.ids[0]->id()!=ParticleID::g ? virtualMasses[0] : ZERO;
      double z = fb.kinematics->z();
      Energy2 pt2 = z*(1.-z)*(z*(1.-z)*sqr(fb.kinematics->scale()) + sqr(masses[0]))
   	- sqr(masses[1])*(1.-z) - sqr(masses[2])*z;
      if(pt2>=ZERO) {
  	break;
      }
      // if only the hard emission have to accept it
      else if ((fc[0].hard && !fc[1].kinematics) ||
	       (fc[1].hard && !fc[0].kinematics) ) {
	break;
      }
      else {
  	// reset the scales for the children
   	for(unsigned int ix=0;ix<2;++ix) {
	  if(fc[ix].hard) continue;
   	  if(fc[ix].kinematics && ! fc[ix].hard )
   	    children[ix]->vetoEmission(fc[ix].type,fc[ix].kinematics->scale());
  	  else
  	    children[ix]->vetoEmission(ShowerPartnerType::QCDColourLine,ZERO);
  	  children[ix]->virtualMass(ZERO);
  	} 
      }
    }
  };
  if(_reconOpt>=2) {
    // shower the first  particle
    if(fc[0].kinematics) {
      // the parent has truncated emission and following line
      if(!fb.hard && fb.iout == 1)
	truncatedTimeLikeShower(children[0],branch,type,fc[0],false);
      // hard emission and subsquent hard emissions
      else if(fb.hard && !branch->children()[0]->children().empty() )
	truncatedTimeLikeShower(children[0],branch->children()[0],type,fc[0],false);
      // normal shower
      else
	timeLikeShower(children[0],type,fc[0],false);
    }
    if(children[0]->spinInfo()) children[0]->spinInfo()->develop();
    // shower the second particle
    if(fc[1].kinematics) {
      // the parent has truncated emission and following line
      if(!fb.hard && fb.iout == 2)
	truncatedTimeLikeShower(children[1],branch,type,fc[1],false);
      // hard emission and subsquent hard emissions
      else if(fb.hard && !branch->children()[1]->children().empty() )
	truncatedTimeLikeShower(children[1],branch->children()[1],type,fc[1],false);
      else
	timeLikeShower(children[1],type,fc[1],false);
    }
    if(children[1]->spinInfo()) children[1]->spinInfo()->develop();
    // branching has happened
    particle->showerKinematics()->updateParent(particle, children,fb.type);
  }
  if(first&&!children.empty())
    particle->showerKinematics()->resetChildren(particle,children);
  if(particle->spinInfo()) particle->spinInfo()->develop();
  return true;
}

bool QTildeShowerHandler::truncatedSpaceLikeShower(tShowerParticlePtr particle, PPtr beam,
				       HardBranchingPtr branch,
				       ShowerInteraction type) {
  tcPDFPtr pdf;
  if(firstPDF().particle()  == beamParticle())
    pdf = firstPDF().pdf();
  if(secondPDF().particle() == beamParticle())
    pdf = secondPDF().pdf();
  Energy freeze = pdfFreezingScale();
  Branching bb;
  // parameters of the force branching
  double z(0.);
  HardBranchingPtr timelike;
  for( unsigned int ix = 0; ix < branch->children().size(); ++ix ) {
    if( branch->children()[ix]->status() ==HardBranching::Outgoing) {
      timelike = branch->children()[ix];
    }
    if( branch->children()[ix]->status() ==HardBranching::Incoming )
      z = branch->children()[ix]->z();
  }
  // generate truncated branching
  tcPDPtr part[2];
  if(z>=0.&&z<=1.) {
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
      part[0] = bb.ids[0];
      part[1] = bb.ids[2];
      double zsplit = bb.kinematics->z();
      // apply the vetos for the truncated shower
      // if doesn't carry most of momentum
      ShowerInteraction type2 = convertInteraction(bb.type);
      if(type2==branch->sudakov()->interactionType() &&
	 zsplit < 0.5) {
	particle->vetoEmission(bb.type,bb.kinematics->scale());
	continue;
      }
      // others
      if( part[0]->id() != particle->id() || // if particle changes type
	  bb.kinematics->pT() > progenitor()->maximumpT(type2) ||   // pt veto
	  bb.kinematics->scale() < branch->scale()) { // angular ordering veto
	particle->vetoEmission(bb.type,bb.kinematics->scale());
	continue;
      }
      // and those from the base class
      if(spaceLikeVetoed(bb,particle)) {
	particle->vetoEmission(bb.type,bb.kinematics->scale());
	continue;
      }
      break;
    }
  }
  if( !bb.kinematics ) {
    //do the hard emission
    ShoKinPtr kinematics =
      branch->sudakov()->createInitialStateBranching( branch->scale(), z, branch->phi(),
    						      branch->children()[0]->pT() );
    // assign the splitting function and shower kinematics
    particle->showerKinematics( kinematics );
    if(kinematics->pT()>progenitor()->highestpT())
      progenitor()->highestpT(kinematics->pT());
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
      updateParent( newParent, theChildren, branch->type());
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
      updateChildren( newParent, theChildren,bb.type,false);
    if(hardOnly()) return true;
    // perform the shower of the final-state particle
    if( timelike->children().empty() ) {
      timeLikeShower( otherChild , type,Branching(),true);
    }
    else {
      truncatedTimeLikeShower( otherChild, timelike , type,Branching(), true);
    }
    updateHistory(otherChild);
    // return the emitted
    return true;
  }
  // assign the splitting function and shower kinematics
  particle->showerKinematics( bb.kinematics );
  if(bb.kinematics->pT()>progenitor()->highestpT())
    progenitor()->highestpT(bb.kinematics->pT());
  // For the time being we are considering only 1->2 branching
  // Now create the actual particles, make the otherChild a final state
  // particle, while the newParent is not
  ShowerParticlePtr newParent = new_ptr( ShowerParticle( part[0], false ) );
  ShowerParticlePtr otherChild = new_ptr( ShowerParticle( part[1], true, true ) );
  ShowerParticleVector theChildren; 
  theChildren.push_back( particle ); 
  theChildren.push_back( otherChild );
  particle->showerKinematics()->
    updateParent( newParent, theChildren, bb.type);
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
  particle->showerKinematics()->
    updateChildren( newParent, theChildren, bb.type,false);
  // perform the shower of the final-state particle
  timeLikeShower( otherChild , type,Branching(),true);
  updateHistory(otherChild);
  // return the emitted
  return true;
}

bool QTildeShowerHandler::
truncatedSpaceLikeDecayShower(tShowerParticlePtr particle, 
			      const ShowerParticle::EvolutionScales & maxScales,
			      Energy minmass, HardBranchingPtr branch,
			      ShowerInteraction type, Branching fb) {
  // select a branching if we don't have one
  if(!fb.kinematics)
    fb = selectSpaceLikeDecayBranching(particle,maxScales,minmass,type,branch);
  // must be an emission, the forced one it not a truncated one
  assert(fb.kinematics);
  ShowerParticleVector children;
  int ntry=0;
  Branching fc[2];
  bool setupChildren = true;
  while (ntry<50) {
    if(!fc[0].hard) fc[0] = Branching();
    if(!fc[1].hard) fc[1] = Branching();
    ++ntry;
    if(setupChildren) {
      ++_nFSR;
      // Assign the shower kinematics to the emitting particle.
      particle->showerKinematics(fb.kinematics);
      if(fb.kinematics->pT()>progenitor()->highestpT())
	progenitor()->highestpT(fb.kinematics->pT());
      // create the ShowerParticle objects for the two children
      children = createTimeLikeChildren(particle,fb.ids);
      // updateChildren the children
      particle->showerKinematics()->
	updateChildren(particle, children, fb.type,_reconOpt>=3);
      setupChildren = false;
    }
    // select branchings for children
    if(!fc[0].kinematics) {
      if(children[0]->id()==particle->id()) {
	// select branching for first particle
	if(!fb.hard)
	  fc[0] = selectSpaceLikeDecayBranching(children[0],maxScales,minmass,type,branch);
	else if(fb.hard && ! branch->children()[0]->children().empty() )
	  fc[0] = selectSpaceLikeDecayBranching(children[0],maxScales,minmass,type,
						branch->children()[0]);
	else
	  fc[0] = selectSpaceLikeDecayBranching(children[0],maxScales,minmass,type,
						HardBranchingPtr());
      }
      else {
	// select branching for first particle
	if(fb.hard && !branch->children()[0]->children().empty() )
	  fc[0] = selectTimeLikeBranching(children[0],type,branch->children()[0]);
	else
	  fc[0] = selectTimeLikeBranching(children[0],type,HardBranchingPtr());
      }
    }
    // select branching for the second particle
     if(!fc[1].kinematics) {
       if(children[1]->id()==particle->id()) {
	 // select branching for first particle
	 if(!fb.hard)
	   fc[1] = selectSpaceLikeDecayBranching(children[1],maxScales,minmass,type,branch);
	 else if(fb.hard && ! branch->children()[1]->children().empty() )
	   fc[1] = selectSpaceLikeDecayBranching(children[1],maxScales,minmass,type,
						 branch->children()[1]);
	 else
	   fc[1] = selectSpaceLikeDecayBranching(children[1],maxScales,minmass,type,
						 HardBranchingPtr());
       }
       else {
	 if(fb.hard && !branch->children()[1]->children().empty() )
	   fc[1] = selectTimeLikeBranching(children[1],type,branch->children()[1]);
	 else
	   fc[1] = selectTimeLikeBranching(children[1],type,HardBranchingPtr());
       }
     }
    // old default
    if(_reconOpt==0 || (_reconOpt==1 && fb.hard) ) {
      // update the history if needed
      currentTree()->updateInitialStateShowerProduct(progenitor(),children[0]);
      currentTree()->addInitialStateBranching(particle,children[0],children[1]);
      // shower the first  particle
      if(fc[0].kinematics) {
	if(children[0]->id()==particle->id()) {
	  if(!fb.hard)
	    truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					   branch,type,fc[0]);
	  else if(fb.hard && ! branch->children()[0]->children().empty() )
	    truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					   branch->children()[0],type,fc[0]);
	  else
	    spaceLikeDecayShower( children[0],maxScales,minmass,type,fc[0]);
	}
	else {
	  if(fb.hard && !branch->children()[0]->children().empty() )
	    truncatedTimeLikeShower(children[0],branch->children()[0],type,fc[0],false);
	  // normal shower
	  else
	    timeLikeShower(children[0],type,fc[0],false);
	}
      }
      // shower the second  particle
      if(fc[1].kinematics) {
	if(children[0]->id()==particle->id()) {
	  if(!fb.hard)
	    truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					   branch,type,fc[1]);
	  else if(fb.hard && ! branch->children()[0]->children().empty() )
	    truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					   branch->children()[0],type,fc[1]);
	  else
	    spaceLikeDecayShower( children[0],maxScales,minmass,type,fc[1]);
	}
	else {
	  if(fb.hard && !branch->children()[0]->children().empty() )
	    truncatedTimeLikeShower(children[0],branch->children()[0],type,fc[1],false);
	  // normal shower
	  else
	    timeLikeShower(children[0],type,fc[1],false);
	}
      }
      updateHistory(children[1]);
      // branching has happened
      break;
    }
    // H7 default
    else if(_reconOpt==1) {
      // update the history if needed
      currentTree()->updateInitialStateShowerProduct(progenitor(),children[0]);
      currentTree()->addInitialStateBranching(particle,children[0],children[1]);
      // shower the first  particle
      if(fc[0].kinematics) {
	if(children[0]->id()==particle->id()) {
	  if(!fb.hard)
	    truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					   branch,type,fc[0]);
	  else if(fb.hard && ! branch->children()[0]->children().empty() )
	    truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					   branch->children()[0],type,fc[0]);
	  else
	    spaceLikeDecayShower( children[0],maxScales,minmass,type,fc[0]);
	}
	else {
	  if(fb.hard && !branch->children()[0]->children().empty() )
	    truncatedTimeLikeShower(children[0],branch->children()[0],type,fc[0],false);
	  // normal shower
	  else
	    timeLikeShower(children[0],type,fc[0],false);
	}
      }
      // shower the second  particle
      if(fc[1].kinematics) {
	if(children[0]->id()==particle->id()) {
	  if(!fb.hard)
	    truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					   branch,type,fc[1]);
	  else if(fb.hard && ! branch->children()[0]->children().empty() )
	    truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					   branch->children()[0],type,fc[1]);
	  else
	    spaceLikeDecayShower( children[0],maxScales,minmass,type,fc[1]);
	}
	else {
	  if(fb.hard && !branch->children()[0]->children().empty() )
	    truncatedTimeLikeShower(children[0],branch->children()[0],type,fc[1],false);
	  // normal shower
	  else
	    timeLikeShower(children[0],type,fc[1],false);
	}
      }
      // clean up the vetoed emission
      if(particle->virtualMass()==ZERO) {
   	particle->showerKinematics(ShoKinPtr());
    	for(unsigned int ix=0;ix<children.size();++ix)
   	  particle->abandonChild(children[ix]);
    	children.clear();
   	particle->vetoEmission(fb.type,fb.kinematics->scale());
    	// generate the new emission
   	fb = selectSpaceLikeDecayBranching(particle,maxScales,minmass,type,branch);
   	// must be at least hard emission
  	assert(fb.kinematics);
   	setupChildren = true;
 	continue;
      }
      else {
	updateHistory(children[1]);
  	break;
      }
    }
    else if(_reconOpt>=2) {
      // cut-off masses for the branching
      const vector<Energy> & virtualMasses = fb.sudakov->virtualMasses(fb.ids);
      // compute the masses of the children
      Energy masses[3];
      // space-like children
      masses[1] = children[0]->virtualMass();
      // time-like child
      if(fc[1].kinematics) {
	const vector<Energy> & vm = fc[1].sudakov->virtualMasses(fc[1].ids);
	Energy2 q2 = 
	  fc[1].kinematics->z()*(1.-fc[1].kinematics->z())*sqr(fc[1].kinematics->scale());
	if(fc[1].ids[0]->id()!=ParticleID::g) q2 += sqr(vm[0]);
	masses[2] = sqrt(q2);
      }
      else {
	masses[2] = virtualMasses[2];
      } 
      masses[0]=particle->virtualMass();
      double z = fb.kinematics->z();
      Energy2 pt2 = (1.-z)*(z*sqr(masses[0])-sqr(masses[1])-z/(1.-z)*sqr(masses[2]));
      if(pt2>=ZERO) {
  	break;
      }
      else {
  	// reset the scales for the children
  	for(unsigned int ix=0;ix<2;++ix) {
  	  if(fc[ix].kinematics)
  	    children[ix]->vetoEmission(fc[ix].type,fc[ix].kinematics->scale());
  	  else {
	    if(ix==0) 
	      children[ix]->vetoEmission(ShowerPartnerType::QCDColourLine,Constants::MaxEnergy);
	    else
	      children[ix]->vetoEmission(ShowerPartnerType::QCDColourLine,ZERO);
	  }
   	} 
	children[0]->virtualMass(_progenitor->progenitor()->mass());
	children[1]->virtualMass(ZERO);
      }
    }
  };
  if(_reconOpt>=2) {
    // update the history if needed
    currentTree()->updateInitialStateShowerProduct(progenitor(),children[0]);
    currentTree()->addInitialStateBranching(particle,children[0],children[1]);
    // shower the first  particle
    if(fc[0].kinematics) {
      if(children[0]->id()==particle->id()) {
	if(!fb.hard)
	  truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					 branch,type,fc[0]);
	else if(fb.hard && ! branch->children()[0]->children().empty() )
	  truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					 branch->children()[0],type,fc[0]);
	else
	  spaceLikeDecayShower( children[0],maxScales,minmass,type,fc[0]);
      }
      else {
	if(fb.hard && !branch->children()[0]->children().empty() )
	  truncatedTimeLikeShower(children[0],branch->children()[0],type,fc[0],false);
	// normal shower
	else
	  timeLikeShower(children[0],type,fc[0],false);
      }
    }
    // shower the second  particle
    if(fc[1].kinematics) {
      if(children[0]->id()==particle->id()) {
	if(!fb.hard)
	  truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					 branch,type,fc[1]);
	else if(fb.hard && ! branch->children()[0]->children().empty() )
	  truncatedSpaceLikeDecayShower( children[0],maxScales,minmass,
					 branch->children()[0],type,fc[1]);
	else
	  spaceLikeDecayShower( children[0],maxScales,minmass,type,fc[1]);
      }
      else {
	if(fb.hard && !branch->children()[0]->children().empty() )
	  truncatedTimeLikeShower(children[0],branch->children()[0],type,fc[1],false);
	// normal shower
	else
	  timeLikeShower(children[0],type,fc[1],false);
      }
    }
    updateHistory(children[1]);
  }
  return true;
}

void QTildeShowerHandler::connectTrees(ShowerTreePtr showerTree, 
				       HardTreePtr hardTree, bool hard ) {
  ShowerParticleVector particles;
  // find the Sudakovs
  for(set<HardBranchingPtr>::iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    // Sudakovs for ISR
    if((**cit).parent()&&(**cit).status()==HardBranching::Incoming) {
      ++_nis;
      array<long,3> br;
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
	IdList ids = cjt->second.particles;
	if(ids[0]->id()==br[0]&&ids[1]->id()==br[1]&&ids[2]->id()==br[2]) {
	  sudakov=cjt->second.sudakov;
	  break;
	}
      }
      if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				     << "QTildeShowerHandler::connectTrees() for ISR" 
				     << Exception::runerror;
      (**cit).parent()->sudakov(sudakov);
    }
    // Sudakovs for FSR
    else if(!(**cit).children().empty()) {
      ++_nfs;
      array<long,3> br;
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
	IdList ids = cjt->second.particles;
	if(ids[0]->id()==br[0]&&ids[1]->id()==br[1]&&ids[2]->id()==br[2]) {
	  sudakov=cjt->second.sudakov;
	  break;
	}
      }
      if(!sudakov) {
	throw Exception() << "Can't find Sudakov for the hard emission in "
			  << "QTildeShowerHandler::connectTrees()" 
			  << Exception::runerror;
      }
      (**cit).sudakov(sudakov);
    }
  }
  // calculate the evolution scale
  for(set<HardBranchingPtr>::iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    particles.push_back((*cit)->branchingParticle());
  }
  showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,!hard,interaction_,true);
  hardTree->partnersSet(true);
  // inverse reconstruction
  if(hard) {
    showerModel()->kinematicsReconstructor()->
      deconstructHardJets(hardTree,interaction_);
  }
  else
    showerModel()->kinematicsReconstructor()->
      deconstructDecayJets(hardTree,interaction_);
  // now reset the momenta of the showering particles
  vector<ShowerProgenitorPtr> particlesToShower=showerTree->extractProgenitors();
  // match them
  map<ShowerProgenitorPtr,HardBranchingPtr> partners;
  for(set<HardBranchingPtr>::const_iterator bit=hardTree->branchings().begin();
      bit!=hardTree->branchings().end();++bit) {
    Energy2 dmin( 1e30*GeV2 );
    ShowerProgenitorPtr partner;
    for(vector<ShowerProgenitorPtr>::const_iterator pit=particlesToShower.begin();
	pit!=particlesToShower.end();++pit) {
      if(partners.find(*pit)!=partners.end()) continue;
      if( (**bit).branchingParticle()->id() !=  (**pit).progenitor()->id() ) continue;
      if( (**bit).branchingParticle()->isFinalState() !=
	  (**pit).progenitor()->isFinalState() ) continue;
      if( (**pit).progenitor()->isFinalState() ) {
	Energy2 dtest =
	  sqr( (**pit).progenitor()->momentum().x() - (**bit).showerMomentum().x() ) +
	  sqr( (**pit).progenitor()->momentum().y() - (**bit).showerMomentum().y() ) +
	  sqr( (**pit).progenitor()->momentum().z() - (**bit).showerMomentum().z() ) +
	  sqr( (**pit).progenitor()->momentum().t() - (**bit).showerMomentum().t() );
	// add mass difference for identical particles (e.g. Z0 Z0 production)
	dtest += 1e10*sqr((**pit).progenitor()->momentum().m()-(**bit).showerMomentum().m());
	if( dtest < dmin ) {
	  partner = *pit;
	  dmin = dtest;
	}
      }
      else {
	// ensure directions are right
	if((**pit).progenitor()->momentum().z()/(**bit).showerMomentum().z()>ZERO) {
	  partner = *pit;
	  break;
	}
      }
    }
    if(!partner) throw Exception() << "Failed to match shower and hard trees in QTildeShowerHandler::hardestEmission"
				   << Exception::eventerror;
    partners[partner] = *bit;
  }
  for(vector<ShowerProgenitorPtr>::const_iterator pit=particlesToShower.begin();
      pit!=particlesToShower.end();++pit) {
    HardBranchingPtr partner = partners[*pit];
    if((**pit).progenitor()->dataPtr()->stable()) {
      (**pit).progenitor()->set5Momentum(partner->showerMomentum());
      (**pit).copy()->set5Momentum(partner->showerMomentum());
    }
    else {
      Lorentz5Momentum oldMomentum = (**pit).progenitor()->momentum();
      Lorentz5Momentum newMomentum = partner->showerMomentum();
      LorentzRotation boost( oldMomentum.findBoostToCM(),oldMomentum.e()/oldMomentum.mass());
      (**pit).progenitor()->transform(boost);
      (**pit).copy()      ->transform(boost);
      boost = LorentzRotation(-newMomentum.findBoostToCM(),newMomentum.e()/newMomentum.mass());
      (**pit).progenitor()->transform(boost);
      (**pit).copy()      ->transform(boost);
    }
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
    decayTree->transform(boost,true);
    boost = LorentzRotation(-newMomentum.findBoostToCM(),newMomentum.e()/newMomentum.mass());
    decayTree->transform(boost,true);
  }
}

void QTildeShowerHandler::doShowering(bool hard,XCPtr xcomb) {
  // zero number of emissions
  _nis = _nfs = 0;
  // if MC@NLO H event and limited emissions
  // indicate both final and initial state emission
  if ( currentTree()->isMCatNLOHEvent() && _limitEmissions != 0 ) {
    _nis = _nfs = 1;
  }
  // extract particles to shower
  vector<ShowerProgenitorPtr> particlesToShower(setupShower(hard));
  // check if we should shower
  bool colCharge = false;
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if(particlesToShower[ix]->progenitor()->dataPtr()->coloured() ||
       particlesToShower[ix]->progenitor()->dataPtr()->charged()) {
      colCharge = true;
      break;
    }
  }
  if(!colCharge) {
    _currenttree->hasShowered(true);
    return;
  }
  // setup the maximum scales for the shower
  if (restrictPhasespace()) setupMaximumScales(particlesToShower,xcomb);
  // set the hard scales for the profiles
  setupHardScales(particlesToShower,xcomb);
  // specific stuff for hard processes and decays
  Energy minmass(ZERO), mIn(ZERO);
  // hard process generate the intrinsic p_T once and for all
  if(hard) {
    generateIntrinsicpT(particlesToShower);
  }
  // decay compute the minimum mass of the final-state
  else {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      if(particlesToShower[ix]->progenitor()->isFinalState()) {
        if(particlesToShower[ix]->progenitor()->dataPtr()->stable()){
          auto dm=  ShowerHandler::currentHandler()->retConstituentMasses()?
          particlesToShower[ix]->progenitor()->dataPtr()->constituentMass():
          particlesToShower[ix]->progenitor()->dataPtr()->mass();
          minmass += dm;
        }else
	  minmass += particlesToShower[ix]->progenitor()->mass();
      }
      else {
	mIn = particlesToShower[ix]->progenitor()->mass();
      }
    }
    // throw exception if decay can't happen
    if ( minmass > mIn ) {
      throw Exception() << "QTildeShowerHandler.cc: Mass of decaying particle is "
			<< "below constituent masses of decay products."
			<< Exception::eventerror;
    }
  }
  // setup for reweighted
  bool reWeighting = _reWeight && hard && ShowerHandler::currentHandler()->firstInteraction();
  double eventWeight=0.;
  unsigned int nTryReWeight(0);
  // create random particle vector (only need to do once)
  vector<ShowerProgenitorPtr> tmp;
  unsigned int nColouredIncoming = 0;
  while(particlesToShower.size()>0){
    unsigned int xx=UseRandom::irnd(particlesToShower.size());
    tmp.push_back(particlesToShower[xx]);
    particlesToShower.erase(particlesToShower.begin()+xx);
  }
  particlesToShower=tmp;
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if(!particlesToShower[ix]->progenitor()->isFinalState() &&
       particlesToShower[ix]->progenitor()->coloured()) ++nColouredIncoming;
  }
  bool switchRecon = hard && nColouredIncoming !=1;
  // main shower loop
  unsigned int ntry(0);
  bool reconstructed = false;
  do {
    // clear results of last attempt if needed
    if(ntry!=0) {
      currentTree()->clear();
      setEvolutionPartners(hard,interaction_,true);
      _nis = _nfs = 0;
      // if MC@NLO H event and limited emissions
      // indicate both final and initial state emission
      if ( currentTree()->isMCatNLOHEvent() && _limitEmissions != 0 ) {
	_nis = _nfs = 1;
      }
      for(unsigned int ix=0; ix<particlesToShower.size();++ix) {
	SpinPtr spin = particlesToShower[ix]->progenitor()->spinInfo();
	if(spin && spin->decayVertex() &&
	   dynamic_ptr_cast<tcSVertexPtr>(spin->decayVertex())) {
	  spin->decayVertex(VertexPtr());
	}
      }
    }
    // loop over particles
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      // extract the progenitor
      progenitor(particlesToShower[ix]);
      // final-state radiation
      if(progenitor()->progenitor()->isFinalState()) {
	if(!doFSR()) continue;
	// perform shower
	progenitor()->hasEmitted(startTimeLikeShower(interaction_));
      }
      // initial-state radiation
      else {
	if(!doISR()) continue;
	// hard process
	if(hard) {
	  // get the PDF
	  setBeamParticle(_progenitor->beam());
	  if(!beamParticle()) {
	    throw Exception() << "Incorrect type of beam particle in "
			      << "QTildeShowerHandler::doShowering(). "
			      << "This should not happen for conventional choices but may happen if you have used a"
			      << " non-default choice and have not changed the create ParticleData line in the input files"
			      << " for this particle to create BeamParticleData."
			      << Exception::runerror;
	  }
	  // perform the shower
	  // set the beam particle
	  tPPtr beamparticle=progenitor()->original();
	  if(!beamparticle->parents().empty()) 
	    beamparticle=beamparticle->parents()[0];
	  // generate the shower
	  progenitor()->hasEmitted(startSpaceLikeShower(beamparticle,
							interaction_));
	}
	// decay
	else {
	  // skip colour and electrically neutral particles
	  if(!progenitor()->progenitor()->dataPtr()->coloured() &&
	     !progenitor()->progenitor()->dataPtr()->charged()) {
	    progenitor()->hasEmitted(false);
	    continue;
	  }
	  // perform shower
	  // set the scales correctly. The current scale is the maximum scale for
	  // emission not the starting scale
	  ShowerParticle::EvolutionScales maxScales(progenitor()->progenitor()->scales());
	  progenitor()->progenitor()->scales() = ShowerParticle::EvolutionScales();
	  if(progenitor()->progenitor()->dataPtr()->charged()) {
	    progenitor()->progenitor()->scales().QED      = progenitor()->progenitor()->mass();
	    progenitor()->progenitor()->scales().QED_noAO = progenitor()->progenitor()->mass();
	  }
	  if(progenitor()->progenitor()->hasColour()) {
	    progenitor()->progenitor()->scales().QCD_c       = progenitor()->progenitor()->mass();
	    progenitor()->progenitor()->scales().QCD_c_noAO  = progenitor()->progenitor()->mass();
	  }
	  if(progenitor()->progenitor()->hasAntiColour()) {
	    progenitor()->progenitor()->scales().QCD_ac      = progenitor()->progenitor()->mass();
	    progenitor()->progenitor()->scales().QCD_ac_noAO = progenitor()->progenitor()->mass();
	  }
	  // perform the shower
	  progenitor()->hasEmitted(startSpaceLikeDecayShower(maxScales,minmass,
							     interaction_));
	}
      }
    }
    // do the kinematic reconstruction, checking if it worked
    reconstructed = hard ?
      showerModel()->kinematicsReconstructor()->
      reconstructHardJets (currentTree(),intrinsicpT(),interaction_,
			   switchRecon && ntry>maximumTries()/2) :
      showerModel()->kinematicsReconstructor()->
      reconstructDecayJets(currentTree(),interaction_);
    if(!reconstructed) continue;
    // apply vetos on the full shower
    for(vector<FullShowerVetoPtr>::const_iterator it=_fullShowerVetoes.begin();
	it!=_fullShowerVetoes.end();++it) {
      int veto = (**it).applyVeto(currentTree());
      if(veto<0) continue;
      // veto the shower
      if(veto==0) {
	reconstructed = false;
	break;
      }
      // veto the shower and reweight
      else if(veto==1) {
	reconstructed = false;
	break;
      }
      // veto the event
      else if(veto==2) {
	throw Veto();
      }
    }
    if(reWeighting) {
      if(reconstructed) eventWeight += 1.;
      reconstructed=false;
      ++nTryReWeight;
      if(nTryReWeight==_nReWeight) {
	reWeighting = false;
	if(eventWeight==0.) throw Veto();
      }
    }
  }
  while(!reconstructed&&maximumTries()>++ntry);
  // check if failed to generate the shower
  if(ntry==maximumTries()) {
    if(hard)
      throw ShowerHandler::ShowerTriesVeto(ntry);
    else
      throw Exception() << "Failed to generate the shower after "
			<< ntry << " attempts in QTildeShowerHandler::showerDecay()"
			<< Exception::eventerror;
  }
  // handle the weights and apply any reweighting required
  if(nTryReWeight>0) {
    tStdEHPtr seh = dynamic_ptr_cast<tStdEHPtr>(generator()->currentEventHandler());
    static bool first = true;
    if(seh) {
      seh->reweight(eventWeight/double(nTryReWeight));
    }
    else if(first) {
      generator()->log() << "Reweighting the shower only works with internal Herwig7 processes"
    			 << "Presumably you are showering Les Houches Events. These will not be"
    			 << "reweighted\n";
      first = false;
    }
  }
  // tree has now showered
  _currenttree->hasShowered(true);
  hardTree(HardTreePtr());
}

void QTildeShowerHandler:: convertHardTree(bool hard,ShowerInteraction type) {
  map<ColinePtr,ColinePtr> cmap;
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=currentTree()->incomingLines().begin();cit!=currentTree()->incomingLines().end();++cit) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      mit = hardTree()->particles().find(cit->first->progenitor());
    // put the colour lines in the map
    ShowerParticlePtr oldParticle = cit->first->progenitor();
    ShowerParticlePtr newParticle = mit->second->branchingParticle();
    ColinePtr cLine = oldParticle->    colourLine();
    ColinePtr aLine = oldParticle->antiColourLine();
    if(newParticle->colourLine() &&
       cmap.find(newParticle->    colourLine())==cmap.end())
      cmap[newParticle->    colourLine()] = cLine;
    if(newParticle->antiColourLine() &&
       cmap.find(newParticle->antiColourLine())==cmap.end())
      cmap[newParticle->antiColourLine()] = aLine;
    // check whether or not particle emits
    bool emission = mit->second->parent();
    if(emission) {
      if(newParticle->colourLine()) {
	ColinePtr ctemp = newParticle->    colourLine();
	ctemp->removeColoured(newParticle);
      }
      if(newParticle->antiColourLine()) {
	ColinePtr ctemp = newParticle->antiColourLine();
	ctemp->removeAntiColoured(newParticle);
      }
      newParticle = mit->second->parent()->branchingParticle();
    }
    // get the new colour lines
    ColinePtr newCLine,newALine;
    // sort out colour lines
    if(newParticle->colourLine()) {
      ColinePtr ctemp = newParticle->    colourLine();
      ctemp->removeColoured(newParticle);
      if(cmap.find(ctemp)!=cmap.end()) {
	newCLine = cmap[ctemp];
      }
      else {
	newCLine = new_ptr(ColourLine());
	cmap[ctemp] = newCLine;
      }
    }
    // and anticolour lines
    if(newParticle->antiColourLine()) {
      ColinePtr ctemp = newParticle->antiColourLine();
      ctemp->removeAntiColoured(newParticle);
      if(cmap.find(ctemp)!=cmap.end()) {
	newALine = cmap[ctemp];
      }
      else {
	newALine = new_ptr(ColourLine());
	cmap[ctemp] = newALine;
      }
    }
    // remove colour lines from old particle
    if(aLine) {
      aLine->removeAntiColoured(cit->first->copy());
      aLine->removeAntiColoured(cit->first->progenitor());
    }
    if(cLine) {
      cLine->removeColoured(cit->first->copy());
      cLine->removeColoured(cit->first->progenitor());
    }
    // add particle to colour lines
    if(newCLine) newCLine->addColoured    (newParticle);
    if(newALine) newALine->addAntiColoured(newParticle);
    // insert new particles
    cit->first->copy(newParticle);
    ShowerParticlePtr sp(new_ptr(ShowerParticle(*newParticle,1,false)));
    cit->first->progenitor(sp);
    currentTree()->incomingLines()[cit->first]=sp;
    cit->first->perturbative(!emission);
    // and the emitted particle if needed
    if(emission) {
      ShowerParticlePtr newOut = mit->second->parent()->children()[1]->branchingParticle();
      if(newOut->colourLine()) {
	ColinePtr ctemp = newOut->    colourLine();
	ctemp->removeColoured(newOut);
	assert(cmap.find(ctemp)!=cmap.end());
	cmap[ctemp]->addColoured    (newOut);
      }
      if(newOut->antiColourLine()) {
	ColinePtr ctemp = newOut->antiColourLine();
	ctemp->removeAntiColoured(newOut);
	assert(cmap.find(ctemp)!=cmap.end());
	cmap[ctemp]->addAntiColoured(newOut);
      }
      ShowerParticlePtr sout=new_ptr(ShowerParticle(*newOut,1,true));
      ShowerProgenitorPtr out=new_ptr(ShowerProgenitor(cit->first->original(),newOut,sout));
      out->perturbative(false);
      currentTree()->outgoingLines().insert(make_pair(out,sout));
    }
    if(hard) {
      // sort out the value of x
      if(mit->second->beam()->momentum().z()>ZERO) {
	sp->x(newParticle->momentum(). plus()/mit->second->beam()->momentum(). plus());
      }
      else {
	sp->x(newParticle->momentum().minus()/mit->second->beam()->momentum().minus());
      }
    }
  }
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=currentTree()->outgoingLines().begin();cit!=currentTree()->outgoingLines().end();++cit) {
    map<tShowerTreePtr,pair<tShowerProgenitorPtr,
			    tShowerParticlePtr> >::const_iterator tit;
    for(tit  = currentTree()->treelinks().begin();
	tit != currentTree()->treelinks().end();++tit) {
      if(tit->second.first && tit->second.second==cit->first->progenitor())
	break;
    }
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      mit = hardTree()->particles().find(cit->first->progenitor());
    if(mit==hardTree()->particles().end()) continue;
    // put the colour lines in the map
    ShowerParticlePtr oldParticle = cit->first->progenitor();
    ShowerParticlePtr newParticle = mit->second->branchingParticle();
    ShowerParticlePtr newOut;
    ColinePtr cLine = oldParticle->    colourLine();
    ColinePtr aLine = oldParticle->antiColourLine();
    if(newParticle->colourLine() &&
       cmap.find(newParticle->    colourLine())==cmap.end())
      cmap[newParticle->    colourLine()] = cLine;
    if(newParticle->antiColourLine() &&
       cmap.find(newParticle->antiColourLine())==cmap.end())
      cmap[newParticle->antiColourLine()] = aLine;
    // check whether or not particle emits
    bool emission = !mit->second->children().empty();
    if(emission) {
      if(newParticle->colourLine()) {
	ColinePtr ctemp = newParticle->    colourLine();
	ctemp->removeColoured(newParticle);
      }
      if(newParticle->antiColourLine()) {
	ColinePtr ctemp = newParticle->antiColourLine();
	ctemp->removeAntiColoured(newParticle);
      }
      newParticle = mit->second->children()[0]->branchingParticle();
      newOut      = mit->second->children()[1]->branchingParticle();
      if(newParticle->id()!=oldParticle->id()&&newParticle->id()==newOut->id())
	swap(newParticle,newOut);
    }
    // get the new colour lines
    ColinePtr newCLine,newALine;
    // sort out colour lines
    if(newParticle->colourLine()) {
      ColinePtr ctemp = newParticle->    colourLine();
      ctemp->removeColoured(newParticle);
      if(cmap.find(ctemp)!=cmap.end()) {
	newCLine = cmap[ctemp];
      }
      else {
	newCLine = new_ptr(ColourLine());
	cmap[ctemp] = newCLine;
      }
    }
    // and anticolour lines
    if(newParticle->antiColourLine()) {
      ColinePtr ctemp = newParticle->antiColourLine();
      ctemp->removeAntiColoured(newParticle);
      if(cmap.find(ctemp)!=cmap.end()) {
	newALine = cmap[ctemp];
      }
      else {
	newALine = new_ptr(ColourLine());
	cmap[ctemp] = newALine;
      }
    }
    // remove colour lines from old particle
    if(aLine) {
      aLine->removeAntiColoured(cit->first->copy());
      aLine->removeAntiColoured(cit->first->progenitor());
    }
    if(cLine) {
      cLine->removeColoured(cit->first->copy());
      cLine->removeColoured(cit->first->progenitor());
    }
    // special for unstable particles
    if(newParticle->id()==oldParticle->id() &&
       (tit!=currentTree()->treelinks().end()||!oldParticle->dataPtr()->stable())) {
      Lorentz5Momentum oldMomentum = oldParticle->momentum();
      Lorentz5Momentum newMomentum = newParticle->momentum();
      LorentzRotation boost( oldMomentum.findBoostToCM(),oldMomentum.e()/oldMomentum.mass());
      if(tit!=currentTree()->treelinks().end()) tit->first->transform(boost,false);
      oldParticle->transform(boost);
      boost = LorentzRotation(-newMomentum.findBoostToCM(),newMomentum.e()/newMomentum.mass());
      oldParticle->transform(boost);
      if(tit!=currentTree()->treelinks().end()) tit->first->transform(boost,false);
      newParticle=oldParticle;
    }
    // add particle to colour lines
    if(newCLine) newCLine->addColoured    (newParticle);
    if(newALine) newALine->addAntiColoured(newParticle);
    // insert new particles
    cit->first->copy(newParticle);
    ShowerParticlePtr sp(new_ptr(ShowerParticle(*newParticle,1,true)));
    cit->first->progenitor(sp);
    currentTree()->outgoingLines()[cit->first]=sp;
    cit->first->perturbative(!emission);
    // and the emitted particle if needed
    if(emission) {
      if(newOut->colourLine()) {
	ColinePtr ctemp = newOut->    colourLine();
	ctemp->removeColoured(newOut);
	assert(cmap.find(ctemp)!=cmap.end());
	cmap[ctemp]->addColoured    (newOut);
      }
      if(newOut->antiColourLine()) {
	ColinePtr ctemp = newOut->antiColourLine();
	ctemp->removeAntiColoured(newOut);
	assert(cmap.find(ctemp)!=cmap.end());
	cmap[ctemp]->addAntiColoured(newOut);
      }
      ShowerParticlePtr sout=new_ptr(ShowerParticle(*newOut,1,true));
      ShowerProgenitorPtr out=new_ptr(ShowerProgenitor(cit->first->original(),newOut,sout));
      out->perturbative(false);
      currentTree()->outgoingLines().insert(make_pair(out,sout));
    }
    // update any decay products
    if(tit!=currentTree()->treelinks().end())
      currentTree()->updateLink(tit->first,make_pair(cit->first,sp));
  }
  // reset the tree
  currentTree()->resetShowerProducts();
  // reextract the particles and set the colour partners
  vector<ShowerParticlePtr> particles = 
    currentTree()->extractProgenitorParticles();
  // clear the partners
  for(unsigned int ix=0;ix<particles.size();++ix) {
    particles[ix]->partner(ShowerParticlePtr());
    particles[ix]->clearPartners();
  }
  // clear the tree
  hardTree(HardTreePtr());
  // Set the initial evolution scales
  showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,!hard,type,!_hardtree);
}

Branching QTildeShowerHandler::selectTimeLikeBranching(tShowerParticlePtr particle,
						       ShowerInteraction type,
						       HardBranchingPtr branch) {
  Branching fb;
  unsigned int iout=0;
  while (true) {
    // break if doing truncated shower and no truncated shower needed
    if(branch && (!isTruncatedShowerON()||hardOnly())) break;
    fb=_splittingGenerator->chooseForwardBranching(*particle,_finalenhance,type);
    // no emission break
    if(!fb.kinematics) break;
    // special for truncated shower
    if(branch) {
      // check haven't evolved too far
      if(fb.kinematics->scale() < branch->scale()) {
	fb=Branching();
	break;
      }
      // find the truncated line
      iout=0;
      if(fb.ids[1]->id()!=fb.ids[2]->id()) {
	if(fb.ids[1]->id()==particle->id())       iout=1;
	else if (fb.ids[2]->id()==particle->id()) iout=2;
      }
      else if(fb.ids[1]->id()==particle->id()) {
	if(fb.kinematics->z()>0.5) iout=1;
	else                       iout=2;
      }
      // apply the vetos for the truncated shower
      // no flavour changing branchings
      if(iout==0) {
	particle->vetoEmission(fb.type,fb.kinematics->scale());
	continue;
      }
      double zsplit = iout==1 ? fb.kinematics->z() : 1-fb.kinematics->z();
      // only if same interaction for forced branching
      ShowerInteraction type2 = convertInteraction(fb.type);
      // and evolution
      if(type2==branch->sudakov()->interactionType()) {
	if(zsplit < 0.5 || // hardest line veto
	   fb.kinematics->scale()*zsplit < branch->scale() ) { // angular ordering veto
	  particle->vetoEmission(fb.type,fb.kinematics->scale());
	  continue;
	}
      }
      // pt veto
      if(fb.kinematics->pT() > progenitor()->maximumpT(type2)) {
	particle->vetoEmission(fb.type,fb.kinematics->scale());
	continue;
      }
    }
    // standard vetos for all emissions
    if(timeLikeVetoed(fb,particle)) {
      particle->vetoEmission(fb.type,fb.kinematics->scale());
      if(particle->spinInfo()) particle->spinInfo()->decayVertex(VertexPtr());
      continue;
    }
    // special for already decayed particles
    // don't allow flavour changing branchings
    bool vetoDecay = false;
    for(map<tShowerTreePtr,pair<tShowerProgenitorPtr,
	  tShowerParticlePtr> >::const_iterator tit  = currentTree()->treelinks().begin();
	tit != currentTree()->treelinks().end();++tit) {
      if(tit->second.first == progenitor()) {
	map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  it = currentTree()->outgoingLines().find(progenitor());
	if(it!=currentTree()->outgoingLines().end() && particle == it->second &&
	   fb.ids[0]!=fb.ids[1] && fb.ids[1]!=fb.ids[2]) {
	  vetoDecay = true;
	  break;
	}
      }
    }
    if(vetoDecay) {
      particle->vetoEmission(fb.type,fb.kinematics->scale());
      if(particle->spinInfo()) particle->spinInfo()->decayVertex(VertexPtr());
      continue;
    }
    break;
  }
  // normal case
  if(!branch) {
    if(fb.kinematics) fb.hard = false;
    return fb;
  }
  // truncated emission
  if(fb.kinematics) {
    fb.hard = false;
    fb.iout = iout;
    return fb;
  }
  // otherwise need to return the hard emission
  // construct the kinematics for the hard emission
  ShoKinPtr showerKin=
    branch->sudakov()->createFinalStateBranching(branch->scale(),
	                                         branch->children()[0]->z(),
	                                         branch->phi(),
						 branch->children()[0]->pT());
  IdList idlist(3);
  idlist[0] = particle->dataPtr();
  idlist[1] = branch->children()[0]->branchingParticle()->dataPtr();
  idlist[2] = branch->children()[1]->branchingParticle()->dataPtr();
  fb = Branching( showerKin, idlist, branch->sudakov(),branch->type() );
  fb.hard = true;
  fb.iout=0;
  // return it
  return fb;
}

Branching QTildeShowerHandler::selectSpaceLikeDecayBranching(tShowerParticlePtr particle,
						 const ShowerParticle::EvolutionScales & maxScales,
						 Energy minmass,ShowerInteraction type,
						 HardBranchingPtr branch) {
  Branching fb;
  unsigned int iout=0;
  while (true) {
    // break if doing truncated shower and no truncated shower needed
    if(branch && (!isTruncatedShowerON()||hardOnly())) break;
    // select branching
    fb=_splittingGenerator->chooseDecayBranching(*particle,maxScales,minmass,
						 _initialenhance,type);
    // return if no radiation
    if(!fb.kinematics) break;
    // special for truncated shower
    if(branch) {
      // check haven't evolved too far
      if(fb.kinematics->scale() < branch->scale()) {
	fb=Branching();
	break;
      }
      // find the truncated line
      iout=0;
      if(fb.ids[1]->id()!=fb.ids[2]->id()) {
	if(fb.ids[1]->id()==particle->id())       iout=1;
	else if (fb.ids[2]->id()==particle->id()) iout=2;
      }
      else if(fb.ids[1]->id()==particle->id()) {
	if(fb.kinematics->z()>0.5) iout=1;
	else                       iout=2;
      }
      // apply the vetos for the truncated shower
      // no flavour changing branchings
      if(iout==0) {
	particle->vetoEmission(fb.type,fb.kinematics->scale());
	continue;
      }
      ShowerInteraction type2 = convertInteraction(fb.type);
      double zsplit = iout==1 ? fb.kinematics->z() : 1-fb.kinematics->z();
      if(type2==branch->sudakov()->interactionType()) {
	if(zsplit < 0.5 || // hardest line veto
	   fb.kinematics->scale()*zsplit < branch->scale() ) { // angular ordering veto
	  particle->vetoEmission(fb.type,fb.kinematics->scale());
	  continue;
	}
      }
      // pt veto
      if(fb.kinematics->pT() > progenitor()->maximumpT(type2)) {
	particle->vetoEmission(fb.type,fb.kinematics->scale());
	continue;
      }
    }
    // if not vetoed break
    if(spaceLikeDecayVetoed(fb,particle)) {
      // otherwise reset scale and continue
      particle->vetoEmission(fb.type,fb.kinematics->scale());
      continue;
    }
    break;
  }
  // normal case
  if(!branch) {
    if(fb.kinematics) fb.hard = false;
    return fb;
  }
  // truncated emission
  if(fb.kinematics) {
    fb.hard = false;
    fb.iout = iout;
    return fb;
  }
  // otherwise need to return the hard emission
  // construct the kinematics for the hard emission
  ShoKinPtr showerKin=
    branch->sudakov()->createDecayBranching(branch->scale(),
					    branch->children()[0]->z(),
					    branch->phi(),
					    branch->children()[0]->pT());
   IdList idlist(3);
   idlist[0] = particle->dataPtr();
   idlist[1] = branch->children()[0]->branchingParticle()->dataPtr();
   idlist[2] = branch->children()[1]->branchingParticle()->dataPtr();
   // create the branching
   fb = Branching( showerKin, idlist, branch->sudakov(),ShowerPartnerType::QCDColourLine  );
   fb.hard=true;
   fb.iout=0;
   // return it
  return fb;
}

void QTildeShowerHandler::checkFlags() {
  string error = "Inconsistent hard emission set-up in QTildeShowerHandler::showerHardProcess(). "; 
  if ( ( currentTree()->isMCatNLOSEvent() || currentTree()->isMCatNLOHEvent() ) ) {
    if (_hardEmission ==2 )
      throw Exception() << error
			<< "Cannot generate POWHEG matching with MC@NLO shower "
			<< "approximation.  Add 'set QTildeShowerHandler:HardEmission 0' to input file."
			<< Exception::runerror;
    if ( canHandleMatchboxTrunc() )
      throw Exception() << error
			<< "Cannot use truncated qtilde shower with MC@NLO shower "
			<< "approximation.  Set LHCGenerator:EventHandler"
			<< ":CascadeHandler to '/Herwig/Shower/ShowerHandler' or "
			<< "'/Herwig/Shower/Dipole/DipoleShowerHandler'."
			<< Exception::runerror;
  }
  else if ( ((currentTree()->isPowhegSEvent() || currentTree()->isPowhegHEvent()) ) &&
	    _hardEmission != 2){
    if ( canHandleMatchboxTrunc())
      throw Exception() << error
			<< "Unmatched events requested for POWHEG shower "
			<< "approximation.  Set QTildeShowerHandler:HardEmission to "
			<< "'POWHEG'."
			<< Exception::runerror;
    else if (_hardEmissionWarn) {
      _hardEmissionWarn = false;
      _hardEmission=2;
      throw Exception() << error
			<< "Unmatched events requested for POWHEG shower "
			<< "approximation. Changing QTildeShowerHandler:HardEmission from "
			<< _hardEmission << " to 2"
			<< Exception::warning;
    }
  }

  if ( currentTree()->isPowhegSEvent() || currentTree()->isPowhegHEvent()) {
    if (currentTree()->showerApproximation()->needsTruncatedShower() &&
	!canHandleMatchboxTrunc() )
      throw Exception() << error
			<< "Current shower handler cannot generate truncated shower.  "
			<< "Set Generator:EventHandler:CascadeHandler to "
			<< "'/Herwig/Shower/PowhegShowerHandler'."
			<< Exception::runerror;
  }
  else if ( currentTree()->truncatedShower() && _missingTruncWarn) {
    _missingTruncWarn=false;
    throw Exception() << "Warning: POWHEG shower approximation used without "
		      << "truncated shower.  Set Generator:EventHandler:"
		      << "CascadeHandler to '/Herwig/Shower/PowhegShowerHandler' and "
		      << "'MEMatching:TruncatedShower Yes'."
		      << Exception::warning;   
  }
  // else if ( !dipme && _hardEmissionMode > 1 && 
  // 	    firstInteraction())
  //   throw Exception() << error
  // 		      << "POWHEG matching requested for LO events.  Include "
  // 		      << "'set Factory:ShowerApproximation MEMatching' in input file."
  // 		      << Exception::runerror;
}


tPPair QTildeShowerHandler::remakeRemnant(tPPair oldp){
  // get the parton extractor
  PartonExtractor & pex = *lastExtractor();
  // get the new partons
  tPPair newp = make_pair(findFirstParton(oldp.first ), 
			  findFirstParton(oldp.second));
  // if the same do nothing
  if(newp == oldp) return oldp;
  // Creates the new remnants and returns the new PartonBinInstances
  // ATTENTION Broken here for very strange configuration
  PBIPair newbins = pex.newRemnants(oldp, newp, newStep());
  newStep()->addIntermediate(newp.first);
  newStep()->addIntermediate(newp.second);
  // return the new partons
  return newp;
}

PPtr QTildeShowerHandler::findFirstParton(tPPtr seed) const{
  if(seed->parents().empty()) return seed;
  tPPtr parent = seed->parents()[0];
  //if no parent there this is a loose end which will 
  //be connected to the remnant soon.
  if(!parent || parent == incomingBeams().first || 
     parent == incomingBeams().second ) return seed;
  else return findFirstParton(parent);
}

void QTildeShowerHandler::decay(ShowerTreePtr tree, ShowerDecayMap & decay) {
  // must be one incoming particle
  assert(tree->incomingLines().size()==1);
  // apply any transforms
  tree->applyTransforms();
  // if already decayed return
  if(!tree->outgoingLines().empty()) return;
  // now we need to replace the particle with a new copy after the shower
  // find particle after the shower
  map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator
    tit = tree->parent()->treelinks().find(tree);
  assert(tit!=tree->parent()->treelinks().end());
  ShowerParticlePtr newparent=tit->second.second;
  PerturbativeProcessPtr newProcess =  new_ptr(PerturbativeProcess());
  newProcess->incoming().push_back(make_pair(newparent,PerturbativeProcessPtr()));
  DecayProcessMap decayMap;
  ShowerHandler::decay(newProcess,decayMap);
  ShowerTree::constructTrees(tree,decay,newProcess,decayMap);
}

namespace {

  ShowerProgenitorPtr
  findFinalStateLine(ShowerTreePtr tree, long id, Lorentz5Momentum momentum) { 
    map<ShowerProgenitorPtr,tShowerParticlePtr>::iterator partner;
    Energy2 dmin(1e30*GeV2);
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::iterator 
	  cit =tree->outgoingLines().begin(); cit!=tree->outgoingLines().end(); ++cit) {
      if(cit->second->id()!=id) continue;
      Energy2 test = 
	sqr(cit->second->momentum().x()-momentum.x())+
	sqr(cit->second->momentum().y()-momentum.y())+
	sqr(cit->second->momentum().z()-momentum.z())+
	sqr(cit->second->momentum().t()-momentum.t());
      if(test<dmin) {
	dmin    = test;
	partner = cit;
      }
    }
    return partner->first;
  }

  ShowerProgenitorPtr 
  findInitialStateLine(ShowerTreePtr tree, long id, Lorentz5Momentum momentum) { 
    map<ShowerProgenitorPtr,ShowerParticlePtr>::iterator partner;
    Energy2 dmin(1e30*GeV2);
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::iterator 
	  cit =tree->incomingLines().begin(); cit!=tree->incomingLines().end(); ++cit) {
      if(cit->second->id()!=id) continue;
      Energy2 test = 
	sqr(cit->second->momentum().x()-momentum.x())+
	sqr(cit->second->momentum().y()-momentum.y())+
	sqr(cit->second->momentum().z()-momentum.z())+
	sqr(cit->second->momentum().t()-momentum.t());
      if(test<dmin) {
	dmin    = test;
	partner = cit;
      }
    }
    return partner->first;
  }

  void fixSpectatorColours(PPtr newSpect,ShowerProgenitorPtr oldSpect,
			   ColinePair & cline,ColinePair & aline, bool reconnect) {
    cline.first  = oldSpect->progenitor()->colourLine();
    cline.second = newSpect->colourLine();
    aline.first  = oldSpect->progenitor()->antiColourLine();
    aline.second = newSpect->antiColourLine();
    if(!reconnect) return;
    if(cline.first) {
      cline.first ->removeColoured(oldSpect->copy());
      cline.first ->removeColoured(oldSpect->progenitor());
      cline.second->removeColoured(newSpect);
      cline.first ->addColoured(newSpect);
    }
    if(aline.first) {
      aline.first ->removeAntiColoured(oldSpect->copy());
      aline.first ->removeAntiColoured(oldSpect->progenitor());
      aline.second->removeAntiColoured(newSpect);
      aline.first ->addAntiColoured(newSpect);
    }
  }

  void fixInitialStateEmitter(ShowerTreePtr tree, PPtr newEmit,PPtr emitted, ShowerProgenitorPtr emitter,
			      ColinePair cline,ColinePair aline,double x) {
    // sort out the colours
    if(emitted->dataPtr()->iColour()==PDT::Colour8) {
      // emitter
      if(cline.first && cline.first == emitter->progenitor()->antiColourLine() &&
	 cline.second !=newEmit->antiColourLine()) {
	// sort out not radiating line
	ColinePtr col = emitter->progenitor()->colourLine();
	if(col) {
	  col->removeColoured(emitter->copy());
	  col->removeColoured(emitter->progenitor());
	  newEmit->colourLine()->removeColoured(newEmit);
	  col->addColoured(newEmit);
	}
      }
      else if(aline.first && aline.first == emitter->progenitor()->colourLine() &&
	      aline.second !=newEmit->colourLine()) {
	// sort out not radiating line
	ColinePtr anti = emitter->progenitor()->antiColourLine();
	if(anti) {
	  anti->removeAntiColoured(emitter->copy());
	  anti->removeAntiColoured(emitter->progenitor());
	  newEmit->colourLine()->removeAntiColoured(newEmit);
	  anti->addAntiColoured(newEmit);
	}
      }
      else
	assert(false);
      // emitted
      if(cline.first && cline.second==emitted->colourLine()) {
	cline.second->removeColoured(emitted);
	cline.first->addColoured(emitted);
      }
      else if(aline.first && aline.second==emitted->antiColourLine()) {
	aline.second->removeAntiColoured(emitted);
	aline.first->addAntiColoured(emitted);
      }
      else
	assert(false);
    }
    else {
      if(emitter->progenitor()->antiColourLine() ) {
	ColinePtr col = emitter->progenitor()->antiColourLine();
	col->removeAntiColoured(emitter->copy());
	col->removeAntiColoured(emitter->progenitor());
 	if(newEmit->antiColourLine()) {
	  newEmit->antiColourLine()->removeAntiColoured(newEmit);
	  col->addAntiColoured(newEmit);
	}
	else if (emitted->colourLine()) {
	  emitted->colourLine()->removeColoured(emitted);
	  col->addColoured(emitted);
	}
	else
	  assert(false);
      }
      if(emitter->progenitor()->colourLine() ) {
	ColinePtr col = emitter->progenitor()->colourLine();
	col->removeColoured(emitter->copy());
	col->removeColoured(emitter->progenitor());
 	if(newEmit->colourLine()) {
	  newEmit->colourLine()->removeColoured(newEmit);
	  col->addColoured(newEmit);
	}
	else if (emitted->antiColourLine()) {
	  emitted->antiColourLine()->removeAntiColoured(emitted);
	  col->addAntiColoured(emitted);
	}
	else
	  assert(false);
      }
    }
    // update the emitter
    emitter->copy(newEmit);
    ShowerParticlePtr sp = new_ptr(ShowerParticle(*newEmit,1,false));
    sp->x(x);
    emitter->progenitor(sp);
    tree->incomingLines()[emitter]=sp;
    emitter->perturbative(false);
    // add emitted
    sp=new_ptr(ShowerParticle(*emitted,1,true));
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(emitter->original(),emitted,sp));
    gluon->perturbative(false);
    tree->outgoingLines().insert(make_pair(gluon,sp));
  }

  void fixFinalStateEmitter(ShowerTreePtr tree, PPtr newEmit,PPtr emitted, ShowerProgenitorPtr emitter,
			    ColinePair cline,ColinePair aline) {
    map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
    // special case if decayed
    for(tit  = tree->treelinks().begin(); tit != tree->treelinks().end();++tit) {
      if(tit->second.first && tit->second.second==emitter->progenitor())
	break;
    }
    // sort out the colour lines
    if(cline.first && cline.first == emitter->progenitor()->antiColourLine() &&
       cline.second !=newEmit->antiColourLine()) {
      // sort out not radiating line
      ColinePtr col = emitter->progenitor()->colourLine();
      if(col) {
	col->removeColoured(emitter->copy());
	col->removeColoured(emitter->progenitor());
	newEmit->colourLine()->removeColoured(newEmit);
	col->addColoured(newEmit);
      }
    }
    else if(aline.first && aline.first == emitter->progenitor()->colourLine() &&
	    aline.second !=newEmit->colourLine()) {
      // sort out not radiating line
      ColinePtr anti = emitter->progenitor()->antiColourLine();
      if(anti) {
	anti->removeAntiColoured(emitter->copy());
	anti->removeAntiColoured(emitter->progenitor());
	newEmit->colourLine()->removeAntiColoured(newEmit);
	anti->addAntiColoured(newEmit);
      }
    }
    else
      assert(false);
    // update the emitter
    emitter->copy(newEmit);
    ShowerParticlePtr sp = new_ptr(ShowerParticle(*newEmit,1,true));
    emitter->progenitor(sp);
    tree->outgoingLines()[emitter]=sp;
    emitter->perturbative(false);
    // update for decaying particles
    if(tit!=tree->treelinks().end())
      tree->updateLink(tit->first,make_pair(emitter,sp));
    // add the emitted particle
    // sort out the colour
    if(cline.first && cline.second==emitted->antiColourLine()) {
      cline.second->removeAntiColoured(emitted);
      cline.first->addAntiColoured(emitted);
    }
    else if(aline.first && aline.second==emitted->colourLine()) {
      aline.second->removeColoured(emitted);
      aline.first->addColoured(emitted);
    }
    else
      assert(false);
    sp=new_ptr(ShowerParticle(*emitted,1,true));
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(emitter->original(),
						       emitted,sp));
    gluon->perturbative(false);
    tree->outgoingLines().insert(make_pair(gluon,sp));
  }
  
}

void QTildeShowerHandler::setupMECorrection(RealEmissionProcessPtr real) {
  assert(real);
  currentTree()->hardMatrixElementCorrection(true);
  // II emission
  if(real->emitter()   < real->incoming().size() &&
     real->spectator() < real->incoming().size()) {
    // recoiling system
    for( map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	   cjt= currentTree()->outgoingLines().begin();
	 cjt != currentTree()->outgoingLines().end();++cjt ) {
      cjt->first->progenitor()->transform(real->transformation());
      cjt->first->copy()->transform(real->transformation());
    }
    // the the radiating system
    ShowerProgenitorPtr emitter,spectator;
    unsigned int iemit  = real->emitter();
    unsigned int ispect = real->spectator();
    int ig = int(real->emitted())-int(real->incoming().size());
    emitter = findInitialStateLine(currentTree(),
				   real->bornIncoming()[iemit]->id(),
				   real->bornIncoming()[iemit]->momentum());
    spectator = findInitialStateLine(currentTree(),
				     real->bornIncoming()[ispect]->id(),
				     real->bornIncoming()[ispect]->momentum());
    // sort out the colours
    ColinePair cline,aline;
    fixSpectatorColours(real->incoming()[ispect],spectator,cline,aline,true);
    // update the spectator
    spectator->copy(real->incoming()[ispect]);
    ShowerParticlePtr sp(new_ptr(ShowerParticle(*real->incoming()[ispect],1,false)));
    sp->x(ispect ==0 ? real->x().first :real->x().second);
    spectator->progenitor(sp);
    currentTree()->incomingLines()[spectator]=sp;
    spectator->perturbative(true);
    // now for the emitter
    fixInitialStateEmitter(currentTree(),real->incoming()[iemit],real->outgoing()[ig],
			   emitter,cline,aline,iemit ==0 ? real->x().first :real->x().second);
  }
  // FF emission
  else if(real->emitter()   >= real->incoming().size() &&
	  real->spectator() >= real->incoming().size()) {
    assert(real->outgoing()[real->emitted()-real->incoming().size()]->id()==ParticleID::g);
    // find the emitter and spectator in the shower tree
    ShowerProgenitorPtr emitter,spectator;
    int iemit = int(real->emitter())-int(real->incoming().size());
    emitter = findFinalStateLine(currentTree(),
				 real->bornOutgoing()[iemit]->id(),
				 real->bornOutgoing()[iemit]->momentum());
    int ispect = int(real->spectator())-int(real->incoming().size());
    spectator = findFinalStateLine(currentTree(),
				   real->bornOutgoing()[ispect]->id(),
				   real->bornOutgoing()[ispect]->momentum());
    map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
    // first the spectator
    // special case if decayed
    for(tit  = currentTree()->treelinks().begin(); tit != currentTree()->treelinks().end();++tit) {
      if(tit->second.first && tit->second.second==spectator->progenitor())
    	break;
    }
    // sort out the colours
    ColinePair cline,aline;
    fixSpectatorColours(real->outgoing()[ispect],spectator,cline,aline,true);
    // update the spectator
    spectator->copy(real->outgoing()[ispect]);
    ShowerParticlePtr sp(new_ptr(ShowerParticle(*real->outgoing()[ispect],1,true)));
    spectator->progenitor(sp);
    currentTree()->outgoingLines()[spectator]=sp;
    spectator->perturbative(true);
    // update for decaying particles
    if(tit!=currentTree()->treelinks().end())
      currentTree()->updateLink(tit->first,make_pair(spectator,sp));
    // now the emitting particle
    int ig = int(real->emitted())-int(real->incoming().size());
    fixFinalStateEmitter(currentTree(),real->outgoing()[iemit],
			 real->outgoing()[ig],
			 emitter,cline,aline);
  }
  // IF emission
  else {
    // scattering process
    if(real->incoming().size()==2) {
      ShowerProgenitorPtr emitter,spectator;
      unsigned int iemit  = real->emitter();
      unsigned int ispect = real->spectator();
      int ig = int(real->emitted())-int(real->incoming().size());
      ColinePair cline,aline;
      // incoming spectator
      if(ispect<2) {
	spectator = findInitialStateLine(currentTree(),
					 real->bornIncoming()[ispect]->id(),
					 real->bornIncoming()[ispect]->momentum());
	fixSpectatorColours(real->incoming()[ispect],spectator,cline,aline,true);
	// update the spectator
	spectator->copy(real->incoming()[ispect]);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*real->incoming()[ispect],1,false)));
	sp->x(ispect ==0 ? real->x().first :real->x().second);
	spectator->progenitor(sp);
	currentTree()->incomingLines()[spectator]=sp;
	spectator->perturbative(true);
      }
      // outgoing spectator
      else {
	spectator = findFinalStateLine(currentTree(),
				       real->bornOutgoing()[ispect-real->incoming().size()]->id(),
				       real->bornOutgoing()[ispect-real->incoming().size()]->momentum());
	// special case if decayed
	map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
	for(tit  = currentTree()->treelinks().begin(); tit != currentTree()->treelinks().end();++tit) {
	  if(tit->second.first && tit->second.second==spectator->progenitor())
	    break;
	}
	fixSpectatorColours(real->outgoing()[ispect-real->incoming().size()],spectator,cline,aline,true);
	// update the spectator
	spectator->copy(real->outgoing()[ispect-real->incoming().size()]);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*real->outgoing()[ispect-real->incoming().size()],1,true)));
	spectator->progenitor(sp);
	currentTree()->outgoingLines()[spectator]=sp;
	spectator->perturbative(true);
	// update for decaying particles
	if(tit!=currentTree()->treelinks().end())
	  currentTree()->updateLink(tit->first,make_pair(spectator,sp));
      }
      // incoming emitter
      if(iemit<2) {
	emitter = findInitialStateLine(currentTree(),
				       real->bornIncoming()[iemit]->id(),
				       real->bornIncoming()[iemit]->momentum());
	fixInitialStateEmitter(currentTree(),real->incoming()[iemit],real->outgoing()[ig],
			       emitter,aline,cline,iemit ==0 ? real->x().first :real->x().second);
      }
      // outgoing emitter
      else {
	emitter = findFinalStateLine(currentTree(),
				     real->bornOutgoing()[iemit-real->incoming().size()]->id(),
				     real->bornOutgoing()[iemit-real->incoming().size()]->momentum());
	fixFinalStateEmitter(currentTree(),real->outgoing()[iemit-real->incoming().size()],
			     real->outgoing()[ig],emitter,aline,cline);
      }
    }
    // decay process
    else {
      assert(real->spectator()==0);
      unsigned int iemit  = real->emitter()-real->incoming().size();
      int ig = int(real->emitted())-int(real->incoming().size());
      ColinePair cline,aline;
      // incoming spectator
      ShowerProgenitorPtr spectator = findInitialStateLine(currentTree(),
							   real->bornIncoming()[0]->id(),
							   real->bornIncoming()[0]->momentum());
      fixSpectatorColours(real->incoming()[0],spectator,cline,aline,false);
      // find the emitter
      ShowerProgenitorPtr emitter = 
	findFinalStateLine(currentTree(),
			   real->bornOutgoing()[iemit]->id(),
			   real->bornOutgoing()[iemit]->momentum());
      // recoiling system
      for( map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	     cjt= currentTree()->outgoingLines().begin();
	   cjt != currentTree()->outgoingLines().end();++cjt ) {
	if(cjt->first==emitter) continue;
	cjt->first->progenitor()->transform(real->transformation());
	cjt->first->copy()->transform(real->transformation());
      }
      // sort out the emitter
      fixFinalStateEmitter(currentTree(),real->outgoing()[iemit],
			   real->outgoing()[ig],emitter,aline,cline);
    }
  }
  // clean up the shower tree
  _currenttree->resetShowerProducts();
}
