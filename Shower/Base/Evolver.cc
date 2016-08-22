// -*- C++ -*-
//
// Evolver.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ShowerKinematics.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Utilities/Throw.h"
#include "ShowerTree.h"
#include "ShowerProgenitor.h"
#include "KinematicsReconstructor.h"
#include "PartnerFinder.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Shower/ShowerHandler.h" 
#include "ThePEG/Utilities/DescribeClass.h"
#include "ShowerVertex.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Base/SubtractedME.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "ThePEG/Handlers/StandardXComb.h"

using namespace Herwig;

namespace {

  /**
   *  A struct to order the particles in the same way as in the DecayMode's
   */
  struct ParticleOrdering {
    /**
     *  Operator for the ordering
     * @param p1 The first ParticleData object
     * @param p2 The second ParticleData object
     */
    bool operator() (tcPDPtr p1, tcPDPtr p2) {
      return abs(p1->id()) > abs(p2->id()) ||
	( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
	( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
    }
  };
  typedef multiset<tcPDPtr,ParticleOrdering> OrderedParticles;

  /**
   * Cached lookup of decay modes.
   * Generator::findDecayMode() is not efficient.
   */
  tDMPtr findDecayMode(const string & tag) {
    static map<string,DMPtr> cache;
    map<string,DMPtr>::const_iterator pos = cache.find(tag);

    if ( pos != cache.end() ) 
    	return pos->second;

    tDMPtr dm = CurrentGenerator::current().findDecayMode(tag);
    cache[tag] = dm;
    return dm;
  }
}

DescribeClass<Evolver,Interfaced>
describeEvolver ("Herwig::Evolver","HwShower.so");

bool Evolver::_hardEmissionModeWarn = true;
bool Evolver::_missingTruncWarn = true;

IBPtr Evolver::clone() const {
  return new_ptr(*this);
}

IBPtr Evolver::fullclone() const {
  return new_ptr(*this);
}

void Evolver::persistentOutput(PersistentOStream & os) const {
  os << _model << _splittingGenerator << _maxtry 
     << _meCorrMode << _hardVetoMode << _hardVetoRead << _hardVetoReadOption
     << _limitEmissions << _spinOpt << _softOpt << _hardPOWHEG
     << ounit(_iptrms,GeV) << _beta << ounit(_gamma,GeV) << ounit(_iptmax,GeV) 
     << _vetoes << _fullShowerVetoes << _nReWeight << _reWeight
     << _trunc_Mode << _hardEmissionMode << _reconOpt
     << _massVetoOption << isMCatNLOSEvent << isMCatNLOHEvent
     << isPowhegSEvent << isPowhegHEvent
     << theFactorizationScaleFactor << theRenormalizationScaleFactor << ounit(muPt,GeV)
     << interaction_ << _maxTryFSR << _maxFailFSR << _fracFSR << interactions_.size();
  for(unsigned int ix=0;ix<interactions_.size();++ix) 
    os << oenum(interactions_[ix]);
}

void Evolver::persistentInput(PersistentIStream & is, int) {
  unsigned int isize;
  is >> _model >> _splittingGenerator >> _maxtry 
     >> _meCorrMode >> _hardVetoMode >> _hardVetoRead >> _hardVetoReadOption
     >> _limitEmissions >> _spinOpt >> _softOpt >> _hardPOWHEG
     >> iunit(_iptrms,GeV) >> _beta >> iunit(_gamma,GeV) >> iunit(_iptmax,GeV)
     >> _vetoes >> _fullShowerVetoes >> _nReWeight >> _reWeight
     >> _trunc_Mode >> _hardEmissionMode >> _reconOpt
     >> _massVetoOption >> isMCatNLOSEvent >> isMCatNLOHEvent
     >> isPowhegSEvent >> isPowhegHEvent
     >> theFactorizationScaleFactor >> theRenormalizationScaleFactor >> iunit(muPt,GeV)
     >> interaction_ >> _maxTryFSR >> _maxFailFSR >> _fracFSR >> isize;
  interactions_.resize(isize);
  for(unsigned int ix=0;ix<interactions_.size();++ix) 
    is >> ienum(interactions_[ix]);
}

void Evolver::doinit() {
  Interfaced::doinit();
  // interactions may have been changed through a setup file so we
  // clear it up here
  interactions_.clear();
  if(interaction_==0) {
    interactions_.push_back(ShowerInteraction::QCD);
    interactions_.push_back(ShowerInteraction::QED);
  }
  else if(interaction_==1) {
    interactions_.push_back(ShowerInteraction::QCD);
  }
  else if(interaction_==2) {
    interactions_.push_back(ShowerInteraction::QED);
    interactions_.push_back(ShowerInteraction::QCD);
  }
  else if(interaction_==3) {
    interactions_.push_back(ShowerInteraction::QED);
  }
  else if(interaction_==4) {
    interactions_.push_back(ShowerInteraction::Both);
  }
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
     &Evolver::_maxtry, 100, 1, 100000,
     false, false, Interface::limited);

  static Parameter<Evolver,unsigned int> interfaceNReWeight
    ("NReWeight",
     "The number of attempts for the shower when reweighting",
     &Evolver::_nReWeight, 100, 10, 10000,
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

  static Switch<Evolver, bool> ifaceHardVetoReadOption
    ("HardVetoReadOption",
     "Apply read-in scale veto to all collisions or just the primary one?",
     &Evolver::_hardVetoReadOption, false, false, false);
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

  static RefVector<Evolver,FullShowerVeto> interfaceFullShowerVetoes
    ("FullShowerVetoes",
     "The vetos to be appliede on the full final state of the shower",
     &Evolver::_fullShowerVetoes, -1, false, false, true, false, false);

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
  static SwitchOption interfaceLimitEmissionsOneEmission
    (interfaceLimitEmissions,
     "OneEmission",
     "Allow one emission in either the final state or initial state, but not both",
     4);

  static Switch<Evolver,bool> interfaceTruncMode
    ("TruncatedShower", "Include the truncated shower?", 
     &Evolver::_trunc_Mode, 1, false, false);
  static SwitchOption interfaceTruncMode0
    (interfaceTruncMode,"No","Truncated Shower is OFF", 0);
  static SwitchOption interfaceTruncMode1
    (interfaceTruncMode,"Yes","Truncated Shower is ON", 1);

  static Switch<Evolver,int> interfaceHardEmissionMode
    ("HardEmissionMode",
     "Whether to use ME corrections or POWHEG for the hardest emission",
     &Evolver::_hardEmissionMode, 0, false, false);
  static SwitchOption interfaceHardEmissionModeDecayMECorrection
    (interfaceHardEmissionMode,
     "DecayMECorrection",
     "Old fashioned ME correction for decays only",
     -1);
  static SwitchOption interfaceHardEmissionModeMECorrection
    (interfaceHardEmissionMode,
     "MECorrection",
     "Old fashioned ME correction",
     0);
  static SwitchOption interfaceHardEmissionModePOWHEG
    (interfaceHardEmissionMode,
     "POWHEG",
     "Powheg style hard emission using internal matrix elements",
     1);
  static SwitchOption interfaceHardEmissionModeMatchboxPOWHEG
    (interfaceHardEmissionMode,
     "MatchboxPOWHEG",
     "Powheg style emission for the hard process using Matchbox",
     2);
  static SwitchOption interfaceHardEmissionModeFullPOWHEG
    (interfaceHardEmissionMode,
     "FullPOWHEG",
     "Powheg style emission for the hard process using Matchbox "
     "and decays using internal matrix elements",
     3);

  static Switch<Evolver,unsigned int > interfaceInteractions
    ("Interactions",
     "The interactions to be used in the shower",
     &Evolver::interaction_, 1, false, false);
  static SwitchOption interfaceInteractionsQCDFirst
    (interfaceInteractions,
     "QCDFirst",
     "QCD first then QED",
     0);
  static SwitchOption interfaceInteractionsQCDOnly
    (interfaceInteractions,
     "QCDOnly",
     "Only QCD",
     1);
  static SwitchOption interfaceInteractionsQEDFirst
    (interfaceInteractions,
     "QEDFirst",
     "QED first then QCD",
     2);
  static SwitchOption interfaceInteractionsQEDOnly
    (interfaceInteractions,
     "QEDOnly",
     "Only QED",
     3);
  static SwitchOption interfaceInteractionsBothAtOnce
    (interfaceInteractions,
     "BothAtOnce",
     "Generate both at the same time",
     4);

  static Switch<Evolver,unsigned int> interfaceReconstructionOption
    ("ReconstructionOption",
     "Treatment of the reconstruction of the transverse momentum of "
     "a branching from the evolution scale.",
     &Evolver::_reconOpt, 0, false, false);
  static SwitchOption interfaceReconstructionOptionCutOff
    (interfaceReconstructionOption,
     "CutOff",
     "Use the cut-off masses in the calculation",
     0);
  static SwitchOption interfaceReconstructionOptionOffShell
    (interfaceReconstructionOption,
     "OffShell",
     "Use the off-shell masses in the calculation",
     1);
  static SwitchOption interfaceReconstructionOptionOffShell2
    (interfaceReconstructionOption,
     "OffShell2",
     "Use the off-shell masses in the calculation but only locally for each branching",
     2);

  static Switch<Evolver,unsigned int> interfaceMassVetoOption
    ("MassVetoOption",
     "Option for the handling of the mass vetos",
     &Evolver::_massVetoOption, 1, false, false);
  static SwitchOption interfaceMassVetoOptionReset
    (interfaceMassVetoOption,
     "Reset",
     "Try another branching without resetting the starting scale",
     0);
  static SwitchOption interfaceMassVetoOptionInclude
    (interfaceMassVetoOption,
     "Include",
     "Include the veto in the scale generation via the veto algorithm",
     1);


  static Switch<Evolver,unsigned int> interfaceSpinCorrelations
    ("SpinCorrelations",
     "Treatment of spin correlations in the parton shower",
     &Evolver::_spinOpt, 1, false, false);
  static SwitchOption interfaceSpinCorrelationsOff
    (interfaceSpinCorrelations,
     "No",
     "No spin correlations",
     0);
  static SwitchOption interfaceSpinCorrelationsSpin
    (interfaceSpinCorrelations,
     "Yes",
     "Include the azimuthal spin correlations only",
     1);

  static Switch<Evolver,unsigned int> interfaceSoftCorrelations
    ("SoftCorrelations",
     "Option for the treatment of soft correlations in the parton shower",
     &Evolver::_softOpt, 2, false, false);
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

  static Switch<Evolver,bool> interfaceHardPOWHEG
    ("HardPOWHEG",
     "Treatment of powheg emissions which are too hard to have a shower interpretation",
     &Evolver::_hardPOWHEG, false, false, false);
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

  static Parameter<Evolver,unsigned int> interfaceMaxTryFSR
    ("MaxTryFSR",
     "The maximum number of attempted FSR emissions in"
     " the generation of the FSR",
     &Evolver::_maxTryFSR, 100000, 10, 100000000,
     false, false, Interface::limited);

  static Parameter<Evolver,unsigned int> interfaceMaxFailFSR
    ("MaxFailFSR",
     "Maximum number of failures generating the FSR",
     &Evolver::_maxFailFSR, 100, 1, 100000000,
     false, false, Interface::limited);


  static Parameter<Evolver,double> interfaceFSRFailureFraction
    ("FSRFailureFraction",
     "Maximum fraction of events allowed to fail due to too many FSR emissions",
     &Evolver::_fracFSR, 0.001, 1e-10, 1,
     false, false, Interface::limited);
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

void Evolver::setupMaximumScales(const vector<ShowerProgenitorPtr> & p,
				 XCPtr xcomb) {
  // let POWHEG events radiate freely
  if(_hardEmissionMode>0&&hardTree()) {
    vector<ShowerProgenitorPtr>::const_iterator ckt = p.begin();
    for (; ckt != p.end(); ckt++) (*ckt)->maxHardPt(Constants::MaxEnergy);
    return;
  }
  // return if no vetos
  if (!hardVetoOn()) return; 
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
  if (!hardVetoXComb()||
      (hardVetoReadOption()&&
       !ShowerHandler::currentHandler()->firstInteraction())) {
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
      if(hardVetoXComb()&&hardVetoReadOption()&&
	 !ShowerHandler::currentHandler()->firstInteraction()) {
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
  ptmax *= ShowerHandler::currentHandler()->hardScaleFactor();
  // set maxHardPt for all progenitors.  For partonic processes this
  // is now the max pt in the FS, for non-partonic processes or
  // processes with no coloured FS the invariant mass of the IS
  vector<ShowerProgenitorPtr>::const_iterator ckt = p.begin();
  for (; ckt != p.end(); ckt++) (*ckt)->maxHardPt(ptmax);
}

void Evolver::setupHardScales(const vector<ShowerProgenitorPtr> & p,
			      XCPtr xcomb) {
  if ( hardVetoXComb() &&
       (!hardVetoReadOption() ||
	ShowerHandler::currentHandler()->firstInteraction()) ) {
    Energy hardScale = ZERO;
    if(currentTree()->isHard()) {
      assert(xcomb);
      hardScale = sqrt( xcomb->lastShowerScale() );
    }
    else {
      hardScale = currentTree()->incomingLines().begin()->first
	->progenitor()->momentum().mass(); 
    }
    hardScale *= ShowerHandler::currentHandler()->hardScaleFactor();
    vector<ShowerProgenitorPtr>::const_iterator ckt = p.begin();
    for (; ckt != p.end(); ckt++) (*ckt)->hardScale(hardScale);
    muPt = hardScale;
  }
}

void Evolver::showerHardProcess(ShowerTreePtr hard, XCPtr xcomb) {


  isMCatNLOSEvent = false;
  isMCatNLOHEvent = false;
  isPowhegSEvent  = false;
  isPowhegHEvent  = false;

  Ptr<SubtractedME>::tptr subme;
  Ptr<MatchboxMEBase>::tptr me;
  Ptr<SubtractionDipole>::tptr dipme;

  Ptr<StandardXComb>::ptr sxc = dynamic_ptr_cast<Ptr<StandardXComb>::ptr>(xcomb);

  if ( sxc ) {
    subme = dynamic_ptr_cast<Ptr<SubtractedME>::tptr>(sxc->matrixElement());
    me = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(sxc->matrixElement());
    dipme = dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(sxc->matrixElement());
  }

  if ( subme ) {
    if ( subme->showerApproximation() ) {
      theShowerApproximation = subme->showerApproximation();
      // separate MCatNLO and POWHEG-type corrections
      if ( !subme->showerApproximation()->needsSplittingGenerator() ) {
	if ( subme->realShowerSubtraction() )
	  isMCatNLOHEvent = true;
	else if ( subme->virtualShowerSubtraction() )
	  isMCatNLOSEvent = true;
      }
      else {
  	if ( subme->realShowerSubtraction() )
  	  isPowhegHEvent = true;
  	else if ( subme->virtualShowerSubtraction() ||  subme->loopSimSubtraction() )
  	  isPowhegSEvent = true;
      }
    }
  } else if ( me ) {
    if ( me->factory()->showerApproximation() ) {
      theShowerApproximation = me->factory()->showerApproximation();
      if ( !me->factory()->showerApproximation()->needsSplittingGenerator() ) 
	isMCatNLOSEvent = true;
      else
  	isPowhegSEvent = true;
    }
  }

  string error = "Inconsistent hard emission set-up in Evolver::showerHardProcess(). "; 
  if ( ( isMCatNLOSEvent || isMCatNLOHEvent ) ){
    if (_hardEmissionMode > 1)
      throw Exception() << error
			<< "Cannot generate POWHEG matching with MC@NLO shower "
			<< "approximation.  Add 'set Evolver:HardEmissionMode 0' to input file."
			<< Exception::runerror;
    if ( ShowerHandler::currentHandler()->canHandleMatchboxTrunc())
      throw Exception() << error
			<< "Cannot use truncated qtilde shower with MC@NLO shower "
			<< "approximation.  Set LHCGenerator:EventHandler"
			<< ":CascadeHandler to '/Herwig/Shower/ShowerHandler' or "
			<< "'/Herwig/DipoleShower/DipoleShowerHandler'."
			<< Exception::runerror;
  }
  else if ( ((isPowhegSEvent || isPowhegHEvent) || dipme) &&
	    _hardEmissionMode < 2){
    if ( ShowerHandler::currentHandler()->canHandleMatchboxTrunc())
      throw Exception() << error
			<< "Unmatched events requested for POWHEG shower "
			<< "approximation.  Set Evolver:HardEmissionMode to "
			<< "'MatchboxPOWHEG' or 'FullPOWHEG'."
			<< Exception::runerror;
    else if (_hardEmissionModeWarn){
      _hardEmissionModeWarn = false;
      _hardEmissionMode+=2;
      throw Exception() << error
			<< "Unmatched events requested for POWHEG shower "
			<< "approximation. Changing Evolver:HardEmissionMode from "
			<< _hardEmissionMode-2 << " to " << _hardEmissionMode
			<< Exception::warning;
    }
  }

  if ( isPowhegSEvent || isPowhegHEvent) {
    if (theShowerApproximation->needsTruncatedShower() &&
	!ShowerHandler::currentHandler()->canHandleMatchboxTrunc() )
      throw Exception() << error
			<< "Current shower handler cannot generate truncated shower.  "
			<< "Set Generator:EventHandler:CascadeHandler to "
			<< "'/Herwig/Shower/PowhegShowerHandler'."
			<< Exception::runerror;
  }
  else if ( dipme && _missingTruncWarn){
    _missingTruncWarn=false;
    throw Exception() << "Warning: POWHEG shower approximation used without "
		      << "truncated shower.  Set Generator:EventHandler:"
		      << "CascadeHandler to '/Herwig/Shower/PowhegShowerHandler' and "
		      << "'MEMatching:TruncatedShower Yes'."
		      << Exception::warning;   
  }
  else if ( !dipme && _hardEmissionMode > 1 && 
	    ShowerHandler::currentHandler()->firstInteraction())
    throw Exception() << error
		      << "POWHEG matching requested for LO events.  Include "
		      << "'set Factory:ShowerApproximation MEMatching' in input file."
		      << Exception::runerror;

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
  // number of attempts if more than one interaction switched on
  unsigned int interactionTry=0;
  do {
    try {
      // generate the showering
      doShowering(true,xcomb);
      // if no vetos return
      return;
    }
    catch (InteractionVeto) {
      currentTree()->clear();
      ++interactionTry;
    }
  }
  while(interactionTry<=5);
  throw Exception() << "Too many tries for shower in "
		    << "Evolver::showerHardProcess()"
		    << Exception::eventerror;
}

void Evolver::hardMatrixElementCorrection(bool hard) {
  // set the initial enhancement factors for the soft correction
  _initialenhance = 1.;
  _finalenhance   = 1.;
  // if hard matrix element switched off return
  if(!MECOn(hard)) return;
  // see if we can get the correction from the matrix element
  // or decayer
  if(hard) {
    if(_hardme&&_hardme->hasMECorrection()) {
      _hardme->initializeMECorrection(_currenttree,
				      _initialenhance,_finalenhance);
      if(hardMEC(hard))
	_hardme->applyHardMatrixElementCorrection(_currenttree);
    }
  }
  else {
    if(_decayme&&_decayme->hasMECorrection()) {
      _decayme->initializeMECorrection(_currenttree,
				       _initialenhance,_finalenhance);
      if(hardMEC(hard))
	_decayme->applyHardMatrixElementCorrection(_currenttree);
    }
  }
}

Branching Evolver::selectTimeLikeBranching(tShowerParticlePtr particle,
					   ShowerInteraction::Type type) {
  Branching fb;
  while (true) {
    fb=_splittingGenerator->chooseForwardBranching(*particle,_finalenhance,type);
    // no emission return
    if(!fb.kinematics) return fb;
    // if emission OK break
    if(!timeLikeVetoed(fb,particle)) break;
    // otherwise reset scale and continue - SO IS involved in veto algorithm
    particle->vetoEmission(fb.type,fb.kinematics->scale());
    if(particle->spinInfo()) particle->spinInfo()->decayVertex(VertexPtr());
  }
  return fb;
}

ShowerParticleVector Evolver::createTimeLikeChildren(tShowerParticlePtr particle, IdList ids) {
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  // if same as definition create particles, otherwise create cc
  tcPDPtr pdata[2];
  for(unsigned int ix=0;ix<2;++ix) pdata[ix]=getParticleData(ids[ix+1]);
  if(particle->id()!=ids[0]) {
    for(unsigned int ix=0;ix<2;++ix) {
      tPDPtr cc(pdata[ix]->CC());
      if(cc) pdata[ix]=cc;
    }
  }
  ShowerParticleVector children;
  for(unsigned int ix=0;ix<2;++ix) {
    children.push_back(new_ptr(ShowerParticle(pdata[ix],true)));
    if(children[ix]->id()==_progenitor->id()&&!pdata[ix]->stable())
      children[ix]->set5Momentum(Lorentz5Momentum(_progenitor->progenitor()->mass()));
    else
      children[ix]->set5Momentum(Lorentz5Momentum(pdata[ix]->mass()));
  }
  return children;
}

bool Evolver::timeLikeShower(tShowerParticlePtr particle, 
			     ShowerInteraction::Type type,
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
			<< "Evolver::timeLikeShower(). Terminating run\n"
			<< Exception::runerror;
    throw Exception() << "Too many attempted emissions in Evolver::timeLikeShower()\n"
		      << Exception::eventerror;
  }
  // generate the emission
  ShowerParticleVector children;
  int ntry=0;
  while (ntry<50) {
    ++ntry;
    // generate the emission
    if(!fb.kinematics) 
      fb = selectTimeLikeBranching(particle,type);
    // no emission, return
    if(!fb.kinematics) {
      if(particle->spinInfo()) particle->spinInfo()->develop();
      return false;
    }
    // has emitted
    // Assign the shower kinematics to the emitting particle.
    ++_nFSR;
    particle->showerKinematics(fb.kinematics);
    // generate phi
    fb.kinematics->phi(fb.sudakov->generatePhiForward(*particle,fb.ids,fb.kinematics));
    // check highest pT
    if(fb.kinematics->pT()>progenitor()->highestpT())
      progenitor()->highestpT(fb.kinematics->pT());
    // create the children
    children = createTimeLikeChildren(particle,fb.ids);
    // update the children
    particle->showerKinematics()->
      updateChildren(particle, children,fb.type);
    // update number of emissions
    ++_nfs;
    if(_limitEmissions!=0) {
      if(children[0]->spinInfo()) children[0]->spinInfo()->develop();
      if(children[1]->spinInfo()) children[1]->spinInfo()->develop();
      if(particle->spinInfo()) particle->spinInfo()->develop();
      return true;
    }
    // select branchings for children
    Branching fc[2] = {selectTimeLikeBranching(children[0],type),
		       selectTimeLikeBranching(children[1],type)};
    // old recon option
    if(_reconOpt==0) {
      // shower the first  particle
      if(fc[0].kinematics) timeLikeShower(children[0],type,fc[0],false);
      if(children[0]->spinInfo()) children[0]->spinInfo()->develop();
      // shower the second particle
      if(fc[1].kinematics) timeLikeShower(children[1],type,fc[1],false);
      if(children[1]->spinInfo()) children[1]->spinInfo()->develop();
    }
    else if(_reconOpt==1) {
      // shower the first  particle
      if(fc[0].kinematics) timeLikeShower(children[0],type,fc[0],false);
      if(children[0]->spinInfo()) children[0]->spinInfo()->develop();
      // shower the second particle
      if(fc[1].kinematics) timeLikeShower(children[1],type,fc[1],false);
      if(children[1]->spinInfo()) children[1]->spinInfo()->develop();
      // branching has happened
      particle->showerKinematics()->
	updateParent(particle, children,fb.type);
      // clean up the vetoed emission
      if(particle->virtualMass()==ZERO) {
	particle->showerKinematics(ShoKinPtr());
	for(unsigned int ix=0;ix<children.size();++ix)
	  particle->abandonChild(children[ix]);
	children.clear();
	if(particle->spinInfo()) particle->spinInfo()->decayVertex(VertexPtr());
	if(_massVetoOption==1) particle->vetoEmission(fb.type,fb.kinematics->scale());
	fb = Branching();
	continue;
      }
    }
    else if(_reconOpt==2) {
      // cut-off masses for the branching
      const vector<Energy> & virtualMasses = fb.sudakov->virtualMasses(fb.ids);
      // compute the masses of the children
      Energy masses[3];
      for(unsigned int ix=0;ix<2;++ix) {
	if(fc[ix].kinematics) {
	  const vector<Energy> & vm = fc[ix].sudakov->virtualMasses(fc[ix].ids);
	  Energy2 q2 = 
	    fc[ix].kinematics->z()*(1.-fc[ix].kinematics->z())*sqr(fc[ix].kinematics->scale());
	  if(fc[ix].ids[0]!=ParticleID::g) q2 += sqr(vm[0]);
	  masses[ix+1] = sqrt(q2);
	}
	else {
	  masses[ix+1] = virtualMasses[ix+1];
	}
      }
      masses[0] = fb.ids[0]!=ParticleID::g ? virtualMasses[0] : ZERO;
      double z = fb.kinematics->z();
      Energy2 pt2 = z*(1.-z)*(z*(1.-z)*sqr(fb.kinematics->scale()) 
			      +sqr(masses[0]))
	- sqr(masses[1])*(1.-z) - sqr(masses[2])*z;
      if(pt2>=ZERO) {
	// branching has happened
	particle->showerKinematics()->
	  updateParent(particle, children,fb.type);
	// shower the first  particle
	if(fc[0].kinematics) timeLikeShower(children[0],type,fc[0],false);
	if(children[0]->spinInfo()) children[0]->spinInfo()->develop();
	// shower the second particle
	if(fc[1].kinematics) timeLikeShower(children[1],type,fc[1],false);
	if(children[1]->spinInfo()) children[1]->spinInfo()->develop();
      }
      else {
	particle->showerKinematics(ShoKinPtr());
	for(unsigned int ix=0;ix<children.size();++ix)
	  particle->abandonChild(children[ix]);
	children.clear();
	if(_massVetoOption==1) particle->vetoEmission(fb.type,fb.kinematics->scale());
	if(particle->spinInfo()) particle->spinInfo()->decayVertex(VertexPtr());
	fb = Branching();
	continue;
      }
    }
    break;
  };
  if(first&&!children.empty())
    particle->showerKinematics()->resetChildren(particle,children);
  if(particle->spinInfo()) particle->spinInfo()->develop();
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
    updateChildren(newParent, theChildren,bb.type);
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
  unsigned int interactionTry=0;
  do {
    try {
      // generate the showering
      doShowering(false,XCPtr());
      // if no vetos 
      // force calculation of spin correlations
      SpinPtr spInfo = decay->incomingLines().begin()->first->progenitor()->spinInfo();
      if(spInfo) {
	if(!spInfo->developed()) spInfo->needsUpdate();
	spInfo->develop();
      }
      // and then return
      return;
    }
    catch (InteractionVeto) {
      currentTree()->clear();
      ++interactionTry;
    }
  }
  while(interactionTry<=5);
  throw Exception() << "Too many tries for QED shower in Evolver::showerDecay()"
		    << Exception::eventerror;
}

bool Evolver::spaceLikeDecayShower(tShowerParticlePtr particle,
				   const ShowerParticle::EvolutionScales & maxScales,
				   Energy minmass,ShowerInteraction::Type type) {
  Branching fb;
  while (true) {
    fb=_splittingGenerator->chooseDecayBranching(*particle,maxScales,minmass,
						 _initialenhance,type);
    // return if no radiation
    if(!fb.kinematics) return false;
    // if not vetoed break
    if(!spaceLikeDecayVetoed(fb,particle)) break;
    // otherwise reset scale and continue
    particle->vetoEmission(fb.type,fb.kinematics->scale());
  }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->showerKinematics(fb.kinematics);
  if(fb.kinematics->pT()>progenitor()->highestpT())
    progenitor()->highestpT(fb.kinematics->pT());
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
  for(unsigned int ix=0;ix<2;++ix) {
    theChildren.push_back(new_ptr(ShowerParticle(pdata[ix],true)));
    if(theChildren[ix]->id()==_progenitor->id()&&!pdata[ix]->stable())
      theChildren[ix]->set5Momentum(Lorentz5Momentum(_progenitor->progenitor()->mass()));
    else
      theChildren[ix]->set5Momentum(Lorentz5Momentum(pdata[ix]->mass()));
  }
  // some code moved to updateChildren
  particle->showerKinematics()->
    updateChildren(particle, theChildren, fb.type);
  // In the case of splittings which involves coloured particles,
  // set properly the colour flow of the branching.
  // update the history if needed
  _currenttree->updateInitialStateShowerProduct(_progenitor,theChildren[0]);
  _currenttree->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
  // shower the first  particle
  spaceLikeDecayShower(theChildren[0],maxScales,minmass,type);
  // shower the second particle
  timeLikeShower(theChildren[1],type,Branching(),true);
  updateHistory(theChildren[1]);
  // branching has happened
  return true;
}

vector<ShowerProgenitorPtr> Evolver::setupShower(bool hard) {
  // generate POWHEG hard emission if needed
  if(_hardEmissionMode>0) hardestEmission(hard);
  ShowerInteraction::Type inter = interactions_[0];
  if(_hardtree&&inter!=ShowerInteraction::Both) {
    inter = _hardtree->interaction();
  }
  // set the initial colour partners
  setEvolutionPartners(hard,inter,false);
  // generate hard me if needed
  if(_hardEmissionMode==0 ||
     (!hard && _hardEmissionMode==-1)) hardMatrixElementCorrection(hard);
  // get the particles to be showered
  vector<ShowerProgenitorPtr> particlesToShower = 
    currentTree()->extractProgenitors();
  // remake the colour partners if needed
  if(_currenttree->hardMatrixElementCorrection()) {
    setEvolutionPartners(hard,interactions_[0],true);
    _currenttree->resetShowerProducts();
  }
  // return the answer
  return particlesToShower;
}

void Evolver::setEvolutionPartners(bool hard,ShowerInteraction::Type type,
				   bool clear) {
  // match the particles in the ShowerTree and hardTree
  if(hardTree() && !hardTree()->connect(currentTree()))
    throw Exception() << "Can't match trees in "
		      << "Evolver::setEvolutionPartners()"
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
      tHardBranchingPtr partner = 
	hardTree()->particles()[particles[ix]]->colourPartner();
      if(!partner) continue;
      for(map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator
	    it=hardTree()->particles().begin();
	  it!=hardTree()->particles().end();++it) {
	if(it->second==partner) particles[ix]->partner(it->first);
      }
      if(!particles[ix]->partner()) 
	throw Exception() << "Can't match partners in "
			  << "Evolver::setEvolutionPartners()"
			  << Exception::eventerror;
    }
  }
  // Set the initial evolution scales
  showerModel()->partnerFinder()->
    setInitialEvolutionScales(particles,!hard,type,!_hardtree);
  if(hardTree() && _hardPOWHEG) {
    bool tooHard=false;
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit=hardTree()->particles().end();
    for(unsigned int ix=0;ix<particles.size();++ix) {
      map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
	mit = hardTree()->particles().find(particles[ix]);
      Energy hardScale(ZERO);
      ShowerPartnerType::Type type(ShowerPartnerType::Undefined);
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

void Evolver::updateHistory(tShowerParticlePtr particle) {
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

bool Evolver::startTimeLikeShower(ShowerInteraction::Type type) {
  _nFSR = 0;
  if(hardTree()) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit=hardTree()->particles().end(),
      mit = hardTree()->particles().find(progenitor()->progenitor());
    if( mit != eit && !mit->second->children().empty() ) {
      bool output=truncatedTimeLikeShower(progenitor()->progenitor(),
					  mit->second ,type,true);
      if(output) updateHistory(progenitor()->progenitor());
      return output;
    }
  }
  bool output = hardOnly() ? false :
    timeLikeShower(progenitor()->progenitor() ,type,Branching(),true) ;
  if(output) updateHistory(progenitor()->progenitor());
  return output;
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

bool Evolver::
startSpaceLikeDecayShower(const ShowerParticle::EvolutionScales & maxScales,
			  Energy minimumMass,ShowerInteraction::Type type) {
  if(hardTree()) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit =hardTree()->particles().end(),
      mit = hardTree()->particles().find(progenitor()->progenitor());
    if( mit != eit && mit->second->parent() ) {
      HardBranchingPtr branch=mit->second;
      while(branch->parent()) branch=branch->parent();
      return truncatedSpaceLikeDecayShower(progenitor()->progenitor(),maxScales,
					   minimumMass, branch ,type);
    }
  }
  return  hardOnly() ? false :
    spaceLikeDecayShower(progenitor()->progenitor(),maxScales,minimumMass,type);
}

bool Evolver::timeLikeVetoed(const Branching & fb,
			     ShowerParticlePtr particle) {
  // work out type of interaction
  ShowerInteraction::Type type = fb.type==ShowerPartnerType::QED ? 
    ShowerInteraction::QED : ShowerInteraction::QCD;
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
  if(fb.kinematics->pT()>_progenitor->maximumpT(type)) return true;
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
  if ( ShowerHandler::currentHandler()->firstInteraction() &&
       ShowerHandler::currentHandler()->profileScales() ) {
    double weight = 
      ShowerHandler::currentHandler()->profileScales()->
      hardScaleProfile(_progenitor->hardScale(),fb.kinematics->pT());
    if ( UseRandom::rnd() > weight )
      return true;
  }
  return false;
}

bool Evolver::spaceLikeVetoed(const Branching & bb,
			      ShowerParticlePtr particle) {
  // work out type of interaction
  ShowerInteraction::Type type = bb.type==ShowerPartnerType::QED ? 
    ShowerInteraction::QED : ShowerInteraction::QCD;
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
  if(bb.kinematics->pT()>_progenitor->maximumpT(type)) return true;

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
  if ( ShowerHandler::currentHandler()->firstInteraction() &&
       ShowerHandler::currentHandler()->profileScales() ) {
    double weight = 
      ShowerHandler::currentHandler()->profileScales()->
      hardScaleProfile(_progenitor->hardScale(),bb.kinematics->pT());
    if ( UseRandom::rnd() > weight )
      return true;
  }
  return false;
}

bool Evolver::spaceLikeDecayVetoed( const Branching & fb,
				    ShowerParticlePtr particle) {
  // work out type of interaction
  ShowerInteraction::Type type = fb.type==ShowerPartnerType::QED ? 
    ShowerInteraction::QED : ShowerInteraction::QCD;
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
  HardTreePtr ISRTree;
  if(  ( _hardme  &&  _hardme->hasPOWHEGCorrection()!=0  && _hardEmissionMode< 2) ||
       ( _decayme && _decayme->hasPOWHEGCorrection()!=0  && _hardEmissionMode!=2) ) {
    if(_hardme) {
      assert(hard);
      if(interaction_==4) {
	vector<ShowerInteraction::Type> inter(2);
	inter[0] = ShowerInteraction::QCD;
	inter[1] = ShowerInteraction::QED;
	_hardtree =  _hardme->generateHardest( currentTree(),inter         );
      }
      else {
	_hardtree =  _hardme->generateHardest( currentTree(),interactions_ );
      }
    }
    else {
      assert(!hard);
      _hardtree = _decayme->generateHardest( currentTree() );
    }
    // store initial state POWHEG radiation
    if(_hardtree && _hardme && _hardme->hasPOWHEGCorrection()==1) 
      ISRTree=_hardtree;
  }

  else if (_hardEmissionMode>1 && hard) {
    // Get minimum pT cutoff used in shower approximation
    Energy maxpt = 1.*GeV;
    int colouredIn  = 0;
    int colouredOut = 0;
    for( map< ShowerProgenitorPtr, tShowerParticlePtr >::iterator it
	   = currentTree()->outgoingLines().begin(); 
	 it != currentTree()->outgoingLines().end(); ++it ) {
      if( it->second->coloured() ) colouredOut+=1;
    }  
    for( map< ShowerProgenitorPtr, ShowerParticlePtr >::iterator it
	   = currentTree()->incomingLines().begin(); 
	 it != currentTree()->incomingLines().end(); ++it ) {
      if( ! it->second->coloured() ) colouredIn+=1;
    }

    if ( theShowerApproximation ){
      if ( theShowerApproximation->ffPtCut() == theShowerApproximation->fiPtCut() &&
	   theShowerApproximation->ffPtCut() == theShowerApproximation->iiPtCut() ) 
	maxpt = theShowerApproximation->ffPtCut();
      else if ( colouredIn == 2 && colouredOut == 0 )
	maxpt = theShowerApproximation->iiPtCut();
      else if ( colouredIn == 0 && colouredOut > 1 )
	maxpt = theShowerApproximation->ffPtCut();
      else if ( colouredIn == 2 && colouredOut == 1 )
	maxpt = min(theShowerApproximation->iiPtCut(), theShowerApproximation->fiPtCut());
      else if ( colouredIn == 1 && colouredOut > 1 )
	maxpt = min(theShowerApproximation->ffPtCut(), theShowerApproximation->fiPtCut());
      else 
	maxpt = min(min(theShowerApproximation->iiPtCut(), theShowerApproximation->fiPtCut()), 
		    theShowerApproximation->ffPtCut());
    }

    // Generate hardtree from born and real emission subprocesses
    _hardtree = ShowerHandler::currentHandler()->generateCKKW(currentTree());

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
    maxpt=min(sqrt(ShowerHandler::currentHandler()->lastXCombPtr()->lastShowerScale()),maxpt);

    // Set maxpt to pT of emission when showering POWHEG real-emission subprocesses
    if (!isPowhegSEvent && !isPowhegHEvent){
      vector<int> outGluon;
      vector<int> outQuark;
      map< ShowerProgenitorPtr, tShowerParticlePtr >::iterator it;
      for( it = currentTree()->outgoingLines().begin(); 
	   it != currentTree()->outgoingLines().end(); ++it ) {
	if ( abs(it->second->id())< 6) outQuark.push_back(it->second->id());
	if ( it->second->id()==21 )    outGluon.push_back(it->second->id());
      } 
      if (outGluon.size() + outQuark.size() == 1){
	for( it = currentTree()->outgoingLines().begin(); 
	     it != currentTree()->outgoingLines().end(); ++it ) {
	  if ( abs(it->second->id())< 6 || it->second->id()==21 )
	    maxpt = it->second->momentum().perp();
	}
      }
      else if (outGluon.size() + outQuark.size() > 1){
	// assume qqbar pair from a Z/gamma
	if (outGluon.size()==1 && outQuark.size() == 2 && outQuark[0]==-outQuark[1]){
	  for( it = currentTree()->outgoingLines().begin(); 
	       it != currentTree()->outgoingLines().end(); ++it ) {
	    if ( it->second->id()==21 )
	      maxpt = it->second->momentum().perp();
	  }
	}
	// otherwise take the lowest pT avoiding born DY events
	else {
	  maxpt = generator()->maximumCMEnergy();
	  for( it = currentTree()->outgoingLines().begin(); 
	       it != currentTree()->outgoingLines().end(); ++it ) {
	    if ( abs(it->second->id())< 6 || it->second->id()==21 )
	      maxpt = min(maxpt,it->second->momentum().perp());
	  }
	}
      }
    } 

    // set maximum pT for subsequent emissions from S events
    if ( isPowhegSEvent  || (!isPowhegSEvent && !isPowhegHEvent)){
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
    _hardtree = ShowerHandler::currentHandler()->generateCKKW(currentTree());

  // if hard me doesn't have a FSR powheg 
  // correction use decay powheg correction
  if (_hardme && _hardme->hasPOWHEGCorrection()<2) {      
    // check for intermediate colour singlet resonance
    const ParticleVector inter =  _hardme->subProcess()->intermediates();
    if (inter.size()!=1 || 
	inter[0]->momentum().m2()/GeV2 < 0 || 
	inter[0]->dataPtr()->iColour()!=PDT::Colour0){
      if(_hardtree) connectTrees(currentTree(),_hardtree,hard);
      return;
    }
   
    map<ShowerProgenitorPtr, tShowerParticlePtr > out = currentTree()->outgoingLines();
    // ignore cases where outgoing particles are not coloured
    if (out.size()!=2 ||
	out. begin()->second->dataPtr()->iColour()==PDT::Colour0 ||
    	out.rbegin()->second->dataPtr()->iColour()==PDT::Colour0) {
      if(_hardtree) connectTrees(currentTree(),_hardtree,hard);
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
    if (!decayer || decayer->hasPOWHEGCorrection()<2){
      if(_hardtree) connectTrees(currentTree(),_hardtree,hard);
      return;
    }
    
    // generate the hardest emission
    ShowerDecayMap decay;
    PPtr in = new_ptr(*inter[0]);
    ShowerTreePtr decayTree = new_ptr(ShowerTree(in, decay));
    HardTreePtr     FSRTree = decayer->generateHardest(decayTree); 
    if (!FSRTree) {
      if(_hardtree) connectTrees(currentTree(),_hardtree,hard);
      return;
    }

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
  if(_hardtree){
    connectTrees(currentTree(),_hardtree,hard); 
  }
}

bool Evolver::truncatedTimeLikeShower(tShowerParticlePtr particle,
				      HardBranchingPtr branch,
				      ShowerInteraction::Type type,bool first) {
  int ntry=0;
  do {
    ++ntry;
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
	particle->vetoEmission(fb.type,fb.kinematics->scale());
	continue;
      }
      double zsplit = iout==1 ? fb.kinematics->z() : 1-fb.kinematics->z();
      // only if same interaction for forced branching 
      ShowerInteraction::Type type2 = fb.type==ShowerPartnerType::QED ? 
	ShowerInteraction::QED : ShowerInteraction::QCD;
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
      // should do base class vetos as well
      if(timeLikeVetoed(fb,particle)) {
	particle->vetoEmission(fb.type,fb.kinematics->scale());
	continue;
      }
      break;
    }
    // if no branching force truncated emission
    if(!fb.kinematics) {
      // construct the kinematics for the hard emission
      ShoKinPtr showerKin=
	branch->sudakov()->createFinalStateBranching(branch->scale(),
						     branch->children()[0]->z(),
						     branch->phi(),
						     branch->children()[0]->pT());
      showerKin->initialize( *particle,PPtr() );
      IdList idlist(3);
      idlist[0] = particle->id();
      idlist[1] = branch->children()[0]->branchingParticle()->id();
      idlist[2] = branch->children()[1]->branchingParticle()->id();
      fb = Branching( showerKin, idlist, branch->sudakov(),branch->type() );
      // Assign the shower kinematics to the emitting particle.
      ++_nFSR;
      particle->showerKinematics( fb.kinematics );
      if(fb.kinematics->pT()>progenitor()->highestpT())
	progenitor()->highestpT(fb.kinematics->pT());
      // Assign the splitting function to the emitting particle. 
      // For the time being we are considering only 1->2 branching
      // Create the ShowerParticle objects for the two children of
      // the emitting particle; set the parent/child relationship
      // if same as definition create particles, otherwise create cc
      ShowerParticleVector theChildren;
      for(unsigned int ix=0;ix<2;++ix) {
	theChildren.push_back(new_ptr(ShowerParticle(branch->children()[ix]->
						     branchingParticle()->dataPtr(),true)));
	if(theChildren[ix]->id()==_progenitor->id()&&!theChildren[ix]->dataPtr()->stable())
	  theChildren[ix]->set5Momentum(Lorentz5Momentum(_progenitor->progenitor()->mass()));
	else
	  theChildren[ix]->set5Momentum(Lorentz5Momentum(theChildren[ix]->dataPtr()->mass()));
      }
      particle->showerKinematics()->
	updateChildren(particle, theChildren,fb.type);
      for(unsigned int ix=0;ix<2;++ix) {
	theChildren[ix]->scales().QED         = min(theChildren[ix]->scales().QED        ,particle->scales().QED        );
	theChildren[ix]->scales().QED_noAO    = min(theChildren[ix]->scales().QED_noAO   ,particle->scales().QED_noAO   );
	theChildren[ix]->scales().QCD_c       = min(theChildren[ix]->scales().QCD_c      ,particle->scales().QCD_c      );
	theChildren[ix]->scales().QCD_c_noAO  = min(theChildren[ix]->scales().QCD_c_noAO ,particle->scales().QCD_c_noAO );
	theChildren[ix]->scales().QCD_ac      = min(theChildren[ix]->scales().QCD_ac     ,particle->scales().QCD_ac     );
	theChildren[ix]->scales().QCD_ac_noAO = min(theChildren[ix]->scales().QCD_ac_noAO,particle->scales().QCD_ac_noAO);
      }
      // shower the first  particle
      if( branch->children()[0]->children().empty() ) {
	if( ! hardOnly() )
	  timeLikeShower(theChildren[0],type,Branching(),false);
      }
      else {
	truncatedTimeLikeShower( theChildren[0],branch->children()[0],type,false);
      } 
      // shower the second particle
      if( branch->children()[1]->children().empty() ) {
	if( ! hardOnly() )
	  timeLikeShower( theChildren[1] , type,Branching(),false);
      }
      else {
	truncatedTimeLikeShower( theChildren[1],branch->children()[1] ,type,false);
      }
      // that's if for old approach
      if(_reconOpt==0) return true;
      // branching has happened
      particle->showerKinematics()->updateParent(particle, theChildren,fb.type);
      // clean up the vetoed emission
      if(particle->virtualMass()==ZERO) {
	particle->showerKinematics(ShoKinPtr());
	for(unsigned int ix=0;ix<theChildren.size();++ix)
	  particle->abandonChild(theChildren[ix]);
	theChildren.clear();
	continue;
      }
      else {
	if(first&&!theChildren.empty())
	  particle->showerKinematics()->resetChildren(particle,theChildren);
	if(particle->spinInfo()) particle->spinInfo()->develop();
	return true;
      }
    }
    // has emitted
    // Assign the shower kinematics to the emitting particle.
    ++_nFSR;
    particle->showerKinematics(fb.kinematics);
    if(fb.kinematics->pT()>progenitor()->highestpT())
      progenitor()->highestpT(fb.kinematics->pT());
    // Assign the splitting function to the emitting particle. 
    // For the time being we are considering only 1->2 branching
    // Create the ShowerParticle objects for the two children of
    // the emitting particle; set the parent/child relationship
    // if same as definition create particles, otherwise create cc
    ShowerParticleVector theChildren; 
    for(unsigned int ix=0;ix<2;++ix) {
      theChildren.push_back( new_ptr( ShowerParticle( pdata[ix], true ) ) );
      if(theChildren[ix]->id()==_progenitor->id()&&!pdata[ix]->stable())
	theChildren[ix]->set5Momentum(Lorentz5Momentum(_progenitor->progenitor()->mass()));
      else
	theChildren[ix]->set5Momentum(Lorentz5Momentum(pdata[ix]->mass()));
    }
    particle->showerKinematics()->
      updateChildren( particle, theChildren , fb.type);
    // shower the first  particle
    if( iout == 1 ) truncatedTimeLikeShower( theChildren[0], branch , type ,false);
    else            timeLikeShower( theChildren[0]  , type,Branching(),false);
    // shower the second particle
    if( iout == 2 ) truncatedTimeLikeShower( theChildren[1], branch , type ,false);
    else            timeLikeShower( theChildren[1]  , type,Branching(),false);
    // that's if for old approach
    if(_reconOpt==0) return true;
    // branching has happened
    particle->showerKinematics()->updateParent(particle, theChildren,fb.type);
    // clean up the vetoed emission
    if(particle->virtualMass()==ZERO) {
      particle->showerKinematics(ShoKinPtr());
      for(unsigned int ix=0;ix<theChildren.size();++ix)
	particle->abandonChild(theChildren[ix]);
      theChildren.clear();
    }
    else {
      if(first&&!theChildren.empty())
	particle->showerKinematics()->resetChildren(particle,theChildren);
      if(particle->spinInfo()) particle->spinInfo()->develop();
      return true;
    }
  }
  while(ntry<50);
  return false;
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
      ShowerInteraction::Type type2 = bb.type==ShowerPartnerType::QED ? 
	ShowerInteraction::QED : ShowerInteraction::QCD;
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
    kinematics->initialize( *particle, beam );
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
      updateChildren( newParent, theChildren,bb.type);
    if(hardOnly()) return true;
    // perform the shower of the final-state particle
    if( timelike->children().empty() ) {
      timeLikeShower( otherChild , type,Branching(),true);
    }
    else {
      truncatedTimeLikeShower( otherChild, timelike , type,true);
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
    updateChildren( newParent, theChildren, bb.type);
  // perform the shower of the final-state particle
  timeLikeShower( otherChild , type,Branching(),true);
  updateHistory(otherChild);
  // return the emitted
  return true;
}

bool Evolver::
truncatedSpaceLikeDecayShower(tShowerParticlePtr particle, 
			      const ShowerParticle::EvolutionScales & maxScales,
			      Energy minmass, HardBranchingPtr branch,
			      ShowerInteraction::Type type) {
  Branching fb;
  unsigned int iout=0;
  tcPDPtr pdata[2];
  while (true) {
    // no truncated shower break
    if(!isTruncatedShowerON()||hardOnly()) break;
    fb=splittingGenerator()->chooseDecayBranching(*particle,maxScales,minmass,1.,type);
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
      particle->vetoEmission(fb.type,fb.kinematics->scale());
      continue;
    }
    ShowerInteraction::Type type2 = fb.type==ShowerPartnerType::QED ? 
      ShowerInteraction::QED : ShowerInteraction::QCD;
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
    // should do base class vetos as well
    // if not vetoed break
    if(!spaceLikeDecayVetoed(fb,particle)) break;
    // otherwise reset scale and continue
    particle->vetoEmission(fb.type,fb.kinematics->scale());
  }
  // this may not be currently used but in principle could be
  // and should be included
 if (!fb.kinematics) {
   // construct the kinematics for the hard emission
   ShoKinPtr showerKin=
     branch->sudakov()->createDecayBranching(branch->scale(),
					      branch->children()[0]->z(),
					      branch->phi(),
					      branch->children()[0]->pT());
   showerKin->initialize( *particle,PPtr() );
   IdList idlist(3);
   idlist[0] = particle->id();
   idlist[1] = branch->children()[0]->branchingParticle()->id();
   idlist[2] = branch->children()[1]->branchingParticle()->id();
   // create the branching
   fb = Branching( showerKin, idlist, branch->sudakov(),ShowerPartnerType::QCDColourLine  );
   // Assign the shower kinematics to the emitting particle.
   particle->showerKinematics( fb.kinematics );
   if(fb.kinematics->pT()>progenitor()->highestpT())
     progenitor()->highestpT(fb.kinematics->pT());
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
     updateChildren(particle, theChildren,fb.type);
   if(theChildren[0]->id()==particle->id()) {
     // update the history if needed
     currentTree()->updateInitialStateShowerProduct(progenitor(),theChildren[0]);
     currentTree()->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
     // shower the space-like particle
     if( branch->children()[0]->children().empty() ) {
	if( ! hardOnly() ) spaceLikeDecayShower(theChildren[0],maxScales,minmass,type);
     }
     else {
	truncatedSpaceLikeDecayShower( theChildren[0],maxScales,minmass,
				       branch->children()[0],type);
     }
     // shower the second particle
     if( branch->children()[1]->children().empty() ) {
       if( ! hardOnly() ) timeLikeShower( theChildren[1] , type,Branching(), true);
     }
     else {
       truncatedTimeLikeShower( theChildren[1],branch->children()[1] ,type,true);
     }
     updateHistory(theChildren[1]);
   }
   else {
     // update the history if needed
     currentTree()->updateInitialStateShowerProduct(progenitor(),theChildren[1]);
     currentTree()->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
     // shower the space-like particle
     if( branch->children()[1]->children().empty() ) {
	if( ! hardOnly() ) spaceLikeDecayShower(theChildren[1],maxScales,minmass,type);
     }
     else {
	truncatedSpaceLikeDecayShower( theChildren[1],maxScales,minmass,
				       branch->children()[1],type);
     }
     // shower the second particle
     if( branch->children()[0]->children().empty() ) {
       if( ! hardOnly() ) timeLikeShower( theChildren[0] , type, Branching(),true);
     }
     else {
       truncatedTimeLikeShower( theChildren[0],branch->children()[0] ,type,true);
     }
     updateHistory(theChildren[0]);
   }
   return true;
 }
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->showerKinematics(fb.kinematics);
  if(fb.kinematics->pT()>progenitor()->highestpT())
    progenitor()->highestpT(fb.kinematics->pT());
  // For the time being we are considering only 1->2 branching
  // Create the ShowerParticle objects for the two children of
  // the emitting particle; set the parent/child relationship
  // if same as definition create particles, otherwise create cc
  ShowerParticleVector theChildren; 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[0],true))); 
  theChildren.push_back(new_ptr(ShowerParticle(pdata[1],true)));
  particle->showerKinematics()->updateChildren(particle, theChildren,fb.type);
  // In the case of splittings which involves coloured particles,
  // set properly the colour flow of the branching.
  // update the history if needed
  currentTree()->updateInitialStateShowerProduct(progenitor(),theChildren[0]);
  currentTree()->addInitialStateBranching(particle,theChildren[0],theChildren[1]);
  // shower the first  particle
  truncatedSpaceLikeDecayShower(theChildren[0],maxScales,minmass,branch,type);
  // shower the second particle
  timeLikeShower(theChildren[1],type,Branching(),true);
  updateHistory(theChildren[1]);
  // branching has happened
  return true;
}

bool Evolver::constructDecayTree(vector<ShowerProgenitorPtr> & particlesToShower,
				 ShowerInteraction::Type inter) {
  Energy ptmax(-GeV);
  // get the maximum pt is all ready a hard tree
  if(hardTree()) {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      if(particlesToShower[ix]->maximumpT(inter)>ptmax&&
	 particlesToShower[ix]->progenitor()->isFinalState()) 
	ptmax = particlesToShower[ix]->maximumpT(inter);
    }
  }
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
  HardTreePtr QCDTree = new_ptr(HardTree(allBranchings,spaceBranchings,inter));
  // set the charge partners
  ShowerParticleVector particles;
  particles.push_back(spaceBranchings.back()->branchingParticle());
  for(set<HardBranchingPtr>::iterator cit=QCDTree->branchings().begin();
      cit!=QCDTree->branchings().end();++cit) {
    if((*cit)->status()==HardBranching::Outgoing)
      particles.push_back((*cit)->branchingParticle());
  }
  // get the partners
  showerModel()->partnerFinder()->setInitialEvolutionScales(particles,true,inter,true);
  // do the inverse recon
  if(!showerModel()->kinematicsReconstructor()->
     deconstructDecayJets(QCDTree,this,inter)) {
    return false;
  }
  // clear the old shower
  currentTree()->clear();
  // set the hard tree
  hardTree(QCDTree);
  // set the charge partners
  setEvolutionPartners(false,inter,false);
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
      cjt!=currentTree()->outgoingLines().end();++cjt) {
    particlesToShower.push_back(((*cjt).first));
    if(ptmax>ZERO) particlesToShower.back()->maximumpT(ptmax,inter);
  }
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
      eit=hardTree()->particles().end(),
      mit = hardTree()->particles().find(particlesToShower[ix]->progenitor());
    if( mit != eit) {
      if(mit->second->status()==HardBranching::Outgoing)
	particlesToShower[ix]->progenitor()->set5Momentum(mit->second->pVector());
    }
  }
  return true;
}

bool Evolver::constructHardTree(vector<ShowerProgenitorPtr> & particlesToShower,
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
	if(particlesToShower[ix]->original()->parents().empty())
	  first->beam(particlesToShower[ix]->original());
	else
	  first->beam(particlesToShower[ix]->original()->parents()[0]);
	last = first;
      }
      spaceBranchings.push_back(first);
      allBranchings.push_back(last);
    }
  }
  if(!noEmission) {
    HardTreePtr QCDTree = new_ptr(HardTree(allBranchings,spaceBranchings,
					   inter));
    // set the charge partners
    ShowerParticleVector particles;
    for(set<HardBranchingPtr>::iterator cit=QCDTree->branchings().begin();
	cit!=QCDTree->branchings().end();++cit) {
      particles.push_back((*cit)->branchingParticle());
    }
    // get the partners
    showerModel()->partnerFinder()->setInitialEvolutionScales(particles,false,
							      inter,true);
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
  setEvolutionPartners(true,inter,false);
  // get the particles to be showered
  particlesToShower = currentTree()->extractProgenitors();
  // reset momenta
  if(hardTree()) {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      map<ShowerParticlePtr,tHardBranchingPtr>::const_iterator 
	eit=hardTree()->particles().end(),
	mit = hardTree()->particles().find(particlesToShower[ix]->progenitor());
      if( mit != eit) {
	particlesToShower[ix]->progenitor()->set5Momentum(mit->second->showerMomentum());
      }
    }
  }
  return true;
}

void Evolver::constructTimeLikeLine(tHardBranchingPtr branch,
				       tShowerParticlePtr particle) {
  for(unsigned int ix=0;ix<particle->children().size();++ix) {
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
  // sort out the type of interaction
  if(!branch->children().empty()) {
    if(branch->branchingParticle()->id()==ParticleID::gamma ||
       branch->children()[0]->branchingParticle()->id()==ParticleID::gamma ||
       branch->children()[1]->branchingParticle()->id()==ParticleID::gamma)
      branch->type(ShowerPartnerType::QED);
    else {
      if(branch->branchingParticle()->id()==
	 branch->children()[0]->branchingParticle()->id()) {
	if(branch->branchingParticle()->dataPtr()->iColour()==PDT::Colour8) {
	  tShowerParticlePtr emittor = 
	    branch->branchingParticle()->showerKinematics()->z()>0.5 ?
	    branch->children()[0]->branchingParticle() : 
	    branch->children()[1]->branchingParticle();
	  if(branch->branchingParticle()->colourLine()==emittor->colourLine())
	    branch->type(ShowerPartnerType::QCDAntiColourLine);
	  else if(branch->branchingParticle()->antiColourLine()==emittor->antiColourLine())
	    branch->type(ShowerPartnerType::QCDColourLine);
	  else
	    assert(false);
	}
	else if(branch->branchingParticle()->colourLine()) {
	  branch->type(ShowerPartnerType::QCDColourLine);
	}
	else if(branch->branchingParticle()->antiColourLine()) {
	  branch->type(ShowerPartnerType::QCDAntiColourLine);
	}
	else
	  assert(false);
      }
      else if(branch->branchingParticle()->id()==ParticleID::g &&
	      branch->children()[0]->branchingParticle()->id()== 
	      -branch->children()[1]->branchingParticle()->id()) {
	if(branch->branchingParticle()->showerKinematics()->z()>0.5)
	  branch->type(ShowerPartnerType::QCDAntiColourLine);
	else
	  branch->type(ShowerPartnerType::QCDColourLine);
	
      }
      else
	assert(false);
    }
  }
}

void Evolver::constructSpaceLikeLine(tShowerParticlePtr particle,
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
  // sort out the type
  if(last->branchingParticle()      ->id() == ParticleID::gamma ||
     newBranch->branchingParticle() ->id() == ParticleID::gamma ||
     timeBranch->branchingParticle()->id() == ParticleID::gamma) {
    last->type(ShowerPartnerType::QED);
  }
  else if(last->branchingParticle()->id()==newBranch->branchingParticle()->id()) {
    if(last->branchingParticle()->id()==ParticleID::g) {
      if(last->branchingParticle()->colourLine()==
	 newBranch->branchingParticle()->colourLine()) {
	last->type(ShowerPartnerType::QCDAntiColourLine);
      }
      else {
	last->type(ShowerPartnerType::QCDColourLine);
      }
    }
    else if(last->branchingParticle()->hasColour()) {
      last->type(ShowerPartnerType::QCDColourLine);
    }
    else if(last->branchingParticle()->hasAntiColour()) {
      last->type(ShowerPartnerType::QCDAntiColourLine);
    }
    else
      assert(false);
  }
  else if(newBranch->branchingParticle()->id()==ParticleID::g) { 
    if(last->branchingParticle()->hasColour()) {
      last->type(ShowerPartnerType::QCDAntiColourLine);
    }
    else if(last->branchingParticle()->hasAntiColour()) {
      last->type(ShowerPartnerType::QCDColourLine);
    }
    else
      assert(false);
  }
  else if(newBranch->branchingParticle()->hasColour()) {
    last->type(ShowerPartnerType::QCDColourLine);
  }
  else if(newBranch->branchingParticle()->hasAntiColour()) {
    last->type(ShowerPartnerType::QCDAntiColourLine);
  }
  else {
    assert(false);
  }
  last=newBranch;
}

void Evolver::connectTrees(ShowerTreePtr showerTree, 
			   HardTreePtr hardTree, bool hard ) {
  ShowerParticleVector particles;
  // find the Sudakovs
  for(set<HardBranchingPtr>::iterator cit=hardTree->branchings().begin();
      cit!=hardTree->branchings().end();++cit) {
    // Sudakovs for ISR
    if((**cit).parent()&&(**cit).status()==HardBranching::Incoming) {
      ++_nis;
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
      ++_nfs;
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
				     << "Evolver::connectTrees()" 
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
    setInitialEvolutionScales(particles,!hard,hardTree->interaction(),
			      !hardTree->partnersSet());
  hardTree->partnersSet(true);
  // inverse reconstruction
  if(hard) {
    showerModel()->kinematicsReconstructor()->
      deconstructHardJets(hardTree,ShowerHandler::currentHandler()->evolver(),
			  hardTree->interaction());
  }
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
    if(!partner) throw Exception() << "Failed to match shower and hard trees in Evolver::hardestEmission"
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

void Evolver::doShowering(bool hard,XCPtr xcomb) {
  // order of the interactions
  bool showerOrder(true);
  // zero number of emissions
  _nis = _nfs = 0;
  // if MC@NLO H event and limited emissions
  // indicate both final and initial state emission
  if ( isMCatNLOHEvent && _limitEmissions != 0 ) {
    _nis = _nfs = 1;
  }
  // extract particles to shower
  vector<ShowerProgenitorPtr> particlesToShower(setupShower(hard));
  // setup the maximum scales for the shower
  if (hardVetoOn()) setupMaximumScales(particlesToShower,xcomb);
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
	if(particlesToShower[ix]->progenitor()->dataPtr()->stable()) 
	  minmass += particlesToShower[ix]->progenitor()->dataPtr()->constituentMass();
	else
	  minmass += particlesToShower[ix]->progenitor()->mass();
      }
      else {
	mIn = particlesToShower[ix]->progenitor()->mass();
      }
    }
    // throw exception if decay can't happen
    if ( minmass > mIn ) {
      throw Exception() << "Evolver.cc: Mass of decaying particle is "
			<< "below constituent masses of decay products."
			<< Exception::eventerror;
    }
  }
  // check if interactions in right order
  if(hardTree() && interaction_!=4 && 
     hardTree()->interaction()!=interactions_[0]) {
    assert(interactions_.size()==2);
    showerOrder = false;
    swap(interactions_[0],interactions_[1]);
  }
  // loop over possible interactions
  bool reWeighting = _reWeight && hard && ShowerHandler::currentHandler()->firstInteraction();
  double eventWeight=0.;
  unsigned int nTryReWeight(0);
  for(unsigned int inter=0;inter<interactions_.size();++inter) {
    // set up for second pass if required
    if(inter!=0) {
      // zero intrinsic pt so only added first time round
      intrinsicpT().clear();
      // construct the tree and throw veto if not possible
      if(!(hard ? 
	   constructHardTree (particlesToShower,interactions_[inter]) :
	   constructDecayTree(particlesToShower,interactions_[inter]))) 
	throw InteractionVeto();
    }
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
	setEvolutionPartners(hard,interactions_[inter],true);
	_nis = _nfs = 0;
	// if MC@NLO H event and limited emissions
	// indicate both final and initial state emission
	if ( isMCatNLOHEvent && _limitEmissions != 0 ) {
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
	  if(!isFSRadiationON()) continue;
	  // perform shower
	  progenitor()->hasEmitted(startTimeLikeShower(interactions_[inter]));
	}
	// initial-state radiation
	else {
	  if(!isISRadiationON()) continue;
	  // hard process
	  if(hard) {
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
							  interactions_[inter]));
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
							       interactions_[inter]));
	  }
	}
      }
      // do the kinematic reconstruction, checking if it worked
      reconstructed = hard ?
	showerModel()->kinematicsReconstructor()->
	reconstructHardJets (currentTree(),intrinsicpT(),interactions_[inter],
			     switchRecon && ntry>maximumTries()/2) :
	showerModel()->kinematicsReconstructor()->
	reconstructDecayJets(currentTree(),interactions_[inter]);
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
			  << ntry << " attempts in Evolver::showerDecay()"
			  << Exception::eventerror;
    }
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
  if(!showerOrder) swap(interactions_[0],interactions_[1]);
  hardTree(HardTreePtr());
}

void Evolver:: convertHardTree(bool hard,ShowerInteraction::Type type) {
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
