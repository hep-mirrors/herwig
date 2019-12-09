// -*- C++ -*-
//
// ShowerHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerHandler class.
//

#include "ShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/PartonBinInstance.h"
#include "Herwig/PDT/StandardMatchers.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Utilities/EnumParticles.h"
#include "Herwig/PDF/MPIPDF.h"
#include "Herwig/PDF/MinBiasPDF.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "Herwig/PDF/HwRemDecayer.h"
#include <cassert>
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"

using namespace Herwig;

DescribeClass<ShowerHandler,CascadeHandler>
describeShowerHandler ("Herwig::ShowerHandler","HwShower.so");

ShowerHandler::~ShowerHandler() {}

tShowerHandlerPtr ShowerHandler::currentHandler_ = tShowerHandlerPtr();

void ShowerHandler::doinit() {
  CascadeHandler::doinit();
  // copy particles to decay before showering from input vector to the 
  // set used in the simulation
  if ( particlesDecayInShower_.empty() ) {
    for(unsigned int ix=0;ix<inputparticlesDecayInShower_.size();++ix)
      particlesDecayInShower_.insert(abs(inputparticlesDecayInShower_[ix]));
  }
  if ( profileScales() ) {
    if ( profileScales()->unrestrictedPhasespace() &&
	 restrictPhasespace() ) {
      generator()->log()
	<< "ShowerApproximation warning: The scale profile chosen requires an unrestricted phase space,\n"
	<< "however, the phase space was set to be restricted. Will switch to unrestricted phase space.\n"
	<< flush;
      restrictPhasespace_ = false;
    }
  }
}

IBPtr ShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr ShowerHandler::fullclone() const {
  return new_ptr(*this);
}

ShowerHandler::ShowerHandler() : 
  maxtry_(10),maxtryMPI_(10),maxtryDP_(10),maxtryDecay_(100),
  factorizationScaleFactor_(1.0),
  renormalizationScaleFactor_(1.0),
  hardScaleFactor_(1.0),
  restrictPhasespace_(true), maxPtIsMuF_(false), 
  spinOpt_(1), pdfFreezingScale_(2.5*GeV),
  doFSR_(true), doISR_(true),
  splitHardProcess_(true),
  includeSpaceTime_(false), vMin_(0.1*GeV2),
  reweight_(1.0) {
  inputparticlesDecayInShower_.push_back( 6  ); //  top 
  inputparticlesDecayInShower_.push_back( 23 ); // Z0
  inputparticlesDecayInShower_.push_back( 24 ); // W+/-
  inputparticlesDecayInShower_.push_back( 25 ); // h0
}

void ShowerHandler::doinitrun(){
  CascadeHandler::doinitrun();
  //can't use isMPIOn here, because the EventHandler is not set at that stage
  if(MPIHandler_) { 
    MPIHandler_->initialize();
    if(MPIHandler_->softInt())
      remDec_->initSoftInteractions(MPIHandler_->Ptmin(), MPIHandler_->beta());
  }
}

void ShowerHandler::dofinish() {
  CascadeHandler::dofinish();
  if(MPIHandler_) MPIHandler_->finalize();
}

void ShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << remDec_ << ounit(pdfFreezingScale_,GeV) << maxtry_ 
     << maxtryMPI_ << maxtryDP_ << maxtryDecay_ 
     << inputparticlesDecayInShower_
     << particlesDecayInShower_ << MPIHandler_ << PDFA_ << PDFB_
     << PDFARemnant_ << PDFBRemnant_
     << includeSpaceTime_ << ounit(vMin_,GeV2)
     << factorizationScaleFactor_ << renormalizationScaleFactor_
     << hardScaleFactor_
     << restrictPhasespace_ << maxPtIsMuF_ << hardScaleProfile_
     << showerVariations_ << doFSR_ << doISR_ << splitHardProcess_
     << spinOpt_ << useConstituentMasses_;
}

void ShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> remDec_ >> iunit(pdfFreezingScale_,GeV) >> maxtry_ 
     >> maxtryMPI_ >> maxtryDP_ >> maxtryDecay_
     >> inputparticlesDecayInShower_
     >> particlesDecayInShower_ >> MPIHandler_ >> PDFA_ >> PDFB_
     >> PDFARemnant_ >> PDFBRemnant_
     >> includeSpaceTime_ >> iunit(vMin_,GeV2)
     >> factorizationScaleFactor_ >> renormalizationScaleFactor_
     >> hardScaleFactor_
     >> restrictPhasespace_ >> maxPtIsMuF_ >> hardScaleProfile_
     >> showerVariations_ >> doFSR_ >> doISR_ >> splitHardProcess_
     >> spinOpt_ >> useConstituentMasses_;
}

void ShowerHandler::Init() {

  static ClassDocumentation<ShowerHandler> documentation
    ("Main driver class for the showering.");

  static Reference<ShowerHandler,HwRemDecayer> 
    interfaceRemDecayer("RemDecayer", 
			"A reference to the Remnant Decayer object", 
			&Herwig::ShowerHandler::remDec_,
			false, false, true, false);

  static Parameter<ShowerHandler,Energy> interfacePDFFreezingScale
    ("PDFFreezingScale",
     "The PDF freezing scale",
     &ShowerHandler::pdfFreezingScale_, GeV, 2.5*GeV, 2.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<ShowerHandler,unsigned int> interfaceMaxTry
    ("MaxTry",
     "The maximum number of attempts for the main showering loop",
     &ShowerHandler::maxtry_, 10, 1, 100,
     false, false, Interface::limited);

  static Parameter<ShowerHandler,unsigned int> interfaceMaxTryMPI
    ("MaxTryMPI",
     "The maximum number of regeneration attempts for an additional scattering",
     &ShowerHandler::maxtryMPI_, 10, 0, 100,
     false, false, Interface::limited);

  static Parameter<ShowerHandler,unsigned int> interfaceMaxTryDP
    ("MaxTryDP",
     "The maximum number of regeneration attempts for an additional hard scattering",
     &ShowerHandler::maxtryDP_, 10, 0, 100,
     false, false, Interface::limited);

  static ParVector<ShowerHandler,long> interfaceDecayInShower
    ("DecayInShower",
     "PDG codes of the particles to be decayed in the shower",
     &ShowerHandler::inputparticlesDecayInShower_, -1, 0l, -10000000l, 10000000l,
     false, false, Interface::limited);

  static Reference<ShowerHandler,UEBase> interfaceMPIHandler
    ("MPIHandler",
     "The object that administers all additional scatterings.",
     &ShowerHandler::MPIHandler_, false, false, true, true);

  static Reference<ShowerHandler,PDFBase> interfacePDFA
    ("PDFA",
     "The PDF for beam particle A. Overrides the particle's own PDF setting."
     "By default used for both the shower and forced splitting in the remnant",
     &ShowerHandler::PDFA_, false, false, true, true, false);

  static Reference<ShowerHandler,PDFBase> interfacePDFB
    ("PDFB",
     "The PDF for beam particle B. Overrides the particle's own PDF setting."
     "By default used for both the shower and forced splitting in the remnant",
     &ShowerHandler::PDFB_, false, false, true, true, false);

  static Reference<ShowerHandler,PDFBase> interfacePDFARemnant
    ("PDFARemnant",
     "The PDF for beam particle A used to generate forced splittings of the remnant."
     " This overrides both the particle's own PDF setting and the value set by PDFA if used.",
     &ShowerHandler::PDFARemnant_, false, false, true, true, false);

  static Reference<ShowerHandler,PDFBase> interfacePDFBRemnant
    ("PDFBRemnant",
     "The PDF for beam particle B used to generate forced splittings of the remnant."
     " This overrides both the particle's own PDF setting and the value set by PDFB if used.",
     &ShowerHandler::PDFBRemnant_, false, false, true, true, false);

  static Switch<ShowerHandler,bool> interfaceIncludeSpaceTime
    ("IncludeSpaceTime",
     "Whether to include the model for the calculation of space-time distances",
     &ShowerHandler::includeSpaceTime_, false, false, false);
  static SwitchOption interfaceIncludeSpaceTimeYes
    (interfaceIncludeSpaceTime,
     "Yes",
     "Include the model",
     true);
  static SwitchOption interfaceIncludeSpaceTimeNo
    (interfaceIncludeSpaceTime,
     "No",
     "Only include the displacement from the particle-s lifetime for decaying particles",
     false);
  
  static Parameter<ShowerHandler,Energy2> interfaceMinimumVirtuality
    ("MinimumVirtuality",
     "The minimum virtuality for the space-time model",
     &ShowerHandler::vMin_, GeV2, 0.1*GeV2, 0.0*GeV2, 1000.0*GeV2,
     false, false, Interface::limited);

  static Parameter<ShowerHandler,double> interfaceFactorizationScaleFactor
    ("FactorizationScaleFactor",
     "The factorization scale factor.",
     &ShowerHandler::factorizationScaleFactor_, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerHandler,double> interfaceRenormalizationScaleFactor
    ("RenormalizationScaleFactor",
     "The renormalization scale factor.",
     &ShowerHandler::renormalizationScaleFactor_, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerHandler,double> interfaceHardScaleFactor
    ("HardScaleFactor",
     "The hard scale factor.",
     &ShowerHandler::hardScaleFactor_, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerHandler,unsigned int> interfaceMaxTryDecay
    ("MaxTryDecay",
     "The maximum number of attempts to generate a decay",
     &ShowerHandler::maxtryDecay_, 200, 10, 0,
     false, false, Interface::lowerlim);

  static Reference<ShowerHandler,HardScaleProfile> interfaceHardScaleProfile
    ("HardScaleProfile",
     "The hard scale profile to use.",
     &ShowerHandler::hardScaleProfile_, false, false, true, true, false);

  static Switch<ShowerHandler,bool> interfaceMaxPtIsMuF
    ("MaxPtIsMuF",
     "",
     &ShowerHandler::maxPtIsMuF_, false, false, false);
  static SwitchOption interfaceMaxPtIsMuFYes
    (interfaceMaxPtIsMuF,
     "Yes",
     "",
     true);
  static SwitchOption interfaceMaxPtIsMuFNo
    (interfaceMaxPtIsMuF,
     "No",
     "",
     false);

  static Switch<ShowerHandler,bool> interfaceRestrictPhasespace
    ("RestrictPhasespace",
     "Switch on or off phasespace restrictions",
     &ShowerHandler::restrictPhasespace_, true, false, false);
  static SwitchOption interfaceRestrictPhasespaceYes
    (interfaceRestrictPhasespace,
     "Yes",
     "Perform phasespace restrictions",
     true);
  static SwitchOption interfaceRestrictPhasespaceNo
    (interfaceRestrictPhasespace,
     "No",
     "Do not perform phasespace restrictions",
     false);

  static Command<ShowerHandler> interfaceAddVariation
    ("AddVariation",
     "Add a shower variation.",
     &ShowerHandler::doAddVariation, false);
 
  static Switch<ShowerHandler,bool> interfaceDoFSR
    ("DoFSR",
     "Switch on or off final state radiation.",
     &ShowerHandler::doFSR_, true, false, false);
  static SwitchOption interfaceDoFSRYes
    (interfaceDoFSR,
     "Yes",
     "Switch on final state radiation.",
     true);
  static SwitchOption interfaceDoFSRNo
    (interfaceDoFSR,
     "No",
     "Switch off final state radiation.",
     false);

  static Switch<ShowerHandler,bool> interfaceDoISR
    ("DoISR",
     "Switch on or off initial state radiation.",
     &ShowerHandler::doISR_, true, false, false);
  static SwitchOption interfaceDoISRYes
    (interfaceDoISR,
     "Yes",
     "Switch on initial state radiation.",
     true);
  static SwitchOption interfaceDoISRNo
    (interfaceDoISR,
     "No",
     "Switch off initial state radiation.",
     false);
 
  static Switch<ShowerHandler,bool> interfaceSplitHardProcess
    ("SplitHardProcess",
     "Whether or not to try and split the hard process into production and decay processes",
     &ShowerHandler::splitHardProcess_, true, false, false);
  static SwitchOption interfaceSplitHardProcessYes
    (interfaceSplitHardProcess,
     "Yes",
     "Split the hard process",
     true);
  static SwitchOption interfaceSplitHardProcessNo
    (interfaceSplitHardProcess,
     "No",
     "Don't split the hard process",
     false);

  static Switch<ShowerHandler,unsigned int> interfaceSpinCorrelations
    ("SpinCorrelations",
     "Treatment of spin correlations in the parton shower",
     &ShowerHandler::spinOpt_, 1, false, false);
  static SwitchOption interfaceSpinCorrelationsNo
    (interfaceSpinCorrelations,
     "No",
     "No spin correlations",
     0);
  static SwitchOption interfaceSpinCorrelationsSpin
    (interfaceSpinCorrelations,
     "Yes",
     "Include the azimuthal spin correlations",
     1);
  
  static Switch<ShowerHandler,bool> interfaceUseConstituentMasses
  ("UseConstituentMasses",
   "Whether or not to use constituent masses for the reconstruction of the particle after showering.",
   &ShowerHandler::useConstituentMasses_, true, false, false);
  static SwitchOption interfaceUseConstituentMassesYes
  (interfaceUseConstituentMasses,
   "Yes",
   "Use constituent masses.",
   true);
  static SwitchOption interfaceUseConstituentMassesNo
  (interfaceUseConstituentMasses,
   "No",
   "Don't use constituent masses.",
   false);

}

Energy ShowerHandler::hardScale() const {
  assert(false);
  return ZERO;
}

void ShowerHandler::cascade() {
  useMe();
  // Initialise the weights in the event object
  // so that any variations are output regardless of
  // whether showering occurs for the given event
  initializeWeights();
  // get the PDF's from ThePEG (if locally overridden use the local versions)
  tcPDFPtr first  = PDFA_ ? tcPDFPtr(PDFA_) : firstPDF().pdf();
  tcPDFPtr second = PDFB_ ? tcPDFPtr(PDFB_) : secondPDF().pdf();
  resetPDFs(make_pair(first,second));
  // set the PDFs for the remnant
  if( ! rempdfs_.first)
    rempdfs_.first  = PDFARemnant_ ? PDFPtr(PDFARemnant_) : const_ptr_cast<PDFPtr>(first);
  if( ! rempdfs_.second)
    rempdfs_.second = PDFBRemnant_ ? PDFPtr(PDFBRemnant_) : const_ptr_cast<PDFPtr>(second);
  // get the incoming partons
  tPPair  incomingPartons = 
    eventHandler()->currentCollision()->primarySubProcess()->incoming();
  // and the parton bins
  PBIPair incomingBins    = 
    make_pair(lastExtractor()->partonBinInstance(incomingPartons.first),
	      lastExtractor()->partonBinInstance(incomingPartons.second));
  // and the incoming hadrons
  tPPair incomingHadrons = 
    eventHandler()->currentCollision()->incoming();
  remnantDecayer()->setHadronContent(incomingHadrons);
  // check if incoming hadron == incoming parton
  // and get the incoming hadron if exists or parton otherwise
  incoming_ = make_pair(incomingBins.first  ? 
			incomingBins.first ->particle() : incomingPartons.first,
			incomingBins.second ? 
			incomingBins.second->particle() : incomingPartons.second);
  // check the collision is of the beam particles
  // and if not boost collision to the right frame
  // i.e. the hadron-hadron CMF of the collision
  bool btotal(false);
  LorentzRotation rtotal;
  if(incoming_.first  != incomingHadrons.first ||
     incoming_.second != incomingHadrons.second ) {
    btotal = true;
    boostCollision(false);
  }
  // set the current ShowerHandler
  setCurrentHandler();
  // first shower the hard process
  try {
    SubProPtr sub = eventHandler()->currentCollision()->primarySubProcess();
    incomingPartons = cascade(sub,lastXCombPtr());
  }
  catch(ShowerTriesVeto &veto){
    throw Exception() << "Failed to generate the shower after "
                      << veto.tries
                      << " attempts in ShowerHandler::cascade()"
                      << Exception::eventerror;
  }
  if(showerHardProcessVeto()) throw Veto();
  // if a non-hadron collision return (both incoming non-hadronic)
  if( ( !incomingBins.first||
        !isResolvedHadron(incomingBins.first ->particle()))&&
      ( !incomingBins.second||
        !isResolvedHadron(incomingBins.second->particle()))) {
    // boost back to lab if needed
    if(btotal) boostCollision(true);
    // perform the reweighting for the hard process shower
    combineWeights();
    // unset the current ShowerHandler
    unSetCurrentHandler();
    return;
  }
  // get the remnants for hadronic collision
  pair<tRemPPtr,tRemPPtr> remnants(getRemnants(incomingBins));
  // set the starting scale of the forced splitting to the PDF freezing scale
  remnantDecayer()->initialize(remnants, incoming_, *currentStep(), pdfFreezingScale());
  // do the first forcedSplitting
  try {
    remnantDecayer()->doSplit(incomingPartons, make_pair(rempdfs_.first,rempdfs_.second), true);
  }
  catch (ExtraScatterVeto) {
    throw Exception() << "Remnant extraction failed in "
                      << "ShowerHandler::cascade() from primary interaction" 
                      << Exception::eventerror;   
  }
  // perform the reweighting for the hard process shower
  combineWeights();
  // if no MPI return
  if( !isMPIOn() ) {
    remnantDecayer()->finalize();
    // boost back to lab if needed
    if(btotal) boostCollision(true);
    // unset the current ShowerHandler
    unSetCurrentHandler();
    return;
  }
  // generate the multiple scatters use modified pdf's now:
  setMPIPDFs();
  // additional "hard" processes
  unsigned int tries(0);
  // This is the loop over additional hard scatters (most of the time
  // only one, but who knows...)
  for(unsigned int i=1; i <= getMPIHandler()->additionalHardProcs(); i++){
    //counter for regeneration
	unsigned int multSecond = 0;
    // generate the additional scatters
    while( multSecond < getMPIHandler()->multiplicity(i) ) {
      // generate the hard scatter 
      tStdXCombPtr lastXC = getMPIHandler()->generate(i);
      SubProPtr sub = lastXC->construct();
      // add to the Step
      newStep()->addSubProcess(sub);
      // increment the counters
      tries++;
      multSecond++;
      if(tries == maxtryDP_)
	throw Exception() << "Failed to establish the requested number " 
			  << "of additional hard processes. If this error "
			  << "occurs often, your selection of additional "
			  << "scatter is probably unphysical"
			  << Exception::eventerror;
      // Generate the shower. If not possible veto the event
      try {
	incomingPartons = cascade(sub,lastXC);
      } 
      catch(ShowerTriesVeto &veto){
	throw Exception() << "Failed to generate the shower of " 
			  << "a secondary hard process after "
			  << veto.tries
			  << " attempts in Evolver::showerHardProcess()"
			  << Exception::eventerror;
      }
      try {
	// do the forcedSplitting
	remnantDecayer()->doSplit(incomingPartons, make_pair(remmpipdfs_.first,remmpipdfs_.second), false);
      } 
      catch(ExtraScatterVeto){
	//remove all particles associated with the subprocess
	newStep()->removeParticle(incomingPartons.first);
	newStep()->removeParticle(incomingPartons.second);
	//remove the subprocess from the list
	newStep()->removeSubProcess(sub);
	//regenerate the scattering
	multSecond--;
	continue;
      }
      // connect with the remnants but don't set Remnant colour,
      // because that causes problems due to the multiple colour lines.
      if ( !remnants.first ->extract(incomingPartons.first , false) ||
	   !remnants.second->extract(incomingPartons.second, false) )
	throw Exception() << "Remnant extraction failed in "
			  << "ShowerHandler::cascade() for additional scatter" 
			  << Exception::runerror;
    }
    // perform the reweighting for the additional hard scatter shower
    combineWeights();
  }
  // the underlying event processes
  unsigned int ptveto(1), veto(0);
  unsigned int max(getMPIHandler()->multiplicity());
  for(unsigned int i=0; i<max; i++) {
    // check how often this scattering has been regenerated
    if(veto > maxtryMPI_) break;
    //generate PSpoint
    tStdXCombPtr lastXC = getMPIHandler()->generate();
    SubProPtr sub = lastXC->construct();
    //If Algorithm=1 additional scatters of the signal type
    // with pt > ptmin have to be vetoed
    //with probability 1/(m+1), where m is the number of occurances in this event
    if( getMPIHandler()->Algorithm() == 1 ){
      //get the pT
      Energy pt = sub->outgoing().front()->momentum().perp();
      if(pt > getMPIHandler()->PtForVeto() && UseRandom::rnd() < 1./(ptveto+1) ){
        ptveto++;
        i--;
        continue;
      } 
    }
    // add to the SubProcess to the step
    newStep()->addSubProcess(sub);
    // Run the Shower. If not possible veto the scattering
    try {
      incomingPartons = cascade(sub,lastXC);
    } 
    // discard this extra scattering, but try the next one
    catch(ShowerTriesVeto) {
      newStep()->removeSubProcess(sub);
      //regenerate the scattering
      veto++;
      i--;
      continue;      
    }
    try{
      //do the forcedSplitting
      remnantDecayer()->doSplit(incomingPartons, make_pair(remmpipdfs_.first,remmpipdfs_.second), false);
    }
    catch (ExtraScatterVeto) {
      //remove all particles associated with the subprocess
      newStep()->removeParticle(incomingPartons.first);
      newStep()->removeParticle(incomingPartons.second);
      //remove the subprocess from the list
      newStep()->removeSubProcess(sub);
      //regenerate the scattering
      veto++;
      i--;
      continue;      
    }
    //connect with the remnants but don't set Remnant colour,
    //because that causes problems due to the multiple colour lines.
    if ( !remnants.first ->extract(incomingPartons.first , false) ||
	 !remnants.second->extract(incomingPartons.second, false) )
      throw Exception() << "Remnant extraction failed in "
			<< "ShowerHandler::cascade() for MPI hard scattering" 
			<< Exception::runerror;
    //reset veto counter
    veto = 0;
    // perform the reweighting for the MPI process shower
    combineWeights();
  }
  // finalize the remnants
  remnantDecayer()->finalize(getMPIHandler()->colourDisrupt(), 
		    getMPIHandler()->softMultiplicity());
  // boost back to lab if needed
  if(btotal) boostCollision(true);
  // unset the current ShowerHandler
  unSetCurrentHandler();
  getMPIHandler()->clean();
  resetPDFs(make_pair(first,second));
}

void ShowerHandler::initializeWeights() {
  if ( !showerVariations().empty() ) {

    tEventPtr event = eventHandler()->currentEvent();

    for ( map<string,ShowerVariation>::const_iterator var =
	    showerVariations().begin();
	  var != showerVariations().end(); ++var ) {

      // Check that this is behaving as intended
      //map<string,double>::iterator wi = event->optionalWeights().find(var->first);
      //assert(wi == event->optionalWeights().end() ); 

      event->optionalWeights()[var->first] = 1.0;
      currentWeights_[var->first] = 1.0;
    }
  }
  reweight_ = 1.0;
}

void ShowerHandler::resetWeights() {
  for ( map<string,double>::iterator w = currentWeights_.begin();
	w != currentWeights_.end(); ++w ) {
    w->second = 1.0;
  }
  reweight_ = 1.0;
}

void ShowerHandler::combineWeights() {
  tEventPtr event = eventHandler()->currentEvent();
  for ( map<string,double>::const_iterator w = 
	  currentWeights_.begin(); w != currentWeights_.end(); ++w ) {
    map<string,double>::iterator ew = event->optionalWeights().find(w->first);
    if ( ew != event->optionalWeights().end() )
      ew->second *= w->second;
    else {
      assert(false && "Weight name unknown.");
      //event->optionalWeights()[w->first] = w->second;
    }
  }
  if ( reweight_ != 1.0 ) {
    Ptr<StandardEventHandler>::tptr eh = 
      dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(eventHandler());
    if ( !eh ) {
      throw Exception() << "ShowerHandler::combineWeights() : Cross section reweighting "
			<< "through the shower is currently only available with standard "
			<< "event generators" << Exception::runerror;
    }
    eh->reweight(reweight_);
  }
}

string ShowerHandler::doAddVariation(string in) {
  if ( in.empty() )
    return "expecting a name and a variation specification";
  string name = StringUtils::car(in);
  ShowerVariation var;
  string res = var.fromInFile(StringUtils::cdr(in));
  if ( res.empty() ) {
    if ( !var.firstInteraction && !var.secondaryInteractions ) {
      // TODO what about decay showers?
      return "variation does not apply to any shower";
    }
    if ( var.renormalizationScaleFactor == 1.0 && 
	 var.factorizationScaleFactor == 1.0 ) {
      return "variation does not vary anything";
    }
    /*
    Repository::clog() << "adding a variation with tag '" << name << "' using\nxir = "
		       << var.renormalizationScaleFactor
		       << " xif = "
		       << var.factorizationScaleFactor
		       << "\napplying to:\n"
		       << "first interaction = " << var.firstInteraction << " "
		       << "secondary interactions = " << var.secondaryInteractions << "\n"
		       << flush;
    */
    showerVariations()[name] = var;
  }
  return res;
}

tPPair ShowerHandler::cascade(tSubProPtr, XCPtr) {
  assert(false);
  return tPPair();
}

ShowerHandler::RemPair 
ShowerHandler::getRemnants(PBIPair incomingBins) {
  RemPair remnants;
  // first beam particle
  if(incomingBins.first&&!incomingBins.first->remnants().empty()) {
    remnants.first  =
      dynamic_ptr_cast<tRemPPtr>(incomingBins.first->remnants()[0] );
    if(remnants.first) {
      ParticleVector children=remnants.first->children();
      for(unsigned int ix=0;ix<children.size();++ix) {
	if(children[ix]->dataPtr()==remnants.first->dataPtr()) 
	  remnants.first = dynamic_ptr_cast<RemPPtr>(children[ix]);
      } 
      //remove existing colour lines from the remnants
      if(remnants.first->colourLine()) 
	remnants.first->colourLine()->removeColoured(remnants.first);
      if(remnants.first->antiColourLine()) 
	remnants.first->antiColourLine()->removeAntiColoured(remnants.first);
    }
  }
  // seconnd beam particle
  if(incomingBins.second&&!incomingBins. second->remnants().empty()) {
    remnants.second = 
      dynamic_ptr_cast<tRemPPtr>(incomingBins.second->remnants()[0] );
    if(remnants.second) {
      ParticleVector children=remnants.second->children();
      for(unsigned int ix=0;ix<children.size();++ix) {
	if(children[ix]->dataPtr()==remnants.second->dataPtr()) 
	  remnants.second = dynamic_ptr_cast<RemPPtr>(children[ix]);
      } 
      //remove existing colour lines from the remnants
      if(remnants.second->colourLine()) 
	remnants.second->colourLine()->removeColoured(remnants.second);
      if(remnants.second->antiColourLine()) 
	remnants.second->antiColourLine()->removeAntiColoured(remnants.second);
    }
  }
  assert(remnants.first || remnants.second);
  return remnants;
}

namespace {

void addChildren(tPPtr in,set<tPPtr> & particles) {
  particles.insert(in);
  for(unsigned int ix=0;ix<in->children().size();++ix)
    addChildren(in->children()[ix],particles);
}
}

void ShowerHandler::boostCollision(bool boost) {
  // calculate boost from lab to rest
  if(!boost) {
    Lorentz5Momentum ptotal=incoming_.first ->momentum()+incoming_.second->momentum();
    boost_ = LorentzRotation(-ptotal.boostVector());
    Axis axis((boost_*incoming_.first ->momentum()).vect().unit());
    if(axis.perp2()>0.) {
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      boost_.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    }
  }
  // first call performs the boost and second inverse
  // get the particles to be boosted 
  set<tPPtr> particles;
  addChildren(incoming_.first,particles);
  addChildren(incoming_.second,particles);
  // apply the boost
  for(set<tPPtr>::const_iterator cit=particles.begin();
      cit!=particles.end();++cit) {
    (*cit)->transform(boost_);
  }
  if(!boost) boost_.invert();
}

void ShowerHandler::setMPIPDFs() {

  if ( !mpipdfs_.first ) {
    // first have to check for MinBiasPDF
    tcMinBiasPDFPtr first = dynamic_ptr_cast<tcMinBiasPDFPtr>(firstPDF().pdf());
    if(first)
      mpipdfs_.first = new_ptr(MPIPDF(first->originalPDF()));
    else
      mpipdfs_.first = new_ptr(MPIPDF(firstPDF().pdf()));
  }

  if ( !mpipdfs_.second ) {
    tcMinBiasPDFPtr second = dynamic_ptr_cast<tcMinBiasPDFPtr>(secondPDF().pdf());
    if(second)
      mpipdfs_.second = new_ptr(MPIPDF(second->originalPDF()));
    else
      mpipdfs_.second = new_ptr(MPIPDF(secondPDF().pdf()));
  }

  if( !remmpipdfs_.first ) {
    tcMinBiasPDFPtr first = dynamic_ptr_cast<tcMinBiasPDFPtr>(rempdfs_.first);
    if(first)
      remmpipdfs_.first = new_ptr(MPIPDF(first->originalPDF()));
    else
      remmpipdfs_.first = new_ptr(MPIPDF(rempdfs_.first));
  }

  if( !remmpipdfs_.second ) {
    tcMinBiasPDFPtr second = dynamic_ptr_cast<tcMinBiasPDFPtr>(rempdfs_.second);
    if(second)
      remmpipdfs_.second = new_ptr(MPIPDF(second->originalPDF()));
    else
      remmpipdfs_.second = new_ptr(MPIPDF(rempdfs_.second));
  }
  // reset the PDFs stored in the base class
  resetPDFs(mpipdfs_);

}

bool ShowerHandler::isResolvedHadron(tPPtr particle) {
  if(!HadronMatcher::Check(particle->data())) return false;
  for(unsigned int ix=0;ix<particle->children().size();++ix) {
    if(particle->children()[ix]->id()==ParticleID::Remnant) return true;
  }
  return false;
}

namespace {

bool decayProduct(tSubProPtr subProcess,
		  tPPtr particle) {
  // must be time-like and not incoming
  if(particle->momentum().m2()<=ZERO||
     particle == subProcess->incoming().first||
     particle == subProcess->incoming().second) return false;
  // if only 1 outgoing and this is it
  if(subProcess->outgoing().size()==1 &&
     subProcess->outgoing()[0]==particle) return true;
  // must not be the s-channel intermediate otherwise
  if(find(subProcess->incoming().first->children().begin(),
	  subProcess->incoming().first->children().end(),particle)!=
     subProcess->incoming().first->children().end()&&
     find(subProcess->incoming().second->children().begin(),
	  subProcess->incoming().second->children().end(),particle)!=
     subProcess->incoming().second->children().end()&&
     subProcess->incoming().first ->children().size()==1&&
     subProcess->incoming().second->children().size()==1)
    return false;
  // if non-coloured this is enough
  if(!particle->dataPtr()->coloured()) return true;
  // if coloured must be unstable
  if(particle->dataPtr()->stable()) return false;
  // must not have same particle type as a child
  int id = particle->id();
  for(unsigned int ix=0;ix<particle->children().size();++ix)
    if(particle->children()[ix]->id()==id) return false;
  // otherwise its a decaying particle
  return true;
}

PPtr findParent(PPtr original, bool & isHard, 
			       set<PPtr> outgoingset,
			       tSubProPtr subProcess) {
  PPtr parent=original;
  isHard |=(outgoingset.find(original) != outgoingset.end());
  if(!original->parents().empty()) {
    PPtr orig=original->parents()[0];
    if(decayProduct(subProcess,orig))
      parent=findParent(orig,isHard,outgoingset,subProcess);
  }
  return parent;
}
}

void ShowerHandler::findDecayProducts(PPtr in,PerturbativeProcessPtr hard,
				      DecayProcessMap & decay) const {
  ParticleVector children=in->children();
  for(ParticleVector::const_iterator it=children.begin(); it!=children.end();++it) {
    // if decayed or should be decayed in shower make the PerturbaitveProcess
    bool radiates = false;
    if(!(**it).children().empty()) {
      // remove d,u,s,c,b quarks and leptons other than on-shell taus
      if( StandardQCDPartonMatcher::Check((**it).id()) ||
	  ( LeptonMatcher::Check((**it).id()) && !(abs((**it).id())==ParticleID::tauminus &&
						   abs((**it).mass()-(**it).dataPtr()->mass())<MeV))) {
	radiates = true;
      }
      else {
 	bool foundParticle(false),foundGauge(false);
 	for(unsigned int iy=0;iy<(**it).children().size();++iy) {
 	  if((**it).children()[iy]->id()==(**it).id()) {
	    foundParticle = true;
	  }
	  else if((**it).children()[iy]->id()==ParticleID::g ||
		  (**it).children()[iy]->id()==ParticleID::gamma) {
	    foundGauge = true;
	  }
 	}
 	radiates = foundParticle && foundGauge;
      }
    }
    if(radiates) {
      findDecayProducts(*it,hard,decay);
    }
    else if(!(**it).children().empty()||
 	    (decaysInShower((**it).id())&&!(**it).dataPtr()->stable())) {
      createDecayProcess(*it,hard,decay);
    }
    else {
      hard->outgoing().push_back(make_pair(*it,PerturbativeProcessPtr()));
    }
  }
}

void ShowerHandler::splitHardProcess(tPVector tagged, PerturbativeProcessPtr & hard,
				     DecayProcessMap & decay) const {
  // temporary storage of the particles
  set<PPtr> hardParticles;
  // tagged particles in a set
  set<PPtr> outgoingset(tagged.begin(),tagged.end());
  bool isHard=false;
  // loop over the tagged particles
  for (tParticleVector::const_iterator taggedP = tagged.begin();
       taggedP != tagged.end(); ++taggedP) {
    // skip remnants
    if (eventHandler()->currentCollision()&&
        eventHandler()->currentCollision()->isRemnant(*taggedP)) continue;
    // find the parent and whether its a decaying particle
    bool isDecayProd=false;
    // check if hard
    isHard |=(outgoingset.find(*taggedP) != outgoingset.end());
    if(splitHardProcess_) {
      tPPtr parent = *taggedP;
      // check if from s channel decaying colourless particle
      while(parent&&!parent->parents().empty()&&!isDecayProd) {
	parent = parent->parents()[0];
 	if(parent == subProcess_->incoming().first ||
 	   parent == subProcess_->incoming().second ) break;
 	isDecayProd = decayProduct(subProcess_,parent);
      }
      if (isDecayProd)
	hardParticles.insert(findParent(parent,isHard,outgoingset,subProcess_));
    }
    if (!isDecayProd) 
      hardParticles.insert(*taggedP);
  }
  // there must be something to shower
  if(hardParticles.empty()) 
    throw Exception() << "No particles to shower in "
		      << "ShowerHandler::splitHardProcess()" 
		      << Exception::eventerror;
  // must be a hard process
  if(!isHard)
    throw Exception() << "Starting on decay not yet implemented in "
		      << "ShowerHandler::splitHardProcess()" 
		      << Exception::runerror;
  // create the hard process
  hard = new_ptr(PerturbativeProcess());
  // incoming particles
  hard->incoming().push_back(make_pair(subProcess_->incoming().first ,PerturbativeProcessPtr()));
  hard->incoming().push_back(make_pair(subProcess_->incoming().second,PerturbativeProcessPtr()));
  // outgoing particles
  for(set<PPtr>::const_iterator it=hardParticles.begin();it!=hardParticles.end();++it) {
    // if decayed or should be decayed in shower make the tree
    PPtr orig = *it;
    bool radiates = false;
    if(!orig->children().empty()) {
      // remove d,u,s,c,b quarks and leptons other than on-shell taus
      if( StandardQCDPartonMatcher::Check(orig->id()) ||
	  ( LeptonMatcher::Check(orig->id()) && 
	    !(abs(orig->id())==ParticleID::tauminus && abs(orig->mass()-orig->dataPtr()->mass())<MeV))) {
	radiates = true;
      }
      else {
	bool foundParticle(false),foundGauge(false);
	for(unsigned int iy=0;iy<orig->children().size();++iy) {
	  if(orig->children()[iy]->id()==orig->id()) {
	    foundParticle = true;
	  }
	  else if(orig->children()[iy]->id()==ParticleID::g ||
		  orig->children()[iy]->id()==ParticleID::gamma) {
	    foundGauge = true;
	  }
	}
	radiates = foundParticle && foundGauge;
      }
    }
    if(radiates) {
      findDecayProducts(orig,hard,decay);
    }
    else if(!(**it).children().empty()||
 	    (decaysInShower((**it).id())&&!(**it).dataPtr()->stable())) {
      createDecayProcess(*it,hard,decay);
    }
    else {
      hard->outgoing().push_back(make_pair(*it,PerturbativeProcessPtr()));
    }
  }
}

void ShowerHandler::createDecayProcess(PPtr in,PerturbativeProcessPtr hard, DecayProcessMap & decay) const {
  // there must be an incoming particle
  assert(in);
  // create the new process and connect with the parent
  PerturbativeProcessPtr newDecay=new_ptr(PerturbativeProcess());
  newDecay->incoming().push_back(make_pair(in,hard));
  Energy width=in->dataPtr()->generateWidth(in->mass());
  decay.insert(make_pair(width,newDecay));
  hard->outgoing().push_back(make_pair(in,newDecay));
  // we need to deal with the decay products if decayed
  ParticleVector children = in->children();
  if(!children.empty()) {
    for(ParticleVector::const_iterator it = children.begin();
	it!= children.end(); ++it) {
      // if decayed or should be decayed in shower make the tree
      in->abandonChild(*it);
      bool radiates = false;
      if(!(**it).children().empty()) {
 	if(StandardQCDPartonMatcher::Check((**it).id())||
 	   (LeptonMatcher::Check((**it).id())&& !(abs((**it).id())==ParticleID::tauminus &&
						  abs((**it).mass()-(**it).dataPtr()->mass())<MeV))) {
 	  radiates = true;
 	}
 	else {
 	  bool foundParticle(false),foundGauge(false);
 	  for(unsigned int iy=0;iy<(**it).children().size();++iy) {
 	    if((**it).children()[iy]->id()==(**it).id()) {
 	      foundParticle = true;
 	    }
 	    else if((**it).children()[iy]->id()==ParticleID::g ||
 		    (**it).children()[iy]->id()==ParticleID::gamma) {
 	      foundGauge = true;
 	    }
 	  }
 	  radiates = foundParticle && foundGauge;
 	}
  	// finally assume all non-decaying particles are in this class
	// pr 27/11/15 not sure about this bit
   	// if(!radiates) {
   	//   radiates = !decaysInShower((**it).id());
   	// }
      }
      if(radiates) {
 	findDecayProducts(*it,newDecay,decay);
      }
      else if(!(**it).children().empty()||
	      (decaysInShower((**it).id())&&!(**it).dataPtr()->stable())) {
	createDecayProcess(*it,newDecay,decay);
      }
      else {
	newDecay->outgoing().push_back(make_pair(*it,PerturbativeProcessPtr()));
      }
    }
  }
}

tDMPtr ShowerHandler::decay(PerturbativeProcessPtr process,
			    DecayProcessMap & decayMap,
			    bool radPhotons ) const {
  PPtr parent = process->incoming()[0].first;
  assert(parent);
  if(parent->spinInfo()) parent->spinInfo()->decay(true);
  unsigned int ntry = 0;
  ParticleVector children;
  tDMPtr dm = DMPtr();
  while (true) {
    // exit if fails
    if (++ntry>=maxtryDecay_)
      throw Exception() << "Failed to perform decay in ShowerHandler::decay()"
 			<< " after " << maxtryDecay_
 			<< " attempts for " << parent->PDGName() 
 			<< Exception::eventerror;
    // select decay mode
    dm = parent->data().selectMode(*parent);

    if(!dm) 
      throw Exception() << "Failed to select decay  mode in ShowerHandler::decay()"
			<< "for " << parent->PDGName()
			<< Exception::eventerror;
    if(!dm->decayer()) 
      throw Exception() << "No Decayer for selected decay mode "
 			<< " in ShowerHandler::decay()"
 			<< Exception::runerror;
    // start of try block
    try {
      children = dm->decayer()->decay(*dm, *parent);
      // if no children have another go
      if(children.empty()) continue;

      if(radPhotons){
	// generate radiation in the decay
	tDecayIntegratorPtr hwdec=dynamic_ptr_cast<tDecayIntegratorPtr>(dm->decayer());
	if (hwdec && hwdec->canGeneratePhotons())
	  children = hwdec->generatePhotons(*parent,children);
      }

      // set up parent
      parent->decayMode(dm);
      // add children
      for (unsigned int i = 0, N = children.size(); i < N; ++i ) {
  	children[i]->setLabVertex(parent->labDecayVertex());
	//parent->addChild(children[i]);
      }
      // if succeeded break out of loop
      break;
    }
    catch(Veto) {
    }
  }
  assert(!children.empty());
  for(ParticleVector::const_iterator it = children.begin();
      it!= children.end(); ++it) {
    if(!(**it).children().empty()||
       (decaysInShower((**it).id())&&!(**it).dataPtr()->stable())) {
      createDecayProcess(*it,process,decayMap);
    }
    else {
      process->outgoing().push_back(make_pair(*it,PerturbativeProcessPtr()));
    }
  }
  return dm;
}


// Note: The tag must be constructed from an ordered particle container.
tDMPtr ShowerHandler::findDecayMode(const string & tag) const {
  static map<string,DMPtr> cache;
  map<string,DMPtr>::const_iterator pos = cache.find(tag);

  if ( pos != cache.end() ) 
    return pos->second;

  tDMPtr dm = CurrentGenerator::current().findDecayMode(tag);
  cache[tag] = dm;
  return dm;
}


/**
 *  Operator for the particle ordering
 * @param p1 The first ParticleData object
 * @param p2 The second ParticleData object
 */

bool ShowerHandler::ParticleOrdering::operator() (tcPDPtr p1, tcPDPtr p2) const {
  return abs(p1->id()) > abs(p2->id()) ||
    ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
    ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
}

