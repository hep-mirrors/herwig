// -*- C++ -*-
//
// ShowerHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig/Shower/Base/Evolver.h"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Utilities/EnumParticles.h"
#include "Herwig/PDF/MPIPDF.h"
#include "Herwig/PDF/MinBiasPDF.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "Herwig/Shower/Base/ShowerTree.h"
#include "Herwig/Shower/Base/HardTree.h"
#include "Herwig/Shower/Base/KinematicsReconstructor.h"
#include "Herwig/Shower/Base/PartnerFinder.h"
#include "Herwig/PDF/HwRemDecayer.h"
#include <cassert>
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeClass<ShowerHandler,CascadeHandler>
describeShowerHandler ("Herwig::ShowerHandler","HwShower.so");

ShowerHandler::~ShowerHandler() {}

ShowerHandler * ShowerHandler::currentHandler_ = 0;

void ShowerHandler::doinit() {
  CascadeHandler::doinit();
  // copy particles to decay before showering from input vector to the 
  // set used in the simulation
  if ( particlesDecayInShower_.empty() )
    particlesDecayInShower_.insert(inputparticlesDecayInShower_.begin(),
				   inputparticlesDecayInShower_.end());
  ShowerTree::_decayInShower = particlesDecayInShower_;
  ShowerTree::_vmin2 = vMin_;
  ShowerTree::_spaceTime = includeSpaceTime_;
}

IBPtr ShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr ShowerHandler::fullclone() const {
  return new_ptr(*this);
}

ShowerHandler::ShowerHandler() : 
  reweight_(1.0),
  pdfFreezingScale_(2.5*GeV),
  maxtry_(10),maxtryMPI_(10),maxtryDP_(10),
  includeSpaceTime_(false), vMin_(0.1*GeV2), subProcess_(),
  factorizationScaleFactor_(1.0),
  renormalizationScaleFactor_(1.0),
  hardScaleFactor_(1.0),
  restrictPhasespace_(true), maxPtIsMuF_(false),
  splitHardProcess_(true) {
  inputparticlesDecayInShower_.push_back( 6  ); //  top 
  inputparticlesDecayInShower_.push_back( 23 ); // Z0
  inputparticlesDecayInShower_.push_back( 24 ); // W+/-
  inputparticlesDecayInShower_.push_back( 25 ); // h0
}

void ShowerHandler::doinitrun(){
  CascadeHandler::doinitrun();
  //can't use isMPIOn here, because the EventHandler is not set at that stage
  if(MPIHandler_){ 
    MPIHandler_->initialize();
    if(MPIHandler_->softInt())
      remDec_->initSoftInteractions(MPIHandler_->Ptmin(), MPIHandler_->beta());
  }
  ShowerTree::_decayInShower = particlesDecayInShower_;
  ShowerTree::_vmin2 = vMin_;
  ShowerTree::_spaceTime = includeSpaceTime_;
}

void ShowerHandler::dofinish(){
  CascadeHandler::dofinish();
  if(MPIHandler_) MPIHandler_->finalize();
}

void ShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << evolver_ << remDec_ << ounit(pdfFreezingScale_,GeV) << maxtry_ 
     << maxtryMPI_ << maxtryDP_ << inputparticlesDecayInShower_
     << particlesDecayInShower_ << MPIHandler_ << PDFA_ << PDFB_
     << PDFARemnant_ << PDFBRemnant_
     << includeSpaceTime_ << ounit(vMin_,GeV2)
     << factorizationScaleFactor_ << renormalizationScaleFactor_
     << hardScaleFactor_
     << restrictPhasespace_ << maxPtIsMuF_ << hardScaleProfile_
     << splitHardProcess_ << showerVariations_;
}

void ShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> evolver_ >> remDec_ >> iunit(pdfFreezingScale_,GeV) >> maxtry_ 
     >> maxtryMPI_ >> maxtryDP_ >> inputparticlesDecayInShower_
     >> particlesDecayInShower_ >> MPIHandler_ >> PDFA_ >> PDFB_
     >> PDFARemnant_ >> PDFBRemnant_
     >> includeSpaceTime_ >> iunit(vMin_,GeV2)
     >> factorizationScaleFactor_ >> renormalizationScaleFactor_
     >> hardScaleFactor_
     >> restrictPhasespace_ >> maxPtIsMuF_ >> hardScaleProfile_
     >> splitHardProcess_ >> showerVariations_;
}

void ShowerHandler::Init() {

  static ClassDocumentation<ShowerHandler> documentation
    ("Main driver class for the showering.",
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

  static Reference<ShowerHandler,Evolver> 
    interfaceEvolver("Evolver", 
		     "A reference to the Evolver object", 
		     &Herwig::ShowerHandler::evolver_,
		     false, false, true, false);

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
  static SwitchOption interfaceRestrictPhasespaceOn
    (interfaceRestrictPhasespace,
     "On",
     "Perform phasespace restrictions",
     true);
  static SwitchOption interfaceRestrictPhasespaceOff
    (interfaceRestrictPhasespace,
     "Off",
     "Do not perform phasespace restrictions",
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

  static Command<ShowerHandler> interfaceAddVariation
    ("AddVariation",
     "Add a shower variation.",
     &ShowerHandler::doAddVariation, false);
 
}

Energy ShowerHandler::hardScale() const {
  return evolver_->hardScale();
}

void ShowerHandler::cascade() {
  
  // Initialise the weights in the event object
  // so that any variations are output regardless of
  // whether showering occurs for the given event
  initializeWeights();

  tcPDFPtr first  = firstPDF().pdf();
  tcPDFPtr second = secondPDF().pdf();

  if ( PDFA_ ) first  = PDFA_;
  if ( PDFB_ ) second = PDFB_;

  resetPDFs(make_pair(first,second));

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
  remDec_->setHadronContent(incomingHadrons);
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
  currentHandler_ = this;
  // first shower the hard process
  useMe();
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
    currentHandler_ = 0;
    return;
  }
  // get the remnants for hadronic collision
  pair<tRemPPtr,tRemPPtr> remnants(getRemnants(incomingBins));
  // set the starting scale of the forced splitting to the PDF freezing scale
  remDec_->initialize(remnants, incoming_, *currentStep(), pdfFreezingScale());
  // do the first forcedSplitting
  try {
    remDec_->doSplit(incomingPartons, make_pair(rempdfs_.first,rempdfs_.second), true);
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
    remDec_->finalize();
    // boost back to lab if needed
    if(btotal) boostCollision(true);
    // unset the current ShowerHandler
    currentHandler_ = 0;
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
	remDec_->doSplit(incomingPartons, make_pair(remmpipdfs_.first,remmpipdfs_.second), false);
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
      remDec_->doSplit(incomingPartons, make_pair(remmpipdfs_.first,remmpipdfs_.second), false);
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
  remDec_->finalize(getMPIHandler()->colourDisrupt(), 
		    getMPIHandler()->softMultiplicity());
  // boost back to lab if needed
  if(btotal) boostCollision(true);
  // unset the current ShowerHandler
  currentHandler_ = 0;
  getMPIHandler()->clean();
}

void ShowerHandler::fillEventRecord() {
  // create a new step 
  StepPtr pstep = newStep();
  assert(!done_.empty());
  assert(done_[0]->isHard());
  // insert the steps
  for(unsigned int ix=0;ix<done_.size();++ix) {
    done_[ix]->fillEventRecord(pstep,
			       evolver_->isISRadiationON(),
			       evolver_->isFSRadiationON());
  }
}

void ShowerHandler::prepareCascade(tSubProPtr sub) { 
  current_ = currentStep(); 
  subProcess_ = sub;
} 

void ShowerHandler::initializeWeights() {
  if ( !showerVariations().empty() ) {

    tEventPtr event = eventHandler()->currentEvent();

    for ( map<string,ShowerHandler::ShowerVariation>::const_iterator var =
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

string ShowerHandler::ShowerVariation::fromInFile(const string& in) {
  // pretty simple for the moment, just to try
  // TODO make this better
  istringstream read(in);
  string where;
  read >> renormalizationScaleFactor >> factorizationScaleFactor >> where;
  if ( !read )
    return "something went wrong with: " + in;
  if ( where != "Hard" && where != "All" && where!= "Secondary" )
    return "The specified process for reweighting does not exist.\nOptions are: Hard, Secondary, All.";
  if ( where == "Hard" || where == "All" )
    firstInteraction = true;
  else
    firstInteraction = false;
  if ( where == "Secondary" || where == "All" )
    secondaryInteractions = true;
  else
    secondaryInteractions = false;
  return "";
}

void ShowerHandler::ShowerVariation::put(PersistentOStream& os) const {
  os << renormalizationScaleFactor << factorizationScaleFactor
     << firstInteraction << secondaryInteractions;
}

void ShowerHandler::ShowerVariation::get(PersistentIStream& is) {
  is >> renormalizationScaleFactor >> factorizationScaleFactor
     >> firstInteraction >> secondaryInteractions;
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

tPPair ShowerHandler::cascade(tSubProPtr sub,
			      XCPtr xcomb) {
  prepareCascade(sub);
  resetWeights();
  // set the scale variation factors; needs to go after prepareCascade
  // to trigger possible different variations for hard and secondary
  // scatters
  evolver()->renormalizationScaleFactor(renormalizationScaleFactor());
  evolver()->factorizationScaleFactor(factorizationScaleFactor());
  evolver()->restrictPhasespace(restrictPhasespace());
  evolver()->hardScaleIsMuF(hardScaleIsMuF());
  // start of the try block for the whole showering process
  unsigned int countFailures=0;
  while (countFailures<maxtry_) {
    try {
      decay_.clear();
      done_.clear();
      ShowerTree::constructTrees(currentSubProcess(),hard_,decay_,
				 firstInteraction() ? tagged() :
				 tPVector(currentSubProcess()->outgoing().begin(),
					  currentSubProcess()->outgoing().end()),
				 splitHardProcess_);
      // if no hard process
      if(!hard_)  throw Exception() << "Shower starting with a decay"
				    << "is not implemented" 
				    << Exception::runerror;
      // perform the shower for the hard process
      evolver_->showerHardProcess(hard_,xcomb);
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
	decayingTree->decay(decay_);
	// now shower the decay
	evolver_->showerDecay(decayingTree);
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
  }
  // if loop exited because of too many tries, throw event away
  if (countFailures >= maxtry_) {
    resetWeights();
    hard_=ShowerTreePtr();
    decay_.clear();
    done_.clear();
    throw Exception() << "Too many tries for main while loop "
		      << "in ShowerHandler::cascade()." 
		      << Exception::eventerror; 	
  }
  //enter the particles in the event record
  fillEventRecord();
  // clear storage
  hard_=ShowerTreePtr();
  decay_.clear();
  done_.clear();
  // non hadronic case return
  if (!isResolvedHadron(incoming_.first ) && 
      !isResolvedHadron(incoming_.second) )
    return incoming_;
  // remake the remnants (needs to be after the colours are sorted
  //                       out in the insertion into the event record)
  if ( firstInteraction() ) return remakeRemnant(sub->incoming());
  //Return the new pair of incoming partons. remakeRemnant is not
  //necessary here, because the secondary interactions are not yet
  //connected to the remnants.
  return make_pair(findFirstParton(sub->incoming().first ),
		   findFirstParton(sub->incoming().second));
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

tPPair ShowerHandler::remakeRemnant(tPPair oldp){
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

PPtr ShowerHandler::findFirstParton(tPPtr seed) const{
  if(seed->parents().empty()) return seed;
  tPPtr parent = seed->parents()[0];
  //if no parent there this is a loose end which will 
  //be connected to the remnant soon.
  if(!parent || parent == incoming_.first || 
     parent == incoming_.second ) return seed;
  else return findFirstParton(parent);
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

HardTreePtr ShowerHandler::generateCKKW(ShowerTreePtr ) const {
  return HardTreePtr();
}

