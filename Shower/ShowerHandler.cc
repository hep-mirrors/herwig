// -*- C++ -*-
//
// ShowerHandler.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/PartonBinInstance.h"
#include "Herwig++/PDT/StandardMatchers.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Utilities/Throw.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include "Herwig++/PDF/MPIPDF.h"
#include "Herwig++/PDF/MinBiasPDF.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/PDF/HwRemDecayer.h"
#include <cassert>

using namespace Herwig;

ShowerHandler::~ShowerHandler() {}

ShowerHandler * ShowerHandler::currentHandler_ = 0;

void ShowerHandler::doinit() {
  CascadeHandler::doinit();
  // copy particles to decay before showering from input vector to the 
  // set used in the simulation
  particlesDecayInShower_.insert(inputparticlesDecayInShower_.begin(),
				 inputparticlesDecayInShower_.end());
  ShowerTree::_decayInShower = particlesDecayInShower_;
}

IBPtr ShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr ShowerHandler::fullclone() const {
  return new_ptr(*this);
}

ShowerHandler::ShowerHandler() : 
  pdfFreezingScale_(2.5*GeV),
  maxtry_(10),maxtryMPI_(10),maxtryDP_(10), subProcess_() {
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
}

void ShowerHandler::dofinish(){
  CascadeHandler::dofinish();
  if(MPIHandler_) MPIHandler_->finalize();
}

void ShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << evolver_ << remDec_ << ounit(pdfFreezingScale_,GeV) << maxtry_ 
     << maxtryMPI_ << maxtryDP_ << inputparticlesDecayInShower_
     << particlesDecayInShower_ << MPIHandler_;
}

void ShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> evolver_ >> remDec_ >> iunit(pdfFreezingScale_,GeV) >> maxtry_ 
     >> maxtryMPI_ >> maxtryDP_ >> inputparticlesDecayInShower_
     >> particlesDecayInShower_ >> MPIHandler_;  
}

ClassDescription<ShowerHandler> ShowerHandler::initShowerHandler;
// Definition of the static class description member.

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
     "The PDF for beam particle A. Overrides the particle's own PDF setting.",
     &ShowerHandler::PDFA_, false, false, true, true, false);

  static Reference<ShowerHandler,PDFBase> interfacePDFB
    ("PDFB",
     "The PDF for beam particle B. Overrides the particle's own PDF setting.",
     &ShowerHandler::PDFB_, false, false, true, true, false);
  
}

void ShowerHandler::cascade() {
  tcPDFPtr first  = firstPDF().pdf();
  tcPDFPtr second = secondPDF().pdf();

  if ( PDFA_ ) first  = PDFA_;
  if ( PDFB_ ) second = PDFB_;

  resetPDFs(make_pair(first,second));

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
                      << " attempts in Evolver::showerHardProcess()"
                      << Exception::eventerror;
  }
  // if a non-hadron collision return (both incoming non-hadronic)
  if( ( !incomingBins.first||
        !HadronMatcher::Check(incomingBins.first ->particle()->data()))&&
      ( !incomingBins.second||
        !HadronMatcher::Check(incomingBins.second->particle()->data()))) {
    // boost back to lab if needed
    if(btotal) boostCollision(true);
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
    remDec_->doSplit(incomingPartons, make_pair(firstPDF() .pdf(), 
						secondPDF().pdf()), true);
  }
  catch (ExtraScatterVeto) {
    throw Exception() << "Remnant extraction failed in "
                      << "ShowerHandler::cascade() from primary interaction" 
                      << Exception::eventerror;   
  }
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
  // We need newpdf to be in scope through the rest of this function.
  pair <PDFPtr, PDFPtr> newpdf;
  setMPIPDFs(newpdf);
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
	remDec_->doSplit(incomingPartons, make_pair(firstPDF().pdf(), 
						    secondPDF().pdf()), false);
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
      remDec_->doSplit(incomingPartons, make_pair(firstPDF().pdf(), 
						  secondPDF().pdf()), false);
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
  }
  // finalize the remnants
  remDec_->finalize(getMPIHandler()->colourDisrupt(), 
		    getMPIHandler()->softMultiplicity());
  // boost back to lab if needed
  if(btotal) boostCollision(true);
  // unset the current ShowerHandler
  currentHandler_ = 0;
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

void ShowerHandler::findShoweringParticles() {
  // clear the storage
  hard_=ShowerTreePtr();
  decay_.clear();
  done_.clear();
  // temporary storage of the particles
  set<PPtr> hardParticles;
  // outgoing particles from the hard process
  PVector outgoing = currentSubProcess()->outgoing();
  set<PPtr> outgoingset(outgoing.begin(),outgoing.end());
  // loop over the tagged particles
  tPVector thetagged;
  if( firstInteraction() ){
    thetagged = tagged();
  }
  else{
    thetagged.insert(thetagged.end(),
		     outgoing.begin(),outgoing.end());
  }
  bool isHard=false;
  for (tParticleVector::const_iterator 
	 taggedP = thetagged.begin();
       taggedP != thetagged.end(); ++taggedP) {
    // if a remnant don't consider
    if(eventHandler()->currentCollision()->isRemnant(*taggedP))
      continue;
    // find the parent and whether its a colourless s-channel resonance
    bool isDecayProd=false;
    tPPtr parent;
    if(!(*taggedP)->parents().empty()) {
      parent = (*taggedP)->parents()[0];
      // check if from s channel decaying colourless particle
      isDecayProd = decayProduct(parent);
    }
    // add to list of outgoing hard particles if needed
    isHard |=(outgoingset.find(*taggedP) != outgoingset.end());
    if(isDecayProd) hardParticles.insert(findParent(parent,isHard,outgoingset));
    else            hardParticles.insert(*taggedP);
  }
  // there must be something to shower
  if(hardParticles.empty()) 
    throw Exception() << "No particles to shower in "
		      << "ShowerHandler::fillShoweringParticles" 
		      << Exception::eventerror;
  if(!isHard)
    throw Exception() << "Starting on decay not yet implemented in "
		      << "ShowerHandler::findShoweringParticles()" 
		      << Exception::runerror;
  // create the hard process ShowerTree
  ParticleVector out(hardParticles.begin(),hardParticles.end());
  hard_=new_ptr(ShowerTree(currentSubProcess()->incoming(),out, decay_));
  hard_->setParents();
}

tPPair ShowerHandler::cascade(tSubProPtr sub,
			      XCPtr xcomb) {
  // get the current step
  current_ = currentStep();
  // get the current subprocess
  subProcess_ = sub;
  // start of the try block for the whole showering process
  unsigned int countFailures=0;
  while (countFailures<maxtry_) {
    try {
      // find the particles in the hard process and the decayed particles to shower
      findShoweringParticles();
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
      ++countFailures;
    }
  }
  // if loop exited because of too many tries, throw event away
  if (countFailures >= maxtry_) {
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
  if (!HadronMatcher::Check(incoming_.first ->data()) && 
      !HadronMatcher::Check(incoming_.second->data()) )
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

PPtr ShowerHandler::findParent(PPtr original, bool & isHard, 
			       set<PPtr> outgoingset) const {
  PPtr parent=original;
  isHard |=(outgoingset.find(original) != outgoingset.end());
  if(!original->parents().empty()) {
    PPtr orig=original->parents()[0];
    if(current_->find(orig)&&decayProduct(orig)) {
      parent=findParent(orig,isHard,outgoingset);
    }
  }
  return parent;
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

bool ShowerHandler::decayProduct(tPPtr particle) const {
  // must be time-like and not incoming
  if(particle->momentum().m2()<=ZERO||
     particle == currentSubProcess()->incoming().first||
     particle == currentSubProcess()->incoming().second) return false;
  // must not be the s-channel intermediate
  if(find(currentSubProcess()->incoming().first->children().begin(),
	  currentSubProcess()->incoming().first->children().end(),particle)!=
     currentSubProcess()->incoming().first->children().end()&&
     find(currentSubProcess()->incoming().second->children().begin(),
	  currentSubProcess()->incoming().second->children().end(),particle)!=
     currentSubProcess()->incoming().second->children().end()&&
     currentSubProcess()->incoming().first ->children().size()==1&&
     currentSubProcess()->incoming().second->children().size()==1)
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

namespace {

void addChildren(tPPtr in,set<tPPtr> particles) {
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
      boost_.rotate(acos(-axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
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

// DO NOT CHANGE THIS SIGNATURE to return the PDFPtr pair. They go out of scope!
void ShowerHandler::setMPIPDFs(pair <PDFPtr, PDFPtr> & newpdf) {

  // first have to check for MinBiasPDF
  tcMinBiasPDFPtr first = dynamic_ptr_cast<tcMinBiasPDFPtr>(firstPDF().pdf());
  if(first)
    newpdf.first = new_ptr(MPIPDF(first->originalPDF()));
  else
    newpdf.first = new_ptr(MPIPDF(firstPDF().pdf()));


  tcMinBiasPDFPtr second = dynamic_ptr_cast<tcMinBiasPDFPtr>(secondPDF().pdf());
  if(second)
    newpdf.second = new_ptr(MPIPDF(second->originalPDF()));
  else
    newpdf.second = new_ptr(MPIPDF(secondPDF().pdf()));


  // reset the PDFs stored in the base class
  resetPDFs(newpdf);
}
