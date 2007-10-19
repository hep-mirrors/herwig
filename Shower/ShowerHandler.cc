// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerHandler class.
//

#include "ShowerHandler.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/PartonBinInstance.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Utilities/Throw.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Utilities/EnumParticles.h"
#include "Herwig++/PDF/MPIPDF.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/CKKW/Clustering/CascadeReconstructor.h"
#include "Herwig++/Shower/CKKW/Reweighting/Reweighter.h"
#include <cassert>

using namespace Herwig;

ShowerHandler::~ShowerHandler() {}

ShowerHandler * ShowerHandler::theHandler = 0;

void ShowerHandler::doinit() throw(InitException) {
  CascadeHandler::doinit();
  // copy particles to decay before showering from input vector to the 
  // set used in the simulation
  _particlesDecayInShower.insert(_inputparticlesDecayInShower.begin(),
				 _inputparticlesDecayInShower.end());


  // check for CKKW and setup if present

  if (_reconstructor && _reweighter) {

    _reconstructor->setup();
    _evolver->useCKKW(_reconstructor,_reweighter);
    _useCKKW = true;


  } else {
    _useCKKW = false;
  }

}

IBPtr ShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr ShowerHandler::fullclone() const {
  return new_ptr(*this);
}

ShowerHandler::ShowerHandler() : 
  theOrderSecondaries(true), theMPIOnOff(true), 
  _maxtry(10), theSubProcess(tSubProPtr()) {
  _inputparticlesDecayInShower.push_back( 6 ); //  top
  _inputparticlesDecayInShower.push_back( 1000001 ); //  SUSY_d_L 
  _inputparticlesDecayInShower.push_back( 1000002 ); //  SUSY_u_L 
  _inputparticlesDecayInShower.push_back( 1000003 ); //  SUSY_s_L 
  _inputparticlesDecayInShower.push_back( 1000004 ); //  SUSY_c_L 
  _inputparticlesDecayInShower.push_back( 1000005 ); //  SUSY_b_1 
  _inputparticlesDecayInShower.push_back( 1000006 ); //  SUSY_t_1 
  _inputparticlesDecayInShower.push_back( 1000011 ); //  SUSY_e_Lminus 
  _inputparticlesDecayInShower.push_back( 1000012 ); //  SUSY_nu_eL 
  _inputparticlesDecayInShower.push_back( 1000013 ); //  SUSY_mu_Lminus 
  _inputparticlesDecayInShower.push_back( 1000014 ); //  SUSY_nu_muL 
  _inputparticlesDecayInShower.push_back( 1000015 ); //  SUSY_tau_1minus 
  _inputparticlesDecayInShower.push_back( 1000016 ); //  SUSY_nu_tauL 
  _inputparticlesDecayInShower.push_back( 1000021 ); //  SUSY_g 
  _inputparticlesDecayInShower.push_back( 1000022 ); //  SUSY_chi_10 
  _inputparticlesDecayInShower.push_back( 1000023 ); //  SUSY_chi_20 
  _inputparticlesDecayInShower.push_back( 1000024 ); //  SUSY_chi_1plus 
  _inputparticlesDecayInShower.push_back( 1000025 ); //  SUSY_chi_30 
  _inputparticlesDecayInShower.push_back( 1000035 ); //  SUSY_chi_40 
  _inputparticlesDecayInShower.push_back( 1000037 ); //  SUSY_chi_2plus 
  _inputparticlesDecayInShower.push_back( 1000039 ); //  SUSY_gravitino 
  _inputparticlesDecayInShower.push_back( 2000001 ); //  SUSY_d_R 
  _inputparticlesDecayInShower.push_back( 2000002 ); //  SUSY_u_R 
  _inputparticlesDecayInShower.push_back( 2000003 ); //  SUSY_s_R 
  _inputparticlesDecayInShower.push_back( 2000004 ); //  SUSY_c_R 
  _inputparticlesDecayInShower.push_back( 2000005 ); //  SUSY_b_2 
  _inputparticlesDecayInShower.push_back( 2000006 ); //  SUSY_t_2 
  _inputparticlesDecayInShower.push_back( 2000011 ); //  SUSY_e_Rminus 
  _inputparticlesDecayInShower.push_back( 2000012 ); //  SUSY_nu_eR 
  _inputparticlesDecayInShower.push_back( 2000013 ); //  SUSY_mu_Rminus 
  _inputparticlesDecayInShower.push_back( 2000014 ); //  SUSY_nu_muR 
  _inputparticlesDecayInShower.push_back( 2000015 ); //  SUSY_tau_2minus 
  _inputparticlesDecayInShower.push_back( 2000016 ); //  SUSY_nu_tauR 
  _inputparticlesDecayInShower.push_back( 25      ); //  h0
  _inputparticlesDecayInShower.push_back( 35      ); //  H0
  _inputparticlesDecayInShower.push_back( 36      ); //  A0
  _inputparticlesDecayInShower.push_back( 37      ); //  H+
  _inputparticlesDecayInShower.push_back( 23      ); // Z0
  _inputparticlesDecayInShower.push_back( 24      ); // W+/-
}

void ShowerHandler::doinitrun(){
  CascadeHandler::doinitrun();
  theMPIHandler->initrun();
  if(IsMPIOn() && abs(generator()->eventHandler()->incoming().first->id()) > 99 &&
    abs(generator()->eventHandler()->incoming().second->id()) > 99)
    theMPIHandler->initialize();

  if (_useCKKW) {
    _reweighter->initialize();
  }
}

void ShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << _evolver << theRemDec << _maxtry << _inputparticlesDecayInShower
     << _particlesDecayInShower << theOrderSecondaries 
     << theMPIOnOff << theMPIHandler << theSubProcess 
     << _useCKKW << _reconstructor << _reweighter;
}

void ShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> _evolver >> theRemDec >> _maxtry
     >> _inputparticlesDecayInShower
     >> _particlesDecayInShower >> theOrderSecondaries 
     >> theMPIOnOff >> theMPIHandler >> theSubProcess
     >> _useCKKW >> _reconstructor >> _reweighter;  
}

ClassDescription<ShowerHandler> ShowerHandler::initShowerHandler;
// Definition of the static class description member.

void ShowerHandler::Init() {

  static ClassDocumentation<ShowerHandler> documentation
    ("Main driver class for the showering.");

  static Reference<ShowerHandler,Evolver> 
    interfaceEvolver("Evolver", 
		     "A reference to the Evolver object", 
		     &Herwig::ShowerHandler::_evolver,
		     false, false, true, false);

  static Reference<ShowerHandler,HwRemDecayer> 
    interfaceRemDecayer("RemDecayer", 
		     "A reference to the Remnant Decayer object", 
		     &Herwig::ShowerHandler::theRemDec,
		     false, false, true, false);

  static Parameter<ShowerHandler,unsigned int> interfaceMaxTry
    ("MaxTry",
     "The maximum number of attempts for the main showering loop",
     &ShowerHandler::_maxtry, 10, 1, 100,
     false, false, Interface::limited);

  static ParVector<ShowerHandler,long> interfaceDecayInShower
    ("DecayInShower",
     "PDG codes of the particles to be decayed in the shower",
     &ShowerHandler::_inputparticlesDecayInShower, -1, 0l, -10000000l, 10000000l,
     false, false, Interface::limited);

  static Reference<ShowerHandler,MPIHandler> interfaceMPIHandler
    ("MPIHandler",
     "The object that admisinsters all additional semihard partonic scatterings.",
     &ShowerHandler::theMPIHandler, false, false, true, false);
  
  
  static Switch<ShowerHandler,bool> interfaceMPIOnOff
    ("MPIOnOff", "flag to switch MPI on or off",
     &ShowerHandler::theMPIOnOff, 1, false, false);

  static SwitchOption interfaceMPIOnOff0                             
    (interfaceMPIOnOff,"MPI-OFF","Multiple parton interactions are OFF", 0);
  static SwitchOption interfaceMPIOnOff1                            
    (interfaceMPIOnOff,"MPI-ON","Multiple parton interactions are ON", 1);

  static Switch<ShowerHandler,bool> interfaceOrderSecondaries
    ("OrderSecondaries", 
     "flag to switch the ordering of the additional interactions on or off",
     &ShowerHandler::theOrderSecondaries, 1, false, false);

  static SwitchOption interfaceOrderSecondaries0                             
    (interfaceOrderSecondaries,"Order-OFF","Multiple parton interactions aren't ordered", 0);
  static SwitchOption interfaceOrderSecondaries1                            
    (interfaceOrderSecondaries,"Order-ON",
     "Multiple parton interactions are ordered according to their scale", 1);

  static Reference<ShowerHandler,CascadeReconstructor> interfaceCascadeReconstructor
    ("CascadeReconstructor",
     "Casacde reconstructor used for ME/PS merging.",
     &ShowerHandler::_reconstructor, false, false, true, true, false);

  static Reference<ShowerHandler,Reweighter> interfaceReweighter
    ("Reweighter",
     "Reweighter used for ME/PS merging.",
     &ShowerHandler::_reweighter, false, false, true, true, false);
}

void ShowerHandler::cascade() {
  theHandler = this;
  tStdXCombPtr lastXC;
  SubProPtr sub;
  tPPair incs;
  Energy2 scale;
  multimap <Energy2, SubProPtr> procs;
  multimap <Energy2, SubProPtr> :: const_iterator pit;
  unsigned int i(0);

  lastXC = dynamic_ptr_cast<StdXCombPtr>(lastXCombPtr());
  sub = eventHandler()->currentCollision()->primarySubProcess();

  //first shower the hard process
  incs = cascade(sub);

  PBIPair incbins = make_pair(lastExtractor()->partonBinInstance(incs.first),
			      lastExtractor()->partonBinInstance(incs.second));

  if (abs(eventHandler()->incoming().first->id()) < 99 ||
    abs(eventHandler()->incoming().second->id()) < 99 )
    return;
 
  pair<tRemPPtr,tRemPPtr> remnants(getRemnants(incbins));

  theRemDec->initialize(remnants, *currentStep());

  //do the first forcedSplitting
  try{
    theRemDec->doSplit(incs, true);
  }catch(ExtraScatterVeto){
    throw Exception() << "Remnant extraction failed in "
                      << "ShowerHandler::cascade()" 
                      << Exception::eventerror;   
  }

  if( !IsMPIOn() ){
    theRemDec->finalize();
    return;
  }

  //use modified pdf's now:
  //const pair <tcPDFPtr, tcPDFPtr> newpdf = 
  //  make_pair(new_ptr(MPIPDF(firstPDF().pdf())), 
  //	      new_ptr(MPIPDF(secondPDF().pdf())));
  //resetPDFs(newpdf);

  int veto(1);
  unsigned int max(theMPIHandler->multiplicity());
  for(i=0; i<max; i++){      
    //generate PSpoint
    lastXC = theMPIHandler->generate();
    sub = lastXC->construct();

    //If Jmueo=1 additional scatters of the signal type with pt > ptmin have to be vetoed
    //with probability 1/(m+1), where m is the number of occurances in this event

    //check if the same process is used for the signal and UE
    //For LesHouches event files the MEBasePtr should be 0
    //That leads to the correct behaviour as long as no QCD2->2 event is read in
    if(sub->handler() == subProcess()->handler() && theMPIHandler->Jmueo() ){
      //get the pT
      Energy pt = sub->outgoing().front()->momentum().perp();
      Energy ptmin = lastCutsPtr()->minKT(sub->outgoing().front()->dataPtr());

      if(pt > ptmin && UseRandom::rnd() < 1./(veto+1) ){
        veto++;
        i--;
        continue;
      } 
    }
    //sort in -scale, because reverse iterator doesn't work with gcc3.x.x
    if( IsOrdered() ) scale = -lastXC->lastScale();
    else scale = 1.*GeV2;

    procs.insert(make_pair(scale, sub));
  }
  
  for( pit=procs.begin(); pit!=procs.end(); ++pit ){
    // cerr << "scale: " << sqrt(-pit->first/GeV2) << endl;
    //add to the EventHandler's list
    newStep()->addSubProcess(pit->second);
    //start the Shower
    incs = cascade(pit->second);

    try{
      //cerr << "do extra scatter forced splitting\n";
      //do the forcedSplitting
      theRemDec->doSplit(incs, false);

      //check if there is enough energy to extract
      if( (remnants.first->momentum() - incs.first->momentum()).e() < 0*MeV ||
	  (remnants.second->momentum() - incs.second->momentum()).e() < 0*MeV )
	throw ExtraScatterVeto();
    }catch(ExtraScatterVeto){
      //cerr << "remove scatter\n";
      //remove all particles associated with the subprocess
      newStep()->removeParticle(incs.first);
      newStep()->removeParticle(incs.second);
      //remove the subprocess from the list
      newStep()->removeSubProcess(pit->second);
      //discard this extra scattering, but try the next one
      continue;
    }
    //connect with the remnants but don't set Remnant colour,
    //because that causes problems due to the multiple colour lines.
    if ( !remnants.first->extract(incs.first, false) ||
	 !remnants.second->extract(incs.second, false) )
      throw Exception() << "Remnant extraction failed in "
			<< "ShowerHandler::cascade()" 
			<< Exception::runerror;   
  } 

  theRemDec->finalize();

  theHandler = 0;
}

void ShowerHandler::fillEventRecord() {
  // create a new step 
  StepPtr pstep = newStep();
  if(_done.empty()) throw Exception() << "Must have some showers to insert in "
				      << "ShowerHandler::fillEventRecord()" 
				      << Exception::runerror;
  if(!_done[0]->isHard()) throw Exception() << "Must start filling with hard process"
					    << " in ShowerHandler::fillEventRecord()" 
					    << Exception::runerror;
  // insert the steps
  for(unsigned int ix=0;ix<_done.size();++ix) {
    _done[ix]->fillEventRecord(pstep,
			       _evolver->isISRadiationON(),
			       _evolver->isFSRadiationON());
  }
} 

void ShowerHandler::findShoweringParticles() {
  // clear the storage
  _hard=ShowerTreePtr();
  _decay.clear();
  _done.clear();
  // temporary storage of the particles
  set<PPtr> hardParticles;
  // outgoing particles from the hard process
  PVector outgoing = currentSubProcess()->outgoing();

  set<PPtr> outgoingset(outgoing.begin(),outgoing.end());
  // loop over the tagged particles
  tPVector thetagged;
  if( FirstInt() ){
    thetagged = tagged();
  }
  else{
    //get the "tagged" particles 
    for(PVector::const_iterator pit = currentSubProcess()->outgoing().begin(); 
	pit != currentSubProcess()->outgoing().end(); ++pit)
      thetagged.push_back(*pit);
  }
  tParticleVector::const_iterator taggedP = thetagged.begin();
  bool isHard=false;
  for (;taggedP != thetagged.end(); ++taggedP) {
    // if a remnant don't consider
    if(eventHandler()->currentCollision()->isRemnant(*taggedP))
      continue;
    // find the parent and if colourless s-channel resonance
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
  _hard=new_ptr(ShowerTree(out, _decay));
  _hard->setParents();
}

tPPair ShowerHandler::
cascade(tSubProPtr sub) {
  // set the current step
  _current=currentStep();
  // set the current subprocess
  theSubProcess = sub;
  //  start of the try block for the whole showering process
  unsigned int countFailures=0;
  ShowerTreePtr hard;
  vector<ShowerTreePtr> decay;
  while (countFailures<_maxtry) {
    try {
      // find the particles in the hard process and the decayed particles to shower
      findShoweringParticles();
      // check if a hard process or decay
      bool isHard = _hard;
      // if a hard process perform the shower for the hard process
      if(isHard) {
	_evolver->showerHardProcess(_hard);
	_done.push_back(_hard);
	_hard->updateAfterShower(_decay);
      }
      // if no decaying particles to shower break out of the loop
      if(_decay.empty()) break;
      // if no hard process
      if(!isHard) 
	throw Exception() << "Shower starting with a decay is not yet implemented" 
			  << Exception::runerror;
      // shower the decay products
      while(!_decay.empty()) {
	multimap<Energy,ShowerTreePtr>::iterator dit=--_decay.end();
	while(!dit->second->parent()->hasShowered() && dit!=_decay.begin()) --dit;
	// get the particle and the width
	ShowerTreePtr decayingTree = dit->second;
	// 	    Energy largestWidthDecayingSystem=(*_decay.rbegin()).first;
	// remove it from the multimap
	_decay.erase(dit);
	// make sure the particle has been decayed
	decayingTree->decay(_decay);
	// now shower the decay
	_evolver->showerDecay(decayingTree);
	_done.push_back(decayingTree);
	decayingTree->updateAfterShower(_decay);
      }
      // suceeded break out of the loop
      break;
    }
    catch (KinematicsReconstructionVeto) {
      ++countFailures;
    }
  }
  // if loop exited because of too many tries, throw event away
  if (countFailures >= _maxtry) {
    throw Exception() << "Too many tries for main while loop "
		      << "in ShowerHandler::cascade()." 
		      << Exception::eventerror; 	
  }
  //enter the particles in the event record
  fillEventRecord();
  //non hadronic case:
  if ( abs(eventHandler()->incoming().first->id()) < 99 ||
    abs(eventHandler()->incoming().second->id()) < 99 ) 
    return eventHandler()->currentCollision()->incoming();

  // remake the remnants (needs to be after the colours are sorted
  //                       out in the insertion into the event record)
  if ( FirstInt() ) return remakeRemnant(sub->incoming());

  //Return the new pair of incoming partons. remakeRemnant is not
  //necessary here, because the secondary interactions are not yet
  //connected to the remnants.
  tPPair inc = generator()->currentEvent()->incoming();
  return make_pair(findFirstParton(sub->incoming().first, inc),
		   findFirstParton(sub->incoming().second, inc));
}

PPtr ShowerHandler::findParent(PPtr original, bool & isHard, 
			       set<PPtr> outgoingset) const {
  PPtr parent=original;
  isHard |=(outgoingset.find(original) != outgoingset.end());
  if(!original->parents().empty()) {
    PPtr orig=original->parents()[0];
    if(_current->find(orig)&&decayProduct(orig)) {
      parent=findParent(orig,isHard,outgoingset);
    }
  }
  return parent;
}

ShowerHandler::RemPair 
ShowerHandler::getRemnants(PBIPair incbins){
  RemPair remnants;
  if ( incbins.first->remnants().size() != 1 ||
       incbins.second->remnants().size() != 1 )
    throw Exception() << "Wrong number of Remnants "
		      << "in ShowerHandler::getRemnants()." 
		      << Exception::runerror; 	
    
  remnants.first = dynamic_ptr_cast<tRemPPtr>(incbins.first->remnants()[0]);
  remnants.second = dynamic_ptr_cast<tRemPPtr>(incbins.second->remnants()[0]);
  
  if(remnants.first && remnants.second){
      //remove existing colour lines from the remnants
    if(remnants.first->colourLine()) 
      remnants.first->colourLine()->removeColoured(remnants.first);
    if(remnants.first->antiColourLine()) 
      remnants.first->antiColourLine()->removeAntiColoured(remnants.first);
    if(remnants.second->colourLine()) 
      remnants.second->colourLine()->removeColoured(remnants.second);
    if(remnants.second->antiColourLine()) 
      remnants.second->antiColourLine()->removeAntiColoured(remnants.second);

    //copy the remnants to the current step, as they may be changed now
    if ( remnants.first->birthStep() != newStep() ) {
      RemPPtr newrem = new_ptr(*remnants.first);
      newStep()->addDecayProduct(remnants.first, newrem, false);
      remnants.first = newrem;
    }
    if ( remnants.second->birthStep() != newStep() ) {
      RemPPtr newrem = new_ptr(*remnants.second);
      newStep()->addDecayProduct(remnants.second, newrem, false);
      remnants.second = newrem;
    }
    return remnants;
  }else{
    throw Exception() << "Remnants are not accessable"
		      << "in ShowerHandler::getRemnants()." 
		      << Exception::runerror; 	    
  }
}

tPPair ShowerHandler::remakeRemnant(tPPair oldp){                     
  PartonExtractor & pex = *lastExtractor();
  tPPair inc = generator()->currentEvent()->incoming();
  
  tPPair newp = make_pair(findFirstParton(oldp.first, inc), findFirstParton(oldp.second, inc));

  if(newp == oldp) return oldp;
  // Get the momentum of the new partons before remnant extraction
  // For normal remnants this does not change, but in general it may
  Lorentz5Momentum p1 = newp.first->momentum();
  Lorentz5Momentum p2 = newp.second->momentum();
  
  // Creates the new remnants and returns the new PartonBinInstances
  PBIPair newbins = pex.newRemnants(oldp, newp, newStep());
  newStep()->addIntermediate(newp.first);
  newStep()->addIntermediate(newp.second);

  // Boosting the remnants is only necessary if the momentum of the
  // incoming has changed, which is not normally true. The two last
  // flags should be false if the first and/or last side do not have
  // remnants at all. It returns an overall boost which should be
  // applied to all partons in the shower.
  //LorentzRotation tot = pex.boostRemnants(newbins, p1, p2, true, true);

  return newp;

}

PPtr ShowerHandler::findFirstParton(tPPtr seed, tPPair incoming) const{
  tPPtr parent = seed->parents()[0];
  //if no parent there this is a loose end which will 
  //be connected to the remnant soon.
  if(!parent){
    assert(StandardQCDPartonMatcher::Check(seed->data()));
    return seed;
  }
  if( parent == incoming.first || parent == incoming.second ){
    assert(StandardQCDPartonMatcher::Check(seed->data()));
    return seed;
  }else{
    return findFirstParton(parent, incoming);
  }
}

bool ShowerHandler::decayProduct(tPPtr particle) const {
  // must be time-like and not incoming
  if(particle->momentum().m2()<=0.0*GeV2||
     particle == currentSubProcess()->incoming().first||
     particle == currentSubProcess()->incoming().second) return false;
  // if non-coloured this is enough
  if(!particle->dataPtr()->coloured()) return true;
  // if coloured must be unstable
  if(particle->dataPtr()->stable()) return false;
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
  // must not have same particle type as a child
  int id = particle->id();
  for(unsigned int ix=0;ix<particle->children().size();++ix)
    if(particle->children()[ix]->id()==id) return false;
  // otherwise its a decaying particle
  return true;
}

double ShowerHandler::reweightCKKW(int minMult, int maxMult) {
  // return if not doing CKKW
  if(!_useCKKW) return 1.;
  
#ifdef HERWIG_DEBUG_CKKW
  generator()->log() << "== ShowerHandler::reweightCKKW" << endl;
#endif

  // get the hard subprocess particles
  
  PPair in = lastXCombPtr()->subProcess()->incoming();
  ParticleVector out  = lastXCombPtr()->subProcess()->outgoing();
  pair<double,double> x = make_pair(lastXCombPtr()->lastX1(),lastXCombPtr()->lastX2());
  
  bool gotHistory = false;
  
  try {
    
    // check resolution cut
    
    _reweighter->unresolvedCut(in,out);
    
    // set the generation alpha_s
    
    _reweighter->MEalpha(lastXCombPtr()->lastAlphaS());
    
    // reconstruct a history
    
    gotHistory = _reconstructor->reconstruct(in,x,out);
    
  } catch (Veto) {
    
    // as Veto is not handled if the subprocess has not been setup
    // completely, we return weight 0, which, according to Leif,
    // does the same job.
    
    return 0.;
    
  }
  
  if (!gotHistory)
    throw Exception() << "Shower : ShowerHandler::reweightCKKW : no cascade history could be obtained."
		      << Exception::eventerror;
  
  double weight = _reweighter->reweight(_reconstructor->history(),out.size(),minMult);
  
  _evolver->initCKKWShower(out.size(),maxMult);
  
  return weight;

}
