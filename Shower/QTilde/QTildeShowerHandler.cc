// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeShowerHandler class.
//

#include "QTildeShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Shower/QTilde/Base/Evolver.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "Herwig/PDF/MPIPDF.h"
#include "Herwig/PDF/MinBiasPDF.h"
#include "Herwig/Shower/QTilde/Base/ShowerTree.h"
#include "Herwig/Shower/QTilde/Base/KinematicsReconstructor.h"
#include "Herwig/Shower/QTilde/Base/PartnerFinder.h"
#include "Herwig/PDF/HwRemDecayer.h"

using namespace Herwig;

QTildeShowerHandler::QTildeShowerHandler() :
  splitHardProcess_(true)
{}

QTildeShowerHandler::~QTildeShowerHandler() {}

IBPtr QTildeShowerHandler::clone() const {
  return new_ptr(*this);
}

IBPtr QTildeShowerHandler::fullclone() const {
  return new_ptr(*this);
}

void QTildeShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << evolver_ << splitHardProcess_;
}

void QTildeShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> evolver_ >> splitHardProcess_;
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

  static Reference<QTildeShowerHandler,Evolver> 
    interfaceEvolver("Evolver", 
		     "A reference to the Evolver object", 
		     &Herwig::QTildeShowerHandler::evolver_,
		     false, false, true, false);

  static Switch<QTildeShowerHandler,bool> interfaceSplitHardProcess
    ("SplitHardProcess",
     "Whether or not to try and split the hard process into production and decay processes",
     &QTildeShowerHandler::splitHardProcess_, true, false, false);
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

}


tPPair QTildeShowerHandler::cascade(tSubProPtr sub,
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
  while (countFailures<maxtry()) {
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
  if (countFailures >= maxtry()) {
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
    done_[ix]->fillEventRecord(pstep,
			       evolver_->isISRadiationON(),
			       evolver_->isFSRadiationON());
  }
}

Energy QTildeShowerHandler::hardScale() const {
  return evolver_->hardScale();
}
