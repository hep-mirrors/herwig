// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ForwardShowerEvolver class.
//

#include "ForwardShowerEvolver.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
// #include "Pythia7/Interface/Parameter.h" 
#include "Pythia7/Interface/Reference.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "ShowerParticle.h"
#include "ShowerKinematics.h"

using namespace Herwig;


ForwardShowerEvolver::~ForwardShowerEvolver() {}


void ForwardShowerEvolver::persistentOutput(PersistentOStream & os) const {
  os << _pointerSplittingGenerator
     << _pointerRhoDMatrixPropagator;
}


void ForwardShowerEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _pointerSplittingGenerator
     >> _pointerRhoDMatrixPropagator;
}


ClassDescription<ForwardShowerEvolver> ForwardShowerEvolver::initForwardShowerEvolver;
// Definition of the static class description member.


void ForwardShowerEvolver::Init() {

  static ClassDocumentation<ForwardShowerEvolver> documentation
    ("This class is responsible for the forward showering of time-like particles",
     "It does also the special forward showering for decaying particles",
     "in which the angular ordering is reversed");

  static Reference<ForwardShowerEvolver,SplittingGenerator> 
    interfaceSplittingGenerator("SplittingGenerator", 
				"A reference to the SplittingGenerator object", 
                                &Herwig::ForwardShowerEvolver::_pointerSplittingGenerator,
				false, false, true, false);
  static Reference<ForwardShowerEvolver,RhoDMatrixPropagator> 
    interfaceRhoDMatrixPropagator("RhoDMatrixPropagator", 
				  "A reference to the RhoDMatrixPropagator object", 
				  &Herwig::ForwardShowerEvolver::_pointerRhoDMatrixPropagator,
				    false, false, true, false);

}

//------------------------------------------------------------------------------

bool ForwardShowerEvolver::
timeLikeShower( tPartCollHdlPtr ch, 
		const tShoConstrPtr showerConstrainer, 
		const tMECorrectionPtr meCorrectionPtr,
		tShoParPtr particle, 
		CollecShoParPtr & collecShoPar,
		const bool specialDecay )  throw (Veto, Stop, Exception) {

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "ForwardShowerEvolver::timeLikeShower "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl;
  }

  bool hasEmitted = false;
  vector<tShoParPtr> particlesYetToShower;
  particlesYetToShower.push_back( particle );

  do {

    tShoParPtr part = particlesYetToShower.back();
    particlesYetToShower.pop_back();

    //***LOOKHERE***  update rhoD matrix of  part ;

    pair<Energy, tSudakovFormFactorPtr> pairScaleSudakov = 
      _pointerSplittingGenerator->chooseForwardBranching(ch, *part, specialDecay);

    //***LOOKHERE***  accept it according to the  showerConstrainer  and soft correction;

    if ( pairScaleSudakov.first == Energy()  ||  
	 pairScaleSudakov.second == tSudakovFormFactorPtr() ) {

      //***LOOKHERE*** rhoD propagation;

    } else {

      hasEmitted = true;
      ShoKinPtr shoKin = _pointerSplittingGenerator->
	generateBranchingKinematics(ch, *part, pairScaleSudakov.first, 
				    pairScaleSudakov.second);
      //FINISH  pairScaleSudakov.second->splitFun()->idEmitter()
      //FINISH pairScaleSudakov.second->splitFun()->massEmitter()
		       
      // create the new ShowerParticles and then call the  part
      //      method  addChildren( children );
      // store them into both collecShoPar and particlesYetToShower;
      // store also the shoKin;
      
    }

  } while ( ! particlesYetToShower.empty() );

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "ForwardShowerEvolver::timeLikeShower "
		       << " ===> END DEBUGGING <=== " 
		       << endl;
  }

  return hasEmitted;

}



