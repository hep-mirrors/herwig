// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BackwardShowerEvolver class.
//

#include "BackwardShowerEvolver.h"
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


BackwardShowerEvolver::~BackwardShowerEvolver() {}


void BackwardShowerEvolver::persistentOutput(PersistentOStream & os) const {
  os << _pointerSplittingGenerator
     << _pointerRhoDMatrixPropagator
     << _pointerForwardShowerEvolver;
}


void BackwardShowerEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _pointerSplittingGenerator
     >> _pointerRhoDMatrixPropagator
     >> _pointerForwardShowerEvolver;
}


ClassDescription<BackwardShowerEvolver> BackwardShowerEvolver::initBackwardShowerEvolver;
// Definition of the static class description member.

void BackwardShowerEvolver::Init() {

  static ClassDocumentation<BackwardShowerEvolver> documentation
    ("This class is responsible for the backward showering of space-like particles");

  static Reference<BackwardShowerEvolver,SplittingGenerator> 
    interfaceSplittingGenerator("SplittingGenerator", 
				"A reference to the SplittingGenerator object", 
                                &Herwig::BackwardShowerEvolver::_pointerSplittingGenerator,
				false, false, true, false);
  static Reference<BackwardShowerEvolver,RhoDMatrixPropagator> 
    interfaceRhoDMatrixPropagator("RhoDMatrixPropagator", 
				  "A reference to the RhoDMatrixPropagator object", 
				  &Herwig::BackwardShowerEvolver::_pointerRhoDMatrixPropagator,
				  false, false, true, false);
  static Reference<BackwardShowerEvolver,ForwardShowerEvolver> 
    interfaceForwardShowerEvolver("ForwardShowerEvolver", 
				  "A reference to the ForwardShowerEvolver object", 
				  &Herwig::BackwardShowerEvolver::_pointerForwardShowerEvolver,
				  false, false, true, false);

}

//------------------------------------------------------------------------------

bool BackwardShowerEvolver::spaceLikeShower( tPartCollHdlPtr ch, 
					     const tShoConstrPtr showerConstrainer, 
					     const tMECorrectionPtr meCorrectionPtr,
					     tShoParPtr particle, 
					     CollecShoParPtr & collecShoPar ) 
  throw (Veto, Stop, Exception) {
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "BackwardShowerEvolver::spaceLikeShower "
		       << " ===> START DEBUGGING <=== "
		       << "   EventNumber=" << generator()->currentEventNumber() 
		       << endl;
  }

  bool hasEmitted = false;
  tShoParPtr spaceLikePart = particle;
  vector<tShoParPtr> particlesYetToShower;   // only time-like particles

  do {

    //***LOOKHERE***  update rhoD matrix of  spaceLikePart;
                      pair<ShoKinPtr, tSudakovFormFactorPtr> pairShowerKinSudakov = 
                        _pointerSplittingGenerator->chooseBackwardBranching( ch, *spaceLikePart );
    //                accept it according to the  showerConstrainer  and soft correction;
    //                if ( does not branch ) {
    //                  rhoD propagation;
    //                  spaceLikePart = tShoParPtr();
    //                } else {
    //                  hasEmitted = true;
		        _pointerSplittingGenerator->
			  generateBranchingKinematics( ch, *spaceLikePart, pairShowerKinSudakov.first, 
						       pairShowerKinSudakov.second );
    //                  create the new ShowerParticles and then store the
    //                    unique space-like one in spaceLikePart, whereas 
    //                    the others are stored into  particleYetToShower;
    //                  store also the shoKin;
    //                }
    // 
    //                NB) To access the PDF:
    //                      PDF myPDF = ch.pdf( parton );
    //                    where parton is a pointer to a Pythia7 particle
    //                    and the method returns a PDF object for the given
    //                    particle. Such method is defined in LastXCombInfo
    //                    from which PartialCollisionHandler inherits from.
    //***endLOOKHERE***

  } while ( ! spaceLikePart );

  while ( ! particlesYetToShower.empty() ) {

    //***LOOKHERE***  update rhoD of part;

    tShoParPtr part = particlesYetToShower.back();
    particlesYetToShower.pop_back();
    hasEmitted = hasEmitted || 
      _pointerForwardShowerEvolver->timeLikeShower(ch, showerConstrainer, meCorrectionPtr,
						   part, collecShoPar);

    //***LOOKHERE***  update rhoD of the parent of part;
    
  } 

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "BackwardShowerEvolver::spaceLikeShower "
		       << " ===> END DEBUGGING <=== "
		       << endl;

  }

  return hasEmitted;

}


