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
#include "SplitFun1to2.h"

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

    pair<ShoKinPtr, tSudakovFormFactorPtr> pairShowerKinSudakov = 
      _pointerSplittingGenerator->chooseForwardBranching(ch, *part, specialDecay);

    //***LOOKHERE***  accept it according to the  showerConstrainer  and soft correction;

    if ( pairShowerKinSudakov.first == ShoKinPtr()  ||  
	 pairShowerKinSudakov.second == tSudakovFormFactorPtr() ) {

      //***LOOKHERE*** rhoD propagation;

    } else {

      hasEmitted = true;
      if ( _pointerSplittingGenerator->
	   generateBranchingKinematics(ch, *part, pairShowerKinSudakov.first, 
				       pairShowerKinSudakov.second) ) {

	// Assign the splitting function and the shower kinematics
	// to the emitting particle.
	part->showerKinematics( pairShowerKinSudakov.first );
	part->splitFun( pairShowerKinSudakov.second->splitFun() ); 

        // For the time being we are considering only 1->2 branching
	tSplitFun1to2Ptr splitFun = 
	  dynamic_ptr_cast< tSplitFun1to2Ptr >( pairShowerKinSudakov.second->splitFun() );
	if ( splitFun ) {	  

          // Create the ShowerParticle objects for the two children
	  // of the emitting particle; set the parent/child relationship;
	  // add them to the  collecShoPar  and  particlesYetToShower  collections.
	  // Notice that the momenta of the shower products is not set: only
	  // at the end of the showering, during the kinematics reconstruction
	  // such momenta are calculated and set.
	  ShoParPtr showerProduct1 = new_ptr( ShowerParticle() );
	  ShoParPtr showerProduct2 = new_ptr( ShowerParticle() );
	  showerProduct1->dataPtr( getParticleData( splitFun->idFirstProduct() ) );
	  showerProduct2->dataPtr( getParticleData( splitFun->idSecondProduct() ) );

	  if ( splitFun->interactionType() == ShowerIndex::QCD ) {

	    //***LOOKHERE*** In the case the splitting is of QCD type
            //               the colour or anticolour line for the two
	    //               children must be set. A part the case of
            //               a gluon splitting in quark-antiquark, a new
	    //               colour line must be created...
            //               Probably the SplitFun1to2 class should have
	    //               some methods that provide information about
	    //               the colour connection about the emitting
	    //               particle and the children products...
            //  
            //               showerProduct1->setAntiColourLine(...);
            //               showerProduct1->setColourLine(...);
            //               showerProduct2->setAntiColourLine(...);
            //               showerProduct2->setColourLine(...);

	  }

          part->addChild( showerProduct1 );
          part->addChild( showerProduct2 );
	  collecShoPar.insert( collecShoPar.end(), showerProduct1 );
	  collecShoPar.insert( collecShoPar.end(), showerProduct2 );
          particlesYetToShower.push_back( showerProduct1 );
          particlesYetToShower.push_back( showerProduct2 );
         
	}
      } else {
	// Something goes wrong: Skip the event!
        throw Exception("ForwardShowerEvolver::timeLikeShower "
                        "***Skip event: problem in the shower kinematics***",
                        Exception::eventerror);            
      }
    }

  } while ( ! particlesYetToShower.empty() );

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "ForwardShowerEvolver::timeLikeShower "
		       << " ===> END DEBUGGING <=== " 
		       << endl;
  }

  return hasEmitted;

}



