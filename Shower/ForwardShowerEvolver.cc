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
#include "ShowerColourLine.h"

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

      _pointerSplittingGenerator->
	generateBranchingKinematics(ch, *part, pairShowerKinSudakov.first, 
				    pairShowerKinSudakov.second); 

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
	// Remember that we our splitting functions are associated
	// only with particles (id>0) and not with antiparticles (id<0)
	// (because we are assuming CP-conserving vertices): therefore
	// the signs of the decays products must be set by hand explicitly;
	// to simplify this, we assume that the first product must 
	// always have the sign as the parent: q->q+g, qbar->qbar+g,
	// g->q+qbar, etc.
	// Notice that the momenta of the shower products is not set: only
	// at the end of the showering, during the kinematics reconstruction
	// such momenta are calculated and set.
	ShoParPtr showerProduct1 = new_ptr( ShowerParticle() );
	ShoParPtr showerProduct2 = new_ptr( ShowerParticle() );
	if ( part->data().id() > 0 ) {
	  showerProduct1->dataPtr( getParticleData( splitFun->idFirstProduct() ) );
	} else {
	  showerProduct1->dataPtr( getParticleData( - splitFun->idFirstProduct() ) );
	}
	showerProduct2->dataPtr( getParticleData( splitFun->idSecondProduct() ) );

	// Set the Sudakov kinematics variables of the branching products.
        vector<double> sudAlphaProducts;
	vector<Energy> sudPxProducts, sudPyProducts;
        part->showerKinematics()->
	  updateChildren( part->sudAlpha(), part->sudPx(), part->sudPy(),
			  sudAlphaProducts, sudPxProducts, sudPyProducts );  
	if ( sudAlphaProducts.size() == 2  &&
	     sudPxProducts.size() == 2  && sudPyProducts.size() == 2 ) {
	  showerProduct1->sudAlpha( sudAlphaProducts[0] ); 
	  showerProduct1->sudPx( sudPxProducts[0] ); 
	  showerProduct1->sudPy( sudPyProducts[0] ); 
	  showerProduct2->sudAlpha( sudAlphaProducts[1] ); 
	  showerProduct2->sudPx( sudPxProducts[1] ); 
	  showerProduct2->sudPy( sudPyProducts[1] ); 
	}

	// In the case of splittings which involves coloured particles,
        // set properly the colour flow of the branching.
	// Notice that the methods:  ShowerColourLine::addColoured  and
	// ShowerColourLine::addAntiColoured  automatically set also,
	// respectively, the colourLine and antiColourLine of the 
	// ShowerParticle  object they received as argument.
	ShoColinePair parentShoColinePair = ShoColinePair( part->colourLine(), 
							   part->antiColourLine() );
	ShoColinePair showerProduct1ShoColinePair = ShoColinePair();
	ShoColinePair showerProduct2ShoColinePair = ShoColinePair();
	splitFun->colourConnection( parentShoColinePair,
				    showerProduct1ShoColinePair, showerProduct2ShoColinePair );
	if ( showerProduct1ShoColinePair.first ) {
	  showerProduct1ShoColinePair.first->addColoured( showerProduct1 );
	}
	if ( showerProduct1ShoColinePair.second ) {
	  showerProduct1ShoColinePair.second->addAntiColoured( showerProduct1 );
	}
	if ( showerProduct2ShoColinePair.first ) {
	  showerProduct2ShoColinePair.first->addColoured( showerProduct2 );
	}
	if ( showerProduct2ShoColinePair.second ) {
	  showerProduct2ShoColinePair.second->addAntiColoured( showerProduct2 );
	}
	
	part->addChild( showerProduct1 );
	part->addChild( showerProduct2 );
	collecShoPar.insert( collecShoPar.end(), showerProduct1 );
	collecShoPar.insert( collecShoPar.end(), showerProduct2 );
	particlesYetToShower.push_back( showerProduct1 );
	particlesYetToShower.push_back( showerProduct2 );	

	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log() << " Splitting : " << part->data().PDGName()
                             << "  ---> " << showerProduct1->data().PDGName()
                             << "  +  " <<  showerProduct2->data().PDGName() << endl
			     << " \t \t \t \t    colourLine \t antiColourLine " << endl;
	  // Create a map of pointers to ShowerColourLine object,
	  // in order to have a nice printing of colour lines,
	  // numbered from 601 onwards, instead of printing 
	  // physical memory addresses.
	  map<tShoColinePtr,int> mapShoColine;
	  int countShoColine = 601;
	  for ( int i=0; i<3; i++ ) {
	    tShoParPtr theParticle;
	    if ( i == 0 ) {
	      theParticle = part;
	      generator()->log() << "\t   parent:     ";
	    } else if ( i == 1 ) {
	      theParticle = showerProduct1;
	      generator()->log() << "\t   product1:   ";
	    } else if ( i == 2 ) {
	      theParticle = showerProduct2;
	      generator()->log() << "\t   product1:   ";
	    }
	    generator()->log() << theParticle->data().PDGName() << "\t \t";
	    if ( theParticle->colourLine() ) {
	      if ( mapShoColine.find( theParticle->colourLine() ) == mapShoColine.end() ) {
		mapShoColine.insert( pair<tShoColinePtr,int>( theParticle->colourLine(), 
							      countShoColine++ ) );
	      }
	      generator()->log() << mapShoColine.find( theParticle->colourLine() )->second;
	    } else {
	      generator()->log() << "0";
	    }
	    generator()->log() << "\t \t";
	    if ( theParticle->antiColourLine() ) {
	      if ( mapShoColine.find( theParticle->antiColourLine() ) == mapShoColine.end() ) {
		mapShoColine.insert( pair<tShoColinePtr,int>( theParticle->antiColourLine(), 
							      countShoColine++ ) );
	      }
	      generator()->log() << mapShoColine.find( theParticle->antiColourLine() )->second;
	    } else {
	      generator()->log() << "0";
	    }
	    generator()->log() << endl;
	  }
	}
	
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



