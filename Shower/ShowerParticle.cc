// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerParticle class.
//

#include "ShowerParticle.h"
#include "Pythia7/EventRecord/Particle.h"
#include "ShowerIndex.h"

using namespace Herwig;


ShowerParticle::ShowerParticle()
  : _pdptr( Pythia7::tcPDPtr() ),
    _parent( tShoParPtr() ),
    _isFinalState( true ),
    _sudAlpha( 1.0 ), 
    _sudBeta( -1.0 ), 
    _sudPx( Energy() ),
    _sudPy( Energy() ),
    _splitFun( tSplitFunPtr() ),
    _decayer( tDecayerPtr() ),
    _showerKinematics( tShoKinPtr() ),
    _scales( ShowerIndex::NumInteractionTypes ),
    _partners( ShowerIndex::NumInteractionTypes ),
    _rhoDUpdate( false ),
    _rhoD( 2, vector< std::complex<double> >(2) )
{ 
  // Initialize the scales to the lowest Energy value, and null partners.
  for ( int i=0; i < ShowerIndex::NumInteractionTypes; i++ ) {
    _scales[i] = Energy();
    _partners[i] = tShoParPtr();
  }

  // Initialize _rhoD matrix to the 2x2 identity, normalisee to trace equal 1.
  _rhoD[0][0] = 0.5;
  _rhoD[0][1] = 0.0;
  _rhoD[1][0] = 0.0;
  _rhoD[1][1] = 0.5;
}


ShowerParticle::~ShowerParticle() {}


ShowerParticle::ShowerParticle(const Pythia7::Particle & inputP7Particle)
  : _pdptr( Pythia7::tcPDPtr() ),
    _parent( tShoParPtr() ),
    _isFinalState( true ),
    _sudAlpha( 1.0 ),
    _sudBeta( -1.0 ),
    _sudPx( Energy() ),
    _sudPy( Energy() ),
    _splitFun( tSplitFunPtr() ),
    _decayer( tDecayerPtr() ),
    _showerKinematics( tShoKinPtr() ),
    _scales( ShowerIndex::NumInteractionTypes ),
    _partners( ShowerIndex::NumInteractionTypes ),
    _rhoDUpdate( false ),
    _rhoD( 2, vector< std::complex<double> >(2) )  
{
  // Initialize the scales to the lowest Energy value, and null partners.
  for ( int i=0; i < ShowerIndex::NumInteractionTypes; i++ ) {
    _scales[i] = Energy();
    _partners[i] = tShoParPtr();
  }

  // Initialize _rhoD matrix to the 2x2 identity, normalisee to trace equal 1.
  _rhoD[0][0] = 0.5;
  _rhoD[0][1] = 0.0;
  _rhoD[1][0] = 0.0;
  _rhoD[1][1] = 0.5;

  dataPtr( inputP7Particle.dataPtr() );
  momentum( inputP7Particle.momentum() );
  position( inputP7Particle.labVertex() );
  //***LOOKHERE*** Notice that the colour lines should be set not here 
  //               but at the beginning of ShowerHandler::cascade() 
}


Pythia7::PPtr ShowerParticle::createPythia7Particle() const {
  PPtr particleP7 = getParticle( data().id() );
  if ( particleP7 ) {
    particleP7->set5Momentum( momentum() );
    particleP7->setLabVertex( position() );
    //***LOOKHERE*** Notice that the parent/children relationships
    //               and the colour lines should be set not here but 
    //               at the end of ShowerHandler::cascade() 
  }
  return particleP7;
}

 
bool ShowerParticle::addChildren(const tCollecShoParPtr & inputChildren) {
  bool isOK = true;
  if (! children().empty() ) {
    // When the current particle has already children (this would happen,
    // under normal conditions without bugs, if that particle has decayed 
    // before showering) then find between the (new)  inputChildren
    // which one has the same id as the current (decaying) particle:
    // this should be present because it represents the emitting
    // particle after the emission. If you can't find it then some
    // bugs is lurking somewhere; otherwise, transfer the original
    // children (decay products) of the current particle to the
    // one with the same identity in the (new)  inputChildren.
    long idParticle = data().id();
    tShoParPtr partAfterEmission = tShoParPtr();
    for ( tCollecShoParPtr::const_iterator cit = inputChildren.begin();
	  cit != inputChildren.end(); ++cit ) {
      if ( (*cit)->data().id() == idParticle ) {
	partAfterEmission = *cit;
      } 
    }    
    if ( partAfterEmission  &&  partAfterEmission->children().empty() ) {
      for ( CollecShoParPtr::const_iterator cit = children().begin();
	    cit != children().end(); ++cit ) {
	partAfterEmission->addChild( *cit );
      }    
      _children.clear();
    } else {  
      isOK = false;  // something wrong.
    }
  }
  if ( isOK ) {
    // Add normally the (new) children one by one, setting their parent pointers.
    for ( tCollecShoParPtr::const_iterator cit = inputChildren.begin();
	  cit != inputChildren.end(); ++cit ) {
      addChild( *cit );
    }
  }
  return isOK;
}

