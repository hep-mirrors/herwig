// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Component class.
//

#include "Component.h"

// using namespace Pythia7;
using namespace Herwig;
 

Ptr<GlobalParameters>::pointer Component::_pointerGlobalParameters = Ptr<GlobalParameters>::pointer();

    
Component::Component(tPPtr pptr) : _beamRemnant(false) {
  if ( pptr ) {
    _id = pptr->data().id();
    Energy Q0 = Energy();
    if ( _pointerGlobalParameters ) Q0 = _pointerGlobalParameters->effectiveGluonMass();
    if ( pptr->scale() > Q0*Q0 ) {
      _perturbative = true;
    } else {
      _perturbative = false;
    }

    // Notice that the mass (5-th component) of the component is set
    // by hand to the constituent mass (whereas the original mass
    // of the pointed particle could be the current mass). It is not
    // harmful to leave the 5-momentum of the component in an
    // inconsistent state (that is its mass is not the invariant
    // mass obtained from the 4-momentum).
    _momentum = pptr->momentum();
    _momentum.setMass( pptr->data().constituentMass() );

    _position = pptr->labVertex();
  } else {
    _id = 0;
    _perturbative = false;
  }
  _pptr = pptr;
}

