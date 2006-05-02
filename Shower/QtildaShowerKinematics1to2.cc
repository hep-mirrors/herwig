// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtildaShowerKinematics1to2 class.
//

#include "QtildaShowerKinematics1to2.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;


QtildaShowerKinematics1to2::~QtildaShowerKinematics1to2() {}

vector<Lorentz5Momentum> QtildaShowerKinematics1to2::getBasis() {
  vector<Lorentz5Momentum> dum; 
  dum.clear(); 
  dum.push_back( _pVector );
  dum.push_back( _nVector );
  return dum; 
}

// --- protected tools ---

Lorentz5Momentum QtildaShowerKinematics1to2::
sudakov2Momentum(double alpha, double beta, Energy px, Energy py) {

  // gives loads of output to check the transformations... 

  const Hep3Vector beta_bb = -(_pVector + _nVector).boostVector();

  // see also these methods!
  //    Hep3Vector findBoostToCM() const;
  //    // Boost needed to get to center-of-mass  frame:
  //            // w.findBoostToCM() == - w.boostVector()
  //            // w.boost(w.findBoostToCM()) == w.rest4Vector()
  
  //    Hep3Vector findBoostToCM( const HepLorentzVector & w ) const;
  //    // Boost needed to get to combined center-of-mass frame:
  //            // w1.findBoostToCM(w2) == w2.findBoostToCM(w1)
  //            // w.findBoostToCM(w) == w.findBoostToCM()

  Lorentz5Momentum p_bb = _pVector;
  Lorentz5Momentum n_bb = _nVector; 
  p_bb.boost( beta_bb );
  n_bb.boost( beta_bb );


  // set first in b2b frame along z-axis (assuming that p and n are
  // b2b as checked above)
  Lorentz5Momentum dq(0.0, 0.0, (alpha - beta)*p_bb.vect().mag(), 
     alpha*p_bb[3] + beta*n_bb[3] );

  // add transverse components
  dq.setPx(px);
  dq.setPy(py);
  
  // rotate to have z-axis parallel to p
  dq.rotateUz( p_bb.vect()/p_bb.vect().mag() );

  // boost back 
  dq.boost( -beta_bb ); 
  dq.rescaleMass(); 

  // check consistency by getting back the Sudakov components from
  // the constructed momentum in the given basis 

  return dq; 
  
}

