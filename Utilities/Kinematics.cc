// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Kinematics class.
//

#include "Kinematics.h"
#include "Pythia7/CLHEPWrap/Lorentz5Vector.h"
#include "Pythia7/CLHEPWrap/LorentzVector.h"
#include "Pythia7/CLHEPWrap/LorentzRotation.h"

using namespace Herwig;
// using namespace Pythia7;


void Kinematics::twoBodyDecay(const Lorentz5Momentum & p,                       // in  
			      const double m1, const double m2,                 // in
			      const Vector3 & unitDir1,                         // in 
			      Lorentz5Momentum & p1, Lorentz5Momentum & p2 ) {  // out
  if ( p.m() >= m1 + m2  &&  m1 >= 0.0  &&  m2 >= 0.0  ) {
    Momentum3 pstarVector = unitDir1;
    pstarVector *= pstarTwoBodyDecay(p.m(),m1,m2);
    p1 = Lorentz5Momentum(m1,pstarVector);
    p2 = Lorentz5Momentum(m2,-pstarVector);
    p1.boost( p.boostVector() );   // boost from CM to LAB
    p2.boost( p.boostVector() );
  }
}


