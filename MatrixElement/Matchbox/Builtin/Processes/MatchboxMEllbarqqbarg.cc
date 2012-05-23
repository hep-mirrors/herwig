// -*- C++ -*-
//
// MatchboxMEllbarqqbarg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "MatchboxMEllbarqqbarg.h"

using namespace Herwig;
using namespace SpinorHelicity;

double MatchboxMEllbarqqbarg::evaluateME2(bool photon) const {

  double Q2 = (momentum(0)+momentum(1)).m2();
  double M2 = sqr(BMass/amplitudeScale());
  double G2 = sqr(BWidth/amplitudeScale());

  double BB1 =
    (16*(4*al*aq*vl*vq*((invariant(0,4)*invariant(1,3) - invariant(0,3)*invariant(1,4))*invariant(2,4) + 
          2*(-((invariant(0,2) + invariant(0,4))*invariant(1,3)) + invariant(0,3)*(invariant(1,2) + invariant(1,4)))*
           sqr(mass(2))) + sqr(al)*(2*(aq - vq)*(aq + vq)*(invariant(0,1) - 4*mass(0)*mass(1))*mass(2)*mass(3)*
           (invariant(2,4) + 2*sqr(mass(2))) + (sqr(aq) + sqr(vq))*
           (invariant(2,4)*(invariant(0,4)*invariant(1,3) + invariant(0,3)*invariant(1,4) - 
                2*invariant(3,4)*mass(0)*mass(1)) - 
             2*((invariant(0,2) + invariant(0,4))*invariant(1,3) + invariant(0,3)*(invariant(1,2) + invariant(1,4)) - 
                2*(invariant(2,3) + invariant(3,4))*mass(0)*mass(1))*sqr(mass(2)))) + 
       sqr(vl)*(2*(aq - vq)*(aq + vq)*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3)*
           (invariant(2,4) + 2*sqr(mass(2))) + (sqr(aq) + sqr(vq))*
           (invariant(2,4)*(invariant(0,4)*invariant(1,3) + invariant(0,3)*invariant(1,4) + 
                2*invariant(3,4)*mass(0)*mass(1)) - 
             2*((invariant(0,2) + invariant(0,4))*invariant(1,3) + invariant(0,3)*(invariant(1,2) + invariant(1,4)) + 
                2*(invariant(2,3) + invariant(3,4))*mass(0)*mass(1))*sqr(mass(2))))))/sqr(invariant(2,4));

  double BB2 =
    (16*(4*al*aq*vl*vq*((-(invariant(0,4)*invariant(1,2)) + invariant(0,2)*invariant(1,4))*invariant(3,4) + 
          2*((invariant(0,3) + invariant(0,4))*invariant(1,2) - invariant(0,2)*(invariant(1,3) + invariant(1,4)))*
           sqr(mass(3))) + sqr(al)*(2*(aq - vq)*(aq + vq)*invariant(3,4)*(invariant(0,1) - 4*mass(0)*mass(1))*mass(2)*
           mass(3) + invariant(3,4)*(invariant(0,4)*invariant(1,2) + invariant(0,2)*invariant(1,4) - 
             2*invariant(2,4)*mass(0)*mass(1))*(sqr(aq) + sqr(vq)) + 
          4*(aq - vq)*(aq + vq)*(invariant(0,1) - 4*mass(0)*mass(1))*mass(2)*mass(3)*sqr(mass(3)) - 
          2*((invariant(0,3) + invariant(0,4))*invariant(1,2) + invariant(0,2)*(invariant(1,3) + invariant(1,4)) - 
             2*(invariant(2,3) + invariant(2,4))*mass(0)*mass(1))*(sqr(aq) + sqr(vq))*sqr(mass(3))) + 
       sqr(vl)*(2*(aq - vq)*(aq + vq)*invariant(3,4)*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3) + 
          invariant(3,4)*(invariant(0,4)*invariant(1,2) + invariant(0,2)*invariant(1,4) + 
             2*invariant(2,4)*mass(0)*mass(1))*(sqr(aq) + sqr(vq)) + 
          4*(aq - vq)*(aq + vq)*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3)*sqr(mass(3)) - 
          2*((invariant(0,3) + invariant(0,4))*invariant(1,2) + invariant(0,2)*(invariant(1,3) + invariant(1,4)) + 
             2*(invariant(2,3) + invariant(2,4))*mass(0)*mass(1))*(sqr(aq) + sqr(vq))*sqr(mass(3)))))/sqr(invariant(3,4));

  double BB12 = 
    (16*(4*al*aq*vl*vq*(invariant(0,4)*(-invariant(1,2) + invariant(1,3))*invariant(2,3) - 
          invariant(0,3)*(invariant(1,4)*invariant(2,3) + 
             invariant(1,2)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))) + 
          invariant(0,2)*(invariant(1,4)*invariant(2,3) + 
             invariant(1,3)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4)))) + 
       sqr(vl)*(sqr(vq)*(2*invariant(0,2)*invariant(1,3)*invariant(2,3) + 
             invariant(0,2)*invariant(1,4)*invariant(2,3) + invariant(0,2)*invariant(1,3)*invariant(2,4) - 
             2*invariant(0,2)*invariant(1,2)*invariant(3,4) + invariant(0,2)*invariant(1,3)*invariant(3,4) + 
             invariant(0,3)*(invariant(1,4)*invariant(2,3) - 2*invariant(1,3)*invariant(2,4) + 
                invariant(1,2)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))) + 
             4*invariant(2,3)*invariant(2,4)*mass(0)*mass(1) + 4*invariant(2,3)*invariant(3,4)*mass(0)*mass(1) + 
             2*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*
              mass(3) + invariant(0,4)*((invariant(1,2) + invariant(1,3))*invariant(2,3) - 
                4*invariant(1,4)*mass(2)*mass(3)) + 4*mass(0)*mass(1)*sqr(invariant(2,3)) - 
             4*invariant(3,4)*mass(0)*mass(1)*sqr(mass(2)) - 4*invariant(2,4)*mass(0)*mass(1)*sqr(mass(3))) + 
          sqr(aq)*(2*invariant(0,2)*invariant(1,3)*invariant(2,3) + invariant(0,2)*invariant(1,4)*invariant(2,3) + 
             invariant(0,2)*invariant(1,3)*invariant(2,4) - 2*invariant(0,2)*invariant(1,2)*invariant(3,4) + 
             invariant(0,2)*invariant(1,3)*invariant(3,4) + 
             invariant(0,3)*(invariant(1,4)*invariant(2,3) - 2*invariant(1,3)*invariant(2,4) + 
                invariant(1,2)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))) + 
             4*invariant(2,3)*invariant(2,4)*mass(0)*mass(1) + 4*invariant(2,3)*invariant(3,4)*mass(0)*mass(1) - 
             2*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*
              mass(3) + invariant(0,4)*((invariant(1,2) + invariant(1,3))*invariant(2,3) + 
                4*invariant(1,4)*mass(2)*mass(3)) + 4*mass(0)*mass(1)*sqr(invariant(2,3)) - 
             4*invariant(3,4)*mass(0)*mass(1)*sqr(mass(2)) - 4*invariant(2,4)*mass(0)*mass(1)*sqr(mass(3)))) + 
       sqr(al)*(sqr(vq)*(2*invariant(0,2)*invariant(1,3)*invariant(2,3) + 
             invariant(0,2)*invariant(1,4)*invariant(2,3) + invariant(0,2)*invariant(1,3)*invariant(2,4) - 
             2*invariant(0,2)*invariant(1,2)*invariant(3,4) + invariant(0,2)*invariant(1,3)*invariant(3,4) + 
             invariant(0,3)*(invariant(1,4)*invariant(2,3) - 2*invariant(1,3)*invariant(2,4) + 
                invariant(1,2)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))) - 
             4*invariant(2,3)*invariant(2,4)*mass(0)*mass(1) - 4*invariant(2,3)*invariant(3,4)*mass(0)*mass(1) + 
             2*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))*(invariant(0,1) - 4*mass(0)*mass(1))*mass(2)*
              mass(3) + invariant(0,4)*((invariant(1,2) + invariant(1,3))*invariant(2,3) - 
                4*invariant(1,4)*mass(2)*mass(3)) - 4*mass(0)*mass(1)*sqr(invariant(2,3)) + 
             4*invariant(3,4)*mass(0)*mass(1)*sqr(mass(2)) + 4*invariant(2,4)*mass(0)*mass(1)*sqr(mass(3))) + 
          sqr(aq)*(2*invariant(0,2)*invariant(1,3)*invariant(2,3) + invariant(0,2)*invariant(1,4)*invariant(2,3) + 
             invariant(0,2)*invariant(1,3)*invariant(2,4) - 2*invariant(0,2)*invariant(1,2)*invariant(3,4) + 
             invariant(0,2)*invariant(1,3)*invariant(3,4) + 
             invariant(0,3)*(invariant(1,4)*invariant(2,3) - 2*invariant(1,3)*invariant(2,4) + 
                invariant(1,2)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))) - 
             4*invariant(2,3)*invariant(2,4)*mass(0)*mass(1) - 4*invariant(2,3)*invariant(3,4)*mass(0)*mass(1) - 
             2*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))*(invariant(0,1) - 4*mass(0)*mass(1))*mass(2)*
              mass(3) + invariant(0,4)*((invariant(1,2) + invariant(1,3))*invariant(2,3) + 
                4*invariant(1,4)*mass(2)*mass(3)) - 4*mass(0)*mass(1)*sqr(invariant(2,3)) + 
             4*invariant(3,4)*mass(0)*mass(1)*sqr(mass(2)) + 4*invariant(2,4)*mass(0)*mass(1)*sqr(mass(3))))))/
    (invariant(2,4)*invariant(3,4));

  double res =
    ( BB1 + BB2 + BB12 ) / ( sqr(Q2-M2) + M2*G2 );

  if ( photon ) {

  double PP1 = 
    (16*sqr(el)*sqr(eq)*(invariant(0,4)*invariant(1,3)*(invariant(2,4) - 2*sqr(mass(2))) + 
       invariant(0,3)*(invariant(1,4)*invariant(2,4) - 2*(invariant(1,2) + invariant(1,4))*sqr(mass(2))) + 
       2*(-(invariant(0,2)*invariant(1,3)*sqr(mass(2))) - 
          (invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3)*(invariant(2,4) + 2*sqr(mass(2))) + 
          mass(0)*mass(1)*(invariant(2,4)*invariant(3,4) - 2*(invariant(2,3) + invariant(3,4))*sqr(mass(2))))))/
    sqr(invariant(2,4));

  double PP2 = 
    (16*sqr(el)*sqr(eq)*(invariant(3,4)*(invariant(0,4)*invariant(1,2) + invariant(0,2)*invariant(1,4) + 
          2*invariant(2,4)*mass(0)*mass(1)) - 2*invariant(3,4)*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3) - 
       2*((invariant(0,3) + invariant(0,4))*invariant(1,2) + invariant(0,2)*(invariant(1,3) + invariant(1,4)) + 
          2*(invariant(2,3) + invariant(2,4))*mass(0)*mass(1))*sqr(mass(3)) - 
			 4*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3)*sqr(mass(3))))/sqr(invariant(3,4));

  double PP12 = 
    (16*sqr(el)*sqr(eq)*(2*invariant(0,2)*invariant(1,3)*invariant(2,3) + invariant(0,2)*invariant(1,4)*invariant(2,3) + 
       invariant(0,2)*invariant(1,3)*invariant(2,4) - 2*invariant(0,2)*invariant(1,2)*invariant(3,4) + 
       invariant(0,2)*invariant(1,3)*invariant(3,4) + 
       invariant(0,3)*(invariant(1,4)*invariant(2,3) - 2*invariant(1,3)*invariant(2,4) + 
          invariant(1,2)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))) + 
       4*invariant(2,3)*invariant(2,4)*mass(0)*mass(1) + 4*invariant(2,3)*invariant(3,4)*mass(0)*mass(1) + 
       2*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3) + 
       invariant(0,4)*((invariant(1,2) + invariant(1,3))*invariant(2,3) - 4*invariant(1,4)*mass(2)*mass(3)) + 
       4*mass(0)*mass(1)*sqr(invariant(2,3)) - 4*invariant(3,4)*mass(0)*mass(1)*sqr(mass(2)) - 
			 4*invariant(2,4)*mass(0)*mass(1)*sqr(mass(3))))/(invariant(2,4)*invariant(3,4));

  double PB1 = 
    (-16*el*eq*(al*aq*((-(invariant(0,4)*invariant(1,3)) + invariant(0,3)*invariant(1,4))*invariant(2,4) + 
          2*((invariant(0,2) + invariant(0,4))*invariant(1,3) - invariant(0,3)*(invariant(1,2) + invariant(1,4)))*
           sqr(mass(2))) + vl*vq*(-(invariant(0,4)*invariant(1,3)*(invariant(2,4) - 2*sqr(mass(2)))) + 
          invariant(0,3)*(-(invariant(1,4)*invariant(2,4)) + 2*(invariant(1,2) + invariant(1,4))*sqr(mass(2))) + 
          2*(invariant(0,2)*invariant(1,3)*sqr(mass(2)) + 
             (invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3)*(invariant(2,4) + 2*sqr(mass(2))) + 
             mass(0)*mass(1)*(-(invariant(2,4)*invariant(3,4)) + 2*(invariant(2,3) + invariant(3,4))*sqr(mass(2)))))))/
    sqr(invariant(2,4));

  double PB2 = 
    (16*el*eq*(al*aq*((-(invariant(0,4)*invariant(1,2)) + invariant(0,2)*invariant(1,4))*invariant(3,4) + 
          2*((invariant(0,3) + invariant(0,4))*invariant(1,2) - invariant(0,2)*(invariant(1,3) + invariant(1,4)))*
           sqr(mass(3))) + vl*vq*(invariant(3,4)*
           (invariant(0,4)*invariant(1,2) + invariant(0,2)*invariant(1,4) + 2*invariant(2,4)*mass(0)*mass(1)) - 
          2*invariant(3,4)*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3) - 
          2*((invariant(0,3) + invariant(0,4))*invariant(1,2) + invariant(0,2)*(invariant(1,3) + invariant(1,4)) + 
             2*(invariant(2,3) + invariant(2,4))*mass(0)*mass(1))*sqr(mass(3)) - 
				  4*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3)*sqr(mass(3)))))/sqr(invariant(3,4));

  double PB12 = 
    (16*el*eq*(al*aq*(invariant(0,4)*(-invariant(1,2) + invariant(1,3))*invariant(2,3) - 
          invariant(0,3)*(invariant(1,4)*invariant(2,3) + 
             invariant(1,2)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))) + 
          invariant(0,2)*(invariant(1,4)*invariant(2,3) + 
             invariant(1,3)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4)))) + 
       vl*vq*(2*invariant(0,2)*invariant(1,3)*invariant(2,3) + invariant(0,2)*invariant(1,4)*invariant(2,3) + 
          invariant(0,2)*invariant(1,3)*invariant(2,4) - 2*invariant(0,2)*invariant(1,2)*invariant(3,4) + 
          invariant(0,2)*invariant(1,3)*invariant(3,4) + 
          invariant(0,3)*(invariant(1,4)*invariant(2,3) - 2*invariant(1,3)*invariant(2,4) + 
             invariant(1,2)*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))) + 
          4*invariant(2,3)*invariant(2,4)*mass(0)*mass(1) + 4*invariant(2,3)*invariant(3,4)*mass(0)*mass(1) + 
          2*(2*invariant(2,3) + invariant(2,4) + invariant(3,4))*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3) + 
          invariant(0,4)*((invariant(1,2) + invariant(1,3))*invariant(2,3) - 4*invariant(1,4)*mass(2)*mass(3))) + 
       4*vl*vq*mass(0)*mass(1)*(sqr(invariant(2,3)) - invariant(3,4)*sqr(mass(2)) - invariant(2,4)*sqr(mass(3)))))/
    (invariant(2,4)*invariant(3,4));

  res +=
    ( PP1 + PP2 + PP12 ) / sqr(Q2) +
    2.*(1.-M2/Q2)*( PB1 + PB2 + PB12 ) / ( sqr(Q2-M2) + M2*G2 );

  }

  res *= 2.*Constants::pi*alphaS*(sqr(Nc)-1.);

  return res;

}
