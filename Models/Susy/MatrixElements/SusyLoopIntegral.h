// -*- C++ -*-
//
// SusysLoopIntegral.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SusysLoopIntegral_H
#define HERWIG_SusysLoopIntegral_H
//
// This is the declaration of the SusysLoopIntegral class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Utilities/Maths.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The SusysLoopIntegral class is a pure static class which implements the various loop
 * functions required for the NLO cross sections
 */
class SusyLoopIntegral {

public:

  /**
   *  The \fB_0\f$ function
   */
  static double B0(Energy2 s, Energy m1, Energy m2, Energy2 mu2);

  /**
   *  The derivative of the \f$B_0\f$ function
   */
  static InvEnergy2 B0P(Energy2 s, Energy m1, Energy m2, Energy2 mu2);

  static complex<Energy2> kappa(Energy2 a, Energy2 b, Energy2 c) {
    complex<double> arg = double((sqr(a)+sqr(b)+sqr(c)
				  -2.*(a*b+a*c+b*c))*UnitRemoval::InvE4);
    return sqrt(arg)*UnitRemoval::E2;
  }

  /**
   *  The \f$C_0\f$ function
   */
  static complex<InvEnergy2> C0(Energy2 p1, Energy2 p2, Energy2 p3,
				Energy m1, Energy m2, Energy m3);

  /**
   *  The divergent \f%C_0\f$ function
   */
  static double C0div(Energy2 t,Energy m1,Energy,Energy ms,
		      Energy2 mu2) {
    Energy2 ms2=sqr(ms),m1s=sqr(m1);
    complex<Energy2> msc2 = ms2*Complex(1.,-epsi);
    return real(spenceFunction(double(t/ms2)))-real(spenceFunction(double(m1s/ms2)))
      + real(sqr(log(1.-t  /msc2))) 
      - real(sqr(log(1.-m1s/msc2)))
      + log(ms2/mu2)*log(abs(-(t-ms2)/(ms2-m1s)));
  }

  /**
   *  Real version of the \f$D_0\f$ function
   */
  static InvEnergy4 D_fin(Energy2 p1 ,Energy2 p2 ,Energy2 p3 ,Energy2 p4 ,
			  Energy2 p12,Energy2 p23,Energy m1,Energy m2,
			  Energy m3,Energy m4) {
    return real( D0(p1,p2,p3,p4,p12,p23,m1,m2,m3,m4) );
  }

  /**
   * divergent \f$D_0\f$ function
   */
  static double D_div(Energy2 t,Energy m1,Energy m2,Energy ms,
		      Energy2 s,Energy2 mu2) {
    double zeta2 = sqr(Constants::pi)/6.;
    Energy2 m1s=sqr(m1),m2s=sqr(m2),ms2=sqr(ms);
    complex<Energy2> msc2 = ms2*Complex(1.,-epsi);
    complex<Energy2>   sc = s  *Complex(1., epsi); 
    return 
      -2.*real(spenceFunction(1.+(ms2-m1s)/(t-ms2)))
      -2.*real(spenceFunction(1.+(ms2-m2s)/(t-ms2)))
      -   real(spenceFunction(1.+(ms2-m1s)*(ms2-m2s)/s/ms2 ))
      -1.5*zeta2
      - real(log(1.+(msc2-m1s)*(msc2-m2s)/s/msc2)*
	     (log(-(msc2-m1s)*(msc2-m2s)/s/msc2)
	      -log((msc2-m1s)/mu2)-log((msc2-m2s)/mu2)
	      +log(-s*msc2/sqr(mu2))))
      + 0.5*sqr(log(s/mu2)) - 0.5*sqr(log(s/ms2))
      + 2.*real(log(-sc/mu2)*log(-(t-msc2)/msc2))
      - real(log((msc2-m1s)/mu2)*log((msc2-m1s)/msc2))
      - real(log((msc2-m2s)/mu2)*log((msc2-m2s)/msc2));
  }

private:

  static complex<InvEnergy4> D0(Energy2 p1,Energy2 p2,Energy2 p3,Energy2 p4,
				Energy2 p12,Energy2 p23,
				Energy m1,Energy m2,Energy m3,Energy m4);
                       
  static Complex spenceFunction(Complex z);
  
  static Complex eta(Complex c1, Complex c2) {
    double im1  = c1.imag();
    double im2  = c2.imag();
    double im12 = (c1*c2).imag();
    if( im1 < 0. && im2 < 0.  && im12 > 0.) {
      return Complex(0.,Constants::twopi);
    }
    else if ( im1>0. && im2 > 0.  && im12 < 0.) {
      return Complex(0.,-Constants::twopi);
    }
    else {
      return 0.;
    }
  }
                                   
  static Complex etas(Complex y, Complex r, Complex rs) {
    if(r.imag()!=0.) {
      return eta(y,r);
    }
    else {
      if(r.real()>0.) {
	return 0.;
      }
      else {
	double imy = y.imag(), imrs = rs.imag();
	double s1 = imy >=0. ? 1. : -1.;
	double s2 = imrs>=0. ? 1. : -1.;
	return 0.5*Complex(0.,Constants::pi)*
	  ((1.-s1)*(1.-s2)-(1.+s1)*(1.+s2));
      }
    }                   
  }

  static Complex quadraticSolution(Complex a,Complex b, Complex c) {
    Complex x1 = 0.5/a*(-b+sqrt(sqr(b)-4.*a*c));
    Complex x2 = 0.5/a*(-b-sqrt(sqr(b)-4.*a*c));
    return abs(x1)>abs(x2) ?  x1 : x2;
  }

  static const double epsi;
};

}

#endif /* HERWIG_SusysLoopIntegral_H */
