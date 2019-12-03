// -*- C++ -*-
//
// AlphaS.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2018-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UtilAlphaS_H
#define HERWIG_UtilAlphaS_H

namespace Herwig {

namespace Math {

  /**
   * The derivative of \f$\alpha_S\f$ with respect to \f$\ln(Q^2/\Lambda^2)\f$
   * @param q The scale
   * @param lam \f$\Lambda_{\rm QCD}\f$
   * @param nf The number of flavours 
   */
inline double derivativeAlphaS(Energy q, Energy lam, 
                               unsigned int nf, unsigned int nloop) {
  using Constants::pi;
  double lx = log(sqr(q/lam));
  double b0 = 11. - 2./3.*nf;
  double b1 = 51. - 19./3.*nf;
  double b2 = 2857. - 5033./9.*nf + 325./27.*sqr(nf);
  if(nloop==1)
    return -4.*pi/(b0*sqr(lx));
  else if(nloop==2)
    return -4.*pi/(b0*sqr(lx))*(1.+2.*b1/sqr(b0)/lx*(1.-2.*log(lx)));
  else
    return -4.*pi/(b0*sqr(lx))*
      (1.  + 2.*b1/sqr(b0)/lx*(1.-2.*log(lx))
       + 4.*sqr(b1)/(sqr(sqr(b0))*sqr(lx))*(1. - 2.*log(lx)
					    + 3.*(sqr(log(lx) - 0.5)+b2*b0/(8.*sqr(b1))-1.25)));
}

  /**
   * The 1,2,3-loop parametrization of \f$\alpha_S\f$.
   * @param q The scale
   * @param lam \f$\Lambda_{\rm QCD}\f$
   * @param nf The number of flavours 
   */
inline double alphaS(Energy q, Energy lam, 
                     unsigned int nf, unsigned int nloop) {
  using Constants::pi;
  double lx(log(sqr(q/lam)));
  double b0 = 11. - 2./3.*nf;
  double b1 = 51. - 19./3.*nf;
  double b2 = 2857. - 5033./9.*nf + 325./27.*sqr(nf);
  // one loop
  if(nloop==1)
    {return 4.*pi/(b0*lx);}
  // two loop
  else if(nloop==2) {
    return 4.*pi/(b0*lx)*(1.-2.*b1/sqr(b0)*log(lx)/lx);
  }
  // three loop
  else
    {return 4.*pi/(b0*lx)*(1.-2.*b1/sqr(b0)*log(lx)/lx + 
			   4.*sqr(b1)/(sqr(sqr(b0))*sqr(lx))*
			   (sqr(log(lx) - 0.5) + b2*b0/(8.*sqr(b1)) - 5./4.));}
}


}

}

#endif
