// -*- C++ -*-
//
// ResonanceHelpers.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2018 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ResonanceHelpers_H
#define HERWIG_ResonanceHelpers_H

namespace Resonance {

/**
 *   The velocity squared
 */
inline double beta2(const Energy2 & s, const Energy & m1, const Energy & m2) {
  return max(0.,(1.-sqr(m1+m2)/s)*(1.-sqr(m1-m2)/s));
}
  
/**
 *   The velocity
 */
inline double beta(const Energy2 & s, const Energy & m1, const Energy & m2) {
  return sqrt(beta2(s,m1,m2));
}
  
/**
 *  The derivative of the function \f$\hat{H}(s)\f$ for the GS Breit-Wigner evaluated
 *  at the resonance mass
 */
inline double dHhatds(const Energy & mRes, const Energy & gamma,
		      const Energy & m1, const Energy & m2) {
  double v2 = beta2(sqr(mRes),m1,m2);
  double v = sqrt(v2);
  double r = (sqr(m1) + sqr(m2))/sqr(mRes);
  return gamma/Constants::pi/mRes/v2*
    ((3.-2.*v2- 3.*r)*log((1.+v)/(1.-v)) + 2.*v*(1.- r/(1.-v2)));
}
 
/**
 *  The  \f$\hat{H}(s)\f$ function for the GS Breit-Wigner
 */
inline Energy2 Hhat(const Energy2 & s, const Energy & mRes, const Energy & gamma,
		    const Energy & m1, const Energy & m2) {
  double vR = beta(sqr(mRes),m1,m2);
  double v  = beta(    s    ,m1,m2);
  return gamma/mRes/Constants::pi*s*pow(v/vR,3)*log((1.+v)/(1.-v));
}

/**
 *  The \f$H(s)\f$ function for the GS Breit-Wigner (input the derivative term)
 */
inline Energy2 H(const Energy2 & s, const Energy & mRes, const Energy & gamma,
		 const Energy & m1, const Energy & m2,
                 const double & dH, const Energy2 & Hres) {
  if(s!=ZERO) 
    return Hhat(s,mRes,gamma,m1,m2) - Hres - (s-sqr(mRes))*dH;
  else
    return -2.*sqr(m1+m2)/Constants::pi*gamma/mRes/pow(beta(sqr(mRes),m1,m2),3) - Hres + sqr(mRes)*dH;
}

/**
 *  The \f$H(s)\f$ function for the GS Breit-Wigner (with calc of derivative term)
 */
inline Energy2 H(const Energy2 & s, const Energy & mRes, const Energy & gamma,
		 const Energy & m1, const Energy & m2) {
  double dH = dHhatds(mRes,gamma,m1,m2);
  Energy2 Hres = Hhat(sqr(mRes),mRes,gamma,m1,m2);
  return H(s,mRes,gamma,m1,m2,dH,Hres);
}

/**
 *    The \f$p\f$-wave runningwidth
 */
inline Energy gammaP(const Energy2 & s, const Energy & mRes, const Energy & gamma,
		     const Energy & m1, const Energy & m2) {
  double v2 = beta2(s,m1,m2);
  if(v2<=0.) return ZERO;
  double vR2 = beta2(sqr(mRes),m1,m2);
  double rp = sqrt(v2/vR2);
  return sqrt(s)/mRes*pow(rp,3)*gamma;
}

/**
 *    The \f$p\f$-wave runningwidth
 */
inline Energy gammaS(const Energy2 & s, const Energy & mRes, const Energy & gamma,
		     const Energy & m1, const Energy & m2) {
  double v2 = beta2(s,m1,m2);
  if(v2<=0.) return ZERO;
  double vR2 = beta2(sqr(mRes),m1,m2);
  double rp = sqrt(v2/vR2);
  return mRes/sqrt(s)*rp*gamma;
}

/**
 *  The GS form of the Breit-Wigner distribution
 */
inline Complex BreitWignerGS(const Energy2 & s, const Energy & mRes, const Energy & gamma,
			     const Energy & m1, const Energy & m2,
                             const Energy2 & H0, const double &dH, const Energy2 & Hres) {
  Energy2 mR2=sqr(mRes);
  return (mR2+H0)/(mR2-s+H(s,mRes,gamma,m1,m2,dH,Hres)-Complex(0.,1.)*sqrt(s)*gammaP(s,mRes,gamma,m1,m2));
}

/**
 *  The GS form of the Breit-Wigner distribution
 */
inline Complex BreitWignerGS(const Energy2 & s, const Energy & mRes,
                             const Energy & gamma,
			     const Energy & m1, const Energy & m2) {
  double dH = dHhatds(mRes,gamma,m1,m2);
  Energy2 Hres = Hhat(sqr(mRes),mRes,gamma,m1,m2);
  Energy2 H0 = H(ZERO,mRes,gamma,m1,m2,dH,Hres);
  return BreitWignerGS(s,mRes,gamma,m1,m2,H0,dH,Hres);
}

/**
 *  Standard \f$p\f$-wave Breit-Wigner
 */
inline Complex BreitWignerPWave(const Energy2 & s, const Energy & mRes, const Energy & gamma,
				const Energy & m1, const Energy & m2) {
  Energy2 mR2=sqr(mRes);
  return mR2/(mR2-s-Complex(0.,1.)*sqrt(s)*gammaP(s,mRes,gamma,m1,m2));
}
  
/**
 *  Standard \f$s\f$-wave Breit-Wigner
 */
inline Complex BreitWignerSWave(const Energy2 & s, const Energy & mRes, const Energy & gamma,
				const Energy & m1, const Energy & m2) {
  Energy2 mR2=sqr(mRes);
  return mR2/(mR2-s-Complex(0.,1.)*sqrt(s)*gammaS(s,mRes,gamma,m1,m2));
}
  
/**
 *  Standard fixed width Breit-Wigner
 */
inline Complex BreitWignerFW(const Energy2 & s, const Energy & mRes, const Energy & gamma) {
  Energy2 mR2=sqr(mRes);
  return mR2/(mR2-s-Complex(0.,1.)*mRes*gamma);
}

/**
 *   The \f$H\f$ function from 0512180
 */
Complex H(const Energy & mass, const Energy & width, const Energy2 & sp, const Energy2 & sm,
	  const Energy2 & s0, const Energy & mp, const Energy & m0) {
  return
    Resonance::BreitWignerPWave(sp,mass,width,mp,m0)+
    Resonance::BreitWignerPWave(sm,mass,width,mp,m0)+
    Resonance::BreitWignerPWave(s0,mass,width,mp,mp);
}
}

#endif
