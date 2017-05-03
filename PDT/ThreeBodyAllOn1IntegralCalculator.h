// -*- C++ -*-
//
// ThreeBodyAllOn1IntegralCalculator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ThreeBodyAllOn1IntegralCalculator_H
#define HERWIG_ThreeBodyAllOn1IntegralCalculator_H
// This is the declaration of the ThreeBodyAllOn1IntegralCalculator class.

#include "WidthCalculatorBase.h"
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Utilities/GSLIntegrator.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 * The ThreeBodyAllOn1IntegralCalculator class is designed to integrate
 * a function which gives \f$d\Gamma/dm^2_{ij}\f$ to give the partial width.
 *
 * @see WidthCalculatorBase
 * @see ThreeBodyAllOn1IntegralOuter
 */
template<class T>
class ThreeBodyAllOn1IntegralCalculator: public WidthCalculatorBase {

public:

  /**
   * Constructor with the \f$d\Gamma/ds\f$ as a function.
   * @param intype The types of the different integration channels.
   * @param inmass The mass for the Jacobian for the different channels.
   * @param inwidth The width for the Jacobian for the different channels.
   * @param inpow  The power for the power-law smoothing function
   * @param indGamma The pointer to the function which gives \f$d\Gamma/ds\f$.
   * @param mode The mode to be calculated
   * @param m1 The mass of the first particle.
   * @param m2 The mass of the second particle.
   * @param m3 The mass of the third  particle.
   */
  ThreeBodyAllOn1IntegralCalculator(int intype, Energy inmass, Energy inwidth,
				    double inpow, T indGamma,int mode,
				    Energy m1,Energy m2,Energy m3)
    : _variabletype(intype),_intmass(inmass),_intwidth(inwidth),
      _intpower(inpow),_mode(mode), _theDgamma(indGamma) {
    _m .resize(4);
    _m2.resize(4);
    _m[1]=m1;_m[2]=m2;_m[3]=m3;
    for(int ix=1;ix<4;++ix)_m2[ix]=sqr(_m[ix]);
  }
  
  /**
   * calculate the width for a given mass
   * @param q2 The mass squared of the decaying particle.
   * @return The partial width.
   */
  Energy partialWidth(Energy2 q2) const;

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @param mass The new value.
   * @return The mass required.
   */
  void resetMass(int imass,Energy mass) {
    assert(imass<4);
    _m[imass]=mass;
    _m2[imass]=mass*mass;
  }

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @return The mass required.
   */
  Energy getMass(const int imass) const {
    assert(imass<4);
    return _m[imass]; 
  }

  /**
   * Get the masses of all bar the one specified. Used to get the limits
   * for integration.
   * @param imass The particle not needed
   * @return The sum of the other masses.
   */
  Energy otherMass(const int imass) const {
    assert(imass>0&&imass<4);
    if(imass==1)      return _m[2]+_m[3];
    else if(imass==2) return _m[1]+_m[3];
    else              return _m[1]+_m[2];
  }

  /**
   * The integrand for the inner integrand.
   * @param argument The mass squared for the inner integral
   * @return The value of the inner integrand.
   */
  Energy operator ()(double argument) const;
  /** Argument type for the GSLIntegrator */
  typedef double ArgType;
  /** Return type for the GSLIntegrator */
  typedef Energy ValType;

private:

  /**
   * Private and non-existent assignment operator.
   */
  ThreeBodyAllOn1IntegralCalculator & 
  operator=(const ThreeBodyAllOn1IntegralCalculator &);

private:

  /**
   * which scale we are using
   */
  int _variabletype;

  /**
   * The mass for the jacobian
   */
  Energy _intmass;

  /**
   * The width for the jacobian
   */
  Energy _intwidth;

  /**
   * The power for power-law smoothing
   */
  double _intpower;


  /**
   *  The mode to be integrated
   */
  int _mode;

  /**
   * masses of the external particles
   */
  mutable vector<Energy>  _m;

  /**
   * mass squareds of the external particles
   */
  mutable vector<Energy2> _m2;

  /**
   * The function for the differential rate
   */
  T _theDgamma;

  /**
   * the integrator
   */
  GSLIntegrator _integrator;

};
}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
 #include "ThreeBodyAllOn1IntegralCalculator.tcc"
#endif

#endif /* HERWIG_ThreeBodyAllOn1IntegralCalculator_H */
