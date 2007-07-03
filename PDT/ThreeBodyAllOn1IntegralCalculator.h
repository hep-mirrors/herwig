// -*- C++ -*-
#ifndef HERWIG_ThreeBodyAllOn1IntegralCalculator_H
#define HERWIG_ThreeBodyAllOn1IntegralCalculator_H
// This is the declaration of the ThreeBodyAllOn1IntegralCalculator class.

#include "WidthCalculatorBase.h"
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Utilities/GaussianIntegrator.h"

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
  inline ThreeBodyAllOn1IntegralCalculator(int intype, Energy inmass, Energy inwidth,
					   double inpow,
					   T indGamma,int mode,
					   Energy m1,Energy m2,Energy m3);

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
  inline void resetMass(int imass,Energy mass);

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @return The mass required.
   */
  inline Energy getMass(const int imass) const;

  /**
   * Get the masses of all bar the one specified. Used to get the limits
   * for integration.
   * @param imass The particle not needed
   * @return The sum of the other masses.
   */
  inline Energy otherMass(const int imass) const;

  /**
   * The integrand for the inner integrand.
   * @param argument The mass squared for the inner integral
   * @return The value of the inner integrand.
   */
  Energy operator ()(double argument) const;
  typedef double ArgType;
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
  GaussianIntegrator _integrator;

};
}

#include "ThreeBodyAllOn1IntegralCalculator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
 #include "ThreeBodyAllOn1IntegralCalculator.tcc"
#endif

#endif /* HERWIG_ThreeBodyAllOn1IntegralCalculator_H */
