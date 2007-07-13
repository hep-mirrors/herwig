// -*- C++ -*-
#ifndef HERWIG_ThreeBodyAllOnCalculator_H
#define HERWIG_ThreeBodyAllOnCalculator_H
// This is the declaration of the ThreeBodyAllOnCalculator class.

#include "WidthCalculatorBase.h"
#include "Herwig++/Utilities/GaussianIntegrator.h"
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"

namespace Herwig {
using namespace ThePEG;

template <class T>
class ThreeBodyAllOnCalculator;


/** \ingroup PDT
 *
 *  The ThreeBodyAllOnCalculator class is designed to integrate 
 *  a three-body matrix element in which all the outgoing particles are
 *  on-shell to give the partial width. A multi-channel type approach is
 *  used together with gaussian quadrature.
 *
 * @see GaussianIntegral
 * @see ThreeBodyAllOnOuter
 * @see ThreeBodyAllOnIner
 *
 */
template <class T>
class ThreeBodyAllOnCalculator: public WidthCalculatorBase {


/** \ingroup PDT
 * The class for the outer integrand of the integral of a three body decay matrix
 * element. This class is used by the ThreeBodyAllOnCalculator
 * to perform the outer integral.
 *
 * @see ThreeBodyAllOnCalculator
 * @see ThreeBodyAllOnInner
 */
struct Outer {
 
  /**
   * Constructor with a pointer to the ThreeBodyAllOnCalculator
   */
  inline Outer(typename Ptr<Herwig::ThreeBodyAllOnCalculator<T> >::pointer); 
  
  /**
   * Retreive function value
   */
  inline Energy4 operator ()(double argument) const;
  typedef double ArgType;
  typedef Energy4 ValType;

  /**
   * pointer to the decay integrator
   */
  typename Ptr<Herwig::ThreeBodyAllOnCalculator<T> >::pointer _integrand;
  
  /**
   * gaussian integrator
   */
  GaussianIntegrator _integrator;
};

public:

  /**
   * The ThreeBodyAllOnOuter class is a friend so it can access the private
   * members and perform the integral.
   */
  friend struct ThreeBodyAllOnOuter;

public:

  /**
   * Constructor with all the parameters
   * @param inweights The weights for the different integration channels
   * @param intype The types of the different integration channels.
   * @param inmass The mass for the Jacobian for the different channels.
   * @param inwidth The width for the Jacobian for the different channels.
   * @param inpow the power for power-law smoothing for a given channel
   * @param inme The pointer to the function which gives the matrix element.
   * @param mode The mode to be integrated
   * @param m1 The mass of the first particle.
   * @param m2 The mass of the second particle.
   * @param m3 The mass of the third  particle.
   */
  inline ThreeBodyAllOnCalculator(vector<double> inweights,
				  vector<int> intype,
				  vector<Energy> inmass,
				  vector<Energy> inwidth,
				  vector<double> inpow,
				  T inme, int mode,
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
  Energy2 operator ()(Energy2 argument) const;
  typedef Energy2 ArgType;
  typedef Energy2 ValType;


protected:

  /**
   * shift the variables for the outer integrand and give limits for the inner one.
   * This member sets the value of the _souter member for the mass squared of the 
   * outer integral and calculates the limits on the mass squared of the inner 
   * integral.
   * @param x The integration variable
   * @param low The lower limit for the inner integral.
   * @param upp The upper limit for the inner integral.
   */
  void outerVariables(const double & x, Energy2 & low, Energy2 & upp);

private:

  /**
   * Private and non-existent assignment operator.
   */
  ThreeBodyAllOnCalculator & operator=(const ThreeBodyAllOnCalculator &);

private:
  
  /**
   * weights for the different channels
   */
  vector<double> _channelweights;

  /**
   * the types for the different channels
   */
  vector<int> _channeltype;

  /**
   * the mass of the resonance for a given channel
   */
  vector<Energy> _channelmass;

  /**
   * the width of the resonance for a given channel
   */
  vector<Energy> _channelwidth;

  /**
   * the power for power-law smoothing for a given channel
   */
  vector<double> _channelpower;

  /**
   * Function giving the matrix element as a function of s12,s13,s23
   */
  T _theME;

  /**
   *  The mode
   */
  int _mode;
 
  /**
   * the channel currently being integrated
   */
  mutable int _thechannel;

  /**
   * the value of s for the outer integral
   */
  mutable Energy2 _souter;

  /**
   * masses of the external particles
   */
  mutable vector<Energy>  _m;

  /**
   * mass squareds of the external particles
   */
  mutable vector<Energy2> _m2;

  /**
   * the outer integrand
   */
  Outer _outer;

  /**
   * member to do the integration
   */
  GaussianIntegrator _integrator;
};
}

#include "ThreeBodyAllOnCalculator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "ThreeBodyAllOnCalculator.tcc"
#endif

#endif /* HERWIG_ThreeBodyAllOnCalculator_H */
