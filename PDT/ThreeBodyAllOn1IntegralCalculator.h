// -*- C++ -*-
#ifndef HERWIG_ThreeBodyAllOn1IntegralCalculator_H
#define HERWIG_ThreeBodyAllOn1IntegralCalculator_H
// This is the declaration of the ThreeBodyAllOn1IntegralCalculator class.

#include "WidthCalculatorBase.h"
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThreeBodyAllOn1IntegralCalculator.fh"
#include "ThreeBodyDGammaDs.h"
#include "Herwig++/Utilities/GaussianIntegral.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 * The <code>ThreeBodyAllOn1IntegralCalculator</code> class is designed to integrate
 * a function which gives \f$d\Gamma/dm^2_{ij}\f$ to give the partial width.
 *
 * @see WidthCalculatorBase
 * @see ThreeBodyAllOn1IntegralOuter
 */
class ThreeBodyAllOn1IntegralCalculator: public WidthCalculatorBase {

public:

  /**
   * The ThreeBodyAllOn1IntegralOuter class is a friend to keep the integration
   *  members private.
   */
  friend class ThreeBodyAllOn1IntegralOuter;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline ThreeBodyAllOn1IntegralCalculator();

  /**
   * Copy constructor
   */
  inline ThreeBodyAllOn1IntegralCalculator(const ThreeBodyAllOn1IntegralCalculator &);

  /**
   * Destructor
   */
  virtual ~ThreeBodyAllOn1IntegralCalculator();
  //@}

public:

  /**
   * Constructor with the \f$d\Gamma/ds\f$ as a function.
   * @param intype The types of the different integration channels.
   * @param inmass The mass for the Jacobian for the different channels.
   * @param inwidth The width for the Jacobian for the different channels.
   * @param indGamma The pointer to the function which gives \f$d\Gamma/ds\f$.
   * @param m1 The mass of the first particle.
   * @param m2 The mass of the second particle.
   * @param m3 The mass of the third  particle.
   */
  inline ThreeBodyAllOn1IntegralCalculator(int intype, Energy inmass, Energy inwidth,
					   Genfun::AbsFunction * indGamma,
					   Energy m1,Energy m2,Energy m3);

  /**
   * Constructor which constructs the \f$d\Gamma/ds\f$ function from a decayer
   * @param intype The types of the different integration channels.
   * @param inmass The mass for the Jacobian for the different channels.
   * @param inwidth The width for the Jacobian for the different channels.
   * @param decay Pointer to the DecayIntegrator class.
   * @param mode The mode in the DecayIntegrator we are integrating.
   * @param m1 The mass of the first particle.
   * @param m2 The mass of the second particle.
   * @param m3 The mass of the third  particle.
   */
  inline ThreeBodyAllOn1IntegralCalculator(int intype, Energy inmass, Energy inwidth,
					   DecayIntegratorPtr decay,int mode,
					   Energy m1,Energy m2,Energy m3);

public:

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

protected:

  /**
   * The integrand
   */
  Energy integrand(double);

private:

  /**
   * Private and non-existent assignment operator.
   */
  ThreeBodyAllOn1IntegralCalculator & operator=(const ThreeBodyAllOn1IntegralCalculator &);

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
   * masses of the external particles
   */
  mutable Energy  _m[4];

  /**
   * mass squareds of the external particles
   */
  mutable Energy2 _m2[4];

  /**
   * The function for the differential rate
   */
  Genfun::AbsFunction *_theDgamma;

  /**
   * the integrand
   */
  Genfun::AbsFunction *_theIntegrand;

  /**
   * the integrator
   */
  GaussianIntegral *_Integrator;

};
}

namespace Herwig {
using namespace Genfun;
using namespace ThePEG; 

/** \ingroup PDT
 * The class for the outer integrand of the integral of a three body decay matrix
 * element where one of the integrals has been performed analytically.
 * This class is used by the ThreeBodyAllOn1IntegralCalculator
 * to perform the outer integral.
 *
 * @see ThreeBodyAllOnCalculator
 */
class ThreeBodyAllOn1IntegralOuter : public Genfun::AbsFunction {
    
public:
  
  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;
  
  /**
   * Clone method
   */
  ThreeBodyAllOn1IntegralOuter *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;

public:
 
  /**
   * Constructor with a pointer to the ThreeBodyAllOn1IntegralCalculator
   */
  ThreeBodyAllOn1IntegralOuter(ThreeBodyAllOn1IntegralCalculatorPtr);
  
  /**
   * Destructor
   */
  virtual ~ThreeBodyAllOn1IntegralOuter();
  
  /**
   * Copy constructor
   */
  ThreeBodyAllOn1IntegralOuter(const ThreeBodyAllOn1IntegralOuter &right);

  /**
   * Retreive function value
   */
  virtual double operator ()(double argument) const;

  /**
   * Retreive function value
   */
  virtual double operator ()(const Argument & a) const {return operator() (a[0]);}


private:
  
  /**
   * It is illegal to assign a function
   */
  const ThreeBodyAllOn1IntegralOuter & operator=(const ThreeBodyAllOn1IntegralOuter &right);

private:
  
  /**
   * pointer to the decay integrator
   */
  ThreeBodyAllOn1IntegralCalculatorPtr _theIntegrator;

};
}





















#include "ThreeBodyAllOn1IntegralCalculator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreeBodyAllOn1IntegralCalculator.tcc"
#endif

#endif /* HERWIG_ThreeBodyAllOn1IntegralCalculator_H */
