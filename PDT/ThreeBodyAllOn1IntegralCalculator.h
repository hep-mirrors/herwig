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

/**
 * The <code>ThreeBodyAllOn1IntegralCalculator</code> class is designed to integrate
 * a function which gives dGamma/ds to give the partial width.
 *
 *
 * <a href="WidthCalculatorBase.html">WidthCalculatorBase.h</a>.
 * 
 */
class ThreeBodyAllOn1IntegralCalculator: public WidthCalculatorBase {

public:

  /**
   * the friend classes to keep the integration members private
   */

  friend class ThreeBodyAllOn1IntegralOuter;

public:

  /**
   * Standard ctors and dtor.
   */
  inline ThreeBodyAllOn1IntegralCalculator();
  /**
   * Standard ctors and dtor.
   */
  inline ThreeBodyAllOn1IntegralCalculator(const ThreeBodyAllOn1IntegralCalculator &);
  /**
   * Standard ctors and dtor.
   */
  virtual ~ThreeBodyAllOn1IntegralCalculator();

public:

  /**
   * constructor with the dgamma/ds as a function
   */
  inline ThreeBodyAllOn1IntegralCalculator(int intype, Energy inmass, Energy inwidth,
				  Genfun::AbsFunction * indGamma,
				  Energy m1,Energy m2,Energy m3);

  /**
   * constructor which constructs the dgamma/ds function from a decayer
   */
  inline ThreeBodyAllOn1IntegralCalculator(int intype, Energy inmass, Energy inwidth,
				  DecayIntegratorPtr,int,
				  Energy m1,Energy m2,Energy m3);

public:

  /**
   * the partial width
   */
  Energy partialWidth(Energy2) const;

  /**
   * mass reset function
   */
  inline void resetMass(int,Energy);

  /**
   * mass get function 
   */
  inline Energy getMass(const int) const;

  inline Energy otherMass(const int) const;

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * the integrand
   */
  Energy integrand(double);

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<ThreeBodyAllOn1IntegralCalculator> initThreeBodyAllOn1IntegralCalculator;

  /**
   * Private and non-existent assignment operator.
   */
  ThreeBodyAllOn1IntegralCalculator & operator=(const ThreeBodyAllOn1IntegralCalculator &);

private:

  /**
   * which scale we are using
   * the mass and the width for the jacobian
   * masses of the external particles
   * the dgammads
   * the integrand
   * the integrator
   */
  int _variabletype;
  /**
   * which scale we are using
   * the mass and the width for the jacobian
   * masses of the external particles
   * the dgammads
   * the integrand
   * the integrator
   */
  Energy _intmass,_intwidth;
  /**
   * which scale we are using
   * the mass and the width for the jacobian
   * masses of the external particles
   * the dgammads
   * the integrand
   * the integrator
   */
  mutable Energy _m[4]; mutable Energy2 _m2[4];
  /**
   * which scale we are using
   * the mass and the width for the jacobian
   * masses of the external particles
   * the dgammads
   * the integrand
   * the integrator
   */
  Genfun::AbsFunction *_theDgamma;
  /**
   * which scale we are using
   * the mass and the width for the jacobian
   * masses of the external particles
   * the dgammads
   * the integrand
   * the integrator
   */
  Genfun::AbsFunction *_theIntegrand;
  /**
   * which scale we are using
   * the mass and the width for the jacobian
   * masses of the external particles
   * the dgammads
   * the integrand
   * the integrator
   */
  GaussianIntegral *_Integrator;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of ThreeBodyAllOn1IntegralCalculator.
 */
struct BaseClassTrait<Herwig::ThreeBodyAllOn1IntegralCalculator,1> {
  typedef Herwig::WidthCalculatorBase NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::ThreeBodyAllOn1IntegralCalculator>
  /**
   * Return the class name.
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  : public ClassTraitsBase<Herwig::ThreeBodyAllOn1IntegralCalculator> {
  static string className() { return "/Herwig++/ThreeBodyAllOn1IntegralCalculator"; }
  static string library() { return "libHwPDT..so"; }

};

}

namespace Herwig {
using namespace Genfun;
using namespace ThePEG; 

/**
 * the class for the integrand
 */
class ThreeBodyAllOn1IntegralOuter : public Genfun::AbsFunction {
    
FUNCTION_OBJECT_DEF(ThreeBodyAllOn1IntegralOuter)
    
public:
 
  /**
   * Constructor
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
