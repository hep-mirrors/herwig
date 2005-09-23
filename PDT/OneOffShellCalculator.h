// -*- C++ -*-
#ifndef HERWIG_OneOffShellCalculator_H
#define HERWIG_OneOffShellCalculator_H
//
// This is the declaration of the OneOffShellCalculator class.
//
#include "GenericMassGenerator.h"
#include "WidthCalculatorBase.h"
#include "OneOffShellCalculator.fh"
#include "Herwig++/Utilities/GaussianIntegral.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 *  Use another <code>WidthCalculatorBase</code> object to integrate over the
 *  mass of on of the external particles which can be off-shell for running
 *  width calculations.
 *
 * @see WidthCalculatorBase
 * @see OneOffShellIntegrand
 * 
 */
class OneOffShellCalculator: public WidthCalculatorBase {

public:

  /**
   *  The OneOffShellIntegrand is a friend to allow access to the members needed
   *  for the integration without making the members public.
   */
  friend class OneOffShellIntegrand;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline OneOffShellCalculator();

  /**
   * Copy constructor
   */
  inline OneOffShellCalculator(const OneOffShellCalculator &);

  /**
   * Destructor
   */
  virtual ~OneOffShellCalculator();

  /**
   * Constructor which should be used setting all the required members.
   * @param inloc The mass which is off-shell and to be integrated over.
   * @param inwidth Pointer to the WidthGeneratorBase object which calculates
   * the partial width for a given mass of the off-shell particle.
   * @param inmass Pointer to the GenericMassGenerator for the off-shell particle.
   * @param inmin The minimum mass for the off-shell particle.
   */
  inline OneOffShellCalculator(int inloc,WidthCalculatorBasePtr inwidth, 
			       GenericMassGeneratorPtr inmass,
			       Energy inmin);
  //@}

public:

  /**
   * member to calculate the partial width.
   * @param scale The mass squared for the decaying particle.
   * @return The partial width.
   */
  Energy partialWidth(Energy2 scale) const;

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

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * The integrand.
   * @param mass The mass of the off-shell particle,
   * @return The differential rate.
   */
  inline  Energy dGamma(Energy mass);

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<OneOffShellCalculator> initOneOffShellCalculator;

  /**
   * Private and non-existent assignment operator.
   */
  OneOffShellCalculator & operator=(const OneOffShellCalculator &);

private:

  /**
   * which mass is offshell.
   */
  int _themass;

  /**
   * the minimum allowed mass.
   */
  Energy _minmass;

  /**
   * pointer to object calculating the on-shell width.
   */
  WidthCalculatorBasePtr _onshellwidth;

  /**
   * pointer to object calculating the mass of the particle.
   */
  GenericMassGeneratorPtr _massptr;

  /**
   * integrator
   */
  GaussianIntegral * _Integrator;

  /**
   * the integrand
   */
  Genfun::AbsFunction *_integrand;

  /**
   * the mass squared of the decaying particle
   */
  mutable Energy2 _scale;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of OneOffShellCalculator.
 */
template <>
struct BaseClassTrait<Herwig::OneOffShellCalculator,1> {
  /** Typedef of the base class of OneOffShellCalculator. */
  typedef Herwig::WidthCalculatorBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::OneOffShellCalculator>
  : public ClassTraitsBase<Herwig::OneOffShellCalculator> {
  /** Return the class name. */
  static string className() { return "Herwig++::OneOffShellCalculator"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return ""; }

};

}

namespace Herwig {
using namespace Genfun;
using namespace ThePEG; 

/** \ingroup PDT
 * Class for the integrand of a matrix element where one of the outgoing
 * particles is off-shell.This class is used by the OneOffShellCalculator class
 * to perform the integral.
 *
 * @see OneOffShellCalculator
 */
class OneOffShellIntegrand : public Genfun::AbsFunction {
  
public:

/**
 * The OneOffShellCaculator is a friend to allow access to the private members
 * for the integration.
 */
 friend class OneOffShellCaculator;

public:

  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;

  /**
   * Clone method
   */
  OneOffShellIntegrand *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;

public:

  /**
   * Constructor.
   * @param in Pointer to the OneOffShellCalculator class this is doing the 
   * integration for.
   * @param m2 The mass squared of the off-shell particle for the Jacobian 
   * transform.
   * @param mw The mass times width of the off-shell particle for the Jacobian 
   * transform.
   */
  OneOffShellIntegrand(tOneOffShellCalculatorPtr in,Energy2 m2,Energy2 mw);

  /**
   * Destructor
   */
  virtual ~OneOffShellIntegrand();

  /**
   * Copy constructor
   */
  OneOffShellIntegrand(const OneOffShellIntegrand &right);

  /**
   * Retreive the function value
   */
  virtual double operator ()(double argument) const;

  /**
   * Retreive the function value
   */
  virtual double operator ()(const Argument & a) const {return operator() (a[0]);}


private:
  
  /**
   * It is illegal to assign a function
   */
  const OneOffShellIntegrand & operator=(const OneOffShellIntegrand &right);
  
private:

  /**
   * pointer to the decay integrator
   */
  tOneOffShellCalculatorPtr _integrand;

  /**
   * The mass squared for the off-shell particle for the Jacobian transform.
   */
  Energy2 _mass2;
  /**
   * The mass times width for the off-shell particle for the Jacobian transform.
   */
  Energy2 _mwidth;
};
}

#include "OneOffShellCalculator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "OneOffShellCalculator.tcc"
#endif

#endif /* HERWIG_OneOffShellCalculator_H */
