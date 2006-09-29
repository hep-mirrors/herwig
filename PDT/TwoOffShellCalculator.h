// -*- C++ -*-
#ifndef HERWIG_TwoOffShellCalculator_H
#define HERWIG_TwoOffShellCalculator_H
//
// This is the declaration of the TwoOffShellCalculator class.
//
#include "WidthCalculatorBase.h"
#include "GenericMassGenerator.h"
#include "TwoOffShellCalculator.fh"
#include "OneOffShellCalculator.fh"
#include "Herwig++/Utilities/GaussianIntegral.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 *  Use <code>WidthCalculatorBase</code> objects to integrate over the mass of two
 *  external particles which can be off-shell for running width calculations.
 *
 * @see WidthCalculatorBase
 * @see TwoOffShellIntegrand
 */
class TwoOffShellCalculator: public WidthCalculatorBase {

  /**
   *  The TwoOffShellIntegrand class is a friend to allow access to the private
   *  members for the integration.
   */
  friend class TwoOffShellIntegrand;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline TwoOffShellCalculator();

  /**
   * Copy constructor
   */
  inline TwoOffShellCalculator(const TwoOffShellCalculator &);

  /**
   * Destructor
   */
  virtual ~TwoOffShellCalculator();

  /**
   * Constructor which should be used setting all the required members.
   * @param inloc The mass which is off-shell and to be integrated over.
   * @param inwidth Pointer to the WidthGeneratorBase object which calculates
   * the partial width for a given mass of the off-shell particle. This
   * should be a OneOffShellCalculator instance.
   * @param inmass Pointer to the GenericMassGenerator for the off-shell particle.
   * @param inmin1 The minimum mass for the first off-shell particle.
   * @param inmin2 The minimum mass for the second off-shell particle.
   */
  inline TwoOffShellCalculator(int inloc, WidthCalculatorBasePtr inwidth,
			       GenericMassGeneratorPtr inmass,
			       Energy inmin2,Energy inmin1);
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
   * @param mass The mass of the second off-shell particle,
   * @return The differential rate.
   */
  inline  Energy dGamma(Energy mass);

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<TwoOffShellCalculator> initTwoOffShellCalculator;

  /**
   * Private and non-existent assignment operator.
   */
  TwoOffShellCalculator & operator=(const TwoOffShellCalculator &);

private:

  /**
   * The second mass which is offshell
   */

  int _themass;
  /**
   * the minimum allowed mass
   */
  Energy _minmass;

  /**
   * sum of the masses of the other decay products
   */
  Energy _mother;

  /**
   * pointer to object calculating the width for one-off shell particle.
   */
  WidthCalculatorBasePtr _oneoffwidth;

  /**
   * pointer to object calculating the mass of the particle
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
 * base class of TwoOffShellCalculator.
 */
template <>
 struct BaseClassTrait<Herwig::TwoOffShellCalculator,1> {
  /** Typedef of the base class of TwoOffShellCalculator. */
   typedef Herwig::WidthCalculatorBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::TwoOffShellCalculator>
  : public ClassTraitsBase<Herwig::TwoOffShellCalculator> {
  /** Return the class name.*/
  static string className() { return "Herwig++::TwoOffShellCalculator"; }
};

}


namespace Herwig {
using namespace Genfun;
using namespace ThePEG; 

/** \ingroup PDT
 * Class for the integrand of a matrix element where two of the outgoing
 * particles is off-shell. This class is used by the TwoOffShellCalculator class
 * to perform the integral.
 */
class TwoOffShellIntegrand : public Genfun::AbsFunction {

public:    

  /**
   * The TwoOffShellCaculator is a friend to allow access to the private members
   * for the integration.
   */
  friend class TwoOffShellCaculator;
  
public:

  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;
  
  /**
   * Clone method
   */
  TwoOffShellIntegrand *clone() const;

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
  TwoOffShellIntegrand(TwoOffShellCalculatorPtr in,Energy2 m2,Energy2 mw);

  /**
   * Copy constructor
   */
  TwoOffShellIntegrand(const TwoOffShellIntegrand & x) 
    : Genfun::AbsFunction(), _integrand(x._integrand), _mass2(x._mass2), _mwidth(x._mwidth) {}
  // this one's required because CLHEP::Genfun::AbsFunction 
  // doesn't implement its copy constructor. Aaaargh!

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
  const TwoOffShellIntegrand & operator=(const TwoOffShellIntegrand &right);
  
private:

  /**
   * pointer to the decay integrator
   */
  TwoOffShellCalculatorPtr _integrand;

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


#include "TwoOffShellCalculator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TwoOffShellCalculator.tcc"
#endif

#endif /* HERWIG_TwoOffShellCalculator_H */
