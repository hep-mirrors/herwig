// -*- C++ -*-
#ifndef HERWIG_ThreeBodyAllOnCalculator_H
#define HERWIG_ThreeBodyAllOnCalculator_H
// This is the declaration of the ThreeBodyAllOnCalculator class.

#include "WidthCalculatorBase.h"
#include "ThreeBodyOnShellME.h"
#include "ThreeBodyAllOnCalculator.fh"
#include "CLHEP/GenericFunctions/AbsFunction.hh"
#include "Herwig++/Utilities/GaussianIntegral.h"
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 *  The <code>ThreeBodyAllOnCalculator</code> class is designed to integrate 
 *  a three-body matrix element in which all the outgoing particles are
 *  on-shell to give the partial width. A multi-channel type approach is
 *  used together with gaussian quadrature.
 *
 * @see GaussianIntegral
 * @see ThreeBodyAllOnOuter
 * @see ThreeBodyAllOnIner
 *
 */
class ThreeBodyAllOnCalculator: public WidthCalculatorBase {

public:

  /**
   * The ThreeBodyAllOnOuter class is a friend so it can access the private
   * members and perform the integral.
   */
  friend class ThreeBodyAllOnOuter;

  /**
   *  he ThreeBodyAllOnInner class is a friend so it can access the private
   * members and perform the integral.
   */
  friend class ThreeBodyAllOnInner;

public:

  /** @name Standard constructors and destructors. */
  //@{

  /**
   * Default constructor
   */
  inline ThreeBodyAllOnCalculator();

  /**
   * Copy constructor
   */
  inline ThreeBodyAllOnCalculator(const ThreeBodyAllOnCalculator &);

  /**
   * Destructor
   */
  virtual ~ThreeBodyAllOnCalculator();

  /**
   * Constructor with all the parameters
   * @param inweights The weights for the different integration channels
   * @param intype The types of the different integration channels.
   * @param inmass The mass for the Jacobian for the different channels.
   * @param inwidth The width for the Jacobian for the different channels.
   * @param inme The pointer to the function which gives the matrix element.
   * @param m1 The mass of the first particle.
   * @param m2 The mass of the second particle.
   * @param m3 The mass of the third  particle.
   */
  inline ThreeBodyAllOnCalculator(vector<double> inweights,
				  vector<int> intype,
				  vector<Energy> inmass,
				  vector<Energy> inwidth,
				  Genfun::AbsFunction * inme,
				  Energy m1,Energy m2,Energy m3);

  /**
   * constructor which constructs the me function from a decayer.
   * @param inweights The weights for the different integration channels
   * @param intype The types of the different integration channels.
   * @param inmass The mass for the Jacobian for the different channels.
   * @param inwidth The width for the Jacobian for the different channels.
   * @param decay Pointer to the decayer.
   * @param mode The mode in the decayer being integrated.
   * @param m1 The mass of the first particle.
   * @param m2 The mass of the second particle.
   * @param m3 The mass of the third  particle.
   */
  inline ThreeBodyAllOnCalculator(vector<double> inweights,
				  vector<int> intype,
				  vector<Energy> inmass,
				  vector<Energy> inwidth,
				  DecayIntegratorPtr decay,int mode,
				  Energy m1,Energy m2,Energy m3);
  //@}

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

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
   * shift the variables for the outer integrand and give limits for the inner one.
   * This member sets the value of the _souter member for the mass squared of the 
   * outer integral and calculates the limits on the mass squared of the inner 
   * integral.
   * @param x The integration variable
   * @param low The lower limit for the inner integral.
   * @param upp The upper limit for the inner integral.
   */
  void outerVariables(const double & x, double & low, double & upp);

  /**
   * The integrand for the inner integrand.
   * @param y The mass squared for the inner integral
   * @return The value of the inner integrand.
   */
  double innerIntegrand(const double & y);
  
  /**
   * Pointer to the function for the inner integrand
   */
  inline Genfun::AbsFunction * InnerIntegrand();
  
  /**
   * Pointer to the matrix element object
   */
  inline Genfun::AbsFunction * ME();

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<ThreeBodyAllOnCalculator> initThreeBodyAllOnCalculator;

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
   * pointer to a function giving the matrix element as a function of s12,s13,s23
   */
  Genfun::AbsFunction *_theME;

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
  mutable Energy  _m[4];

  /**
   * mass squareds of the external particles
   */
  mutable Energy2 _m2[4];

  /**
   * the inner integrand
   */
  Genfun::AbsFunction *_theInnerIntegrand;

  /**
   * the outer integrand
   */
  Genfun::AbsFunction *_theOuterIntegrand;

  /**
   * member to do the integration
   */
  GaussianIntegral *_Integrator;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of ThreeBodyAllOnCalculator.
 */
 struct BaseClassTrait<Herwig::ThreeBodyAllOnCalculator,1> {
  /** Typedef of the base class of ThreeBodyAllOnCalculator. */
   typedef Herwig::WidthCalculatorBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ThreeBodyAllOnCalculator>
  : public ClassTraitsBase<Herwig::ThreeBodyAllOnCalculator> {
  /** Return the class name.*/
  static string className() { return "/Herwig++/ThreeBodyAllOnCalculator"; }
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
 * The class for the outer integrand of the integral of a three body decay matrix
 * element. This class is used by the ThreeBodyAllOnCalculator
 * to perform the outer integral.
 *
 * @see ThreeBodyAllOnCalculator
 * @see ThreeBodyAllOnInner
 */
class ThreeBodyAllOnOuter : public Genfun::AbsFunction {

public:

  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;

  /**
   * Clone method
   */
  ThreeBodyAllOnOuter *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;

public:

  /**
   * The ThreeBodyAllOnCalculator is a friend so it can access the private
   * members and perform the integral.
   */
 friend class ThreeBodyAllOnCalculator; 

  /**
   * The ThreeBodyAllOnInner is a friend so it can access the private
   * members and perform the integral.
   */
 friend class ThreeBodyAllOnInner; 

public:
 
  /**
   * Constructor with a pointer to the ThreeBodyAllOnCalculator
   */
 ThreeBodyAllOnOuter(ThreeBodyAllOnCalculatorPtr); 

  /**
   * Destructor
   */
  virtual ~ThreeBodyAllOnOuter();

  /**
   * Copy constructor
   */
  ThreeBodyAllOnOuter(const ThreeBodyAllOnOuter &right);

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
  const ThreeBodyAllOnOuter & operator=(const ThreeBodyAllOnOuter &right);

private:
  
  /**
   * pointer to the decay integrator
   */
  ThreeBodyAllOnCalculatorPtr _integrand;

  /**
   * gaussian integrator
   */
  GaussianIntegral *_Integrator;
  
};

}

namespace Herwig {

using namespace Genfun;
using namespace ThePEG; 

/** \ingroup PDT
 * The class for the inner integrand of the integral of a three body decay matrix
 * element. This class is used by the ThreeBodyAllOnCalcuator
 *  to perform the inner integral.
 *
 * @see ThreeBodyAllOnCalculator
 * @see ThreeBodyAllOnOuter
 */
class ThreeBodyAllOnInner : public Genfun::AbsFunction {

  /**
   * The ThreeBodyAllOnCalculator is a friend so it can access the private
   * members and perform the integral.
   */
  friend class ThreeBodyAllOnCalculator; 

  /**
   * The ThreeBodyAllOnOuter is a friend so it can access the private
   * members and perform the integral.
   */
  friend class ThreeBodyAllOnOuter; 
    
public:

  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;

  /**
   * Clone method
   */
  ThreeBodyAllOnInner *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;
  
public:

  /**
   * Constructor with a pointer to the ThreeBodyAllOnCalculator
   */
  ThreeBodyAllOnInner(ThreeBodyAllOnCalculatorPtr);

  /**
   * Destructor
   */
  virtual ~ThreeBodyAllOnInner();

  /**
   * Copy constructor
   */
  ThreeBodyAllOnInner(const ThreeBodyAllOnInner &right);

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
  const ThreeBodyAllOnInner & operator=(const ThreeBodyAllOnInner &right);
  

private:

  /**
   * pointer to the decay integrator
   */
  ThreeBodyAllOnCalculatorPtr _integrand;
};
}

#include "ThreeBodyAllOnCalculator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreeBodyAllOnCalculator.tcc"
#endif

#endif /* HERWIG_ThreeBodyAllOnCalculator_H */
