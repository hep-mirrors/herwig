// -*- C++ -*-
#ifndef HERWIG_ThreeBodyIntegrator_H
#define HERWIG_ThreeBodyIntegrator_H
//
// This is the declaration of the ThreeBodyIntegrator class.
//
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "CLHEP/GenericFunctions/AbsFunction.hh"
#include "Herwig++/Utilities/GaussianIntegral.h"
#include "ThreeBodyIntegrator.fh"


namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  This class is designed to integrate, using non-MonteCarlo
 *  methods, ie trapizium rule/Simpson's/Gaussian Quadrature a three body matrix
 *  element using a smoothing like a multi-channel phase-space integrator
 *
 * @see GaussianIntegral
 * @see ThreeBodyOuterIntegrand
 * @see ThreeBodyInnerIntegrand
 * 
 *  \author Peter Richardson
 *
 */
class ThreeBodyIntegrator {

  /**
   * The ThreeBodyOuterIntegrand class is a friend so it can access the private
   * members and perform the integral.
   */
  friend class ThreeBodyOuterIntegrand;

  /**
   * The ThreeBodyInnerIntegrand class is a friend so it can access the private
   * members and perform the integral.
   */
  friend class ThreeBodyInnerIntegrand;
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
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
  inline ThreeBodyIntegrator(vector<double> inweights,
			     vector<int> intype,
			     vector<Energy> inmass,
			     vector<Energy> inwidth,
			     Genfun::AbsFunction * inme,
			     Energy m1,Energy m2,Energy m3);

  /**
   * copy constructor
   */
  inline ThreeBodyIntegrator(const ThreeBodyIntegrator &);

  /**
   * Destructor
   */
  virtual ~ThreeBodyIntegrator();
  //@}

public:
  
  /**
   * calculate the width for a given mass
   * @param q2 The mass squared of the decaying particle.
   * @return The partial width.
   */
  Energy width(Energy2 q2) const;
  
public:
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
private:
  
  /**
   * Private and non-existent assignment operator.
   */
  ThreeBodyIntegrator & operator=(const ThreeBodyIntegrator &);
  
private:

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
   * pointer to a function giving the matrix element as a function of \f$m_{12}\f$,
   * \f$m_{13}\f$, \f$m_{23}\f$.
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
};
  
}
namespace Herwig {

using namespace Genfun;
using namespace ThePEG; 

/** \ingroup Decay
 * The class for the outer integrand of the integral of a three body decay matrix
 * element. This class is used by the ThreeBodyIntegrator to perform the outer integral.
 *
 * @see ThreeBodyIntegrator
 * @see ThreeBodyInnerIntegrand
 */
class ThreeBodyOuterIntegrand : public Genfun::AbsFunction {
  
public:

  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;

  /**
   * Clone method
   */
  ThreeBodyOuterIntegrand *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;

public:

  /**
   * Constructor with a pointer to the ThreeBodyIntegrator
   */
  ThreeBodyOuterIntegrand(ThreeBodyIntegrator *);
  
  /**
   * Destructor
   */
  virtual ~ThreeBodyOuterIntegrand();
  
  /**
   * Copy constructor
   */
  ThreeBodyOuterIntegrand(const ThreeBodyOuterIntegrand &right);
  
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
  const ThreeBodyOuterIntegrand & operator=(const ThreeBodyOuterIntegrand &right);
  
private:
  
  /**
   * pointer to the decay integrator
   */
  ThreeBodyIntegrator * _theIntegrator;
};
  
}

namespace Herwig {

using namespace Genfun;
using namespace ThePEG; 

/** \ingroup Decay
 * The class for the inner integrand of the integral of a three body decay matrix
 * element. This class is used by the ThreeBodyIntegrator to perform the inner integral.
 *
 * @see ThreeBodyIntegrator
 * @see ThreeBodyOuterIntegrand
 */
class ThreeBodyInnerIntegrand : public Genfun::AbsFunction {
    
public:

  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;

  /**
   * Clone method
   */
  ThreeBodyInnerIntegrand *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;
  
public:

  /**
   * Constructor with a pointer to the ThreeBodyIntegrator
   */
  ThreeBodyInnerIntegrand(ThreeBodyIntegrator *);
  
  /**
   * Destructor
   */
  virtual ~ThreeBodyInnerIntegrand();
  
  /**
   * Copy constructor
   */
  ThreeBodyInnerIntegrand(const ThreeBodyInnerIntegrand &right);
  
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
  const ThreeBodyInnerIntegrand & operator=(const ThreeBodyInnerIntegrand &right);
  
private:
  
  /**
   * pointer to the decay integrator
   */
  ThreeBodyIntegrator * _theIntegrator;
};

}






#include "ThreeBodyIntegrator.icc"

#endif /* HERWIG_ThreeBodyIntegrator_H */
