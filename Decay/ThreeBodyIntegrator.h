// -*- C++ -*-
#ifndef HERWIG_ThreeBodyIntegrator_H
#define HERWIG_ThreeBodyIntegrator_H
//
// This is the declaration of the <!id>ThreeBodyIntegrator<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class is designed to integrate using non-MonteCarlo
//  methods, ie trapizum rule/Simpson's/Gaussian Quadrature a three body matrix
//  element using a smoothing like a multi-channel phase-space integrator
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 
//  Author: Peter Richardson
//

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "CLHEP/GenericFunctions/AbsFunction.hh"
#include "Herwig++/Utilities/GaussianIntegral.h"

namespace Herwig {

using namespace ThePEG;

class ThreeBodyIntegrator {
  
public:
  
  // constructor with all the parameters
  inline ThreeBodyIntegrator(vector<double> inweights,
			     vector<int> intype,
			     vector<Energy> inmass,
			     vector<Energy> inwidth,
			     Genfun::AbsFunction * inme,
			     Energy m1,Energy m2,Energy m3);

  // copy constructor
  inline ThreeBodyIntegrator(const ThreeBodyIntegrator &);

  virtual ~ThreeBodyIntegrator();
  // Standard ctors and dtor.

public:
  
  void outerVariables(const double & x, double & low, double & upp);
  // shift the variables for the outer integrand and give limits for the inner one

  double innerIntegrand(const double & y);
  // the integrand for the inner integral
  
  inline Genfun::AbsFunction * InnerIntegrand();
  // pointer to the function for the inner integrand
  
  inline Genfun::AbsFunction * ME();
  // pointer to the matrix element object
  
  Energy width(Energy2 q2);
  // calculate the width for a given mass
  
public:
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
private:
  
  ThreeBodyIntegrator & operator=(const ThreeBodyIntegrator &);
  // Private and non-existent assignment operator.
  
private:
  
  vector<double> _channelweights;
  // weights for the different channels
  vector<int> _channeltype;
  // the types for the different channels
  vector<Energy> _channelmass;
  // the mass of the resonance for a given channel
  vector<Energy> _channelwidth;
  // the width of the resonance for a given channel
  Genfun::AbsFunction *_theME;
  // pointer to a function giving the matrix element as a function of s12,s13,s23
  int _thechannel;
  // the channel currently being integrated
  Energy2 _souter;
  // the value of s for the outer integral
  Energy _m[4]; Energy2 _m2[4];
  // masses of the external particles
  
  Genfun::AbsFunction *_theInnerIntegrand;
  Genfun::AbsFunction *_theOuterIntegrand;
  // the integrands
};
  
}
// the class for the outer integrand
namespace Herwig {

using namespace Genfun;
using namespace ThePEG; 

class ThreeBodyOuterIntegrand : public Genfun::AbsFunction {
  
FUNCTION_OBJECT_DEF(ThreeBodyOuterIntegrand)
    
public:

// Constructor
ThreeBodyOuterIntegrand(ThreeBodyIntegrator *);
  
  // Destructor
  virtual ~ThreeBodyOuterIntegrand();
  
  // Copy constructor
  ThreeBodyOuterIntegrand(const ThreeBodyOuterIntegrand &right);
  
  // Retreive function value
  virtual double operator ()(double argument) const;
  virtual double operator ()(const Argument & a) const {return operator() (a[0]);}
  
private:
  
  // It is illegal to assign a function
  const ThreeBodyOuterIntegrand & operator=(const ThreeBodyOuterIntegrand &right);
  
private:
  
  ThreeBodyIntegrator * _theIntegrator;
  // pointer to the decay integrator
};
  
}

// the class for the inner integrand
namespace Herwig {

using namespace Genfun;
using namespace ThePEG; 

class ThreeBodyInnerIntegrand : public Genfun::AbsFunction {
    
FUNCTION_OBJECT_DEF(ThreeBodyInnerIntegrand)
  
  public:

// Constructor
ThreeBodyInnerIntegrand(ThreeBodyIntegrator *);
  
  // Destructor
  virtual ~ThreeBodyInnerIntegrand();
  
  // Copy constructor
  ThreeBodyInnerIntegrand(const ThreeBodyInnerIntegrand &right);
  
  // Retreive function value
  virtual double operator ()(double argument) const;
  virtual double operator ()(const Argument & a) const {return operator() (a[0]);}
  
  inline void sets(double);
  
private:
  
  // It is illegal to assign a function
  const ThreeBodyInnerIntegrand & operator=(const ThreeBodyInnerIntegrand &right);
  
private:
  
  ThreeBodyIntegrator * _theIntegrator;
  // pointer to the decay integrator
};

}






#include "ThreeBodyIntegrator.icc"

#endif /* HERWIG_ThreeBodyIntegrator_H */
