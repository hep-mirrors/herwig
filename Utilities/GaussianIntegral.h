// -*- C++ -*-
#ifndef HERWIG_GaussianIntegral_H
#define HERWIG_GaussianIntegral_H
//
// This is the declaration of the <!id>GaussianIntegral<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class is designed to perform the integral of a function
// using Gaussian quadrature and is based on the GenericFunctions
// provided by CLHEP. The method is adaptive based on using 6th,12th,
// 24th,48th, or 96th order Gaussian quadrature combined with
// subdivision of the integral if this is insufficient.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 
//  Author: Peter Richardson
//

#include <vector>
#include "CLHEP/GenericFunctions/AbsFunctional.hh"
#include "CLHEP/GenericFunctions/AbsFunction.hh"
// #include "GaussianIntegral.fh"
// #include "GaussianIntegral.xh"

namespace Herwig {
  
using namespace Genfun;
using std::vector;

class GaussianIntegral: public AbsFunctional {

public:

  // constructors (default with only the limits)
  inline GaussianIntegral(double lower, double upper);

  // specify all the parameters
  inline GaussianIntegral(double lower, double upper, 
			  double abserr, double relerr, double binwidth,
			  int maxint, int maxeval);

  ~GaussianIntegral();
  // Standard ctors and dtor.

  // Take the definite integral of a function between the bounds:
  virtual double operator [] (const AbsFunction & function) const;
  
private:

  // initialise the weights and abscissae
  inline void Init();

private:

  // the weights for the gaussian quadrature
  vector< vector<double> > _weights;
  // the abscissae
  vector< vector <double> > _abscissae;
  // limits of the integral
  double _lower,_upper;
  // the parameters controlling the error
  double _abserr,_relerr;
  // the minimum width of a bin as a fraction of the integration region
  double _binwidth;
  // maximum number of bins
  int _maxint;
  // maximum number of function evaluations
  int _maxeval;
};

}

// CLASSDOC OFF

#include "GaussianIntegral.icc"

#endif /* HERWIG_GaussianIntegral_H */
