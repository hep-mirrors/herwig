// -*- C++ -*-
#ifndef HERWIG_GaussianIntegral_H
#define HERWIG_GaussianIntegral_H
//
// This is the declaration of the GaussianIntegral class.

#include <vector>
#include "CLHEP/GenericFunctions/AbsFunctional.hh"
#include "CLHEP/GenericFunctions/AbsFunction.hh"
// #include "GaussianIntegral.fh"
// #include "GaussianIntegral.xh"

namespace Herwig {
  
using namespace Genfun;
using std::vector;

/** \ingroup Utilities
 *  \author  Peter Richardson  
 *
 *  This class is designed to perform the integral of a function
 *  using Gaussian quadrature and is based on the GenericFunctions
 *  provided by CLHEP. The method is adaptive based on using 6th,12th,
 *  24th,48th, or 96th order Gaussian quadrature combined with
 *  subdivision of the integral if this is insufficient.
 */
class GaussianIntegral: public AbsFunctional {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor (default with only the limits).
   * @param lower The lower limit of integration.
   * @param upper The upper limit of integration.
   */
  inline GaussianIntegral(double lower, double upper);

  /**
   * Specify all the parameters.
   * @param lower The lower limit of integration.
   * @param upper The upper limit of integration.
   * @param abserr Absolute error.
   * @param relerr Relative error.
   * @param binwidth Width of the bin as a fraction of the integration region.
   * @param maxint Maximum number of intervals
   * @param maxeval Maximum number of function evaluations
   */
  inline GaussianIntegral(double lower, double upper, 
			  double abserr, double relerr, double binwidth,
			  int maxint, int maxeval);

  /**
   * Destructor.
   */
  ~GaussianIntegral();
  //@}

  /**
   * Take the definite integral of a function between the bounds.
   */
  virtual double operator [] (const AbsFunction & function) const;
  
  /**
   * Reset the integration limits.
   * @param lower The lower limit of integration.
   * @param upper The upper limit of integration.
   */
  inline void resetLimits(double lower,double upper);
private:

  /**
   * Initialise the weights and abscissae.
   */
  inline void Init();

private:

  /**
   * The weights for the gaussian quadrature.
   */
  vector< vector<double> > _weights;

  /**
   * The abscissae.
   */
  vector< vector <double> > _abscissae;

  /**
   * Lower limit on the integral
   */
  double _lower;

  /**
   * Upper limit on the integral
   */
  double _upper;

  /**
   * The parameters controlling the error.
   */
  double _abserr,_relerr;

  /**
   * The minimum width of a bin as a fraction of the integration region.
   */
  double _binwidth;

  /**
   * Maximum number of bins.
   */
  int _maxint;

  /**
   * Maximum number of function evaluations.
   */
  int _maxeval;

};

}

#include "GaussianIntegral.icc"

#endif /* HERWIG_GaussianIntegral_H */
