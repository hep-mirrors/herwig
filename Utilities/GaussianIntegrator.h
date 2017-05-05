// -*- C++ -*-
//
// GaussianIntegrator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GaussianIntegrator_H
#define HERWIG_GaussianIntegrator_H
//
// This is the declaration of the GaussianIntegrator class.
//

#include "ThePEG/Pointer/ReferenceCounted.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include <vector>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Utilities
 *  \author  Peter Richardson  
 *  This class is designed to perform the integral of a function
 *  using Gaussian quadrature.The method is adaptive based on using 6th,12th,
 *  24th,48th, or 96th order Gaussian quadrature combined with
 *  subdivision of the integral if this is insufficient.
 *
 *  The class is templated on a simple class which should provide a 
 *  T::operator () (double) const which provides the integrand for the function.
 */
class GaussianIntegrator : public Pointer::ReferenceCounted {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default Constructor
   */
  GaussianIntegrator() 
    : _abserr(1.E-35), _relerr(5.E-5), _binwidth(1.E-5), 
      _maxint(100), _maxeval(100000) {
    // setup the weights and abscissae
    Init();
  }
  
  /**
   * Specify all the parameters.
   * @param abserr Absolute error.
   * @param relerr Relative error.
   * @param binwidth Width of the bin as a fraction of the integration region.
   * @param maxint Maximum number of intervals
   * @param maxeval Maximum number of function evaluations
   */
  GaussianIntegrator(double abserr, double relerr, double binwidth,
		     int maxint, int maxeval)
    : _abserr(abserr), _relerr(relerr), _binwidth(binwidth), _maxint(maxint),
      _maxeval(maxeval) {
    // setup the weights and abscissae
    Init();
  }

  /**
   * The value of the integral
   * @param lower The lower limit of integration.
   * @param upper The upper limit of integration.
   */
  template <class T>
  inline typename BinaryOpTraits<typename T::ValType,
				 typename T::ArgType>::MulT
  value(const T &, 
	const typename T::ArgType lower,
	const typename T::ArgType upper) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GaussianIntegrator & operator=(const GaussianIntegrator &);

  /**
   * Initialise the weights and abscissae.
   */
  void Init();

private:

  /**
   * The weights for the gaussian quadrature.
   */
  std::vector< std::vector<double> > _weights;

  /**
   * The abscissae.
   */
  std::vector< std::vector <double> > _abscissae;

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

#include "GaussianIntegrator.tcc"

#endif /* HERWIG_GaussianIntegrator_H */
