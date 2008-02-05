// -*- C++ -*-
//
// GSLIntegrator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GSLIntegrator_H
#define HERWIG_GSLIntegrator_H
//
// This is the declaration of the GSLIntegrator class.
//

#include "ThePEG/Pointer/ReferenceCounted.h"
#include "gsl/gsl_integration.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Utilities
 * This class is designed to integrate a given function between
 * 2 limits using the gsl QAGS integration subroutine.
 * 
 * The function is supplied using a templated class that must define
 * operator(argument). The units of the argument ArgType and return 
 * type ValType must be supplied in the integrand class using a typedef
 * i.e. <br>
 * <code> struct integrand { </code><br>
 * <code> ... </code> <BR>
 * <code>Energy operator(double arg) const;</code><BR>
 * <code>typedef double ArgType</code><BR>
 * <code>typedef Energy ValType</code><BR>
 * <code> ... </code> <BR>
 * <code>}</code> <BR>
 */
class GSLIntegrator : public Pointer::ReferenceCounted {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default Constructor
   */
  inline GSLIntegrator();

  /**
   * Specify all the parameters.
   * @param abserr Absolute error.
   * @param relerr Relative error.
   * @param binwidth Width of the bin as a fraction of the integration region.
   * @param maxint Maximum number of intervals
   * @param maxeval Maximum number of function evaluations
   */
  inline GSLIntegrator(double abserr, double relerr, int nbins);
  //@}

  /**
   * The value of the integral
   * @param function The integrand class that defines operator()
   * @param lower The lower limit of integration.
   * @param upper The upper limit of integration.
   */
  template <class T>
  inline typename BinaryOpTraits<typename T::ValType,
				 typename T::ArgType>::MulT
  value(const T & function, 
	const typename T::ArgType lower,
	const typename T::ArgType upper) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GSLIntegrator & operator=(const GSLIntegrator &);

private:

  /**
   * The parameters controlling the error.
   */
  double _abserr,_relerr;

  /**
   * The maximum number of intervals to use.
   */
  int _nbins;
};

}

#include "GSLIntegrator.icc"
#include "GSLIntegrator.tcc"

#endif /* HERWIG_GSLIntegrator_H */
