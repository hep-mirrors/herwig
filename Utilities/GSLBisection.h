// -*- C++ -*-
//
// GSLBisection.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GSLBisection_H
#define HERWIG_GSLBisection_H
//
// This is the declaration of the GSLBisection class.
//

#include "ThePEG/Pointer/ReferenceCounted.h"
#include "Herwig/Utilities/GSLHelper.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Utilities
 * This class is designed to find the root of a given function between
 * 2 limits using bisection methods.
 *
 * \author Manuel B\"ahr
 * 
 * The function is supplied using a templated class that must define
 * operator(argument). The units of the argument ArgType and return type
 * ValType must be supplied in the integrand class using a typedef. In
 * addition the baseunit should be supplied by static methods vUnit()
 * and aUnit() to avoid numerical problems that arise when the centrally
 * defined baseunit is several orders of magnitude off the one you
 * need. As an example see: <br>
 * <code> struct integrand { </code><br>
 * <code> ... </code> <BR>
 * <code>Energy operator(CrossSection arg) const;</code><BR>
 * <code>typedef CrossSection ArgType</code><BR>
 * <code>typedef Energy ValType</code><BR>
 * <code>static ArgType aUnit(){return 1.*millibarn;} </code> <BR>
 * <code>static ValType vUnit(){return 1.*MeV;} </code> <BR>
 * <code> ... </code> <BR>
 * <code>}</code> <BR>
 * This can be facilitated by deriving from the GSLHelper struct. Which
 * implents the vUnit() and aUnit() methods using the baseunit static
 * method. Also the typedefs are written there.
 */
class GSLBisection : public Pointer::ReferenceCounted {

public:

  /**
   * Struct that is used to throw and catch GSL errors
   */
  struct GSLerror {};

  /**
   * Struct that is used to throw and catch GSL errors
   */
  struct IntervalError {};

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default Constructor
   */
  GSLBisection() : abserr_(0), relerr_(1.E-8), maxPoints_(100) {}

  /**
   * Specify all the parameters.
   * @param abserr Absolute error.
   * @param relerr Relative error.
   * @param max Maximum number of intervals
   */
  inline GSLBisection(double abserr, double relerr, int max) :
    abserr_(abserr), relerr_(relerr), maxPoints_(max) {}

  //@}

  /**
   * Function to overwrite the default GSL error handling
   */
  static void GSLsubstHandler(const char *, const char *, 
			      int, int){
    throw GSLerror();
  }

  /**
   * The result of the root finding.
   * @param function The integrand class that defines operator()
   * @param lower The lower limit of integration.
   * @param upper The upper limit of integration.
   */
  template <class T>
  inline typename T::ArgType value(const T & function, 
				   const typename T::ArgType lower,
				   const typename T::ArgType upper) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GSLBisection & operator=(const GSLBisection &);

private:

  /**
   * The parameters controlling the absolute error.
   */
  double abserr_;

  /**
   * The parameters controlling the relatve error.
   */
  double relerr_;

  /**
   * The maximum number of evaluations to use.
   */
  int maxPoints_;
};

}

#include "GSLBisection.tcc"

#endif /* HERWIG_GSLBisection_H */
