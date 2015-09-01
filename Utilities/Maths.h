// -*- C++ -*-
//
// Maths.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Math_H
#define HERWIG_Math_H

#include <cmath>
#include <vector>
#include "ThePEG/Config/Complex.h"
#include "ThePEG/Utilities/Maths.h"

namespace Herwig {
using ThePEG::Complex;
using std::complex;

/** The Math namespace includes the declaration of some useful
 *  mathematical functions. */
namespace Math {

  /**
   * The dilog function taken from FORTRAN Herwig
   */
  Complex Li2(Complex);

  /**
   * The real part of the dilog function taken from FORTRAN Herwig
   */
  long double ReLi2(long double);
  
  /**
   * Fold angles into the range (0,2 Pi)
   */
  inline double angleZeroTo2Pi(double angle) {
    double ret = fmod(angle, 2 * M_PI);
    if (ret < 0) ret += 2 * M_PI;
    return ret;
  }

  /**
   * Fold angles into the range (-Pi,Pi)
   */
  inline double angleMinusPiToPi(double angle) {
    double ret = angleZeroTo2Pi(angle);
    if (ret > M_PI) ret -= 2 * M_PI;
    return ret;
  }

  /**
   * Calculates the (lower) median of a vector of T objects. T has to be
   * comparable, i.e. T::operator< must be defined.
   */
  template <typename T>
  inline T median(std::vector<T> v) {
    if (v.empty()) return T();
    sort ( v.begin(), v.end() );
    const size_t N = v.size();
    return (N % 2) ? v.at((N+1)/2 - 1) : v.at(N/2 - 1);
  }

}

}

#endif /* HERWIG_Math_H */
