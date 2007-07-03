// -*- C++ -*-
#ifndef HERWIG_Math_H
#define HERWIG_Math_H

#include <cmath>
#include "ThePEG/Config/Complex.h"
#include "ThePEG/Utilities/Math.h"

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
   * Fold angles into the range (-Pi,Pi)
   */
  inline double angleMinusPiToPi(double angle);

  /**
   * Fold angles into the range (0,2 Pi)
   */
  inline double angleZeroTo2Pi(double angle);
}

}


#include "Math.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Math.tcc"
#endif

#endif /* HERWIG_Math_H */
