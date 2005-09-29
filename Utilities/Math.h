// -*- C++ -*-
#ifndef HERWIG_Math_H
#define HERWIG_Math_H

#include <cmath>
#include "ThePEG/Config/Complex.h"
#include "ThePEG/Utilities/Math.h"

namespace Herwig {
using ThePEG::Complex;

/** The Math namespace includes the declaration of some useful
 *  mathematical functions. */
namespace Math {

  /**
   * The dilog function taken from FORTRAN Herwig
   */
  Complex Li2(Complex);

}

}


#include "Math.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Math.tcc"
#endif

#endif /* HERWIG_Math_H */
