// -*- C++ -*-
/*
 * parameters.h
 * This file is indented to contain program wide parameters.
 * Created on: Oct 26, 2010
 * Author: Malin Sjodahl
 */

#ifndef COLORFULL_Parameters_h
#define COLORFULL_Parameters_h

namespace ColorFull {

/// This sets the accuracy for numerical comparisons.
/// It is used for numerical comparison for example
/// when the symmetry of a matrix is checked,
/// when it's checked if a matrix is diagonal
/// when realness is checked
/// and when turning complex numbers into real numbers.
/// It is also used by << for dmatr and
/// the Polynomial member function simplify .
const double accuracy=0.0000000000001;

}

#endif /* COLORFULL_Parameters_h_ */
