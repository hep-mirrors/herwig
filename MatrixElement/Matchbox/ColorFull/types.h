// -*- C++ -*-
/*
 * types.h
 * Contains the definition of a few types used by ColorFull.
 * Created on: Jul 7, 2010
 * Author: Malin Sjodahl
 */

#ifndef COLORFULL_types_h
#define COLORFULL_types_h

#include <complex>
#include <vector>


namespace ColorFull {
// ///////////////  Define a few basic types  //////////////////

/// A complex number, needed as amplitudes in general are complex.
typedef std::complex< double > cnum;

/// A vector of complex numbers.
typedef std::vector<cnum>  cvec;

/// A matrix of complex numbers.
typedef std::vector< std::vector <cnum> > cmatr;

/// A vector of double numbers.
typedef std::vector<double>  dvec;

/// A matrix of double numbers
typedef std::vector< std::vector <double> > dmatr;

/// An unsigned int.
typedef unsigned int uint;

}

#endif /* COLORFULL_types_h */
