// -*- C++ -*-
#ifndef HERWIG_Interpolator_H
#define HERWIG_Interpolator_H
//
// This is the declaration of the <!id>Interpolator<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class implments a polynominal interpolation of a table of values, it is
//  based on the interpolation code in FORTRAN HERWIG and uses the GenericFunctions
//  classes of CLHEP
//  
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 
//  Author: Peter Richardson
//

#include <vector>
#include <cmath>
#include "CLHEP/GenericFunctions/AbsFunction.hh"
// #include "Interpolator.fh"
// #include "Interpolator.xh"

namespace Herwig {

using namespace Genfun;
using std::vector;

class Interpolator: public Genfun:: AbsFunction {

  FUNCTION_OBJECT_DEF(Interpolator)
    
public:
  
  // constructor with data as vectors
  Interpolator(vector<double> f, vector<double> x, int order);

  // Copy constructor
  Interpolator(const Interpolator &right);
  
  virtual ~Interpolator();
  // Standard ctors and dtor.
  
  // Retreive function value
  virtual double operator ()(double argument) const;
  virtual double operator ()(const Argument & a) const {return operator() (a[0]);}
  
private:
  
  // It is illegal to assign a function
  const Interpolator & operator=(const Interpolator &right);
  
private:
  
  vector<double> _xval;
  // the x values
  vector<double> _fun;
  // the function values
  unsigned int _order;
  // the order of interpolation
};
  
}

// CLASSDOC OFF

#endif /* HERWIG_Interpolator_H */
