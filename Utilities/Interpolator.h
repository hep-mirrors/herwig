// -*- C++ -*-
#ifndef HERWIG_Interpolator_H
#define HERWIG_Interpolator_H
//
// This is the declaration of the Interpolator class.

#include <vector>
#include <cmath>
#include "CLHEP/GenericFunctions/AbsFunction.hh"
// #include "Interpolator.fh"
// #include "Interpolator.xh"

namespace Herwig {

using namespace Genfun;
using std::vector;

/** \ingroup Utilities
 *  \author Peter Richardson
 *
 *  This class implments a polynominal interpolation of a table of values, it is
 *  based on the interpolation code in FORTRAN HERWIG and uses the GenericFunctions
 *  classes of CLHEP.
 */
class Interpolator: public Genfun:: AbsFunction {

  FUNCTION_OBJECT_DEF(Interpolator)
    
public:
  
  /**
   * Constructor with data as vectors.
   */
  Interpolator(vector<double> f, vector<double> x, int order);
   
  /**
   * Copy constructor
   */
  Interpolator(const Interpolator & x) 
    : Genfun::AbsFunction(), _xval(x._xval), _fun(x._fun), _order(x._order) {}
  // this one's required because CLHEP::Genfun::AbsFunction 
  // doesn't implement its copy constructor. Aaaargh!
  // Yes, -Wextra will warn about this.

 
  /**
   * Retreive function value.
   */
  virtual double operator ()(double argument) const;
  
  /**
   * Retreive function value.
   */
  virtual double operator ()(const Argument & a) const {return operator() (a[0]);}
  
private:
  
  /**
   * It is illegal to assign a function.
   */
  const Interpolator & operator=(const Interpolator &right);
  
private:
  
  /**
   * The x values.
   */
  vector<double> _xval;

  /**
   * the function values.
   */
  vector<double> _fun;

  /**
   * the order of interpolation.
   */
  unsigned int _order;

};
  
}

#endif /* HERWIG_Interpolator_H */
