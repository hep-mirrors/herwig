// -*- C++ -*-
//
// Statistic.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Statistic_H
#define HERWIG_Statistic_H
#include <cmath>

//
// This is the declaration of the Statistic class.
//

namespace Herwig {

/**
 * The Statistic class is a simple class designed to 
 * store a variable for statistical analysis
 */
class Statistic {

public:

  /**
   * The default constructor.
   */
  Statistic() : _n(0), _xsum(0.), _x2sum(0.),
		_min(-1e100), _max(1e100) {}

  /**
   *  The minimum value
   */
  double minimum() const { return _min; }

  /**
   *  The maximum value
   */
  double maximum() const { return _max; }

  /**
   *  Operator to add another point
   */
  void operator+=(double input) 
  {
    ++_n;
    _xsum  += input;
    _x2sum += input * input;
    if (_min > input || _n == 1) _min = input;
    if (_max < input || _n == 1) _max = input;
  }

  /**
   *  Number of points
   */
  unsigned int numberOfPoints() const { return _n; }
  
  /**
   *  Mean
   */
  double mean() const 
  {
    return _n > 0  ?  _xsum / _n : 0.; 
  }

  /**
   *  Error on the mean estimate. Needed for example for Profile
   *  histograms, where this should be used to compute a chi2
   *  or significance level of deviation to data, rather than stdDeV.
   *  This is obvious because the error on the estimate should go to 
   *  zero for N -> infinity.
   */
  double mean_stdDev() const { return std::sqrt(mean_var()); }

  /**
   *  Variance on the mean estimate. Needed for example for Profile
   *  histograms, where this should be used to compute a chi2
   *  or significance level of deviation to data, rather than stdDeV
   *  This is obvious because the error on the estimate should go to 
   *  zero for N -> infinity.
   */
  double mean_var() const 
  {
    return _n > 1  ?  var() / _n : 0.; 
  }

  /**
   *  Standard Deviation
   */
  double stdDev() const { return std::sqrt(var()); }

  /**
   *  Variance
   */
  double var() const 
  { 
    return _n > 1  ?  ( _x2sum - _xsum*_xsum/_n ) / ( _n - 1 ) : 0.; 
  }

  /**
   *  Total entry
   */
  double total() const { return _xsum; }

private:

  /**
   *   Number of entries
   */
  unsigned int _n;

  /**
   *  Sum of the values
   */ 
  double _xsum;

  /**
   *  Sum of the squares of the values
   */
  double _x2sum;

  /**
   *  The minimum value
   */
  double _min;
  
  /**
   *  The maximum value
   */
  double _max;
};

}

#endif /* HERWIG_Statistic_H */
