// -*- C++ -*-
//
// linear_interpolator.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2017 Simon Platzer -- simon.plaetzer@desy.de, The Herwig Collaboration
//
// ExSample is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_linear_interpolator_h_included
#define EXSAMPLE_linear_interpolator_h_included

#include "utility.h"

namespace exsample {

  /// \brief Exception thrown, if inversion of the interpolation has
  /// no solution
  struct inversion_has_no_solution { };

  /// \brief Exception thrown, if a constant piece of the
  /// interpolation has been hit
  struct constant_interpolation {

    /// standard constructor
    constant_interpolation(double xlow,
			   double xhigh,
			   double h)
      : range(xlow,xhigh), value(h) {}

    /// The range in which the interpolation is constant
    std::pair<double,double> range;

    /// The value of the interpolation in the range
    double value;

  };

  /// \brief A linear interpolation allowing for inversion of the
  /// linear interpolation
  class linear_interpolator {

  public:

    /// default constructor
    linear_interpolator();

    /// construct a linear interpolator given the
    /// map of interpolation points and values
    explicit linear_interpolator(const std::map<double,double>& points);

  public:

    /// return the interpolations value at the given point
    double operator()(double x) const;

    /// return the range of the interpolation
    std::pair<double,double> range() const { return range_; }

  public:

    /// return true, if an inverse exists for
    /// the given value
    bool invertible(double f) const {
      return (range_.first <= f &&
	      range_.second >= f);
    }

    /// return the inverse, assuming the inverse
    /// is unique
    double unique_inverse(double f) const;

    /// access the interpolation map
    std::map<double,double>& interpolation() {
      return interpolation_;
    }

    /// set the value at the given point
    void set_interpolation(double point, double value);

    /// reset after interpolation points have been changed
    void reset();

  public:

    /// put to ostream
    template<class OStream>
    void put(OStream& os) const;

    /// get from istream
    template<class IStream>
    void get(IStream& is);

  private:

    /// map points to values
    std::map<double, double> interpolation_;

    /// the range of the interpolation
    std::pair<double,double> range_;

  };

  /// ostream operator for linear_interpolator objects
  template<class OStream>
  inline OStream& operator<<(OStream& os, const linear_interpolator& ip) {
    ip.put(os);
    return os;
  }

  /// istream operator for linear_interpolator objects
  template<class IStream>
  inline IStream& operator>>(IStream& is, linear_interpolator& ip) {
    ip.get(is);
    return is;
  }

}

#include "linear_interpolator.icc"

#endif // EXSAMPLE_linear_interpolator_h_included
