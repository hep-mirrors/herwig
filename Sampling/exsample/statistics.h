// -*- C++ -*-
//
// statistics.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2017 Simon Platzer -- simon.plaetzer@desy.de, The Herwig Collaboration
//
// ExSample is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_statistics_h_included
#define EXSAMPLE_statistics_h_included

#include "utility.h"

namespace exsample {

  /// \brief statistics is a helper class for keeping track of event
  /// generation statistics.
  class statistics {

  public:

    /// default constructor
    statistics();

    /// update the statistics for a weight encountered during
    /// presampling
    void presampled(double weight);

    /// indicate that a weight has been selected; optionally preven
    /// the weight from entering the caluclation of the integral
    void select(double weight,
		bool calculate_integral = true);

    /// indicate that a point has been accepted
    void accept(double weight) {
      ++accepted_;
      if (weight < 0) ++accepted_negative_;
    }

    /// reject a prviously accepted event
    void reject(double weight) {
      --accepted_;
      if (weight < 0) --accepted_negative_;
    }

    /// reset the statistics object
    void reset();

    /// return the integral's estimate and its uncertainty at the
    /// currently accumulated statistics
    std::pair<double,double> current() const;

    /// the average weight
    double average_weight() const { 
      return (n_iterations_ == 0 ? 0. : average_weight_/n_iterations_);
    }
    
    /// the average absolute weight
    double average_abs_weight() const { 
      return (n_iterations_ == 0 ? 0. : average_abs_weight_/n_iterations_);
    }
    
    /// the variance of the weight
    double average_weight_variance() const { 
      return (n_iterations_ == 0 ? 0. : average_weight_variance_/n_iterations_);
    }

    /// the number of points in this iteration
    unsigned long iteration_points() const { return iteration_points_; }

    /// the number of iterations
    unsigned long n_iterations() const { return n_iterations_; }
    
    /// the total number of attempted in this bin
    unsigned long attempted() const { return attempted_; }
    
    /// the total number of finally accepted events in this bin
    unsigned long accepted() const { return accepted_; }
    
    /// the total number of acceptet events with negative weights
    unsigned long accepted_negative() const { return accepted_negative_; }
    
    /// the sum of weights
    double sum_weights() const { return sum_weights_; }
    
    /// the sum of absolute values of the weights
    double sum_abs_weights() const { return sum_abs_weights_; }
    
    /// the sum of weights squared
    double sum_weights_squared() const { return sum_weights_squared_; }
    
    /// the maximum weight
    double max_weight() const { return max_weight_; }

  public:

    /// put ostream
    template<class OStream>
    void put(OStream& os) const;

    /// get from istream
    template<class IStream>
    void get(IStream& is);

  private:

    /// the average weight
    double average_weight_;
    
    /// the average absolute weight
    double average_abs_weight_;
    
    /// the variance of the weight
    double average_weight_variance_;

    /// the number of points in this iteration
    unsigned long iteration_points_;
    
    /// the total number of attempted in this bin
    unsigned long attempted_;
    
    /// the total number of finally accepted events in this bin
    unsigned long accepted_;
    
    /// the total number of acceptet events with negative weights
    unsigned long accepted_negative_;
    
    /// the sum of weights
    double sum_weights_;
    
    /// the sum of absolute values of the weights
    double sum_abs_weights_;
    
    /// the sum of weights squared
    double sum_weights_squared_;
    
    /// the maximum weight
    double max_weight_;

    /// the number of iterations used to calculate the integral
    unsigned long n_iterations_;
    
  };

}

#include "statistics.icc"

#endif // EXSAMPLE_statistics_h_included
