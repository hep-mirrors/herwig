// -*- C++ -*-
//
// generator.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2019 Simon Platzer -- simon.plaetzer@desy.de, The Herwig Collaboration
//
// ExSample is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_generator_h_included
#define EXSAMPLE_generator_h_included

#include "cell.h"
#include "selectors.h"
#include "statistics.h"
#include "binary_tree.h"

namespace exsample {

  /// \brief Exception thrown, if the generator has just changed its
  /// state. The attempt of generating an event should be repeated.
  struct generator_update{};

  /// \brief A generator for plain sampling and integrating
  template<class Function, class Random>
  class generator {

  public:

    /// default constructor
    generator()
      : function_(0), statistics_(), check_events_(0),
	adaption_info_(), root_cell_(),
	rnd_gen_(), did_split_(false), initialized_(false),
	last_cell_(), last_point_(), last_value_(0.),
	compensating_(false) {}

  public:

    /// initialize this generator
    template<class SlaveStatistics>
    void initialize(SlaveStatistics&);

    /// generate an event, returning
    /// the sign of the weight
    template<class SlaveStatistics>
    double generate(SlaveStatistics&);

    /// return the last sampled phase space point
    const std::vector<double>& last_point() const { return last_point_; }

    /// indicate that the last generated point has been rejected
    void reject() {
      statistics_.reject(last_value_);
      last_cell_->info().reject();
    }

    /// finalize this generator
    void finalize() {
      statistics_.reset();
    }

  public:

    /// return true, if this generator has been initialized
    bool initialized() const { return initialized_; }

    /// return true, if at least one split has been performed
    bool did_split() const { return did_split_; }

    /// access the function
    Function& function() { return *function_; }

    /// set the function
    void function(Function * f) { function_ = f; }

    /// return the statistics object
    const statistics& stats() const { return statistics_; }

    /// return the sampled volume
    double volume() const {
      return exsample::volume(adaption_info_.lower_left, adaption_info_.upper_right);
    }

    /// return the integral
    double integral() const {
      return volume() * statistics_.average_weight();
    }

    /// return the error on the integral
    double integral_uncertainty() const {
      return volume() * std::sqrt(statistics_.average_weight_variance());
    }

    /// return the integral
    double current_integral() const {
      return volume() * statistics_.current().first;
    }

    /// return the error on the integral
    double current_integral_uncertainty() const {
      return volume() * std::sqrt(statistics_.current().second);
    }

    /// return the variance of the integral estimate
    double integral_variance() const {
      return sqr(volume()) * statistics_.average_weight_variance();
    }

    /// access the adaption_info object
    adaption_info& sampling_parameters() { return adaption_info_; }

    /// return true, if still compensating
    bool compensating() const { return compensating_; }

  public:

    /// put to ostream
    template<class OStream>
    void put(OStream& os) const;

    /// get from istream
    template<class IStream>
    void get(IStream& is);

  private:

    /// check for and possibly split
    /// the last selected cell
    bool split();

    /// compensate the last selected cell indicating the
    /// value and position of the new overestimate
    void compensate();

    /// function to be sampled
    Function * function_;

    /// the global statistics object
    statistics statistics_;

    /// the number of events after which
    /// a cell is checked for splits
    unsigned long check_events_;

    /// the adaption info object
    adaption_info adaption_info_;

    /// the root cell
    binary_tree<cell> root_cell_;

    /// the random number generator to be used
    rnd_generator<Random> rnd_gen_;

    /// wether a split has already been performed
    bool did_split_;

    /// wether this generator has been initialized
    bool initialized_;

    /// the last selected cell
    binary_tree<cell>::iterator last_cell_;      

    /// the last sampled phasespace point
    std::vector<double> last_point_;

    /// the last function value
    double last_value_;

    /// wether or not we are compensating
    bool compensating_;

  };

}

#include "generator.icc"

#endif // EXSAMPLE_generator_h_included

