// -*- C++ -*-
//
// exponential_generator.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2011 Simon Platzer -- simon.plaetzer@desy.de
//
// ExSample is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_exponential_generator_h_included
#define EXSAMPLE_exponential_generator_h_included

#include "cell.h"
#include "selectors.h"
#include "statistics.h"
#include "linear_interpolator.h"
#include "binary_tree.h"

namespace exsample {

  /// \brief Exception thrown, if the exponential_generator has just changed its
  /// state. The attempt of generating an event should be repeated.
  struct exponential_regenerate{};

  /// \brief The generator for sudakov-type distributions.
  template<class Function, class Random>
  class exponential_generator {

  public:

    /// default constructor
    exponential_generator()
      : function_(0), check_events_(0), adaption_info_(), root_cell_(),
	rnd_gen_(), did_split_(false), initialized_(false),
	evolution_variable_(0), evolution_cutoff_(0.),
	sample_variables_(), sample_other_variables_(),
	parameter_splits_(),
	last_cell_(), last_point_(), last_value_(0.),
	last_parameter_bin_(), exponents_(),
	last_exponent_integrand_(),
	last_exponent_(), compensating_(false),
        integral_accessor_(), missing_accessor_(),
	parametric_selector_(), exponent_selector_(),
	parametric_sampler_(), attempts_(0), accepts_(0),
      splits_(0), docompensate_(false), detuning_(1.0) {}

  public:

    /// initialize this generator
    void initialize();

    /// finalize this generator
    void finalize() {}

    /// generate an event, returning
    /// the sign of the weight or zero
    /// for an event below the evolution cutoff
    double generate();

    /// generate an event, returning
    /// the sign of the weight or zero
    /// for an event below the evolution cutoff
    double generate(double cutoff) {
      double oldcut = evolution_cutoff_;
      evolution_cutoff_ = cutoff;
      double w = 0.0;
      try {
	w = generate();
      } catch(...) {
	evolution_cutoff_ = oldcut;
	throw;
      }
      evolution_cutoff_ = oldcut;
      return w;
    }

    /// return the last sampled phase space point
    const std::vector<double>& last_point() const { return last_point_; }

    /// return the last evaluated function
    double last_value() const { return last_value_; }

    /// indicate that the last generated point has been rejected
    void reject() {
      last_cell_->info().reject();
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

    /// access the adaption_info object
    adaption_info& sampling_parameters() { return adaption_info_; }

    /// indicate, if compensation should be applied
    void docompensate(bool yes = true) { docompensate_ = yes; }

    /// set the detuning parameter
    void detuning(double val) { detuning_ = val; }

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

    /// get the projection of the density integrating over every
    /// variable to be sampled, except the evolution variable for the
    /// indicated parameter point.  the k'th entry in
    /// last_exponent_integrand_ is the value in the evolution
    /// variable bin from evolution_splits_[k] to
    /// evolution_splits_[k+1]
    void get_exponent();

    /// compensate 
    void compensate();

    /// get all parameter points to build
    /// all possible sub tree hashes
    std::set<std::vector<double> > parameter_points();

    /// get all parameter points to build
    /// all possible sub tree hashes
    void recursive_parameter_points(std::set<std::vector<double> >&,
				    std::vector<double>&,
				    size_t);

    /// function to be sampled
    Function * function_;

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

    /// the position of the evolution variable
    std::size_t evolution_variable_;

    /// the cutoff on the evolution variable
    double evolution_cutoff_;

    /// flags of variables to be sampled
    /// including the evolution variable
    std::vector<bool> sample_variables_;

    /// flags of variables to be sampled
    /// excluding the evolution variable
    std::vector<bool> sample_other_variables_;

    /// the splits in any parameter done so far
    /// (including the evolution variable)
    std::map<std::size_t,std::vector<double> > parameter_splits_;

    /// the last selected cell
    binary_tree<cell>::iterator last_cell_;      

    /// the last sampled phasespace point
    std::vector<double> last_point_;

    /// the last function value
    double last_value_;

    /// the last parameter bin id
    bit_container<parameter_hash_bits> last_parameter_bin_;

    /// map parameter bin ids to exponent interpolations
    std::map<bit_container<parameter_hash_bits>,linear_interpolator > exponents_;

    /// the last exponent integrand
    std::vector<double> last_exponent_integrand_;

    /// the last exponent
    std::map<bit_container<parameter_hash_bits>,linear_interpolator >::iterator last_exponent_;

    /// wether or not we are compensating
    bool compensating_;

    /// the integral accessor to be used
    integral_accessor integral_accessor_;

    /// the missing events accessor to be used
    parametric_missing_accessor missing_accessor_;

    /// the parametric selector to be used
    parametric_selector parametric_selector_;

    /// the parametric selector to be used for parameter bins
    parametric_selector exponent_selector_;

    /// the parametric sampler to be used
    parametric_sampling_selector<rnd_generator<Random> > parametric_sampler_;

    /// the number of trials in the veto loo so far
    unsigned long attempts_;

    /// the number of accepted events so far
    unsigned long accepts_;

    /// number of splits done
    unsigned long splits_;

    /// true, if compensation should be applied
    bool docompensate_;

    /// a detuning factor to be applied to the overestimate
    double detuning_;

  };

}

#include "exponential_generator.icc"

#endif // EXSAMPLE_exponential_generator_h_included
