// -*- C++ -*-
//
// selectors.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2011 Simon Platzer -- simon.plaetzer@desy.de
//
// ExSample is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_selectors_h_included
#define EXSAMPLE_selectors_h_included

#include "cell.h"

namespace exsample {

  /// \brief flat sampling selector
  template<class Random>
  struct sampling_selector {

    /// the default constructor
    sampling_selector()
      : rnd_gen(), compensate(false) {}

    /// the standard constructor
    explicit sampling_selector(const Random& r, bool comp = true)
      : rnd_gen(r), compensate(comp) {}

    /// return which of the children cells should be considered
    std::pair<bool,bool> use(cell& parent,
			     const cell& first_child,
			     const cell& second_child) const;

    /// return true, if the leaf cell should be considered
    bool use(cell& leaf) const;

    /// The random number generator to be used
    Random rnd_gen;

    /// Whether or not compensation is needed
    bool compensate;

  };

  /// \brief selector selecting only bins
  /// which contain the given parameter point
  class parametric_selector {

  public:

    /// the default constructor
    parametric_selector ()
      : point_(), sampled_variables_() {}

    /// construct from reference to point and flags to sample variables
    parametric_selector(std::vector<double> * point,
			const std::vector<bool>& sample)
      : point_(point), sampled_variables_(sample) {}

  public:

    /// return which of the children cells should be considered
    std::pair<bool,bool> use(const cell& parent,
			     const cell&,
			     const cell&) const;

    /// return true, if the leaf cell should be considered
    bool use(const cell&) const { return true; }

  private:

    /// the point chosen
    std::vector<double> * point_;

    /// flags for variables to be sampled
    std::vector<bool> sampled_variables_;

  };

  /// \brief sampling selector selecting only bins
  /// which contain the given parameter point
  template<class Random>
  class parametric_sampling_selector {

  public:

    /// the default constructor
    parametric_sampling_selector()
      : point_(), bin_id_(),
	sampled_variables_(), rnd_gen_(),
	compensate_(false) {}

    /// construct from reference to point, subtree hash, flags of
    /// variables to be sampled, and random number generator
    parametric_sampling_selector(std::vector<double> * p,
				 bit_container<parameter_hash_bits> * bin_id,
				 const std::vector<bool>& sample,
				 const Random& rnd_gen)
      : point_(p), bin_id_(bin_id),
	sampled_variables_(sample), rnd_gen_(rnd_gen),
	compensate_(false) {}

  public:

    /// return which of the children cells should be considered
    std::pair<bool,bool> use(cell& parent,
			     const cell& first_child,
			     const cell& second_child) const;

    /// return true, if the leaf cell should be considered
    bool use(cell& leaf) const;

    /// indicate that compensation is to take place
    void compensate (bool doit = true) { compensate_ = doit; }

  private:

    /// the point chosen
    std::vector<double> * point_;

    /// the corresponding bin id
    bit_container<parameter_hash_bits> * bin_id_;

    /// flags for variables to be sampled
    std::vector<bool> sampled_variables_;

    /// the random number generator
    Random rnd_gen_;

    /// wether or not compensation is needed
    bool compensate_;

  };

  /// \brief accessor returning the integral of a cell
  struct integral_accessor {    

    /// update and return the value
    double& set(cell& node) const {
      return node.integral();
    }

    /// get the value
    double get(const cell& node, bool) const {
      return node.integral();
    }

  };

  /// \brief accessor returning the number of missing events
  struct missing_accessor {    

    /// update and return the value
    int& set(cell& node) const {
      return node.missing_events();
    }

    /// get the value
    int get(const cell& node, bool) const {
      return node.missing_events();
    }

  };

  /// \brief accessor returning the number of missing events
  /// for given parameter bin id
  struct parametric_missing_accessor {    

    /// the default constructor
    parametric_missing_accessor ()
      : id_() {}

    /// construct from subtree hash to consider
    explicit parametric_missing_accessor (bit_container<parameter_hash_bits>* id)
      : id_ (id) {}

    /// update and return the value
    int& set(cell& node) const {
      return node.missing_events();
    }

    /// get the value
    int get(const cell& node, bool isleaf) const {
      if (isleaf)
	return node.info().parametric_missing(*id_);
      return node.missing_events();
    }

  private:

    /// the subtree hash to consider
    bit_container<parameter_hash_bits>* id_;

  };

}

#include "selectors.icc"

#endif // EXSAMPLE_selectors_h_included
