// -*- C++ -*-
//
// adaption_info.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2011 Simon Platzer -- simon.plaetzer@desy.de
//
// ExSample is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_adaption_info_h_included
#define EXSAMPLE_adaption_info_h_included

namespace exsample {

  /// \brief adaption_info is a container for parameters relevant to
  /// sampling and adaption.
  struct adaption_info {

    /// the default constructor
    adaption_info()
      : dimension(0), lower_left(), upper_right(),
	presampling_points(100000),
	histo_depth(2), adapt(), freeze_grid(0),
	maxtry(200000), efficiency_threshold(.9),
	gain_threshold(.1) {}

    /// the phasespace dimension
    std::size_t dimension;

    /// the lower left corner of the function's support
    std::vector<double> lower_left;

    /// the upper right corner of the function's support
    std::vector<double> upper_right;

    /// the number of presampling points
    unsigned long presampling_points;

    /// use 2^histo_depth bins in efficiency histograms
    std::size_t histo_depth;

    /// indicate which dimensions should be adapted
    std::vector<bool> adapt;

    /// the number of accepted events after the grid is frozen
    unsigned long freeze_grid;

    /// the maximum number of misses allowed
    unsigned long maxtry;

    /// the efficiency threshold below which splits are considered
    double efficiency_threshold;

    /// a minimum gain for splits to be performed
    double gain_threshold;

    /// put to ostream
    template<class OStream>
    void put(OStream& os) const;

    /// get from istream
    template<class IStream>
    void get(IStream& is);

  };

}

#include "adaption_info.icc"

#endif // EXSAMPLE_adaption_info_h_included
