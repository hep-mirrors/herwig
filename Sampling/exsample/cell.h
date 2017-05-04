// -*- C++ -*-
//
// cell.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2017 Simon Platzer -- simon.plaetzer@desy.de, The Herwig Collaboration
//
// ExSample is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_cell_h_included
#define EXSAMPLE_cell_h_included

#include "utility.h"
#include "adaption_info.h"
#include "statistics.h"

namespace exsample {

  /// \brief Information contained in a leaf cell
  class cell_info {

  public:

    /// the default constructor
    cell_info();

    /// construct from boundaries and adaption info
    cell_info(const std::vector<double>& ll,
	      const std::vector<double>& ur,
	      const adaption_info& ainfo);

    /// construct from boundaries, flags for variables to be sampled,
    /// and adaption info
    cell_info(const std::vector<double>& ll,
	      const std::vector<double>& ur,
	      const std::vector<bool>& sampled_variables,
	      const adaption_info& ainfo);

  public:

    /// generate a flat trial point in this cell
    template<class Random>
    void select(Random&,
		std::vector<double>&);

    /// generate a flat trial point in this cell
    /// only for the variables falgged as true
    template<class Random>
    void select(Random&,
		std::vector<double>&,
		const std::vector<bool>&);

    /// indicate a function value for the given point
    void selected(const std::vector<double>&,
		  double,
		  const adaption_info&);

    /// indicate that a point has been
    /// accepted in this cell
    void accept() { ++accepted_; }

    /// reject a previously accepted event
    void reject() { --accepted_; }

  public:

    /// return true, if below efficiency threshold
    bool bad(const adaption_info& ainfo) const {
      return ((static_cast<double>(accepted_)/static_cast<double>(attempted_)) < 
	      ainfo.efficiency_threshold);
    }

    /// suggest a split and indicate wether it is worth
    /// to be performed
    std::pair<std::size_t,double> get_split(const adaption_info&,
					    bool&) const;

    /// explore this cell performing a flat sampling,
    /// updating the given statistics object and pre-filling
    /// the efficiency histogram by a trial unweighting
    template<class Random, class Function, class SlaveStatistics>
    void explore(Random&, const adaption_info&, Function*, statistics*,
		 SlaveStatistics& opt, double detuning);

    /// explore this cell in a more refined way, which
    /// is however not suited for already calculating integrals
    /// and stuff
    template<class Random, class Function>
    void explore(Random&, const adaption_info&, Function*,
		 double detuning);

  public:

    /// get the current overestimate
    double overestimate() const { return overestimate_; }

    /// return the position of the last maximum
    const std::vector<double>& last_max_position() const { return last_max_position_; }

    /// set the current overestimate and maximum position
    void overestimate(double v, const std::vector<double>& pos,
		      double detuning) { 
      overestimate_ = detuning * v;
      last_max_position_ = pos;
    }

    /// get the volume
    double volume() const { return volume_; }

    /// get the lower left corner
    const std::vector<double>& lower_left() const { return lower_left_; }

    /// get the upper right corner
    const std::vector<double>& upper_right() const { return upper_right_; }

    /// get the number of attempted events
    unsigned long attempted() const { return attempted_; }

    /// get the number of accepted events
    unsigned long accepted() const { return accepted_; }

  public:

    /// return the number of missing events
    /// for the given parameter bin id
    int parametric_missing(const bit_container<parameter_hash_bits>& id) const;

    /// set the number of missing events
    /// for the given parameter bin id
    void parametric_missing(const bit_container<parameter_hash_bits>& id, int n);

    /// increase to the number of missing events
    /// for the given parameter bin id
    void increase_parametric_missing(const bit_container<parameter_hash_bits>& id);

    /// decrease to the number of missing events
    /// for the given parameter bin id
    void decrease_parametric_missing(const bit_container<parameter_hash_bits>& id);

    /// return true, if the cell is compensating in
    /// at least one parameter bin
    bool parametric_compensating() const {
      return !parametric_missing_map_.empty();
    }

    /// return true, if the cell contains the
    /// indicated parameter point
    bool contains_parameter(const std::vector<double>& point,
			    const std::vector<bool>& sampled) const;

  public:

    /// put to ostream
    template<class OStream>
    void put(OStream& os) const;

    /// get from istream
    template<class IStream>
    void get(IStream& is);

  private:

    /// the value of the overestimate in this cell
    double overestimate_;

    /// the volume of this cell
    double volume_;

    /// the lower left corner of this cell
    std::vector<double> lower_left_;

    /// the upper right corner of this cell
    std::vector<double> upper_right_;

    /// midpoint of this cell
    std::vector<double> mid_point_;

    /// the position of the last encountered
    /// maximum in this cell
    std::vector<double> last_max_position_;

    /// left-right statistics of average weight
    std::vector<std::pair<double,double> > avg_weight_;

    /// the number of attempts in this cell
    unsigned long attempted_;

    /// the number of accepted events in this cell
    unsigned long accepted_;

    /// an optional map of parameter bin ids
    /// to the number of missing events
    std::map<bit_container<parameter_hash_bits>,int> parametric_missing_map_;

  };

  /// \brief the general cell class
  class cell {

  public:

    /// default constructor
    cell();

    /// construct from boundaries and adaption info
    cell(const std::vector<double>& ll,
	 const std::vector<double>& ur,
	 const adaption_info& ainfo);

    /// construct from boundaries, flags for variables to be sampled,
    /// and adaption info
    cell(const std::vector<double>& ll,
	 const std::vector<double>& ur,
	 const std::vector<bool>& sampled_variables,
	 const adaption_info& ainfo);

    /// copy constructor
    cell(const cell& x);

    /// assignment
    cell& operator=(const cell& x);

  public:

    /// split this cell, exploring the
    /// child not containing the current overestimate
    template<class Random, class Function>
    std::pair<cell,cell> split(std::pair<std::size_t,double> split_d,
			       Random& rnd_gen,
			       Function* f,
			       const adaption_info& ainfo,
			       const std::vector<bool>& sampled = 
			       std::vector<bool>(),
			       double detuning = 1.0);

  public:

    /// return the split dimension
    std::size_t split_dimension() const { return split_dimension_; }

    /// return the split value
    double split_point() const { return split_point_; }

    /// return the integral
    double integral() const { return integral_; }

    /// access the integral
    double& integral() { return integral_; }

    /// set the integral
    void integral(double v) { integral_ = v; }

    /// access the number of missing events
    int& missing_events() { return missing_events_; }

    /// return the number of missing events
    int missing_events() const { return missing_events_; }

    /// set the number of missing events
    void missing_events(int n) { missing_events_ = n; }

    /// access the cell_info object
    cell_info& info() { assert(cell_info_); return *cell_info_; }

    /// return the cell_info object
    const cell_info& info() const { assert(cell_info_); return *cell_info_; }

  public:

    /// put to ostream
    template<class OStream>
    void put(OStream& os) const;

    /// get from istream
    template<class IStream>
    void get(IStream& is);

  private:

    /// the dimension along this cell
    /// was split
    std::size_t split_dimension_;

    /// the value, where this cell was split
    double split_point_;

    /// the integral of the absolute value
    /// of the overestimate over all the
    /// children cells
    double integral_;

    /// the number of missing events in this cell
    int missing_events_;

    /// a pointer to the cell info object,
    /// if this is a leaf cell
    std::unique_ptr<cell_info> cell_info_;


  };

}

#include "cell.icc"

#endif // EXSAMPLE_cell_h_included
