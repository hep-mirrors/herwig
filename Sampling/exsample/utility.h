// -*- C++ -*-
//
// utility.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2011 Simon Platzer -- simon.plaetzer@desy.de
//
// ExSample is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_utility_h_included
#define EXSAMPLE_utility_h_included

#include "config.h"
#include <memory>

namespace exsample {

  /// \brief separate quantities written to an ostream
  template<class OStream>
  struct ostream_traits {

    /// put the separator to the ostream
    static void separator(OStream& os) { os << "\n"; }

  };

#ifdef EXSAMPLE_has_ThePEG

  /// \brief separate quantities written to a ThePEG::PersistentOStream
  template<>
  struct ostream_traits<ThePEG::PersistentOStream> {

    /// put the separator to the ostream
    static void separator(ThePEG::PersistentOStream&) { }

  };  

#endif // EXSAMPLE_has_ThePEG

  /// \brief Compile time conversion of unsigned long to bool
  template<unsigned long>
  struct static_binary {
    enum { value = 1 };
  };

  /// \brief Compile time conversion of unsigned long to bool
  template<>
  struct static_binary<0> {
    enum { value = 0 };
  };

  /// \brief Fixed-size, packed vector of bools
  template<unsigned long bits>
  struct bit_container {

    enum {
      /// the number of bits contained
      n_bits = bits,
      /// the number of bits fitting in a unsigned long
      uint_bits = CHAR_BIT * sizeof(unsigned long),
      /// the number of unsigned long segments needed
      n_segments = bits / uint_bits + static_binary<bits % uint_bits>::value
    };


    /// the default constructor
    bit_container() {
      for (std::size_t i = 0; i < n_segments; ++i)
	segments[i] = 0;
    }

    /// put all values to false
    void reset() {
      for (std::size_t i = 0; i < n_segments; ++i)
	segments[i] = 0;
    }

    /// compare for equality
    bool operator==(const bit_container& x) const {
      for (std::size_t i = 0; i < n_segments; ++i)
	if(segments[i] != x.segments[i])
	  return false;
      return true;
    }

    /// compare for ordering
    bool operator<(const bit_container& x) const {
      for (std::size_t i = 0; i < n_segments; ++i)
	if(segments[i] != x.segments[i])
	  return (segments[i] < x.segments[i]);
      return false;
    }

    /// set the k'th bit
    void bit(unsigned long k, bool value) {
      assert(k<n_bits);
      if (value)
	segments[n_segments-k/uint_bits-1] |= (1ul << (k % uint_bits));
      else
	segments[n_segments-k/uint_bits-1] &= ~(1ul << (k % uint_bits));
    }

    /// get the k'th bit
    bool bit(unsigned long k) const {
      assert(k<n_bits);
      return (segments[n_segments-k/uint_bits-1] & (1ul << (k % uint_bits)));
    }

    /// print to ostream
    template<class OStream>
    void dump(OStream& os) const {
      for ( unsigned int k = 0; k < n_segments; ++k )
	os << segments[k] << " ";
    }

    /// put to ostream
    template<class OStream>
    void put(OStream& os) const {
      for ( size_t k = 0; k < n_segments; ++k ) {
	os << segments[k];
	ostream_traits<OStream>::separator(os);
      }
    }

    /// get from istream
    template<class IStream>
    void get(IStream& is) {
      for ( size_t k = 0; k < n_segments; ++k ) {
	is >> segments[k];
      }
    }

  private:

    /// segments needed to keep the hash value
    unsigned long segments[n_segments];

  };

  /// \brief square a number
  template<class T>
  T sqr(T x) {
    return x*x;
  }

  /// \brief cube a number
  template<class T>
  T cube(T x) {
    return x*x*x;
  }


  /// \brief Round a floating point value to an integer value of the
  /// same type.
  template<class T>
  T round(T x) {
    T f = std::floor(x);
    T c = std::ceil(x);
    if (x < (f+c)/2.)
      return f;
    return c;
  }

  /// \brief Calculate fast powers of two.
  inline std::size_t two_to(std::size_t n) {
    assert(n <= sizeof(std::size_t)*CHAR_BIT);
    return (1 << n);
  }

  /// \brief Fast, zero memory-overhead one-dimensional
  /// histogram with 2^n equally spaced bins
  template<class Statistics>
  struct fast_small_histogram {

    /// default constructor
    fast_small_histogram()
      : depth(0), bins(0) {}

    /// copy constructor
    fast_small_histogram(const fast_small_histogram& x)
      : depth(x.depth), bins(0) {
      if (x.bins) {
	bins.reset(new Statistics[two_to(depth)]);
	for(std::size_t k = 0; k < two_to(depth); ++k)
	  bins[k] = x.bins[k];
      }
    }

    /// assignment
    fast_small_histogram& operator=(const fast_small_histogram& x) {
      if (&x == this)
	return *this;
      depth = x.depth;
      bins.reset(0);
      if (x.bins) {
	bins.reset(new Statistics[two_to(depth)]);
	for(std::size_t k = 0; k < two_to(depth); ++k)
	  bins[k] = x.bins[k];
      }
      return *this;
    }

    /// construct from depth d, creating 2^d bins
    explicit fast_small_histogram(std::size_t d)
      : depth(d), bins(0) {
      bins.reset(new Statistics[two_to(d)]);
    }

    /// return the bin from event belongs to given outer boundaries
    Statistics& bin(double lower,
		    double upper,
		    double event) {
      double thelower = lower;
      double theupper = upper;
      std::size_t bindex = 0;
      std::size_t current_depth = 0;
      while (true) {
	double cut
	  = (thelower+theupper)/2.;
	if (event < cut) {
	  theupper = cut;
	} else {
	  thelower = cut;
	  bindex += two_to(depth-current_depth-1);
	}
	if(++current_depth == depth)
	  break;
      }
      return bins[bindex];
    }

    /// the depth, defining a histogram of 2^depth bins
    std::size_t depth;

    /// the contained statistics objects
    boost::scoped_array<Statistics> bins;

    /// put histogram to an ostream
    template<class OStream>
    void put(OStream& os) const {
      os << depth;
      ostream_traits<OStream>::separator(os);
      for (std::size_t k = 0; k < two_to(depth); ++k) {
	bins[k].put(os);
      }
    }

    /// get histogram from an istream
    template<class IStream>
    void get(IStream& is) {
      is >> depth;
      bins.reset(new Statistics[two_to(depth)]);
      for(std::size_t k = 0; k < two_to(depth); ++k) {
	bins[k].get(is);
      }
    }
     
  };
 
  /// \brief Generalize the transform algorithm to only apply
  /// depending on a range of flags accompanying the input range
  template<class FirstInputIterator, 
	   class SecondInputIterator, 
	   class FlagIterator,
	   class OutputIterator,
	   class BinaryOperation>
  OutputIterator conditional_transform(FirstInputIterator first1,
				       FirstInputIterator last1,
				       SecondInputIterator first2,
				       FlagIterator firstf,
				       OutputIterator result,
				       BinaryOperation binary_op) {
    for (; first1 != last1; ++first1, ++first2, ++firstf, ++result)
      if (*firstf)
	*result = binary_op(*first1, *first2);
    return result;
  }

  /// \brief calculate a volume given lower left and upper right
  /// corners
  inline double volume(const std::vector<double>& lower_left,
		       const std::vector<double>& upper_right) {
    std::vector<double> delta;
    std::transform(upper_right.begin(),upper_right.end(),
		   lower_left.begin(),std::back_inserter(delta),
		   std::minus<double>());
    return
      std::accumulate(delta.begin(),delta.end(),1.,std::multiplies<double>());
  }

  /// \brief calculate a volume given lower left and upper right
  /// corners, taking into account only part of the dimensions, which
  /// are flagged with true in the correspponding random access
  /// container
  inline double volume(const std::vector<double>& lower_left,
		       const std::vector<double>& upper_right,
		       const std::vector<bool>& flags) {
    std::vector<double> delta;
    conditional_transform(upper_right.begin(),upper_right.end(),
			  lower_left.begin(),flags.begin(),
			  std::back_inserter(delta),
			  std::minus<double>());
    return
      std::accumulate(delta.begin(),delta.end(),1.,std::multiplies<double>());
  }


  /// \brief Exception thrown if the maximum number of attempts to
  /// select a cell has been reached.
  struct selection_maxtry{};

  /// \brief Exception thrown, if the maximum number of misses has
  /// been reached.
  struct hit_and_miss_maxtry{};

  /// \brief Random generator traits.
  template<class Random>
  struct rnd_generator {

    ///Generate uniform random number on [0,1]
    double operator()() const {
      return Random::rnd();
    }

    ///Generate uniform random number on [0,a]
    double operator()(double a) const {
      return a*Random::rnd();
    }

    ///Generate uniform random number on [a,b]
    double operator()(double a, double b) const {
      return (a + (b-a)*Random::rnd());
    }

  };

}

#endif // EXSAMPLE_utility_h_included
