// -*- C++ -*-
//
// GeneralStatictis.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_BinnedStatistics_H
#define Herwig_BinnedStatistics_H
//
// This is the declaration of the BinnedStatistics class.
//

#include "GeneralStatistics.h"
#include "ThePEG/Repository/UseRandom.h"

#include <boost/utility.hpp>

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Aka histogram, yet not intented for analyses.
 *
 */
class BinnedStatistics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  BinnedStatistics() 
    : lastPoint(0.), lastStatistics(0), theWeightThreshold(0.001) {}

  /**
   * The standard constructor.
   */
  BinnedStatistics(unsigned int bins, double threshold = 0.001) 
    : lastPoint(0.), lastStatistics(0) {
    initialize(bins);
    theWeightThreshold = threshold;
  }

  /**
   * The destructor.
   */
  virtual ~BinnedStatistics();
  //@}

public:

  /**
   * Sample a point and return its weight to be divided out as a bias.
   */
  double sample(double& point) {
    const pair<double,double>& range =
      selectorMap.upper_bound(UseRandom::rnd())->second;
    lastPoint = UseRandom::rnd(range.first,range.second);
    point = lastPoint;
    lastStatistics = &(statisticsMap.upper_bound(lastPoint)->second);
    double weight = weightMap.upper_bound(lastPoint)->second;
    return 1./weight;
  }

  /**
   * Get a bin corresponding to a given point.
   */
  void bin(double point) {
    lastPoint = point;
    lastStatistics = &(statisticsMap.upper_bound(lastPoint)->second);
  }

  /**
   * Select the last sampled point with a given weight.
   */
  void select(double w) {
    lastStatistics->select(w);
  }

  /**
   * Accept the last sampled point.
   */
  void accept() {
    lastStatistics->accept();
  }

  /**
   * Reject the last sampled point.
   */
  void reject() {
    lastStatistics->reject();
  }

  /**
   * Initialize with flat sampling over the complete range,
   * using the given number of bins to accumulate statistics.
   */
  void initialize(unsigned int bins);

  /**
   * Return the bins.
   */
  const map<double,GeneralStatistics>& statistics() const {
    return statisticsMap;
  }

  /**
   * Update the sampling bins to reflect the accumulated statistics and
   * binning used.
   */
  template<class Adaptor>
  void update(const Adaptor& ap) {
    double avgweight = 0.;
    size_t bins = 0;
    for ( map<double,GeneralStatistics>::const_iterator s =
	    statisticsMap.begin(); s != statisticsMap.end(); ++s ) {
      avgweight += ap.importanceMeasure(s->second);
      ++bins;
    }
    avgweight /= bins;
    weightMap.clear();
    double norm = 0.;
    for ( map<double,GeneralStatistics>::const_iterator s =
	    statisticsMap.begin(); s != statisticsMap.end(); ++s ) {
      double weight = ap.importanceMeasure(s->second);
      if ( weight < theWeightThreshold*avgweight )
	weight = theWeightThreshold*avgweight;
      weightMap[s->first] = weight;
      norm += 
	weight *
	(s != statisticsMap.begin() ? (s->first - boost::prior(s)->first) : s->first);
    }
    selectorMap.clear();
    double current = 0.;
    for ( map<double,double>::iterator bw = weightMap.begin(); 
	  bw != weightMap.end(); ++bw ) {
      bw->second /= norm;
      pair<double,double> range = 
	make_pair(bw != weightMap.begin() ? boost::prior(bw)->first : 0.,
		  bw->first);
      current += bw->second*(range.second-range.first);
      selectorMap[current] = range;
    }
  }

  /**
   * Half those bins, which meet the given predicate
   * and update the statistics.
   */
  template<class Adaptor>
  void adapt(const Adaptor& ap) {
    update(ap);
    map<double,GeneralStatistics> newBins;
    for ( map<double,GeneralStatistics>::const_iterator b
	    = statisticsMap.begin(); b != statisticsMap.end(); ++b ) {
      newBins[b->first] = GeneralStatistics();
      if ( ap.adapt(b->second) ) {
	double bound =
	  b != statisticsMap.begin() ? (boost::prior(b)->first + b->first)/2. : b->first/2.;
	newBins[bound] = GeneralStatistics();
      }
    }
    statisticsMap = newBins;
  }

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void put(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void get(PersistentIStream & is);
  //@}

private:

  /**
   * Map upper bin boundaries to statistics contained.
   * The lower bin boundary is always 0.
   */
  map<double,GeneralStatistics> statisticsMap;

  /**
   * Map upper bin boundaries to bin weights currently used.
   */
  map<double,double> weightMap;

  /**
   * Selector map to sample a point
   */
  map<double,pair<double,double> > selectorMap;

  /**
   * The last sampled point.
   */
  double lastPoint;

  /**
   * The statistics relevant to the last sampled point.
   */
  GeneralStatistics* lastStatistics;

  /**
   * The weight threshold which governs the minimum bin weight.
   */
  double theWeightThreshold;

};

inline PersistentOStream& operator<<(PersistentOStream& os, const BinnedStatistics& s) {
  s.put(os); return os;
}

inline PersistentIStream& operator>>(PersistentIStream& is, BinnedStatistics& s) {
  s.get(is); return is;
}

}

#endif /* Herwig_BinnedStatistics_H */
