// -*- C++ -*-
//
// EventContribution.hpp is a part of myStatistics
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//
#ifndef MYSTATISTICS_EventContribution_hpp_included
#define MYSTATISTICS_EventContribution_hpp_included

#include <cfloat>

namespace Statistics {
  
  /**
   * \brief A pointlike or boxlike eventContribution; serves to define the EventContribution concept
   * \author Simon Platzer
   */
  class EventContribution {
    
  public:

    /**
     * Construct an eventContribution with given width
     */
    EventContribution(double newCentralValue,
		      double newWeight,
		      double newWidth = 0.0)
      : theCentralValue(newCentralValue),
	theSupport(newCentralValue - newWidth/2.,
		   newCentralValue + newWidth/2.),
	theWeight(newWeight) {}

  public:

    /**
     * Return the central value
     */
    double centralValue() const { return theCentralValue; }

    /**
     * Return the support
     */
    const std::pair<double,double>& support() const { return theSupport; }

    /**
     * Return the normalized overlap with an interval
     */
    double overlap(const std::pair<double,double>& interval) const {
      return calculateOverlap(interval,support(),support().second-support().first);
    }

    /**
     * Return the eventContribution weight
     */
    double weight() const { return theWeight; }

  public:

    /**
     * Remap the central value and support given a periodicity
     * interval; if the support exceeds the periodicity it is ajusted
     * to the peridocity interval.
     */
    void periodic(const std::pair<double,double>& periodicity) {
      double delta = periodicity.second - periodicity.first;
      if ( support().second - support().first > delta ) {
	theSupport.first = centralValue() - delta/2.;
	theSupport.second = centralValue() + delta/2.;
      }
      double shift = 0.;
      if ( centralValue() >= periodicity.second ) {
	while ( centralValue() >= periodicity.second ) {
	  shift -= delta;
	  theCentralValue -= delta;
	}
      } else if ( centralValue() < periodicity.first ) {
	while ( centralValue() <= periodicity.first ) {
	  shift += delta;
	  theCentralValue += delta;
	}
      }
      theSupport.first += shift;
      theSupport.second += shift;
      if ( theSupport.first < periodicity.first )
	theSupport.first += delta;
      if ( theSupport.second > periodicity.second )
	theSupport.second -= delta;
    }

    /**
     * Adjust to lower boundary
     */
    void noUnderflow(double lower) {
      if ( support().first < lower ) {
	double shift = lower - support().first;
	theSupport.first += shift;
	theSupport.second += shift;
	theCentralValue += shift;
      }
    }

    /**
     * Adjust to upper boundary
     */
    void noOverflow(double upper) {
      if ( support().second > upper ) {
	double shift = upper - support().second - DBL_EPSILON;
	theSupport.first += shift;
	theSupport.second += shift;
	theCentralValue += shift;
      }
    }

    /**
     * Return the normalized overlap with an interval, assuming a
     * periodic quantity
     */
    double overlap(const std::pair<double,double>& interval,
		   const std::pair<double,double>& periodicity) const {
      if ( support().first < support().second )
	return calculateOverlap(interval,support(),support().second-support().first);
      double norm = 
	support().second - periodicity.first +
	periodicity.second - support().first;
      return
	calculateOverlap(interval,std::make_pair(periodicity.first,support().second),norm) +
	calculateOverlap(interval,std::make_pair(support().first,periodicity.second),norm);
    }

  private:

    /**
     * Calculate the normalized overlap with an interval
     */
    double calculateOverlap(const std::pair<double,double>& interval,
			    const std::pair<double,double>& newSupport,
			    double norm) const {
      if ( newSupport.first == newSupport.second ) {
	if ( newSupport.first >= interval.first &&
	     newSupport.second < interval.second )
	  return 1.0;
	return 0.0;
      }
      if ( newSupport.first <= interval.first &&
	   newSupport.second <= interval.first )
	return 0.0;
      if ( newSupport.first >= interval.second &&
	   newSupport.second >= interval.second )
	return 0.0;
      double lower = std::max(newSupport.first,interval.first);
      double upper = std::min(newSupport.second,interval.second);
      return std::max(0.0,(upper-lower)/norm);
    }

    /**
     * The central value
     */
    double theCentralValue;

    /**
     * The support
     */
    std::pair<double,double> theSupport;

    /**
     * The eventContribution weight
     */
    double theWeight;

  };

}

#endif // MYSTATISTICS_EventContribution_hpp_included
