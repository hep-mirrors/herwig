// -*- C++ -*-
//
// Histogram.hpp is a part of myStatistics
// Copyright (C) 2012-2013 Simon Platzer
//
// myStatistics is licenced under version 2 of the GPL, see COPYING for details.
//
#ifndef MYSTATISTICS_Histogram_hpp_included
#define MYSTATISTICS_Histogram_hpp_included

#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "EventContribution.h"
#include "Bin.h"

namespace Statistics {

  /**
   * \brief A (one dimensional) histogram
   * \author Simon Platzer
   */
  class Histogram {

  public:

    /**
     * Default constructor
     */
    Histogram();

    /**
     * Destructor
     */
    virtual ~Histogram();

    /**
     * Constructor giving a name and bin boundaries; bin boundaries
     * are subsequent boundaries.
     */
    Histogram(const std::string& newName,
	      const std::vector<double>& newBoundaries,
	      bool newNoUnderflow = false,
	      bool newNoOverflow = false);

    /**
     * Constructor giving a name, bin boundaries and periodicity
     * interval; bin boundaries are subsequent boundaries.
     */
    Histogram(const std::string& newName,
	      const std::vector<double>& newBoundaries,
	      const std::pair<double,double>& newPeriodicity);

    /**
     * Create equally spaced bin edges
     */
    static std::vector<double> regularBinEdges(double lower, double upper, size_t nBins);

    /**
     * Create logarithmicaly spaced bin edges
     */
    static std::vector<double> logBinEdges(double lower, double upper, size_t nBins);

  public:

    /**
     * Initialize this histogram.
     */
    void initialize();

    /**
     * Reset this histogram.
     */
    void reset();

    /**
     * Finalize this histogram.
     */
    void finalize();

  public:

    /**
     * Book a contribution to the current event.
     */
    bool count(EventContribution event, size_t id);

  public:

    /**
     * Return true, if this histogram is compatible with another one.
     */
    bool isCompatible(const Histogram& other) const;

    /**
     * Add a histogram to this histogram
     */
    Histogram& operator+=(const Histogram& other);

    /**
     * Subtract a histogram from this histogram
     */
    Histogram& operator-=(const Histogram& other);

  public:

    /**
     * Return the id of the histogram
     */
    const std::string& name() const { return theName; }

    /**
     * Return the underflow bin
     */
    const Bin& underflow() const { return theUnderflow; }

    /**
     * Return the bins
     */
    const std::vector<Bin>& bins() const { return theBins; }

    /**
     * Return the overflow bin
     */
    const Bin& overflow() const { return theOverflow; }

    /**
     * True, if there is no underflow
     */
    bool noUnderflow() const { return theNoUnderflow; }

    /**
     * True, if there is no overflow
     */
    bool noOverflow() const { return theNoOverflow; }


    /**
     * Return true, if the quantity considered is periodic
     */
    bool isPeriodic() const { return theIsPeriodic; }

    /**
     * Return the periodicity interval, if appropriate
     */
    const std::pair<double,double>& periodicity() const { return thePeriodicity; }

  public:

    /**
     * Fill histogram data from an XML element
     */
    void fromXML(const XML::Element&);

    /**
     * Return an XML element for the data of this histogram
     */
    XML::Element toXML() const;

  private:

    /**
     * The id of the histogram
     */
    std::string theName;

    /**
     * The underflow bin
     */
    Bin theUnderflow;

    /**
     * The bins
     */
    std::vector<Bin> theBins;

    /**
     * The overflow bin
     */
    Bin theOverflow;

    /**
     * True, if there is no underflow
     */
    bool theNoUnderflow;

    /**
     * True, if there is no overflow
     */
    bool theNoOverflow;

    /**
     * True, if the quantity considered is periodic
     */
    bool theIsPeriodic;

    /**
     * The periodicity of the quantity considered
     */
    std::pair<double,double> thePeriodicity;

    /**
     * Fill the bin map
     */
    void fillBinMap();

    /**
     * Return a bin by index
     */
    Bin& binByIndex(int);

    /**
     * Map bin upper boundaries to bins
     */
    std::map<double,int> binMap;

  };

}

#endif // MYSTATISTICS_Histogram_hpp_included
