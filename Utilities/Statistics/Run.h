// -*- C++ -*-
//
// Run.hpp is a part of myStatistics
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//
#ifndef MYSTATISTICS_Run_hpp_included
#define MYSTATISTICS_Run_hpp_included

#include "Histogram.h"

namespace Statistics {

  /**
   * \brief A simulation run
   * \author Simon Platzer
   */
  class Run {

  public:

    /**
     * Default constructor
     */
    Run();

    /**
     * Construct giving a name and seed
     */
    explicit Run(const std::string& newName);

    /**
     * Destructor
     */
    virtual ~Run();

  public:

    /**
     * Initialize this run.
     */
    void initialize();

    /**
     * Reset this run.
     */
    void reset();

    /**
     * Finalize this run.
     */
    void finalize(size_t newAttemptedPoints);

  public:

    /**
     * Add a point to this run
     */
    void count(double weight) {
      theSumOfWeights += weight;
      theSumOfSquaredWeights += weight*weight;
    }

    /**
     * Add a histogram
     */
    Histogram& addHistogram(const std::string& newName,
			    const std::vector<double>& newBoundaries);

    /**
     * Add a histogram
     */
    Histogram& addHistogram(const std::string& newName,
			    const std::vector<double>& newBoundaries,
			    const std::pair<double,double>& newPeriodicity);

    /**
     * Return a given histogram
     */
    Histogram& histogram(const std::string& histoName);

    /**
     * Return a given histogram
     */
    const Histogram& histogram(const std::string& histoName) const;

    /**
     * Return the histograms
     */
    const std::map<std::string,Histogram>& histograms() const { return theHistograms; }

    /**
     * Return the name of the run
     */
    const std::string& name() const { return theName; }

    /**
     * Set the name of the run
     */
    void name(const std::string& newName) { theName = newName; }

    /**
     * Return the total number of attempted points
     */
    size_t attemptedPoints() const { return theAttemptedPoints; }

    /**
     * The sum of weights
     */
    double sumOfWeights() const { return theSumOfWeights; }

    /**
     * The sum of squared weights
     */
    double sumOfSquaredWeights() const { return theSumOfSquaredWeights; }

  public:

    /**
     * Add a run to this run
     */
    Run& operator+=(const Run& other);

  public:

    /**
     * Fill run data from an XML element
     */
    void fromXML(const XML::Element&);

    /**
     * Return an XML element for the data of this run
     */
    XML::Element toXML() const;

  private:

    /**
     * The name of the run
     */
    std::string theName;

    /**
     * The total number of attempted points
     */
    size_t theAttemptedPoints;

    /**
     * The sum of weights
     */
    double theSumOfWeights;

    /**
     * The sum of squared weights
     */
    double theSumOfSquaredWeights;

    /**
     * The histograms
     */
    std::map<std::string,Histogram> theHistograms;


  };

}

#endif // MYSTATISTICS_Run_hpp_included
