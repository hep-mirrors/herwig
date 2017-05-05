// -*- C++ -*-
//
// Distribution.hpp is a part of myStatistics
// Copyright (C) 2012-2017 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//
#ifndef MYSTATISTICS_Distribution_hpp_included
#define MYSTATISTICS_Distribution_hpp_included

#include "Histogram.h"

#include <iostream>

namespace Statistics {

  /**
   * \brief A (one dimensional) distribution
   * \author Simon Platzer
   */
  class Distribution {

  public:

    /**
     * A bin in a distribution
     */
    struct DistributionBin {

      /**
       * Default constructor
       */
      DistributionBin();

      /**
       * Construct from a histogram bin, given a number of attempted
       * points
       */
      DistributionBin(const Bin& bin,
		      double nPoints);

      /**
       * Add two distributionBins
       */
      DistributionBin& operator+=(const DistributionBin& other);

      /**
       * Subtract two distributionBins
       */
      DistributionBin& operator-=(const DistributionBin& other);

      /**
       * Multiply two distributionBins
       */
      DistributionBin& operator*=(const DistributionBin& other);

      /**
       * Divide two distributionBins
       */
      DistributionBin& operator/=(const DistributionBin& other);

      /**
       * Fill distribution data from an XML element
       */
      void fromXML(const XML::Element&);

      /**
       * Return an XML element for the data of this distribution
       */
      XML::Element toXML() const;

      /**
       * The bin boundaries
       */
      std::pair<double,double> boundaries;

      /**
       * The bin value
       */
      double value;

      /**
       * The bin error squared
       */
      double errorSquared;

    };

  public:

    /**
     * Default constructor
     */
    Distribution();

    /**
     * Construct from a histogram
     */
    explicit Distribution(const Histogram& histo,
			  double nPoints);

    /**
     * Destructor
     */
    virtual ~Distribution();

  public:

    /**
     * Return the name of the distribution
     */
    const std::string& name() const { return theName; }

    /**
     * Return the bins in the distribution
     */
    const std::vector<DistributionBin>& bins() const { return theBins; }

  public:

    /**
     * Return true, if this distribution is compatible with another one.
     */
    bool isCompatible(const Distribution& other) const;

    /**
     * Add two distributions
     */
    Distribution& operator+=(const Distribution& other);

    /**
     * Subtract two distributions
     */
    Distribution& operator-=(const Distribution& other);

    /**
     * Multiply two distributions
     */
    Distribution& operator*=(const Distribution& other);

    /**
     * Divide two distributions
     */
    Distribution& operator/=(const Distribution& other);

  public:

    /**
     * Fill distribution data from an XML element
     */
    void fromXML(const XML::Element&);

    /**
     * Return an XML element for the data of this distribution
     */
    XML::Element toXML() const;

  public:

    /**
     * Write out data ready for make-plots, given an analysis name and
     * plot options
     */
    void toMakePlots(const std::string& analysis,
		     const std::string& options = "") const;

    /**
     * Write out data ready for make-plots, given an analysis name and
     * plot options
     */
    void toMakePlots(const std::string& analysis,
		     const Distribution& lower,
		     const Distribution& upper,
		     const std::string& options = "") const;

    /**
     * Write out data ready for make-plots, given an analysis name and
     * plot options
     */
    void appendToMakePlots(std::ostream& os,
			   const std::string& analysis,
			   const std::string& options = "") const;

    /**
     * Write out data ready for make-plots, given an analysis name and
     * plot options
     */
    void appendToMakePlots(std::ostream& os,
			   const std::string& analysis,
			   const Distribution& lower,
			   const Distribution& upper,
			   const std::string& options = "") const;

  private:

    /**
     * The name of the distribution
     */
    std::string theName;

    /**
     * The bins in the distribution
     */
    std::vector<DistributionBin> theBins;

  };

}

#endif // MYSTATISTICS_Distribution_hpp_included
