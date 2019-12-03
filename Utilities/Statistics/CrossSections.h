// -*- C++ -*-
//
// CrossSections.hpp is a part of myStatistics
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//
#ifndef MYSTATISTICS_CrossSections_hpp_included
#define MYSTATISTICS_CrossSections_hpp_included

#include "Distribution.h"

namespace Statistics {

  /**
   * \brief A simulation run
   * \author Simon Platzer
   */
  class CrossSections {

  public:

    /**
     * Default constructor
     */
    CrossSections();

    /**
     * Construct giving a name and seed
     */
    explicit CrossSections(const std::string& newName);

    /**
     * Destructor
     */
    virtual ~CrossSections();

  public:

    /**
     * Add a distribution
     */
    Distribution& addDistribution(const Distribution&);

    /**
     * Return a given distribution
     */
    Distribution& distribution(const std::string& histoName);

    /**
     * Return a given distribution
     */
    const Distribution& distribution(const std::string& histoName) const;

    /**
     * Return the distributions
     */
    const std::map<std::string,Distribution>& distributions() const { return theDistributions; }

    /**
     * Return the name of the run
     */
    const std::string& name() const { return theName; }

    /**
     * Set the name of the run
     */
    void name(const std::string& newName) { theName = newName; }

    /**
     * The integral
     */
    double integral() const { return theIntegral; }

    /**
     * The variance of the integral
     */
    double varianceOfIntegral() const { return theVarianceOfIntegral; }

  public:

    /**
     * Add 
     */
    CrossSections& operator+=(const CrossSections& other);

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
     * The integral
     */
    double theIntegral;

    /**
     * The variance of the integral
     */
    double theVarianceOfIntegral;

    /**
     * The distributions
     */
    std::map<std::string,Distribution> theDistributions;


  };

}

#endif // MYSTATISTICS_CrossSections_hpp_included
