// -*- C++ -*-
//
// Bin.hpp is a part of myStatistics
// Copyright (C) 2012-2017 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//
#ifndef MYSTATISTICS_Bin_hpp_included
#define MYSTATISTICS_Bin_hpp_included

#include <utility>
#include <cfloat>

#include "Counter.h"

namespace Statistics {

  /**
   * \brief A bin in a (one dimensional) histogram
   * \author Simon Platzer
   */
  class Bin 
    : public Counter {

  public:

    /**
     * Construct a bin given boundaries
     */
    explicit Bin(const std::pair<double,double>& newBoundaries =
		 std::make_pair(-DBL_MAX,DBL_MAX));

    /**
     * Destruct a bin
     */
    virtual ~Bin();

    /**
     * Return the bin boundaries
     */
    const std::pair<double,double>& boundaries() const { return theBoundaries; }

  public:

    /**
     * Add a bin to this bin
     */
    Bin& operator+=(const Bin& other);

    /**
     * Subtract a bin from this bin
     */
    Bin& operator-=(const Bin& other);

  public:

    /**
     * Fill bin data from an XML element
     */
    void fromXML(const XML::Element&);

    /**
     * Return an XML element for the data of this bin
     */
    XML::Element toXML() const;

  private:

    /**
     * The bin boundaries
     */
    std::pair<double,double> theBoundaries;

  };

}

#endif // MYSTATISTICS_Bin_hpp_included
