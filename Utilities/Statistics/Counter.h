// -*- C++ -*-
//
// Counter.hpp is a part of myStatistics
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//
#ifndef MYSTATISTICS_Counter_hpp_included
#define MYSTATISTICS_Counter_hpp_included

#include "Herwig/Utilities/XML/Element.h"

namespace Statistics {

  /**
   * \brief A (weighted) counter
   * \author Simon Platzer
   */
  class Counter {

  public:

    /**
     * Construct a new counter.
     */
    Counter();

    /**
     * Destruct a counter.
     */
    virtual ~Counter();

  public:

    /**
     * Initialize this counter.
     */
    void initialize() {}

    /**
     * Reset this counter.
     */
    void reset() { *this = Counter(); }

    /**
     * Finalize this counter.
     */
    void finalize() { if ( isOpen() ) close(); }

  public:

    /**
     * Open the counter for the next event of the given id.
     */
    void open(size_t id);

    /**
     * Close the counter and update the overall statistics.
     */
    void close();

    /**
     * Book a contribution to the current event. Implies close() and
     * open(id) if the id is different from the currently considered
     * event id.
     */
    void count(double weight, size_t id);

  public:

    /**
     * Add a counter to this counter
     */
    Counter& operator+=(const Counter& other);

    /**
     * Subtract a counter from this counter
     */
    Counter& operator-=(const Counter& other);

  public:

    /**
     * Return the sum of weights
     */
    double sumOfWeights() const { return theSumOfWeights; }

    /**
     * Return the sum of squared weights
     */
    double sumOfSquaredWeights() const { return theSumOfSquaredWeights; }

    /**
     * Return the sum of weights before the next event occured
     */
    double sumOfEventWeights() const { return theSumOfEventWeights; }

    /**
     * Return the current event id
     */
    size_t eventId() const { return theEventId; }

    /**
     * Return true, if this counter is open
     */
    bool isOpen() const { return eventId() != 0; }

    /**
     * Return true, if this counter is closed
     */
    bool isClosed() const { return !isOpen(); }

  public:

    /**
     * Given a number of sampled points, return the average of this
     * counter
     */
    double average(double nPoints) const;

    /**
     * Given a number of sampled points, return the variance of the
     * average of this counter
     */
    double varianceOfAverage(double nPoints) const;

  public:

    /**
     * Fill counter data from an XML element
     */
    void fromXML(const XML::Element&);

    /**
     * Return an XML element for the data of this counter
     */
    XML::Element toXML() const;

  private:

    /**
     * The sum of weights counted
     */
    double theSumOfWeights;

    /**
     * The sum of squared weights counted
     */
    double theSumOfSquaredWeights;

    /**
     * The sum of weights counted before the next event occured;
     * within one event, all weights are considered to be 100%
     * correlated
     */
    double theSumOfEventWeights;

    /**
     * The current event id. The counter is closed if this is zero.
     */
    size_t theEventId;

  };

}

#endif // MYSTATISTICS_Counter_hpp_included
