// -*- C++ -*-
//
// Counter.cpp is a part of myStatistics
// Copyright (C) 2012-2013 Simon Platzer
//
// myStatistics is licenced under version 2 of the GPL, see COPYING for details.
//

#include "Counter.h"

#include <exception>
#include <stdexcept>
#include <cmath>

using namespace Statistics;
using namespace std;

template<class T>
inline T sqr(const T& x) {
  return x*x;
}

Counter::Counter()
  : theSumOfWeights(0.0), theSumOfSquaredWeights(0.0),
    theSumOfEventWeights(0.0), theEventId(0) {}

Counter::~Counter() {}

Counter& Counter::operator+=(const Counter& other) {
  theSumOfWeights += other.sumOfWeights();
  theSumOfSquaredWeights += other.sumOfSquaredWeights();
  return *this;
}

Counter& Counter::operator-=(const Counter& other) {
  theSumOfWeights -= other.sumOfWeights();
  theSumOfSquaredWeights += other.sumOfSquaredWeights();
  return *this;
}

void Counter::open(size_t id) {
  if ( !isClosed() )
    throw runtime_error("[Statistics::Counter] attempt to open an unclosed counter.");
  if ( id == 0 )
    throw runtime_error("[Statistics::Counter] attempt to open a counter with invalid event id.");
  theSumOfEventWeights = 0.0;
  theEventId = id;
}

void Counter::close() {
  if ( !isOpen() )
    throw runtime_error("[Statistics::Counter] attempt to close an unopened counter.");
  theSumOfWeights += theSumOfEventWeights;
  theSumOfSquaredWeights += sqr(theSumOfEventWeights);
  theSumOfEventWeights = 0.;
  theEventId = 0;
}

void Counter::count(double weight, size_t id) {
  if ( isClosed() )
    open(id);
  else if ( isOpen() && id != eventId() ) {
    close();
    open(id);
  }
  theSumOfEventWeights += weight;
}

double Counter::average(double nPoints) const {
  if ( nPoints <= 1 )
    throw runtime_error("[Statistics::Counter] cannot perform a reliable estimate for less than two points.");
  return
    sumOfWeights()/nPoints;
}

double Counter::varianceOfAverage(double nPoints) const {
  if ( nPoints <= 1 )
    throw runtime_error("[Statistics::Counter] cannot perform a reliable estimate for less than two points.");
  return
    abs(sumOfSquaredWeights()/nPoints - sqr(sumOfWeights()/nPoints))/(nPoints-1);
}

void Counter::fromXML(const XML::Element& elem) {

  elem.getFromAttribute("sumOfWeights",theSumOfWeights);
  elem.getFromAttribute("sumOfSquaredWeights",theSumOfSquaredWeights);
  elem.getFromAttribute("sumOfEventWeights",theSumOfEventWeights);
  elem.getFromAttribute("eventId",theEventId);  

}

XML::Element Counter::toXML() const {

  XML::Element elem(XML::ElementTypes::EmptyElement,"Counter");

  elem.appendAttribute("sumOfWeights",theSumOfWeights);
  elem.appendAttribute("sumOfSquaredWeights",theSumOfSquaredWeights);
  elem.appendAttribute("sumOfEventWeights",theSumOfEventWeights);
  elem.appendAttribute("eventId",theEventId);

  return elem;

}

