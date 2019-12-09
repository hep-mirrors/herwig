// -*- C++ -*-
//
// Bin.cpp is a part of myStatistics
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//

#include "Bin.h"

#include <exception>
#include <stdexcept>
#include <algorithm>

using namespace Statistics;
using namespace std;

Bin::Bin(const std::pair<double,double>& newBoundaries)
  : Counter(), theBoundaries(newBoundaries) {}

Bin::~Bin() {}

Bin& Bin::operator+=(const Bin& other) {
  if ( boundaries() != other.boundaries() )
    throw runtime_error("[Statistics::Bin] Incompatible bin sizes.");
  Counter::operator+=(other);
  return *this;
}

Bin& Bin::operator-=(const Bin& other) {
  if ( boundaries() != other.boundaries() )
    throw runtime_error("[Statistics::Bin] Incompatible bin sizes.");
  Counter::operator-=(*this);
  return *this;
}

void Bin::fromXML(const XML::Element& elem) {

  elem.getFromAttribute("lowerBound",theBoundaries.first);
  elem.getFromAttribute("upperBound",theBoundaries.second);

  list<XML::Element>::const_iterator cit =
    elem.findFirst(XML::ElementTypes::EmptyElement,"Counter");

  if ( cit == elem.children().end() )
    throw runtime_error("[Statistics::Bin] Expected a counter element.");
  
  Counter::fromXML(*cit);

}

XML::Element Bin::toXML() const {

  XML::Element elem(XML::ElementTypes::Element,"Bin");

  elem.appendAttribute("lowerBound",theBoundaries.first);
  elem.appendAttribute("upperBound",theBoundaries.second);

  elem.append(Counter::toXML());
  
  return elem;

}

