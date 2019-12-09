// -*- C++ -*-
//
// Histogram.cpp is a part of myStatistics
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//

#include "Histogram.h"

#include <exception>
#include <stdexcept>
#include <cmath>

using namespace Statistics;
using namespace std;

Histogram::Histogram()
  : theNoUnderflow(false),
    theNoOverflow(false),
    theIsPeriodic(false),
    thePeriodicity(-DBL_MAX,DBL_MAX) {}

Histogram::~Histogram() {}

Histogram::Histogram(const string& newName,
		     const vector<double>& newBoundaries,
		     bool newNoUnderflow,
		     bool newNoOverflow)
  : theName(newName),
    theUnderflow(make_pair(-DBL_MAX,newBoundaries.front())),
    theOverflow(make_pair(newBoundaries.back(),DBL_MAX)),
    theNoUnderflow(newNoUnderflow),
    theNoOverflow(newNoOverflow),
    theIsPeriodic(false),
    thePeriodicity(-DBL_MAX,DBL_MAX) {
  for ( vector<double>::const_iterator bb = newBoundaries.begin();
	bb != newBoundaries.end() - 1; ++bb ) {
    theBins.push_back(Bin(make_pair(*bb,*(bb+1))));
  }
  fillBinMap();
}

Histogram::Histogram(const string& newName,
		     const vector<double>& newBoundaries,
		     const pair<double,double>& newPeriodicity)
  : theName(newName),
    theNoUnderflow(true),
    theNoOverflow(true),
    theIsPeriodic(true),
    thePeriodicity(newPeriodicity) {
  for ( vector<double>::const_iterator bb = newBoundaries.begin();
	bb != newBoundaries.end() - 1; ++bb ) {
    theBins.push_back(Bin(make_pair(*bb,*(bb+1))));
  }
  fillBinMap();
}

vector<double> Histogram::regularBinEdges(double lower, double upper, size_t nBins) {
  vector<double> res;
  double step = (upper - lower)/nBins;
  for (size_t k = 0; k <= nBins; ++k )
    res.push_back(lower + k*step);
  return res;
}

vector<double> Histogram::logBinEdges(double lower, double upper, size_t nBins) {
  vector<double> res;
  double step = log10(upper/lower)/nBins;
  for (size_t k = 0; k <= nBins; ++k )
    res.push_back(lower*pow(10.0,k*step));
  return res;
}

void Histogram::initialize() {
  fillBinMap();
  for ( vector<Bin>::iterator b = theBins.begin();
	b != theBins.end(); ++b )
    b->initialize();
}

void Histogram::reset() {
  for ( vector<Bin>::iterator b = theBins.begin();
	b != theBins.end(); ++b )
    b->reset();
}

void Histogram::finalize() {
  for ( vector<Bin>::iterator b = theBins.begin();
	b != theBins.end(); ++b )
    b->finalize();
  if ( !noUnderflow() )
    theUnderflow.finalize();
  if ( !noOverflow() )
    theOverflow.finalize();
}

void Histogram::fillBinMap() {
  if ( !binMap.empty() )
    return;
  for ( size_t k = 0; k < theBins.size(); ++k )
    binMap[theBins[k].boundaries().second] = k;
  if ( !noUnderflow() )
    binMap[underflow().boundaries().second] = -1;
  if ( !noOverflow() )
    binMap[overflow().boundaries().second] = -2;
}

Bin& Histogram::binByIndex(int idx) {
  if ( idx >= 0 )
    return theBins[idx];
  if ( idx == -1 )
    return theUnderflow;
  if ( idx == -2 )
    return theOverflow;
  static Bin empty;
  return empty;
}

bool Histogram::count(EventContribution event, size_t id) {

  if ( isPeriodic() ) {
    event.periodic(periodicity());
  } else {
    if ( noUnderflow() )
      event.noUnderflow(bins().front().boundaries().first);
    if ( noOverflow() )
      event.noOverflow(bins().back().boundaries().second);
  }

  map<double,int>::iterator lbin = binMap.upper_bound(event.support().first);
  map<double,int>::iterator ubin = binMap.upper_bound(event.support().second);

  if ( lbin == binMap.end() ||
       ubin == binMap.end() )
    return false;

  ++ubin;

  bool counted = false;

  if ( event.support().first <= event.support().second ) {
    for ( map<double,int>::iterator b = lbin; b != ubin; ++b ) {
      Bin& bin = binByIndex(b->second);
      double overlap = event.overlap(bin.boundaries());
      bin.count(overlap*event.weight(),id);
      counted = true;
    }
  } else {
    if ( !isPeriodic() )
      return false;
    for ( map<double,int>::iterator b = binMap.begin(); b != ubin; ++b ) {
      Bin& bin = binByIndex(b->second);
      double overlap = event.overlap(bin.boundaries(),periodicity());
      bin.count(overlap*event.weight(),id);
      counted = true;
    }
    for ( map<double,int>::iterator b = lbin; b != binMap.end(); ++b ) {
      Bin& bin = binByIndex(b->second);
      double overlap = event.overlap(bin.boundaries(),periodicity());
      bin.count(overlap*event.weight(),id);
      counted = true;
    }
  }

  return counted;

}

bool Histogram::isCompatible(const Histogram& other) const {
  if ( bins().size() != other.bins().size() )
    return false;
  vector<Bin>::const_iterator b = bins().begin();
  vector<Bin>::const_iterator ob = other.bins().begin();
  for ( ; b != bins().end(); ++b, ++ob ) {
    if ( b->boundaries() != ob->boundaries() )
      return false;
  }
  if ( isPeriodic() && !other.isPeriodic() )
    return false;
  if ( !isPeriodic() && other.isPeriodic() )
    return false;
  if ( isPeriodic() &&
       periodicity() != other.periodicity() )
    return false;
  if ( noUnderflow() && !other.noUnderflow() )
    return false;
  if ( !noUnderflow() && other.noUnderflow() )
    return false;
  if ( noOverflow() && !other.noOverflow() )
    return false;
  if ( !noOverflow() && other.noOverflow() )
    return false;
  return true;
}

Histogram& Histogram::operator+=(const Histogram& other) {
  if ( !isCompatible(other) )
    throw runtime_error("[Statistics::Histogram] Incompatible histograms.");
  theUnderflow += other.underflow();
  theOverflow += other.overflow();
  vector<Bin>::iterator b = theBins.begin();
  vector<Bin>::const_iterator ob = other.bins().begin();
  for ( ; b != theBins.end(); ++b, ++ob )
    *b += *ob;
  return *this;
}

Histogram& Histogram::operator-=(const Histogram& other) {
  if ( !isCompatible(other) )
    throw runtime_error("[Statistics::Histogram] Incompatible histograms.");
  theUnderflow -= other.underflow();
  theOverflow -= other.overflow();
  vector<Bin>::iterator b = theBins.begin();
  vector<Bin>::const_iterator ob = other.bins().begin();
  for ( ; b != theBins.end(); ++b, ++ob )
    *b -= *ob;
  return *this;
}

void Histogram::fromXML(const XML::Element& elem) {

  elem.getFromAttribute("name",theName);
  elem.getFromAttribute("isPeriodic",theIsPeriodic);

  list<XML::Element>::const_iterator cit;

  if ( isPeriodic() ) {
    elem.getFromAttribute("periodicityLowerBound",thePeriodicity.first);
    elem.getFromAttribute("periodicityUpperBound",thePeriodicity.second);
    theNoUnderflow = true;
    theNoOverflow = true;
  } 

  cit = elem.findFirst(XML::ElementTypes::Element,"Underflow");
  if ( cit == elem.children().end() ) {
    theNoUnderflow = true;
  } else {
    list<XML::Element>::const_iterator b = cit->children().begin();
    for ( ; b != cit->children().end(); ++b ) {
      if ( b->type() == XML::ElementTypes::Element &&
	   b->name() == "Bin" ) {
	theUnderflow.fromXML(*b);
	break;
      }
    }
    if ( b == cit->children().end() )
      throw runtime_error("[Statistics::Histogram] Expected bin information in the underflow element.");
    theNoUnderflow = false;
  }

  cit = elem.findFirst(XML::ElementTypes::Element,"Overflow");
  if ( cit == elem.children().end() ) {
    theNoOverflow = true;
  } else {
    list<XML::Element>::const_iterator b = cit->children().begin();
    for ( ; b != cit->children().end(); ++b ) {
      if ( b->type() == XML::ElementTypes::Element &&
	   b->name() == "Bin" ) {
	theOverflow.fromXML(*b);
	break;
      }
    }
    if ( b == cit->children().end() )
      throw runtime_error("[Statistics::Histogram] Expected bin information in the overflow element.");
    theNoOverflow = false;
  }


  cit = elem.findFirst(XML::ElementTypes::Element,"Bins");
  if ( cit == elem.children().end() )
    throw runtime_error("[Statistics::Histogram] Expected a Bins element.");

  for ( list<XML::Element>::const_iterator b = cit->children().begin();
	b != cit->children().end(); ++b ) {
    if ( b->type() == XML::ElementTypes::Element &&
	 b->name() == "Bin" ) {
      Bin bin;
      bin.fromXML(*b);
      theBins.push_back(bin);
    }
  }

  fillBinMap();

}

XML::Element Histogram::toXML() const {

  XML::Element elem(XML::ElementTypes::Element,"Histogram");

  elem.appendAttribute("name",theName);
  elem.appendAttribute("isPeriodic",theIsPeriodic);

  if ( isPeriodic() ) {
    elem.appendAttribute("periodicityLowerBound",thePeriodicity.first);
    elem.appendAttribute("periodicityUpperBound",thePeriodicity.second);
  } 

  if ( !noUnderflow() ) {
    XML::Element uf(XML::ElementTypes::Element,"Underflow");
    uf.append(theUnderflow.toXML());
    elem.append(uf);
  }

  if ( !noOverflow() ) {
    XML::Element of(XML::ElementTypes::Element,"Overflow");
    of.append(theOverflow.toXML());
    elem.append(of);
  }

  XML::Element xbins(XML::ElementTypes::Element,"Bins");

  for ( vector<Bin>::const_iterator b = bins().begin();
	b != bins().end(); ++b ) {
    xbins.append(b->toXML());
  }

  elem.append(xbins);

  return elem;

}

