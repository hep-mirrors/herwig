// -*- C++ -*-
//
// Run.cpp is a part of myStatistics
// Copyright (C) 2012-2017 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//

#include "Run.h"

#include <exception>
#include <stdexcept>

using namespace Statistics;
using namespace std;

Run::Run()
  : theName(""),
    theAttemptedPoints(0),
    theSumOfWeights(0.0), theSumOfSquaredWeights(0.0) {}

Run::Run(const string& newName)
  : theName(newName),
    theAttemptedPoints(0),
    theSumOfWeights(0.0), theSumOfSquaredWeights(0.0) {}

Run::~Run() {}

void Run::initialize() {
  for ( map<string,Histogram>::iterator h = theHistograms.begin();
	h != theHistograms.end(); ++h )
    h->second.initialize();
}

void Run::reset() {
  for ( map<string,Histogram>::iterator h = theHistograms.begin();
	h != theHistograms.end(); ++h )
    h->second.reset();
}

void Run::finalize(size_t newAttemptedPoints) {
  for ( map<string,Histogram>::iterator h = theHistograms.begin();
	h != theHistograms.end(); ++h )
    h->second.finalize();
  theAttemptedPoints = newAttemptedPoints;
}

Histogram& Run::addHistogram(const string& newName,
			     const vector<double>& newBoundaries) {
  theHistograms[newName] = Histogram(newName,newBoundaries);
  return histogram(newName);
}

Histogram& Run::addHistogram(const string& newName,
			     const vector<double>& newBoundaries,
			     const pair<double,double>& newPeriodicity) {
  theHistograms[newName] = Histogram(newName,newBoundaries,newPeriodicity);
  return histogram(newName);
}

Histogram& Run::histogram(const string& histoName) {
  map<string,Histogram>::iterator h = theHistograms.find(histoName);
  if ( h == theHistograms.end() )
    throw runtime_error("[Statistics::Run] No such histogram.");
  return h->second;
}

const Histogram& Run::histogram(const string& histoName) const {
  map<string,Histogram>::const_iterator h = theHistograms.find(histoName);
  if ( h == theHistograms.end() )
    throw runtime_error("[Statistics::Run] No such histogram.");
  return h->second;
}

Run& Run::operator+=(const Run& other) {
  theAttemptedPoints += other.attemptedPoints();
  theSumOfWeights += other.sumOfWeights();
  theSumOfSquaredWeights += other.sumOfSquaredWeights();
  for ( map<string,Histogram>::iterator h = theHistograms.begin();
	h != theHistograms.end(); ++h ) {
    if ( other.histograms().find(h->first) == other.histograms().end() )
      continue;
    h->second += other.histogram(h->first);
  }
  for ( map<string,Histogram>::const_iterator h = other.histograms().begin();
	h != other.histograms().end(); ++h ) {
    if ( theHistograms.find(h->first) == theHistograms.end() )
      theHistograms[h->first] = h->second;
  }
  return *this;
}

void Run::fromXML(const XML::Element& elem) {

  elem.getFromAttribute("name",theName);
  elem.getFromAttribute("attemptedPoints",theAttemptedPoints);
  elem.getFromAttribute("sumOfWeights",theSumOfWeights);
  elem.getFromAttribute("sumOfSquaredWeights",theSumOfSquaredWeights);

  list<XML::Element>::const_iterator cit =
    elem.findFirst(XML::ElementTypes::Element,"Histograms");

  if ( cit == elem.children().end() )
    throw runtime_error("[Statistics::Run] Expected a histograms element.");

  for ( list<XML::Element>::const_iterator h = cit->children().begin();
	h != cit->children().end(); ++h ) {
    if ( h->type() == XML::ElementTypes::Element &&
	 h->name() == "Histogram" ) {
      Histogram histo;
      histo.fromXML(*h);
      theHistograms[histo.name()] = histo;
    }
  }

}

XML::Element Run::toXML() const {

  XML::Element elem(XML::ElementTypes::Element,"Run");

  elem.appendAttribute("name",theName);
  elem.appendAttribute("attemptedPoints",theAttemptedPoints);
  elem.appendAttribute("sumOfWeights",theSumOfWeights);
  elem.appendAttribute("sumOfSquaredWeights",theSumOfSquaredWeights);

  XML::Element xhistos(XML::ElementTypes::Element,"Histograms");

  for ( map<string,Histogram>::const_iterator h = theHistograms.begin();
	h != theHistograms.end(); ++h )
    xhistos.append(h->second.toXML());

  elem.append(xhistos);

  return elem;

}


