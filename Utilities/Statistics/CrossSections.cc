// -*- C++ -*-
//
// CrossSections.cpp is a part of myStatistics
// Copyright (C) 2012-2019 Simon Platzer, The Herwig Collaboration
//
// myStatistics is licenced under version 3 of the GPL, see COPYING for details.
//

#include "CrossSections.h"

#include <exception>
#include <stdexcept>

using namespace Statistics;
using namespace std;

CrossSections::CrossSections()
  : theName(""),
    theIntegral(0.0), theVarianceOfIntegral(0.0) {}

CrossSections::CrossSections(const string& newName)
  : theName(newName),
    theIntegral(0.0), theVarianceOfIntegral(0.0) {}

CrossSections::~CrossSections() {}


Distribution& CrossSections::addDistribution(const Distribution& dist) {
  theDistributions[dist.name()] = dist;
  return distribution(dist.name());
}

Distribution& CrossSections::distribution(const string& histoName) {
  map<string,Distribution>::iterator h = theDistributions.find(histoName);
  if ( h == theDistributions.end() )
    throw runtime_error("[Statistics::CrossSections] No such distribution.");
  return h->second;
}

const Distribution& CrossSections::distribution(const string& histoName) const {
  map<string,Distribution>::const_iterator h = theDistributions.find(histoName);
  if ( h == theDistributions.end() )
    throw runtime_error("[Statistics::CrossSections] No such distribution.");
  return h->second;
}

CrossSections& CrossSections::operator+=(const CrossSections& other) {
  theIntegral += other.integral();
  theVarianceOfIntegral += other.varianceOfIntegral();
  for ( map<string,Distribution>::iterator h = theDistributions.begin();
	h != theDistributions.end(); ++h ) {
    if ( other.distributions().find(h->first) == other.distributions().end() )
      continue;
    h->second += other.distribution(h->first);
  }
  for ( map<string,Distribution>::const_iterator h = other.distributions().begin();
	h != other.distributions().end(); ++h ) {
    if ( theDistributions.find(h->first) == theDistributions.end() )
      theDistributions[h->first] = h->second;
  }
  return *this;
}

void CrossSections::fromXML(const XML::Element& elem) {

  elem.getFromAttribute("name",theName);
  elem.getFromAttribute("integral",theIntegral);
  elem.getFromAttribute("varianceOfIntegral",theVarianceOfIntegral);

  list<XML::Element>::const_iterator cit =
    elem.findFirst(XML::ElementTypes::Element,"Distributions");

  if ( cit == elem.children().end() )
    //throw runtime_error("[Statistics::CrossSections] Expected a distributions element.");
    return;

  for ( list<XML::Element>::const_iterator h = cit->children().begin();
	h != cit->children().end(); ++h ) {
    if ( h->type() == XML::ElementTypes::Element &&
	 h->name() == "Distribution" ) {
      Distribution histo;
      histo.fromXML(*h);
      theDistributions[histo.name()] = histo;
    }
  }

}

XML::Element CrossSections::toXML() const {

  XML::Element elem(XML::ElementTypes::Element,"CrossSections");

  elem.appendAttribute("name",theName);
  elem.appendAttribute("integral",theIntegral);
  elem.appendAttribute("varianceOfIntegral",theVarianceOfIntegral);

  XML::Element xhistos(XML::ElementTypes::Element,"Distributions");

  for ( map<string,Distribution>::const_iterator h = theDistributions.begin();
	h != theDistributions.end(); ++h )
    xhistos.append(h->second.toXML());

  elem.append(xhistos);

  return elem;

}

