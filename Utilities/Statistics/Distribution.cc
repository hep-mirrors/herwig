// -*- C++ -*-
//
// Distribution.cpp is a part of myStatistics
// Copyright (C) 2012-2013 Simon Platzer
//
// myStatistics is licenced under version 2 of the GPL, see COPYING for details.
//

#include "Distribution.h"

#include <exception>
#include <stdexcept>
#include <cmath>

#include <fstream>

using namespace XML;
using namespace Statistics;
using namespace std;

template<class T>
inline T sqr(const T& x) {
  return x*x;
}

Distribution::DistributionBin::DistributionBin()
  : boundaries(-DBL_MAX,DBL_MAX),
    value(0.0), errorSquared(0.0) {}

Distribution::DistributionBin::DistributionBin(const Bin& bin,
					       double nPoints)
  : boundaries(bin.boundaries()),
    value(bin.average(nPoints)/(boundaries.second-boundaries.first)),
    errorSquared(bin.varianceOfAverage(nPoints)/sqr(boundaries.second-boundaries.first)) {}

Distribution::DistributionBin& 
Distribution::DistributionBin::operator+=
(const Distribution::DistributionBin& other) {
  value += other.value;
  errorSquared += other.errorSquared;
  return *this;
}

Distribution::DistributionBin& 
Distribution::DistributionBin::operator-=
(const Distribution::DistributionBin& other) {
  value -= other.value;
  errorSquared += other.errorSquared;
  return *this;
}

Distribution::DistributionBin& 
Distribution::DistributionBin::operator*=
(const Distribution::DistributionBin& other) {
  errorSquared = 
    sqr(value*other.value)*
    (errorSquared/sqr(value) + other.errorSquared/sqr(other.value));
  value *= other.value;
  return *this;
}

Distribution::DistributionBin& 
Distribution::DistributionBin::operator/=
(const Distribution::DistributionBin& other) {
  errorSquared = 
    sqr(value/other.value)*
    (errorSquared/sqr(value) + other.errorSquared/sqr(other.value));
  value /= other.value;
  return *this;
}

void Distribution::DistributionBin::fromXML(const XML::Element& elem) {

  elem.getFromAttribute("lowerBound",boundaries.first);
  elem.getFromAttribute("upperBound",boundaries.second);
  elem.getFromAttribute("value",value);
  elem.getFromAttribute("errorSquared",errorSquared);

}

XML::Element Distribution::DistributionBin::toXML() const {

  XML::Element elem(XML::ElementTypes::EmptyElement,"DistributionBin");

  elem.appendAttribute("lowerBound",boundaries.first);
  elem.appendAttribute("upperBound",boundaries.second);
  elem.appendAttribute("value",value);
  elem.appendAttribute("errorSquared",errorSquared);

  return elem;

}

Distribution::Distribution()
  : theName("") {}

Distribution::Distribution(const Histogram& histo,
			   double nPoints)
  : theName(histo.name()) {
  for ( vector<Bin>::const_iterator b = histo.bins().begin();
	b != histo.bins().end(); ++b ) {
    theBins.push_back(DistributionBin(*b,nPoints));
  }
}

Distribution::~Distribution() {}

bool Distribution::isCompatible(const Distribution& other) const {
  if ( bins().size() != other.bins().size() )
    return false;
  vector<DistributionBin>::const_iterator b = bins().begin();
  vector<DistributionBin>::const_iterator ob = other.bins().begin();
  for ( ; b != bins().end(); ++b, ++ob ) {
    if ( b->boundaries != ob->boundaries )
      return false;
  }
  return true;
}

Distribution& Distribution::operator+=(const Distribution& other) {
  if ( !isCompatible(other) )
    throw runtime_error("[Statistics::Histogram] Incompatible distributions.");
  vector<DistributionBin>::iterator b = theBins.begin();
  vector<DistributionBin>::const_iterator ob = other.bins().begin();
  for ( ; b != theBins.end(); ++b, ++ob )
    *b += *ob;
  return *this;
}

Distribution& Distribution::operator-=(const Distribution& other) {
  if ( !isCompatible(other) )
    throw runtime_error("[Statistics::Histogram] Incompatible distributions.");
  vector<DistributionBin>::iterator b = theBins.begin();
  vector<DistributionBin>::const_iterator ob = other.bins().begin();
  for ( ; b != theBins.end(); ++b, ++ob )
    *b -= *ob;
  return *this;
}

Distribution& Distribution::operator*=(const Distribution& other) {
  if ( !isCompatible(other) )
    throw runtime_error("[Statistics::Histogram] Incompatible distributions.");
  vector<DistributionBin>::iterator b = theBins.begin();
  vector<DistributionBin>::const_iterator ob = other.bins().begin();
  for ( ; b != theBins.end(); ++b, ++ob )
    *b *= *ob;
  return *this;
}

Distribution& Distribution::operator/=(const Distribution& other) {
  if ( !isCompatible(other) )
    throw runtime_error("[Statistics::Histogram] Incompatible distributions.");
  vector<DistributionBin>::iterator b = theBins.begin();
  vector<DistributionBin>::const_iterator ob = other.bins().begin();
  for ( ; b != theBins.end(); ++b, ++ob )
    *b /= *ob;
  return *this;
}

void Distribution::fromXML(const XML::Element& elem) {

  elem.getFromAttribute("name",theName);

  for ( list<XML::Element>::const_iterator b = elem.children().begin();
	b != elem.children().end(); ++b ) {
    if ( b->type() == XML::ElementTypes::EmptyElement &&
	 b->name() == "DistributionBin" ) {
      DistributionBin bin;
      bin.fromXML(*b);
      theBins.push_back(bin);
    }
  }

}

XML::Element Distribution::toXML() const {

  XML::Element elem(XML::ElementTypes::Element,"Distribution");

  elem.appendAttribute("name",theName);

  for ( vector<DistributionBin>::const_iterator b = bins().begin();
	b != bins().end(); ++b ) {
    elem.append(b->toXML());
  }

  return elem;

}


void Distribution::appendToMakePlots(ostream& out,
				     const std::string& analysis,
				     const std::string& options) const {

  string aidaPath = "/" + analysis + "/" + name();

  out << "# BEGIN HISTOGRAM " << aidaPath << "\n";

  if ( options != "" )
    out << options;

  for ( vector<DistributionBin>::const_iterator b = bins().begin();
	b != bins().end(); ++b ) {
    out << b->boundaries.first << " "
	<< b->boundaries.second << " "
	<< b->value << " "
	<< sqrt(b->errorSquared) << " "
	<< sqrt(b->errorSquared) << "\n";
  }

  out << "# END HISTOGRAM\n\n";

}

void Distribution::appendToMakePlots(ostream& out,
				     const std::string& analysis,
				     const Distribution& lower,
				     const Distribution& upper,
				     const std::string& options) const {
  if ( !isCompatible(lower) || !isCompatible(upper) )
    throw runtime_error("[Statistics::Histogram] Incompatible distributions.");

  string aidaPath = "/" + analysis + "/" + name();

  out << "# BEGIN HISTOGRAM " << aidaPath << "\n";

  if ( options != "" )
    out << options;
  
  vector<DistributionBin>::const_iterator b = bins().begin();
  vector<DistributionBin>::const_iterator lb = lower.bins().begin();
  vector<DistributionBin>::const_iterator ub = upper.bins().begin();

  for ( ; b != bins().end(); ++b, ++lb, ++ub ) {
    out << b->boundaries.first << " "
	<< b->boundaries.second << " "
	<< b->value << " "
	<< abs(b->value - lb->value) << " "
	<< abs(b->value - ub->value) << "\n";
  }

  out << "# END HISTOGRAM\n\n";

}


void Distribution::toMakePlots(const std::string& analysis,
			       const std::string& options) const {

  ofstream out((name()+".dat").c_str());

  appendToMakePlots(out,analysis,options);

}

void Distribution::toMakePlots(const std::string& analysis,
			       const Distribution& lower,
			       const Distribution& upper,
			       const std::string& options) const {

  ofstream out((name()+".dat").c_str());

  appendToMakePlots(out,analysis,lower,upper,options);

}

