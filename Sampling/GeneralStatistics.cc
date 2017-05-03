// -*- C++ -*-
//
// GeneralStatictis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralStatistics class.
//

#include "GeneralStatistics.h"

using namespace Herwig;

GeneralStatistics::~GeneralStatistics() {}

void GeneralStatistics::put(PersistentOStream & os) const {
  os << theMaxWeight << theMinWeight << theSumWeights
     << theSumSquaredWeights << theSumAbsWeights 
     << theSelectedPoints << theAcceptedPoints
     << theNanPoints << theAllPoints << theLastWeight;
}

void GeneralStatistics::get(PersistentIStream & is) {
  is >> theMaxWeight >> theMinWeight >> theSumWeights
     >> theSumSquaredWeights >> theSumAbsWeights 
     >> theSelectedPoints >> theAcceptedPoints
     >> theNanPoints >> theAllPoints >> theLastWeight;
}

void GeneralStatistics::fromXML(const XML::Element& elem) {

  elem.getFromAttribute("maxWeight", theMaxWeight);
  elem.getFromAttribute("minWeight", theMinWeight);
  elem.getFromAttribute("sumWeights", theSumWeights);
  elem.getFromAttribute("sumSquaredWeights", theSumSquaredWeights);
  elem.getFromAttribute("sumAbsWeights", theSumAbsWeights);
  elem.getFromAttribute("selectedPoints", theSelectedPoints);
  elem.getFromAttribute("acceptedPoints", theAcceptedPoints);
  elem.getFromAttribute("nanPoints", theNanPoints);
  elem.getFromAttribute("allPoints", theAllPoints);
  elem.getFromAttribute("lastWeight", theLastWeight);

}

XML::Element GeneralStatistics::toXML() const {

  XML::Element elem(XML::ElementTypes::Element,"GeneralStatistics");

  elem.appendAttribute("maxWeight", theMaxWeight);
  elem.appendAttribute("minWeight", theMinWeight);
  elem.appendAttribute("sumWeights", theSumWeights);
  elem.appendAttribute("sumSquaredWeights", theSumSquaredWeights);
  elem.appendAttribute("sumAbsWeights", theSumAbsWeights);
  elem.appendAttribute("selectedPoints", theSelectedPoints);
  elem.appendAttribute("acceptedPoints", theAcceptedPoints);
  elem.appendAttribute("nanPoints", theNanPoints);
  elem.appendAttribute("allPoints", theAllPoints);
  elem.appendAttribute("lastWeight", theLastWeight);

  return elem;

}
