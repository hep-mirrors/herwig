// -*- C++ -*-
//
// MultiIterationStatictis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MultiIterationStatistics class.
//

#include "MultiIterationStatistics.h"

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/InterfacedBase.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include <cassert>

using namespace Herwig;

MultiIterationStatistics::MultiIterationStatistics() 
  : Interfaced(), GeneralStatistics(),
    theMinIterationPoints(100), theUseAllIterations(false) {}

MultiIterationStatistics::~MultiIterationStatistics() {}

IBPtr MultiIterationStatistics::clone() const {
  return new_ptr(*this);
}

IBPtr MultiIterationStatistics::fullclone() const {
  return new_ptr(*this);
}

void MultiIterationStatistics::put(PersistentOStream & os) const {
  GeneralStatistics::put(os);
  os << theIterations << theMinIterationPoints << theUseAllIterations;
}

void MultiIterationStatistics::get(PersistentIStream & is) {
  GeneralStatistics::get(is);
  is >> theIterations >> theMinIterationPoints >> theUseAllIterations;
}

void MultiIterationStatistics::persistentOutput(PersistentOStream & os) const {
  put(os);
}

void MultiIterationStatistics::persistentInput(PersistentIStream & is, int) {
  get(is);
}

void MultiIterationStatistics::fromXML(const XML::Element& elem) {

  size_t nit;
  elem.getFromAttribute("nIterations",nit);
  elem.getFromAttribute("minIterationPoints",theMinIterationPoints);
  elem.getFromAttribute("useAllIterations",theUseAllIterations);

  theIterations.resize(nit);

  list<XML::Element>::const_iterator stats =
    elem.findFirst(XML::ElementTypes::Element,"GeneralStatistics");

  if ( stats == elem.children().end() )
    throw Exception()
      << "MultiIterationStatistics expected a GeneralStatistics element"
      << Exception::abortnow;

  GeneralStatistics::fromXML(*stats);

  stats = elem.findFirst(XML::ElementTypes::Element,"Iterations");

  if ( stats == elem.children().end() )
    throw Exception()
      << "MultiIterationStatistics expected an Iterations element"
      << Exception::abortnow;

  pair<multimap<pair<int,string>,list<XML::Element>::iterator>::const_iterator,
       multimap<pair<int,string>,list<XML::Element>::iterator>::const_iterator>
    range = stats->findAll(XML::ElementTypes::Element,"GeneralStatistics");

  if ( (size_t)(std::distance(range.first,range.second)) != nit )
    throw Exception() << "MultiIterationStatistics expected "
		      << nit << " iterations but only found "
		      << std::distance(range.first,range.second)
		      << Exception::abortnow;

  set<size_t> which;

  for ( multimap<pair<int,string>,list<XML::Element>::iterator>::const_iterator a = range.first;
	a != range.second; ++a ) {
    size_t id;
    a->second->getFromAttribute("number",id);
    which.insert(id);
    theIterations[id].fromXML(*(a->second));
  }

  if ( which.size() != nit )
    throw Exception() << "MultiIterationStatistics expected "
		      << nit << " iterations but only found "
		      << which.size()
		      << Exception::abortnow;

}

XML::Element MultiIterationStatistics::toXML() const {

  XML::Element elem(XML::ElementTypes::Element,"MultiIterationStatistics");

  elem.appendAttribute("nIterations",theIterations.size());
  elem.appendAttribute("minIterationPoints",theMinIterationPoints);
  elem.appendAttribute("useAllIterations",theUseAllIterations);

  XML::Element stats = GeneralStatistics::toXML();
  elem.append(stats);

  XML::Element its(XML::ElementTypes::Element,"Iterations");
  for ( size_t k = 0; k < theIterations.size(); ++k ) {
    XML::Element it = theIterations[k].toXML();
    it.appendAttribute("number",k);
    its.append(it);
  }
  elem.append(its);

  return elem;

}


double MultiIterationStatistics::chi2() const {
  assert(!iterations().empty());
  double current = averageWeight(true);
  double res = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < minIterationPoints() || s->averageWeightVariance() == 0.0 )
      continue;
    res += sqr(s->averageWeight()-current)/s->averageWeightVariance();
  }
  res += 
    selectedPoints() >= minIterationPoints() && GeneralStatistics::averageWeightVariance() != 0.0 ?
    sqr(GeneralStatistics::averageWeight()-current)/
    GeneralStatistics::averageWeightVariance() : 0.;
  res /= iterations().size();
  return res;
}

double MultiIterationStatistics::averageWeight(bool useAll) const {
  double invSigmaBar = 0.;
  double res = 0.;
  if ( useAllIterations() || useAll ) {
    for ( vector<GeneralStatistics>::const_iterator s =
	    iterations().begin(); s != iterations().end(); ++s ) {
      if ( s->selectedPoints() < minIterationPoints() || s->averageWeightVariance() == 0.0 )
	continue;
      invSigmaBar += 1./s->averageWeightVariance();
      res += s->averageWeight()/s->averageWeightVariance();
    }
  }
  invSigmaBar += 
    selectedPoints() >= minIterationPoints() && GeneralStatistics::averageWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageWeightVariance() : 0.;
  res += 
    selectedPoints() >= minIterationPoints() && GeneralStatistics::averageWeightVariance() != 0.0 ?
    GeneralStatistics::averageWeight()/GeneralStatistics::averageWeightVariance() : 0.;
  if ( invSigmaBar != 0.0 )
    res /= invSigmaBar;
  return res;
}

double MultiIterationStatistics::averageWeightVariance(bool useAll) const {
  double invSigmaBar = 0.;
  if ( useAllIterations() || useAll ) {
    for ( vector<GeneralStatistics>::const_iterator s =
	    iterations().begin(); s != iterations().end(); ++s ) {
      if ( s->selectedPoints() < minIterationPoints() || s->averageWeightVariance() == 0.0 )
	continue;
      invSigmaBar += 1./s->averageWeightVariance();
    }
  }
  invSigmaBar += 
    selectedPoints() >= minIterationPoints() && GeneralStatistics::averageWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageWeightVariance() : 0.;
  return invSigmaBar != 0.0 ? 1./invSigmaBar : 0.0;
}

double MultiIterationStatistics::averageAbsWeight(bool useAll) const {
  double invSigmaBar = 0.;
  double res = 0.;
  if ( useAllIterations() || useAll ) {
    for ( vector<GeneralStatistics>::const_iterator s =
	    iterations().begin(); s != iterations().end(); ++s ) {
      if ( s->selectedPoints() < minIterationPoints() || s->averageAbsWeightVariance() == 0.0 )
	continue;
      invSigmaBar += 1./s->averageAbsWeightVariance();
      res += s->averageAbsWeight()/s->averageAbsWeightVariance();
    }
  }
  invSigmaBar += 
    selectedPoints() >= minIterationPoints() && GeneralStatistics::averageAbsWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageAbsWeightVariance() : 0.;
  res += 
    selectedPoints() >= minIterationPoints() && GeneralStatistics::averageAbsWeightVariance() != 0.0 ?
    GeneralStatistics::averageAbsWeight()/GeneralStatistics::averageAbsWeightVariance() : 0.;
  if ( invSigmaBar != 0.0 )
    res /= invSigmaBar;
  return invSigmaBar != 0.0 ? res : 0.0;
}

double MultiIterationStatistics::averageAbsWeightVariance(bool useAll) const {
  double invSigmaBar = 0.;
  if ( useAllIterations() || useAll ) {
    for ( vector<GeneralStatistics>::const_iterator s =
	    iterations().begin(); s != iterations().end(); ++s ) {
      if ( s->selectedPoints() < minIterationPoints() || s->averageAbsWeightVariance() == 0.0 )
	continue;
      invSigmaBar += 1./s->averageAbsWeightVariance();
    }
  }
  invSigmaBar += 
    selectedPoints() >= minIterationPoints() && GeneralStatistics::averageAbsWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageAbsWeightVariance() : 0.;
  return invSigmaBar != 0.0 ? 1./invSigmaBar : 0.0;
}

DescribeClass<MultiIterationStatistics,Interfaced>
  describeHerwigMultiIterationStatistics("Herwig::MultiIterationStatistics", "HwSampling.so");

void MultiIterationStatistics::Init() {

  static ClassDocumentation<MultiIterationStatistics> documentation
    ("MultiIterationStatistics");

  static ThePEG::Parameter<MultiIterationStatistics,unsigned int> interfaceMinIterationPoints
    ("MinIterationPoints",
     "Discard iterations with less than the given number of points.",
     &MultiIterationStatistics::theMinIterationPoints, 100, 2, 0,
     false, false, ThePEG::Interface::lowerlim);

  static Switch<MultiIterationStatistics,bool> interfaceUseAllIterations
    ("UseAllIterations",
     "Combine integrals from all iterations.",
     &MultiIterationStatistics::theUseAllIterations, false, false, false);
  static SwitchOption interfaceUseAllIterationsYes
    (interfaceUseAllIterations,
     "Yes",
     "Combine integrals from all iterations.",
     true);
  static SwitchOption interfaceUseAllIterationsNo
    (interfaceUseAllIterations,
     "No",
     "Calculate integral from last iteration only.",
     true);

}
