// -*- C++ -*-
//
// MultiIterationStatictis.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "ThePEG/Utilities/DescribeClass.h"

#include <cassert>

using namespace Herwig;

MultiIterationStatistics::MultiIterationStatistics() 
  : Interfaced(), GeneralStatistics(), theMinNumEventsPerIteration(100) {}

MultiIterationStatistics::~MultiIterationStatistics() {}

IBPtr MultiIterationStatistics::clone() const {
  return new_ptr(*this);
}

IBPtr MultiIterationStatistics::fullclone() const {
  return new_ptr(*this);
}

void MultiIterationStatistics::put(PersistentOStream & os) const {
  GeneralStatistics::put(os);
  os << theIterations << theMinNumEventsPerIteration;
}

void MultiIterationStatistics::get(PersistentIStream & is) {
  GeneralStatistics::get(is);
  is >> theIterations >> theMinNumEventsPerIteration;
}

void MultiIterationStatistics::persistentOutput(PersistentOStream & os) const {
  put(os);
}

void MultiIterationStatistics::persistentInput(PersistentIStream & is, int) {
  get(is);
}


double MultiIterationStatistics::chi2() const {
  assert(!iterations().empty());
  double current = averageWeight();
  double res = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < theMinNumEventsPerIteration || s->averageWeightVariance() == 0.0 )
      continue;
    res += sqr(s->averageWeight()-current)/s->averageWeightVariance();
  }
  res += 
    selectedPoints() >= theMinNumEventsPerIteration && GeneralStatistics::averageWeightVariance() != 0.0 ?
    sqr(GeneralStatistics::averageWeight()-current)/
    GeneralStatistics::averageWeightVariance() : 0.;
  res /= iterations().size();
  return res;
}

double MultiIterationStatistics::averageWeight() const {
  double invSigmaBar = 0.;
  double res = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < theMinNumEventsPerIteration || s->averageWeightVariance() == 0.0 )
      continue;
    invSigmaBar += 1./s->averageWeightVariance();
    res += s->averageWeight()/s->averageWeightVariance();
  }
  invSigmaBar += 
    selectedPoints() >= theMinNumEventsPerIteration && GeneralStatistics::averageWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageWeightVariance() : 0.;
  res += 
    selectedPoints() >= theMinNumEventsPerIteration && GeneralStatistics::averageWeightVariance() != 0.0 ?
    GeneralStatistics::averageWeight()/GeneralStatistics::averageWeightVariance() : 0.;
  if ( invSigmaBar != 0.0 )
    res /= invSigmaBar;
  return res;
}

double MultiIterationStatistics::averageWeightVariance() const {
  double invSigmaBar = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < theMinNumEventsPerIteration || s->averageWeightVariance() == 0.0 )
      continue;
    invSigmaBar += 1./s->averageWeightVariance();
  }
  invSigmaBar += 
    selectedPoints() >= theMinNumEventsPerIteration && GeneralStatistics::averageWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageWeightVariance() : 0.;
  return invSigmaBar != 0.0 ? 1./invSigmaBar : 0.0;
}

double MultiIterationStatistics::averageAbsWeight() const {
  double invSigmaBar = 0.;
  double res = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < theMinNumEventsPerIteration || s->averageAbsWeightVariance() == 0.0 )
      continue;
    invSigmaBar += 1./s->averageAbsWeightVariance();
    res += s->averageAbsWeight()/s->averageAbsWeightVariance();
  }
  invSigmaBar += 
    selectedPoints() >= theMinNumEventsPerIteration && GeneralStatistics::averageAbsWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageAbsWeightVariance() : 0.;
  res += 
    selectedPoints() >= theMinNumEventsPerIteration && GeneralStatistics::averageAbsWeightVariance() != 0.0 ?
    GeneralStatistics::averageAbsWeight()/GeneralStatistics::averageAbsWeightVariance() : 0.;
  if ( invSigmaBar != 0.0 )
    res /= invSigmaBar;
  return invSigmaBar != 0.0 ? res : 0.0;
}

double MultiIterationStatistics::averageAbsWeightVariance() const {
  double invSigmaBar = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < theMinNumEventsPerIteration || s->averageAbsWeightVariance() == 0.0 )
      continue;
    invSigmaBar += 1./s->averageAbsWeightVariance();
  }
  invSigmaBar += 
    selectedPoints() >= theMinNumEventsPerIteration && GeneralStatistics::averageAbsWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageAbsWeightVariance() : 0.;
  return 1./invSigmaBar;
}



DescribeClass<MultiIterationStatistics,Herwig::GeneralStatistics>
  describeHerwigMultiIterationStatistics("Herwig::MultiIterationStatistics", "HwExsample2.so");

void MultiIterationStatistics::Init() {

  static ClassDocumentation<MultiIterationStatistics> documentation
    ("MultiIterationStatistics");

  static ThePEG::Parameter<MultiIterationStatistics,unsigned int> interfaceMinNumEventsPerIteration
    ("MinNumEventsPerIteration",
     "Set the number of presampling points per cell",
     &MultiIterationStatistics::theMinNumEventsPerIteration, 100, 2, 0,
     false, false, ThePEG::Interface::lowerlim);

}
