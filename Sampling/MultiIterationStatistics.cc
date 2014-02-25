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

DescribeClass<MultiIterationStatistics,Herwig::GeneralStatistics>
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
