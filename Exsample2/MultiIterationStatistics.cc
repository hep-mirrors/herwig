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

#include <cassert>

using namespace Herwig;

MultiIterationStatistics::MultiIterationStatistics() 
  : GeneralStatistics() {}

MultiIterationStatistics::~MultiIterationStatistics() {}

void MultiIterationStatistics::put(PersistentOStream & os) const {
  GeneralStatistics::put(os);
  os << theIterations;
}

void MultiIterationStatistics::get(PersistentIStream & is) {
  GeneralStatistics::get(is);
  is >> theIterations;
}

double MultiIterationStatistics::chi2() const {
  assert(!iterations().empty());
  double current = averageWeight();
  double res = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < 2 || s->averageWeightVariance() == 0.0 )
      continue;
    res += sqr(s->averageWeight()-current)/s->averageWeightVariance();
  }
  res += 
    selectedPoints() > 1 && GeneralStatistics::averageWeightVariance() != 0.0 ?
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
    if ( s->selectedPoints() < 2 || s->averageWeightVariance() == 0.0 )
      continue;
    invSigmaBar += 1./s->averageWeightVariance();
    res += s->averageWeight()/s->averageWeightVariance();
  }
  invSigmaBar += 
    selectedPoints() > 1 && GeneralStatistics::averageWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageWeightVariance() : 0.;
  res += 
    selectedPoints() > 1 && GeneralStatistics::averageWeightVariance() != 0.0 ?
    GeneralStatistics::averageWeight()/GeneralStatistics::averageWeightVariance() : 0.;
  res /= invSigmaBar;
  return res;
}

double MultiIterationStatistics::averageWeightVariance() const {
  double invSigmaBar = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < 2 || s->averageWeightVariance() == 0.0 )
      continue;
    invSigmaBar += 1./s->averageWeightVariance();
  }
  invSigmaBar += 
    selectedPoints() > 1 && GeneralStatistics::averageWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageWeightVariance() : 0.;
  return 1./invSigmaBar;
}

double MultiIterationStatistics::averageAbsWeight() const {
  double invSigmaBar = 0.;
  double res = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < 2 || s->averageAbsWeightVariance() == 0.0 )
      continue;
    invSigmaBar += 1./s->averageAbsWeightVariance();
    res += s->averageAbsWeight()/s->averageAbsWeightVariance();
  }
  invSigmaBar += 
    selectedPoints() > 1 && GeneralStatistics::averageAbsWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageAbsWeightVariance() : 0.;
  res += 
    selectedPoints() > 1 && GeneralStatistics::averageAbsWeightVariance() != 0.0 ?
    GeneralStatistics::averageAbsWeight()/GeneralStatistics::averageAbsWeightVariance() : 0.;
  res /= invSigmaBar;
  return res;
}

double MultiIterationStatistics::averageAbsWeightVariance() const {
  double invSigmaBar = 0.;
  for ( vector<GeneralStatistics>::const_iterator s =
	  iterations().begin(); s != iterations().end(); ++s ) {
    if ( s->selectedPoints() < 2 || s->averageAbsWeightVariance() == 0.0 )
      continue;
    invSigmaBar += 1./s->averageAbsWeightVariance();
  }
  invSigmaBar += 
    selectedPoints() > 1 && GeneralStatistics::averageAbsWeightVariance() != 0.0 ?
    1./GeneralStatistics::averageAbsWeightVariance() : 0.;
  return 1./invSigmaBar;
}
