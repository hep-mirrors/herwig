// -*- C++ -*-
//
// GeneralStatictis.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BinnedStatistics class.
//

#include "BinnedStatistics.h"

using namespace Herwig;

BinnedStatistics::~BinnedStatistics() {}

void BinnedStatistics::put(PersistentOStream & os) const {
  os << statisticsMap << weightMap
     << selectorMap << lastPoint;
}

void BinnedStatistics::get(PersistentIStream & is) {
  is >> statisticsMap >> weightMap
     >> selectorMap >> lastPoint;
  lastStatistics = &(statisticsMap.upper_bound(lastPoint)->second);
}

void BinnedStatistics::initialize(unsigned int bins) {
  weightMap[1.] = 1.;
  selectorMap[1.] = make_pair(0.,1.);
  double step = 1./bins;
  for ( unsigned int i = 1; i <= bins; ++i ) {
    statisticsMap[i*step] = GeneralStatistics();
  }
}

