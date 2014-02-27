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
// functions of the GeneralStatistics class.
//

#include "GeneralStatistics.h"

using namespace Herwig;

GeneralStatistics::~GeneralStatistics() {}

void GeneralStatistics::put(PersistentOStream & os) const {
  os << theMaxWeight << theMinWeight << theSumWeights
     << theSumSquaredWeights << theSumAbsWeights 
     << theSelectedPoints << theAcceptedPoints
     << theNanPoints << theAllPoints << theLastWeight
     << theBias;
}

void GeneralStatistics::get(PersistentIStream & is) {
  is >> theMaxWeight >> theMinWeight >> theSumWeights
     >> theSumSquaredWeights >> theSumAbsWeights 
     >> theSelectedPoints >> theAcceptedPoints
     >> theNanPoints >> theAllPoints >> theLastWeight
     >> theBias;
}

