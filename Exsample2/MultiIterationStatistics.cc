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

