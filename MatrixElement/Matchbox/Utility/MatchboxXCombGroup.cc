// -*- C++ -*-
//
// MatchboxXCombGroup.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxXCombGroup class.
//

#include "MatchboxXCombGroup.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

MatchboxXCombGroup::MatchboxXCombGroup() 
  : StdXCombGroup() {
  flushCaches();
}

MatchboxXCombGroup::~MatchboxXCombGroup() {}

MatchboxXCombGroup::MatchboxXCombGroup(Energy newMaxEnergy, const cPDPair & inc,
				       tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
				       tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
				       const PBPair & newPartonBins, tCutsPtr newCuts, tMEGroupPtr newME,
				       const DiagramVector & newDiagrams, bool mir,
				       tStdXCombPtr newHead)
  : StdXCombGroup(newMaxEnergy, inc,
		  newEventHandler, newSubProcessHandler,
		  newExtractor, newCKKW,
		  newPartonBins, newCuts, newME,
		  newDiagrams, mir,
		  newHead),
  MatchboxXCombData(newME) {
  flushCaches();
}

void MatchboxXCombGroup::clean() {
  StdXCombGroup::clean();
  flushCaches();
}

void MatchboxXCombGroup::persistentOutput(PersistentOStream & os) const {
  MatchboxXCombData::persistentOutput(os);
}

void MatchboxXCombGroup::persistentInput(PersistentIStream & is, int version) {
  MatchboxXCombData::persistentInput(is,version);
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxXCombGroup,StdXCombGroup>
  describeHerwigMatchboxXCombGroup("Herwig::MatchboxXCombGroup", "Herwig.so");

void MatchboxXCombGroup::Init() {

  static ClassDocumentation<MatchboxXCombGroup> documentation
    ("MatchboxXCombGroup extents StdXCombGroup.");

}

