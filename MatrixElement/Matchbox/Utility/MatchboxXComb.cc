// -*- C++ -*-
//
// MatchboxXComb.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxXComb class.
//

#include "MatchboxXComb.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

MatchboxXComb::MatchboxXComb() 
  : StandardXComb() {
  flushCaches();
}

MatchboxXComb::~MatchboxXComb() {}

MatchboxXComb::MatchboxXComb(Energy newMaxEnergy, const cPDPair & inc,
	      tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
	      tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
	      const PBPair & newPartonBins, tCutsPtr newCuts, tMEPtr newME,
	      const DiagramVector & newDiagrams, bool mir,
	      tStdXCombPtr newHead)
  : StandardXComb(newMaxEnergy, inc,
		  newEventHandler,newSubProcessHandler,
		  newExtractor, newCKKW,
		  newPartonBins, newCuts, newME,
		  newDiagrams, mir,
		  newHead),
    MatchboxXCombData(newME) {
  flushCaches();
}

MatchboxXComb::MatchboxXComb(tStdXCombPtr newHead,
			     const PBPair & newPartonBins, tMEPtr newME,
			     const DiagramVector & newDiagrams)
  : StandardXComb(newHead, newPartonBins, newME, newDiagrams),
    MatchboxXCombData(newME) {
  flushCaches();
}

void MatchboxXComb::clean() {
  StandardXComb::clean();
  flushCaches();
}

void MatchboxXComb::persistentOutput(PersistentOStream & os) const {
  MatchboxXCombData::persistentOutput(os);
}

void MatchboxXComb::persistentInput(PersistentIStream & is, int version) {
  MatchboxXCombData::persistentInput(is,version);
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxXComb,StandardXComb>
  describeHerwigMatchboxXComb("Herwig::MatchboxXComb", "Herwig.so");

void MatchboxXComb::Init() {

  static ClassDocumentation<MatchboxXComb> documentation
    ("MatchboxXComb extents StandardXComb.");

}

