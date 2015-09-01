// -*- C++ -*-
//
// MatchboxTopMassScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2014 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxTopMassScale class.
//

#include "MatchboxTopMassScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxTopMassScale::MatchboxTopMassScale() {}

MatchboxTopMassScale::~MatchboxTopMassScale() {}

IBPtr MatchboxTopMassScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxTopMassScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxTopMassScale::renormalizationScale() const {
  
  size_t k = 2;
  int top = -1;
  int antitop = -1;

  while ( (top == -1 || antitop == -1) && k < mePartonData().size() ){
    if ( mePartonData()[k]->id() == 6 ) {
      if ( top < 0 )
	top = k;
      else
	assert(false);
    } else if ( mePartonData()[k]->id() == -6 ) {
      if ( antitop < 0 )
	antitop = k;
      else
	assert(false);
    }
    k++;
  }

  if ( top < 2 || antitop < 2 ){
    throw Exception() << "MatchboxTopMassScale: Could not find a top-antitop-pair in the final state!\n"
		      << Exception::runerror;
  }
  // cerr << " sqrt(TopMassScale)  = "
  //      << sqrt((meMomenta()[top]+meMomenta()[antitop]).m2())/GeV
  //      << "\n" << flush;
  return((meMomenta()[top]+meMomenta()[antitop]).m2());

}

Energy2 MatchboxTopMassScale::factorizationScale() const {
  return(renormalizationScale());
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxTopMassScale::persistentOutput(PersistentOStream &) const {}

void MatchboxTopMassScale::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxTopMassScale,MatchboxScaleChoice>
  describeHerwigMatchboxTopMassScale("Herwig::MatchboxTopMassScale", "HwMatchboxScales.so");

void MatchboxTopMassScale::Init() {

  static ClassDocumentation<MatchboxTopMassScale> documentation
    ("MatchboxTopMassScale implements invariant mass of top-antitop pair as scale choice.");


}

