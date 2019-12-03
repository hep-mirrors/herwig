// -*- C++ -*-
//
// MatchboxLeptonMassScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxLeptonMassScale class.
//

#include "MatchboxLeptonMassScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxLeptonMassScale::MatchboxLeptonMassScale() {}

MatchboxLeptonMassScale::~MatchboxLeptonMassScale() {}

IBPtr MatchboxLeptonMassScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxLeptonMassScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxLeptonMassScale::renormalizationScale() const {

  int firstLepton = -1;
  int secondLepton = -1;

  for ( size_t k = 0; k < mePartonData().size(); ++k ) {
    if ( abs(mePartonData()[k]->id()) > 10 && 
	 abs(mePartonData()[k]->id()) < 17 ) {
      if ( firstLepton < 0 ) {
	firstLepton = k;
      } else if ( secondLepton < 0 ) {
	secondLepton = k;
      } else break;
    }
  }

  if ( firstLepton < 0 || secondLepton < 0 )
    throw Exception() << "MatchboxLeptonMassScale::renormalizationScale(): "
		      << "No lepton pair could be found. Check your setup."
		      << Exception::runerror;

  if ( (firstLepton < 2 && secondLepton > 1) || 
       (firstLepton > 1 && secondLepton < 2) )
    return abs((meMomenta()[firstLepton] -
		meMomenta()[secondLepton]).m2());

  return
    (meMomenta()[firstLepton] +
     meMomenta()[secondLepton]).m2();

}

Energy2 MatchboxLeptonMassScale::factorizationScale() const {
  return renormalizationScale();
}

Energy2 MatchboxLeptonMassScale::renormalizationScaleQED() const {
  return renormalizationScale();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxLeptonMassScale::persistentOutput(PersistentOStream &) const {}

void MatchboxLeptonMassScale::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxLeptonMassScale,MatchboxScaleChoice>
  describeHerwigMatchboxLeptonMassScale("Herwig::MatchboxLeptonMassScale", "HwMatchboxScales.so");

void MatchboxLeptonMassScale::Init() {

  static ClassDocumentation<MatchboxLeptonMassScale> documentation
    ("MatchboxLeptonMassScale implements scale choices related "
     "to lepton pair invariant masses.");


}

