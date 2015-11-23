// -*- C++ -*-
//
// MatchboxLeptonPtScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxLeptonPtScale class.
//

#include "MatchboxLeptonPtScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxLeptonPtScale::MatchboxLeptonPtScale() {}

MatchboxLeptonPtScale::~MatchboxLeptonPtScale() {}

IBPtr MatchboxLeptonPtScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxLeptonPtScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxLeptonPtScale::renormalizationScale() const {

  int firstLepton = -1;
  int secondLepton = -1;

  for ( size_t k = 2; k < mePartonData().size(); ++k ) {
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
    throw Exception() << "MatchboxLeptonPtScale::renormalizationScale(): "
		      << "No lepton pair could be found. Check your setup."
		      << Exception::runerror;

  return
    (meMomenta()[firstLepton] +
     meMomenta()[secondLepton]).perp2();

}

Energy2 MatchboxLeptonPtScale::factorizationScale() const {
  return renormalizationScale();
}

Energy2 MatchboxLeptonPtScale::renormalizationScaleQED() const {

  int firstLepton = -1;
  int secondLepton = -1;

  for ( size_t k = 2; k < mePartonData().size(); ++k ) {
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
    throw Exception() << "MatchboxLeptonPtScale::renormalizationScaleQED(): "
		      << "No lepton pair could be found. Check your setup."
		      << Exception::runerror;

  return
    (meMomenta()[firstLepton] +
     meMomenta()[secondLepton]).m2();

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxLeptonPtScale::persistentOutput(PersistentOStream &) const {}

void MatchboxLeptonPtScale::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxLeptonPtScale,MatchboxScaleChoice>
  describeHerwigMatchboxLeptonPtScale("Herwig::MatchboxLeptonPtScale", "HwMatchboxScales.so");

void MatchboxLeptonPtScale::Init() {

  static ClassDocumentation<MatchboxLeptonPtScale> documentation
    ("MatchboxLeptonPtScale implements scale choices related "
     "to lepton pair transverse momenta.");


}

