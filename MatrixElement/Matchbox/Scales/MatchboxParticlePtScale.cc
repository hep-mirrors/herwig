// -*- C++ -*-
//
// MatchboxParticlePtScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxParticlePtScale class.
//

#include "MatchboxParticlePtScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxParticlePtScale::MatchboxParticlePtScale() {}

MatchboxParticlePtScale::~MatchboxParticlePtScale() {}

IBPtr MatchboxParticlePtScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxParticlePtScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxParticlePtScale::renormalizationScale() const {
  tcPDVector pd (mePartonData().begin() + 2, mePartonData().end());
  vector<LorentzMomentum> p (meMomenta().begin() + 2, meMomenta().end());

  Energy2 pt2 = ZERO;
  int found = 0;
  tcPDVector::const_iterator itpd = pd.begin();
  for (vector<LorentzMomentum>::const_iterator itp = p.begin() ;
       itp != p.end(); ++itp, ++itpd )
    if ( theMatcher->check(**itpd) ) {
      found++;
      pt2 = (*itp).perp2();
    }
  if ( found == 1 )
    return pt2;
  else
    throw Exception() << "MatchboxParticlePtScale: Found "
		      << found << " particles of the requested type "
		      << "where exactly 1 was expected."
		      << Exception::runerror;

  return ZERO;
	    
}

Energy2 MatchboxParticlePtScale::factorizationScale() const {
  return renormalizationScale();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxParticlePtScale::persistentOutput(PersistentOStream & os) const {
  os << theMatcher;
}

void MatchboxParticlePtScale::persistentInput(PersistentIStream & is, int) {
  is >> theMatcher;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxParticlePtScale,MatchboxScaleChoice>
  describeHerwigMatchboxParticlePtScale("Herwig::MatchboxParticlePtScale", "HwMatchboxScales.so");

void MatchboxParticlePtScale::Init() {

  static ClassDocumentation<MatchboxParticlePtScale> documentation
    ("MatchboxParticlePtScale implements scale choices related to transverse momenta.");

  static Reference<MatchboxParticlePtScale,MatcherBase> interfaceMatcher
    ("Matcher",
     "A matcher to determine the particle that this scale is working on",
     &MatchboxParticlePtScale::theMatcher, false, false, true, false, false);

}

