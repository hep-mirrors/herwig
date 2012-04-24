// -*- C++ -*-
//
// MatchboxPtScale.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxPtScale class.
//

#include "MatchboxPtScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxPtScale::MatchboxPtScale() {}

MatchboxPtScale::~MatchboxPtScale() {}

IBPtr MatchboxPtScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxPtScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxPtScale::renormalizationScale() const {
  cPDVector::const_iterator pd = mePartonData().begin() + 2;
  vector<Lorentz5Momentum>::const_iterator p = meMomenta().begin() + 2;
  Energy2 maxpt2 = ZERO;
  for ( ; p != meMomenta().end(); ++p, ++pd )
    if ( (**pd).coloured() )
      maxpt2 = max(maxpt2,(*p).perp2());
  return maxpt2;
}

Energy2 MatchboxPtScale::factorizationScale() const {
  return renormalizationScale();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxPtScale::persistentOutput(PersistentOStream &) const {}

void MatchboxPtScale::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxPtScale,MatchboxScaleChoice>
  describeHerwigMatchboxPtScale("Herwig::MatchboxPtScale", "HwMatchbox.so");

void MatchboxPtScale::Init() {

  static ClassDocumentation<MatchboxPtScale> documentation
    ("MatchboxPtScale implements scale choices related to transverse momenta.");


}

