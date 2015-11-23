// -*- C++ -*-
//
// MatchboxSHatScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxSHatScale class.
//

#include "MatchboxSHatScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxSHatScale::MatchboxSHatScale() {}

MatchboxSHatScale::~MatchboxSHatScale() {}

IBPtr MatchboxSHatScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxSHatScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxSHatScale::renormalizationScale() const {
  return lastSHat();
}

Energy2 MatchboxSHatScale::factorizationScale() const {
  return renormalizationScale();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxSHatScale::persistentOutput(PersistentOStream &) const {}

void MatchboxSHatScale::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxSHatScale,MatchboxScaleChoice>
  describeHerwigMatchboxSHatScale("Herwig::MatchboxSHatScale", "HwMatchboxScales.so");

void MatchboxSHatScale::Init() {

  static ClassDocumentation<MatchboxSHatScale> documentation
    ("MatchboxSHatScale implements lastSHat() as scale choice.");


}

