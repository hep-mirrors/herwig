// -*- C++ -*-
//
// MatchboxScaleChoice.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxScaleChoice class.
//

#include "MatchboxScaleChoice.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxScaleChoice::MatchboxScaleChoice() {}

MatchboxScaleChoice::~MatchboxScaleChoice() {}

IBPtr MatchboxScaleChoice::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxScaleChoice::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxScaleChoice::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb;
}

void MatchboxScaleChoice::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxScaleChoice,HandlerBase>
  describeHerwigMatchboxScaleChoice("Herwig::MatchboxScaleChoice", "HwMatchbox.so");

void MatchboxScaleChoice::Init() {

  static ClassDocumentation<MatchboxScaleChoice> documentation
    ("MatchboxScaleChoice is the base class for scale choices "
     "within Matchbox.");

}

