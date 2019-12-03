// -*- C++ -*-
//
// MatchboxZGammaAmplitude.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxZGammaAmplitude class.
//

#include "MatchboxZGammaAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxZGammaAmplitude::MatchboxZGammaAmplitude() 
  : MatchboxAmplitude(), theIncludeZ(true), theIncludeGamma(true) {}

MatchboxZGammaAmplitude::~MatchboxZGammaAmplitude() {}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxZGammaAmplitude::persistentOutput(PersistentOStream & os) const {
  os << theIncludeZ << theIncludeGamma;
}

void MatchboxZGammaAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> theIncludeZ >> theIncludeGamma;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxZGammaAmplitude,MatchboxAmplitude>
  describeHerwigMatchboxZGammaAmplitude("Herwig::MatchboxZGammaAmplitude", "HwMatchboxBuiltin.so");

void MatchboxZGammaAmplitude::Init() {

  static ClassDocumentation<MatchboxZGammaAmplitude> documentation
    ("There is no documentation for the MatchboxZGammaAmplitude class");

  static Switch<MatchboxZGammaAmplitude,bool> interfaceIncludeZ
    ("IncludeZ",
     "Include the Z contribution.",
     &MatchboxZGammaAmplitude::theIncludeZ, true, false, false);
  static SwitchOption interfaceIncludeZYes
    (interfaceIncludeZ,
     "Yes",
     "",
     true);
  static SwitchOption interfaceIncludeZNo
    (interfaceIncludeZ,
     "No",
     "",
     false);

  static Switch<MatchboxZGammaAmplitude,bool> interfaceIncludeGamma
    ("IncludeGamma",
     "Include the photon contribution.",
     &MatchboxZGammaAmplitude::theIncludeGamma, true, false, false);
  static SwitchOption interfaceIncludeGammaYes
    (interfaceIncludeGamma,
     "Yes",
     "",
     true);
  static SwitchOption interfaceIncludeGammaNo
    (interfaceIncludeGamma,
     "No",
     "",
     false);

}

