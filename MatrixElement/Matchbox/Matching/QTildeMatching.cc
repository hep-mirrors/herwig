// -*- C++ -*-
//
// QTildeMatching.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeMatching class.
//

#include "QTildeMatching.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

QTildeMatching::QTildeMatching() {}

QTildeMatching::~QTildeMatching() {}

IBPtr QTildeMatching::clone() const {
  return new_ptr(*this);
}

IBPtr QTildeMatching::fullclone() const {
  return new_ptr(*this);
}

Energy QTildeMatching::hardScale() const {
  return ZERO;
}

double QTildeMatching::hardScaleProfile(Energy, Energy) const {
  return 1.;
}

bool QTildeMatching::isInShowerPhasespace() const {
  return false;
}

bool QTildeMatching::isAboveCutoff() const {
  return false;
}

CrossSection QTildeMatching::dSigHatDR() const {

  return ZERO;

}

double QTildeMatching::me2() const {
  throw Exception() << "Not intented to use. Disable the ShowerApproximationGenerator."
		    << Exception::abortnow;
  return 0.;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void QTildeMatching::persistentOutput(PersistentOStream & os) const {
  os << theLargeNBasis;
}

void QTildeMatching::persistentInput(PersistentIStream & is, int) {
  is >> theLargeNBasis;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<QTildeMatching,Herwig::ShowerApproximation>
  describeHerwigQTildeMatching("Herwig::QTildeMatching", "HwMatchbox.so");

void QTildeMatching::Init() {

  static ClassDocumentation<QTildeMatching> documentation
    ("QTildeMatching implements NLO matching with the default shower.");

  static Reference<QTildeMatching,ColourBasis> interfaceLargeNBasis
    ("LargeNBasis",
     "Set the large-N colour basis implementation.",
     &QTildeMatching::theLargeNBasis, false, false, true, true, false);

}

