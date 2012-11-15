// -*- C++ -*-
//
// ShowerApproximation.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerApproximation class.
//

#include "ShowerApproximation.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

using namespace Herwig;

ShowerApproximation::ShowerApproximation() 
  : HandlerBase(), theBelowCutoff(false) {}

ShowerApproximation::~ShowerApproximation() {}

void ShowerApproximation::setDipole(Ptr<SubtractionDipole>::tcptr dip) { theDipole = dip; }

Ptr<SubtractionDipole>::tcptr ShowerApproximation::dipole() const { return theDipole; }


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ShowerApproximation::persistentOutput(PersistentOStream & os) const {
  os << theBornXComb << theRealXComb << theDipole << theBelowCutoff;
}

void ShowerApproximation::persistentInput(PersistentIStream & is, int) {
  is >> theBornXComb >> theRealXComb >> theDipole >> theBelowCutoff;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<ShowerApproximation,HandlerBase>
  describeHerwigShowerApproximation("Herwig::ShowerApproximation", "HwMatchbox.so");

void ShowerApproximation::Init() {

  static ClassDocumentation<ShowerApproximation> documentation
    ("ShowerApproximation describes the shower emission to be used "
     "in NLO matching.");

}

