// -*- C++ -*-
//
// DipoleMatching.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleMatching class.
//

#include "DipoleMatching.h"
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

DipoleMatching::DipoleMatching() 
  : theShowerKernels(true) {}

DipoleMatching::~DipoleMatching() {}

IBPtr DipoleMatching::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleMatching::fullclone() const {
  return new_ptr(*this);
}

CrossSection DipoleMatching::dSigHatDR() const {

  double xme2 = 0.;

  if ( theShowerKernels ) {
    xme2 = dipole()->me2();
  } else {
    pair<int,int> ij(dipole()->bornEmitter(),
		     dipole()->bornSpectator());
    double ccme2 = 
      dipole()->underlyingBornME()->largeNColourCorrelatedME2(ij,theLargeNBasis);
    xme2 = dipole()->me2Avg(-ccme2);
  }

  xme2 /= dipole()->underlyingBornME()->lastXComb().lastAlphaS();
  xme2 *= bornPDFWeight(dipole()->underlyingBornME()->lastScale());

  return
    sqr(hbarc) * 
    realXComb()->jacobian() * 
    subtractionScaleWeight() *
    xme2 /
    (2. * realXComb()->lastSHat());

}

double DipoleMatching::me2() const {
  throw Exception() << "Not intented to use. Disable the ShowerApproximationGenerator."
		    << Exception::abortnow;
  return 0.;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleMatching::persistentOutput(PersistentOStream & os) const {
  os << theShowerKernels << theLargeNBasis;
}

void DipoleMatching::persistentInput(PersistentIStream & is, int) {
  is >> theShowerKernels >> theLargeNBasis;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DipoleMatching,Herwig::ShowerApproximation>
  describeHerwigDipoleMatching("Herwig::DipoleMatching", "HwMatchbox.so");

void DipoleMatching::Init() {

  static ClassDocumentation<DipoleMatching> documentation
    ("DipoleMatching implements NLO matching with the dipole shower.");

  static Reference<DipoleMatching,ColourBasis> interfaceLargeNBasis
    ("LargeNBasis",
     "Set the large-N colour basis implementation.",
     &DipoleMatching::theLargeNBasis, false, false, true, true, false);


  static Switch<DipoleMatching,bool> interfaceShowerKernels
    ("ShowerKernels",
     "Switch between exact and shower approximated dipole functions.",
     &DipoleMatching::theShowerKernels, true, false, false);
  static SwitchOption interfaceShowerKernelsOn
    (interfaceShowerKernels,
     "On",
     "Switch to shower approximated dipole functions.",
     true);
  static SwitchOption interfaceShowerKernelsOff
    (interfaceShowerKernels,
     "Off",
     "Switch to full dipole functions.",
     false);

}

