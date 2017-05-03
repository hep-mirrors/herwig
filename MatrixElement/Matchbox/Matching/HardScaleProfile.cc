// -*- C++ -*-
//
// HardScaleProfile.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardScaleProfile class.
//

#include "HardScaleProfile.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HardScaleProfile::HardScaleProfile() 
  : theFixedHardScale(ZERO), theProfileRho(0.3),
    theProfileType(resummation) {}

HardScaleProfile::~HardScaleProfile() {}

IBPtr HardScaleProfile::clone() const {
  return new_ptr(*this);
}

IBPtr HardScaleProfile::fullclone() const {
  return new_ptr(*this);
}

double HardScaleProfile::hardScaleProfile(Energy hard, Energy soft) const {
  if ( theFixedHardScale > ZERO )
    hard = theFixedHardScale;
  double x = soft/hard;
  if ( theProfileType == theta ) {
    return x <= 1. ? 1.0 : 0.0;
  }
  if ( theProfileType == resummation ) {
    if ( x > 1. ) {
      return 0.;
    } else if ( x <= 1. && x > 1. - theProfileRho ) {
      return sqr(1.-x)/(2.*sqr(theProfileRho));
    } else if ( x <= 1. - theProfileRho &&
		x > 1. - 2.*theProfileRho ) {
      return 1. - sqr(1.-2.*theProfileRho-x)/(2.*sqr(theProfileRho));
    } else {
      return 1.;
    }
  }
  if ( theProfileType == hfact ) {
    return 1./(1.+x);
  }
  return 1.;
}

bool HardScaleProfile::unrestrictedPhasespace() const {
  if ( theProfileType == theta ||
       theProfileType == resummation ) {
    return false;
  }
  return true;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HardScaleProfile::persistentOutput(PersistentOStream & os) const {
  os << ounit(theFixedHardScale,GeV) << theProfileRho << theProfileType;
}

void HardScaleProfile::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theFixedHardScale,GeV) >> theProfileRho >> theProfileType;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<HardScaleProfile,Interfaced>
  describeHerwigHardScaleProfile("Herwig::HardScaleProfile", "Herwig.so");

void HardScaleProfile::Init() {

  static ClassDocumentation<HardScaleProfile> documentation
    ("HardScaleProfile implements profile scales.");

  static Parameter<HardScaleProfile,Energy> interfaceFixedHardScale
    ("FixedHardScale",
     "A fixed hard scale to be used instead of the process specific choice.",
     &HardScaleProfile::theFixedHardScale, GeV, ZERO, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<HardScaleProfile,double> interfaceProfileRho
    ("ProfileRho",
     "The profile width parameter",
     &HardScaleProfile::theProfileRho, 0.3, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<HardScaleProfile,int> interfaceProfileType
    ("ProfileType",
     "The type of profile to use.",
     &HardScaleProfile::theProfileType, resummation, false, false);
  static SwitchOption interfaceProfileTypeTheta
    (interfaceProfileType,
     "Theta",
     "Use a hard cutoff.",
     theta);
  static SwitchOption interfaceProfileTypeResummation
    (interfaceProfileType,
     "Resummation",
     "Use the resummation profile with quadratic interpolation.",
     resummation);
  static SwitchOption interfaceProfileTypeHFact
    (interfaceProfileType,
     "HFact",
     "Use the hfact profile.",
     hfact);

}

