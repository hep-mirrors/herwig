// -*- C++ -*-
//
// FFLightTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFLightTildeKinematics class.
//

#include "FFLightTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFLightTildeKinematics::FFLightTildeKinematics() {}

FFLightTildeKinematics::~FFLightTildeKinematics() {}

IBPtr FFLightTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFLightTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FFLightTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double y = emission*emitter / (emission*emitter + emission*spectator + emitter*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = z;

  bornEmitterMomentum() = emitter+emission-(y/(1.-y))*spectator;
  bornSpectatorMomentum() = spectator/(1.-y);

  bornEmitterMomentum().setMass(ZERO);
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(ZERO);
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FFLightTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return scale * sqrt(y*z*(1.-z));

}

double FFLightTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFLightTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FFLightTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FFLightTildeKinematics::Init() {

  static ClassDocumentation<FFLightTildeKinematics> documentation
    ("FFLightTildeKinematics implements the 'tilde' kinematics for "
     "a final-final subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFLightTildeKinematics,TildeKinematics>
describeHerwigFFLightTildeKinematics("Herwig::FFLightTildeKinematics", "Herwig.so");
