// -*- C++ -*-
//
// IFLightTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFLightTildeKinematics class.
//

#include "IFLightTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFLightTildeKinematics::IFLightTildeKinematics() {}

IFLightTildeKinematics::~IFLightTildeKinematics() {}

IBPtr IFLightTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFLightTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IFLightTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double x = 
    (- emission*spectator + emitter*spectator + emitter*emission) / 
    (emitter*emission + emitter*spectator);
  double u = emitter*emission / (emitter*emission + emitter*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = u;

  bornEmitterMomentum() = x*emitter;
  bornSpectatorMomentum() = spectator + emission - (1.-x)*emitter;

  bornEmitterMomentum().setMass(ZERO);
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(ZERO);
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy IFLightTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return scale * sqrt(u*(1.-u)*(1.-x)/x);

}

double IFLightTildeKinematics::lastZ() const {
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return 1. - (1.-x)*(1.-u);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFLightTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void IFLightTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void IFLightTildeKinematics::Init() {

  static ClassDocumentation<IFLightTildeKinematics> documentation
    ("IFLightTildeKinematics implements the 'tilde' kinematics for "
     "a initial-final subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFLightTildeKinematics,TildeKinematics>
describeHerwigIFLightTildeKinematics("Herwig::IFLightTildeKinematics", "Herwig.so");
