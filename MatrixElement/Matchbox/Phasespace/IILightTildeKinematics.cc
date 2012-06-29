// -*- C++ -*-
//
// IILightTildeKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IILightTildeKinematics class.
//

#include "IILightTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IILightTildeKinematics::IILightTildeKinematics() {}

IILightTildeKinematics::~IILightTildeKinematics() {}

IBPtr IILightTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IILightTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

Lorentz5Momentum IILightTildeKinematics::transform(const Lorentz5Momentum& k) const {

  LorentzMomentum res =
    k - 2.*((k*(K+Ktilde)/(K+Ktilde).m2())*(K+Ktilde)-((k*K)/(K.m2()))*Ktilde);

  return res;

}

bool IILightTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double x = (emitter*spectator - emitter*emission - spectator*emission)/(emitter*spectator);
  double v = (emitter*emission)/(emitter*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = v;

  bornEmitterMomentum() = x * emitter;
  bornSpectatorMomentum() = spectator;

  bornEmitterMomentum().setMass(ZERO);
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(ZERO);
  bornSpectatorMomentum().rescaleEnergy();

  K = emitter + spectator - emission;
  Ktilde = x * emitter + spectator;

  return true;

}

Energy IILightTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  return scale * sqrt(v*(1.-x-v));

}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IILightTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void IILightTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void IILightTildeKinematics::Init() {

  static ClassDocumentation<IILightTildeKinematics> documentation
    ("IILightTildeKinematics implements the 'tilde' kinematics for "
     "a initial-initial subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IILightTildeKinematics,TildeKinematics>
describeHerwigIILightTildeKinematics("Herwig::IILightTildeKinematics", "HwMatchbox.so");
