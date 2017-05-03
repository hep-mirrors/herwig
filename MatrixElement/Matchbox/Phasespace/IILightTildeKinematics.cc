// -*- C++ -*-
//
// IILightTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
  return scale * sqrt(v*(1.-x-v)/x);

}


Energy IILightTildeKinematics::lastPt(Lorentz5Momentum ,Lorentz5Momentum emission,Lorentz5Momentum )const {
  return emission.perp();
}

pair<double,double> IILightTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  double root = (1.-emitterX())*sqrt(1.-sqr(pt/hardPt));
  return make_pair(0.5*( 1.+emitterX() - root),0.5*( 1.+emitterX() + root));
}


double IILightTildeKinematics::lastZ() const {
  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  return x + v;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IILightTildeKinematics::persistentOutput(PersistentOStream & os) const {
  os << ounit(K,GeV) << ounit(Ktilde,GeV);
}

void IILightTildeKinematics::persistentInput(PersistentIStream & is, int) {
  is >> iunit(K,GeV) >> iunit(Ktilde,GeV);
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
describeHerwigIILightTildeKinematics("Herwig::IILightTildeKinematics", "Herwig.so");
