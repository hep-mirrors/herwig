// -*- C++ -*-
//
// IILightInvertedTildeKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IILightInvertedTildeKinematics class.
//

#include "IILightInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IILightInvertedTildeKinematics::IILightInvertedTildeKinematics() {}

IILightInvertedTildeKinematics::~IILightInvertedTildeKinematics() {}

IBPtr IILightInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IILightInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IILightInvertedTildeKinematics::doMap(const double * r) {

  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();

  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  double z = ptz.second;

  double ratio = sqr(pt/lastScale());
  double x = z*(1.-z) / ( 1. - z + ratio );
  double v = ratio * z / ( 1. - z + ratio );

  if ( x < emitterX() || x > 1. ||
       v < 0. || v > 1.-x ) {
    jacobian(0.0);
    return false;
  }

  mapping /= z*(1.-z);
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = v;

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt = getKt(emitter,spectator,pt,phi);

  realEmitterMomentum() = (1./x)*emitter;
  realEmissionMomentum() = ((1.-x-v)/x)*emitter+v*spectator+kt;
  realSpectatorMomentum() = spectator;

  K = realEmitterMomentum() + realSpectatorMomentum() - realEmissionMomentum();
  K2 = K.m2();

  Ktilde = emitter + spectator;
  KplusKtilde = K + Ktilde;

  KplusKtilde2 = KplusKtilde.m2();

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy IILightInvertedTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  return scale * sqrt(v*(1.-x-v)/x);

}

Energy IILightInvertedTildeKinematics::ptMax() const {
  double tau = emitterX()*spectatorX();
  return (1.-tau)*lastScale()/(2.*sqrt(tau));
}

pair<double,double> IILightInvertedTildeKinematics::zBounds(Energy pt) const {
  double tau = emitterX()*spectatorX();
  double s = sqrt(1.-sqr(pt/ptMax()));
  return make_pair(0.5*(1.+tau-(1.-tau)*s),0.5*(1.+tau+(1.-tau)*s));
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IILightInvertedTildeKinematics::persistentOutput(PersistentOStream & os) const {
  os << ounit(K,GeV) << ounit(K2,GeV2) << ounit(Ktilde,GeV)
     << ounit(KplusKtilde,GeV) << ounit(KplusKtilde2,GeV2);
}

void IILightInvertedTildeKinematics::persistentInput(PersistentIStream & is, int) {
  is >> iunit(K,GeV) >> iunit(K2,GeV2) >> iunit(Ktilde,GeV)
     >> iunit(KplusKtilde,GeV) >> iunit(KplusKtilde2,GeV2);
}

void IILightInvertedTildeKinematics::Init() {

  static ClassDocumentation<IILightInvertedTildeKinematics> documentation
    ("IILightInvertedTildeKinematics inverts the initial-initial tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IILightInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigIILightInvertedTildeKinematics("Herwig::IILightInvertedTildeKinematics", "HwMatchbox.so");
