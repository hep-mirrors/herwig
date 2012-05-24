// -*- C++ -*-
//
// FFLightInvertedTildeKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFLightInvertedTildeKinematics class.
//

#include "FFLightInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFLightInvertedTildeKinematics::FFLightInvertedTildeKinematics() {}

FFLightInvertedTildeKinematics::~FFLightInvertedTildeKinematics() {}

IBPtr FFLightInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFLightInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FFLightInvertedTildeKinematics::doMap(const double * r) {

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

  double y = sqr(pt/lastScale())/(z*(1.-z));

  mapping /= z*(1.-z);
  jacobian(mapping*(1.-y)*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt
    = getKt(emitter,spectator,pt,phi);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = z;

  realEmitterMomentum() = z*emitter + y*(1.-z)*spectator + kt;
  realEmissionMomentum() = (1.-z)*emitter + y*z*spectator - kt;
  realSpectatorMomentum() = (1.-y)*spectator;

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FFLightInvertedTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return scale * sqrt(y*z*(1.-z));

}

Energy FFLightInvertedTildeKinematics::ptMax() const {
  return lastScale()/2.;
}

pair<double,double> FFLightInvertedTildeKinematics::zBounds(Energy pt) const {
  double s = sqrt(1.-sqr(pt/ptMax()));
  return make_pair(0.5*(1.-s),0.5*(1.+s));
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFLightInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FFLightInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FFLightInvertedTildeKinematics::Init() {

  static ClassDocumentation<FFLightInvertedTildeKinematics> documentation
    ("FFLightInvertedTildeKinematics inverts the final-final tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFLightInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFFLightInvertedTildeKinematics("Herwig::FFLightInvertedTildeKinematics", "HwMatchbox.so");
