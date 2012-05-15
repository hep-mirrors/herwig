// -*- C++ -*-
//
// FILightInvertedTildeKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FILightInvertedTildeKinematics class.
//

#include "FILightInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FILightInvertedTildeKinematics::FILightInvertedTildeKinematics() {}

FILightInvertedTildeKinematics::~FILightInvertedTildeKinematics() {}

IBPtr FILightInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FILightInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FILightInvertedTildeKinematics::doMap(const double * r) {

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
  double x = 1./(1.+y);

  if ( x < spectatorX() ) {
    jacobian(0.0);
    return false;
  }

  mapping /= z*(1.-z);
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt
    = getKt(spectator,emitter,pt,phi,true);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = z;

  realEmitterMomentum() = z*emitter + y*(1.-z)*spectator + kt;
  realEmissionMomentum() = (1.-z)*emitter + y*z*spectator - kt;
  realSpectatorMomentum() = (1.+y)*spectator;

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FILightInvertedTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return scale * sqrt(z*(1.-z)*(1.-x)/x);

}

Energy FILightInvertedTildeKinematics::ptMax() const {
  double x = spectatorX();
  return sqrt((1.-x)/x)*lastScale()/2.;
}

pair<double,double> FILightInvertedTildeKinematics::zBounds(Energy pt) const {
  double s = sqrt(1.-sqr(pt/ptMax()));
  return make_pair(0.5*(1.-s),0.5*(1.+s));
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FILightInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FILightInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FILightInvertedTildeKinematics::Init() {

  static ClassDocumentation<FILightInvertedTildeKinematics> documentation
    ("FILightInvertedTildeKinematics inverts the final-initial tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FILightInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFILightInvertedTildeKinematics("Herwig::FILightInvertedTildeKinematics", "HwMatchbox.so");
