// -*- C++ -*-
//
// FILightTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FILightTildeKinematics class.
//

#include "FILightTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FILightTildeKinematics::FILightTildeKinematics() {}

FILightTildeKinematics::~FILightTildeKinematics() {}

IBPtr FILightTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FILightTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FILightTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double x = 
    (- emission*emitter + emission*spectator + emitter*spectator) / 
    (emitter*spectator + emission*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = z;

  bornEmitterMomentum() = emitter+emission-(1.-x)*spectator;
  bornSpectatorMomentum() = x*spectator;

  bornEmitterMomentum().setMass(ZERO);
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(ZERO);
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FILightTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return scale * sqrt(z*(1.-z)*(1.-x)/x);

}

Energy FILightTildeKinematics::lastPt(Lorentz5Momentum emitter,Lorentz5Momentum emission,Lorentz5Momentum spectator)const {

  
  Energy2 scale =  - (emitter+emission-spectator).m2();
  
  
  double x =
  (- emission*emitter + emission*spectator + emitter*spectator) /
  (emitter*spectator + emission*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);
  
  return sqrt( z*(1.-z)*(1.-x)/x*scale ) ;
}

pair<double,double> FILightTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  return make_pair(0.5*(1.-s),0.5*(1.+s));
}



double FILightTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FILightTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FILightTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FILightTildeKinematics::Init() {

  static ClassDocumentation<FILightTildeKinematics> documentation
    ("FILightTildeKinematics implements the 'tilde' kinematics for "
     "a final-initial subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FILightTildeKinematics,TildeKinematics>
describeHerwigFILightTildeKinematics("Herwig::FILightTildeKinematics", "Herwig.so");
