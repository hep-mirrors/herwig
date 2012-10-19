// -*- C++ -*-
//
// IFMassiveTildeKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMassiveTildeKinematics class.
//

#include "IFMassiveTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMassiveTildeKinematics::IFMassiveTildeKinematics() {}

IFMassiveTildeKinematics::~IFMassiveTildeKinematics() {}

IBPtr IFMassiveTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFMassiveTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IFMassiveTildeKinematics::doMap() {

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
  bornSpectatorMomentum().setMass(bornSpectatorData()->mass());
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy IFMassiveTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(realEmissionMomentum()*realEmitterMomentum()-realEmissionMomentum()*realSpectatorMomentum()+realEmitterMomentum()*realSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return scale * sqrt(u*(1.-u)*(1.-x));

}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMassiveTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void IFMassiveTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void IFMassiveTildeKinematics::Init() {

  static ClassDocumentation<IFMassiveTildeKinematics> documentation
    ("IFMassiveTildeKinematics implements the 'tilde' kinematics for "
     "a initial-final subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMassiveTildeKinematics,TildeKinematics>
describeHerwigIFMassiveTildeKinematics("Herwig::IFMassiveTildeKinematics", "HwMatchbox.so");
