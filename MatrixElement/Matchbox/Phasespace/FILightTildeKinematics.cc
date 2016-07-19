// -*- C++ -*-
//
// FILightTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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

double FILightTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}

double FILightTildeKinematics::lastRealR() const {
  double deta2 = sqr(realEmitterMomentum().eta() - realEmissionMomentum().eta());
  double dphi =  abs(realEmitterMomentum().phi() - realEmissionMomentum().phi());
  if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
  double dr = sqrt(deta2 + sqr(dphi));
  return  dr;
}


double FILightTildeKinematics::lastBornR() const {
  double deta2 = sqr(bornEmitterMomentum().eta() - bornSpectatorMomentum().eta());
  double dphi =0.;//pi??  abs(bornEmitterMomentum().phi() - bornSpectatorMomentum().phi());
  if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
  double dr = sqrt(deta2 + sqr(dphi));
  return  dr;
}



double FILightTildeKinematics::jacobian(Energy2 sB,Energy2 sR, int n) const {
  
  
  assert(false);
  
  return 16.*ThePEG::Constants::pi*ThePEG::Constants::pi/(2.*bornEmitterMomentum()*realSpectatorMomentum())*sB;// *pow(sR/sB,n-4);
  
  return 16.*ThePEG::Constants::pi*ThePEG::Constants::pi/(2.*realEmitterMomentum()*realSpectatorMomentum())*sR*pow(sR/sB,n-4);

  

  return 16.*ThePEG::Constants::pi*ThePEG::Constants::pi*4.*(1-subtractionParameters()[0]);
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
