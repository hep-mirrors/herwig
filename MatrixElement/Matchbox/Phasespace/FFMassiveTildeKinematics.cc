// -*- C++ -*-
//
// FFMassiveTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMassiveTildeKinematics class.
//

#include "FFMassiveTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FFMassiveTildeKinematics::FFMassiveTildeKinematics() {}

FFMassiveTildeKinematics::~FFMassiveTildeKinematics() {}

IBPtr FFMassiveTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFMassiveTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FFMassiveTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double y = emission*emitter / (emission*emitter + emission*spectator + emitter*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = z;
  
  Lorentz5Momentum pTot = emitter+emission+spectator;
  Energy scale = pTot.m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  double Mui2 = sqr( bornEmitterData()->hardProcessMass() / scale );
  double Muj2 = sqr( bornSpectatorData()->hardProcessMass() / scale );
  
  double bar = 1.-mui2-mu2-muj2;
  
  // from Catani,Seymour,Dittmaier,Trocsanyi (CSm!=0) (5.9)
  // has the right massless limit
  bornSpectatorMomentum() =
    rootOfKallen( 1., Mui2, Muj2 ) / rootOfKallen( 1., mui2+mu2+y*bar, muj2 ) *
    ( spectator - (pTot*spectator)/sqr(scale)*pTot ) +
    0.5*( 1.+Muj2-Mui2 ) * pTot;
  bornEmitterMomentum() = pTot - bornSpectatorMomentum();

  bornEmitterMomentum().setMass( sqrt(Mui2)*scale );
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass( sqrt(Muj2)*scale );
  bornSpectatorMomentum().rescaleEnergy();
  
  return true;

}

Energy FFMassiveTildeKinematics::lastPt() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double mu2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muj2 = sqr( realSpectatorData()->hardProcessMass() / scale );
    
  Energy ret = scale * sqrt( y * (1.-mui2-mu2-muj2) * z*(1.-z) - sqr(1.-z)*mui2 - sqr(z)*mu2 );
  
  return ret;
}

double FFMassiveTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMassiveTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FFMassiveTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FFMassiveTildeKinematics::Init() {

  static ClassDocumentation<FFMassiveTildeKinematics> documentation
    ("FFMassiveTildeKinematics implements the 'tilde' kinematics for "
     "a final-final subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMassiveTildeKinematics,TildeKinematics>
describeHerwigFFMassiveTildeKinematics("Herwig::FFMassiveTildeKinematics", "Herwig.so");
