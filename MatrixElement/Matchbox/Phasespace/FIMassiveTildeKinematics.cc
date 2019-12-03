// -*- C++ -*-
//
// FIMassiveTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMassiveTildeKinematics class.
//

#include "FIMassiveTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMassiveTildeKinematics::FIMassiveTildeKinematics() {}

FIMassiveTildeKinematics::~FIMassiveTildeKinematics() {}

IBPtr FIMassiveTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FIMassiveTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FIMassiveTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());
  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());

  double x = 
    (- emission*emitter + emission*spectator + emitter*spectator +
     0.5*(Mi2-mi2-m2)) / 
    (emitter*spectator + emission*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = z;

  bornEmitterMomentum() = emitter+emission-(1.-x)*spectator;
  bornSpectatorMomentum() = x*spectator;

  bornEmitterMomentum().setMass(bornEmitterData()->hardProcessMass());
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(bornSpectatorData()->hardProcessMass());
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FIMassiveTildeKinematics::lastPt() const {

  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());
  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());

  Energy2 scale = Mi2 - (realEmitterMomentum()+realEmissionMomentum()-realSpectatorMomentum()).m2();
  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return sqrt( z*(1.-z)*(1.-x)/x*scale -
	       ((1.-z)*mi2+z*m2-z*(1.-z)*Mi2) );

}



Energy FIMassiveTildeKinematics::lastPt(Lorentz5Momentum emitter,Lorentz5Momentum emission,Lorentz5Momentum spectator)const {
    // g->QQ or Q -> Qg
  Energy2 Mi2 = emitter.m()==emission.m()?0.*GeV2:max(emitter.m2(),emission.m2());
  Energy2 mi2 = emitter.m2();
  Energy2 m2  = emission.m2();
  
  Energy2 scale = Mi2 - (emitter+emission-spectator).m2();
  
  double x =
  (- emission*emitter + emission*spectator + emitter*spectator +
   0.5*(Mi2-mi2-m2)) /
  (emitter*spectator + emission*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);
  
  return sqrt( z*(1.-z)*(1.-x)/x*scale -
              ((1.-z)*mi2+z*m2-z*(1.-z)*Mi2) );
}


pair<double,double> FIMassiveTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());
  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());
    // s^star/x
  Energy2 scale=2.*bornEmitterMomentum()*bornSpectatorMomentum();
  Energy2 s = scale * (1.-spectatorX())/spectatorX() +  Mi2;
  
  double zm = .5*( 1.+(mi2-m2)/s - rootOfKallen(s/s,mi2/s,m2/s) *
                  sqrt( 1.-sqr(pt/hardPt) ) );
  double zp = .5*( 1.+(mi2-m2)/s + rootOfKallen(s/s,mi2/s,m2/s) *
                  sqrt( 1.-sqr(pt/hardPt) ) );
  return make_pair(zm, zp);
  
}



double FIMassiveTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMassiveTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FIMassiveTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FIMassiveTildeKinematics::Init() {

  static ClassDocumentation<FIMassiveTildeKinematics> documentation
    ("FIMassiveTildeKinematics implements the 'tilde' kinematics for "
     "a final-initial subtraction dipole involving a massive particle.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMassiveTildeKinematics,TildeKinematics>
describeHerwigFIMassiveTildeKinematics("Herwig::FIMassiveTildeKinematics", "Herwig.so");
