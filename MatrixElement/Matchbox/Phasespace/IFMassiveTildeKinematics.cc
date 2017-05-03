// -*- C++ -*-
//
// IFMassiveTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
  bornSpectatorMomentum().setMass(bornSpectatorData()->hardProcessMass());
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy IFMassiveTildeKinematics::lastPt() const {
  Energy2 scale = 2.*(realEmissionMomentum()*realEmitterMomentum()
-realEmissionMomentum()*realSpectatorMomentum()
+realEmitterMomentum()*realSpectatorMomentum());
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    return sqrt(scale * ( u*(1.-u)*(1.-x)/x - u*u*muk2 ));
   }

Energy IFMassiveTildeKinematics::lastPt(Lorentz5Momentum emitter,Lorentz5Momentum emission,Lorentz5Momentum spectator)const {
  Energy2 scale = 2.*(emission*emitter-emission*spectator+emitter*spectator);
  double x = 0.5*scale / (emitter*emission + emitter*spectator);
  double u = emitter*emission / (emitter*emission + emitter*spectator);
  
    double muk2 = sqr(spectator.mass())/scale;
    return sqrt(scale * ( u*(1.-u)*(1.-x)/x - u*u*muk2 ));
  }
  
pair<double,double> IFMassiveTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  double xe = emitterX();
  return make_pair(0.5*(1.+xe-(1.-xe)*s),0.5*(1.+xe+(1.-xe)*s));
}

double IFMassiveTildeKinematics::lastZ() const {
Energy2 scale = 2.*(realEmissionMomentum()*realEmitterMomentum()
-realEmissionMomentum()*realSpectatorMomentum()
+realEmitterMomentum()*realSpectatorMomentum());
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    return u + x - u*x*(1.-muk2);
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
     "a initial-final subtraction dipole involving a massive particle.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMassiveTildeKinematics,TildeKinematics>
describeHerwigIFMassiveTildeKinematics("Herwig::IFMassiveTildeKinematics", "Herwig.so");
