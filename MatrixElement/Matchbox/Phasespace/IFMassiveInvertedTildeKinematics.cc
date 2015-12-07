// -*- C++ -*-
//
// IFMassiveInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMassiveInvertedTildeKinematics class.
//

#include "IFMassiveInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFMassiveInvertedTildeKinematics::IFMassiveInvertedTildeKinematics() {}

IFMassiveInvertedTildeKinematics::~IFMassiveInvertedTildeKinematics() {}

IBPtr IFMassiveInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFMassiveInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IFMassiveInvertedTildeKinematics::doMap(const double * r) { 
  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();
  Energy2 scale = 2.*(spectator*emitter);

  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ){
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  double z = ptz.second;
  double ratio = sqr(pt)/scale;

  // x
  double u = ratio/(1.-z);
  double x = (z*(1.-z)-ratio)/(1.-z-ratio);
  double up = (1.-x) / (1.-x+(x*sqr(bornSpectatorData()->hardProcessMass())/scale));

  if ( x < emitterX() || x > 1. || u > up) {
    jacobian(0.0);
    return false;
  }

  pt = sqrt(scale*u*(1.-u)*(1.-x));  
  Energy magKt = sqrt(scale*u*(1.-u)*(1.-x)/x - sqr(u*bornSpectatorData()->hardProcessMass()));

  // TODO: why not mapping /= sqr(z*(1.-z)-ratio)/(1.-z) ? (see phd thesis, (5.74))
  mapping /= sqr(z*(1.-z)-ratio)/(1.-z-ratio);
  // mapping *= (1.-u)/(1.-2.*u+u*u*alpha)/x;
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt = getKt(emitter,spectator,magKt,phi,true);
  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = u;
  
  realEmitterMomentum() = (1./x)*emitter;
  realEmissionMomentum() = (-kt*kt-u*u*sqr(bornSpectatorData()->hardProcessMass()))/(u*scale)*emitter +  
    u*spectator - kt;
  realSpectatorMomentum() = (-kt*kt+sqr(bornSpectatorData()->hardProcessMass())*u*(2.-u))/((1.-u)*scale)*emitter + 
    (1.-u)*spectator + kt;

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(bornSpectatorData()->hardProcessMass());
  realSpectatorMomentum().rescaleEnergy();
  return true;

}

Energy IFMassiveInvertedTildeKinematics::lastPt() const {
  Energy2 scale = 2.*(bornEmitterMomentum()*bornSpectatorMomentum());
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  // TODO: can't make sense of the following comment
  // there was no factor 1/x in massless case >> check
  //  return scale * sqrt(u*(1.-u)*(1.-x)/x);
  return sqrt(scale*u*(1.-u)*(1.-x));
}

double IFMassiveInvertedTildeKinematics::lastZ() const {
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return 1. - (1.-x)*(1.-u);
}

Energy IFMassiveInvertedTildeKinematics::ptMax() const {
  Energy2 scale = 2.*(bornEmitterMomentum()*bornSpectatorMomentum());
  double x = emitterX();
  return sqrt(scale*(1.-x))/2.;
}

pair<double,double> IFMassiveInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  double x = emitterX();
  return make_pair(0.5*(1.+x-(1.-x)*s),0.5*(1.+x+(1.-x)*s));
}



// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMassiveInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void IFMassiveInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void IFMassiveInvertedTildeKinematics::Init() {

  static ClassDocumentation<IFMassiveInvertedTildeKinematics> documentation
    ("IFMassiveInvertedTildeKinematics inverts the initial-final tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMassiveInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigIFMassiveInvertedTildeKinematics("Herwig::IFMassiveInvertedTildeKinematics", "Herwig.so");
