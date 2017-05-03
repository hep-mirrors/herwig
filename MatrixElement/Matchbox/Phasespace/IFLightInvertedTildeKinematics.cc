// -*- C++ -*-
//
// IFLightInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFLightInvertedTildeKinematics class.
//

#include "IFLightInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IFLightInvertedTildeKinematics::IFLightInvertedTildeKinematics() {}

IFLightInvertedTildeKinematics::~IFLightInvertedTildeKinematics() {}

IBPtr IFLightInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFLightInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IFLightInvertedTildeKinematics::doMap(const double * r) {

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
  double rho = 1. - 4.*ratio*z*(1.-z) / sqr(1. - z + ratio);
  if ( rho < 0. ) {
    jacobian(0.0);
    return false;
  }

  double x = 0.5*(1./ratio)*(1.-z+ratio)*(1.-sqrt(rho));
  double u = 0.5*(1./(1.-z))*(1.-z+ratio)*(1.-sqrt(rho));

  if ( x < emitterX() || x > 1. || 
       u < 0. || u > 1. ) {
    jacobian(0.0);
    return false;
  }

  // This jacobian is (1/x^2)*dx*du
  mapping *= (1.-x)/((1.-z)*(z*(1.-z)+sqr(x-z)));
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt = getKt(emitter,spectator,pt,phi,true);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = u;

  realEmitterMomentum() = (1./x)*emitter;
  realEmissionMomentum() = ((1.-x)*(1.-u)/x)*emitter + u*spectator + kt;
  realSpectatorMomentum() = ((1.-x)*u/x)*emitter + (1.-u)*spectator - kt;

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy IFLightInvertedTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return scale * sqrt(u*(1.-u)*(1.-x)/x);

}

Energy IFLightInvertedTildeKinematics::ptMax() const {
  double x = emitterX();
  return sqrt((1.-x)/x)*lastScale()/2.;
}

pair<double,double> IFLightInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  double x = emitterX();
  return make_pair(0.5*(1.+x-(1.-x)*s),0.5*(1.+x+(1.-x)*s));
}

double IFLightInvertedTildeKinematics::lastZ() const {
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return 1. - (1.-x)*(1.-u);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFLightInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void IFLightInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void IFLightInvertedTildeKinematics::Init() {

  static ClassDocumentation<IFLightInvertedTildeKinematics> documentation
    ("IFLightInvertedTildeKinematics inverts the initial-final tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFLightInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigIFLightInvertedTildeKinematics("Herwig::IFLightInvertedTildeKinematics", "Herwig.so");
