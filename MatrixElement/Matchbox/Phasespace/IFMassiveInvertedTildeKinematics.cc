// -*- C++ -*-
//
// IFMassiveInvertedTildeKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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

  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  double z = ptz.second;

  double ratio = sqr(pt/lastScale());
  double alpha = 1. - 2.*sqr(bornSpectatorData()->mass()/lastScale());
  double x = alpha == 1. ? ( z*(1.-z) - ratio ) / ( 1. - z - ratio ) :
    ( sqr(alpha)*ratio + 2.*z - alpha*(1.+z) +
      alpha*sqrt( sqr(1.-z+alpha*ratio) - 4.*ratio*(1.-z) ) ) /
    (2.*(1.-alpha));
  double u = ( 1.-z + alpha*ratio -
	       sqrt( sqr(1.-z+alpha*ratio) - 4.*ratio*(1.-z) ) ) /
    (2.*(1.-z));
  // double x = ( z*(1.-z) - ratio ) / ( 1. - z - ratio );
  // double u = ratio/(1.-z);
  double up = (1.-x) /
    ( 1.-x + x*sqr(bornSpectatorData()->mass()/lastScale()) );

  pt = lastScale()*sqrt(u*(1.-u)*(1.-x)/x);

  if ( x < emitterX() || x > 1. || u > up ) {
    jacobian(0.0);
    return false;
  }

  // TODO: why not mapping /= sqr(z*(1.-z)-ratio)/(1.-z) ? (see phd thesis, (5.74))
  //  mapping /= sqr(z*(1.-z)-ratio)/(1.-z-ratio);
  mapping *= (1.-u)/(1.-2.*u+u*u*alpha)/x;
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt = getKt(emitter,spectator,pt,phi,true);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = u;

  realEmitterMomentum() = (1./x)*emitter;
  //  realEmissionMomentum() = ((1.-x)*(1.-u)/x)*emitter + u*spectator + kt;
  //  realSpectatorMomentum() = ((1.-x)*u/x)*emitter + (1.-u)*spectator - kt;
  realEmissionMomentum() = (pt*pt-u*u*sqr(bornSpectatorData()->mass()))/(u*sqr(lastScale()))*emitter +
      u*spectator + kt;
  realSpectatorMomentum() = (pt*pt+sqr(bornSpectatorData()->mass())-sqr(1.-u)*sqr(bornSpectatorData()->mass()))/((1.-u)*sqr(lastScale()))*emitter +
      (1.-u)*spectator - kt;

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(bornSpectatorData()->mass());
  realSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy IFMassiveInvertedTildeKinematics::lastPt() const {

  Energy scale = (-bornEmitterMomentum()+bornSpectatorMomentum()).m();
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  // TODO: can't make sense of the following comment
  // there was no factor 1/x in massless case >> check
  //  return scale * sqrt(u*(1.-u)*(1.-x)/x);
  return scale * sqrt(u*(1.-u)*(1.-x));

}

Energy IFMassiveInvertedTildeKinematics::ptMax() const {
  double x = emitterX();
  return sqrt(1.-x)*lastScale()/2.;
}

pair<double,double> IFMassiveInvertedTildeKinematics::zBounds(Energy pt) const {
  double alpha = 1. - 2.*sqr(bornSpectatorData()->mass()/lastScale());
  double xe = emitterX();
  double zp = 0.5*( alpha + xe - (alpha-1.)*xe +
		    alpha*(1.-xe)*sqrt(1.-sqr(pt/ptMax()) ) );
  double zm = 0.5*( alpha + xe - (alpha-1.)*xe +
		    alpha*(1.-xe)*sqrt(1.-sqr(pt/ptMax()) ) );
  return make_pair(zm,zp);
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
