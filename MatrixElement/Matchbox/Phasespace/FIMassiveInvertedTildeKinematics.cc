// -*- C++ -*-
//
// FIMassiveInvertedTildeKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMassiveInvertedTildeKinematics class.
//

#include "FIMassiveInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FIMassiveInvertedTildeKinematics::FIMassiveInvertedTildeKinematics() {}

FIMassiveInvertedTildeKinematics::~FIMassiveInvertedTildeKinematics() {}

IBPtr FIMassiveInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FIMassiveInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FIMassiveInvertedTildeKinematics::doMap(const double * r) {

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

  Energy2 mi2 = sqr(realEmitterData()->mass());
  Energy2 m2  = sqr(realEmissionData()->mass());
  Energy2 Mi2 = sqr(bornEmitterData()->mass());

  //  double y = sqr(pt/lastScale())/(z*(1.-z));
  double y = (pt*pt+(1.-z)*mi2+z*m2-z*(1.-z)*Mi2) / 
    (z*(1.-z)*sqr(lastScale()));
  double x = 1./(1.+y);

  if ( x < spectatorX() ) {
    jacobian(0.0);
    return false;
  }

  // no additional massive factors
  mapping /= z*(1.-z);
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt
    = getKt(spectator,emitter,pt,phi,true);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = z;

  realEmitterMomentum() = z*emitter +
    (sqr(pt)+mi2-z*z*Mi2)/(z*sqr(lastScale()))*spectator + kt;
  realEmissionMomentum() = (1.-z)*emitter +
    (pt*pt+m2-sqr(1.-z)*Mi2)/((1.-z)*sqr(lastScale()))*spectator - kt;
  realSpectatorMomentum() = (1.+y)*spectator;

  realEmitterMomentum().setMass(sqrt(mi2));
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(sqrt(m2));
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FIMassiveInvertedTildeKinematics::lastPt() const {

  Energy2 mi2 = sqr(realEmitterData()->mass());
  Energy2 m2  = sqr(realEmissionData()->mass());
  Energy2 Mi2 = sqr(bornEmitterData()->mass());

  Energy2 scale = Mi2 - (realEmitterMomentum()+realEmissionMomentum()-realSpectatorMomentum()).m2();
  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return sqrt( z*(1.-z)*(1.-x)/x*scale -
	       ((1.-z)*mi2+z*m2-z*(1.-z)*Mi2) );

}

Energy FIMassiveInvertedTildeKinematics::ptMax() const {
  Energy2 mi2 = sqr(realEmitterData()->mass());
  Energy2 m2  = sqr(realEmissionData()->mass());
  Energy2 Mi2 = sqr(bornEmitterData()->mass());

  double x = spectatorX();
  // s^star/x
  Energy2 s = sqr(lastScale()) * (1.-x)/x + Mi2;

  Energy ptmax = .5 * sqrt(s) * rootOfKallen( s/s, mi2/s, m2/s );
  return ptmax > 0.*GeV ? ptmax : 0.*GeV;
}

pair<double,double> FIMassiveInvertedTildeKinematics::zBounds(Energy pt) const {
  Energy2 mi2 = sqr(realEmitterData()->mass());
  Energy2 m2  = sqr(realEmissionData()->mass());
  Energy2 Mi2 = sqr(bornEmitterData()->mass());
  double x = subtractionParameters()[0];
  double mui2 = x*mi2/sqr(lastScale());
  double mu2  = x*m2/sqr(lastScale());
  double Mui2 = x*Mi2/sqr(lastScale());

  // s^star/x
  Energy2 s = sqr(lastScale()) * (1.-spectatorX())/spectatorX() +
    Mi2;
  double zm1 = .5*( 1.+(mi2-m2)/s - rootOfKallen(s/s,mi2/s,m2/s) *
		    sqrt( 1.-sqr(pt/ptMax()) ) );
  double zp1 = .5*( 1.+(mi2-m2)/s + rootOfKallen(s/s,mi2/s,m2/s) *
		    sqrt( 1.-sqr(pt/ptMax()) ) );
  double zm = .5*( 1.-x+Mui2+mui2-mui2 -
		   sqrt( sqr(1.-x+Mui2-mui2-mu2)-4.*mui2*mu2 ) ) /
    (1.-x+Mui2);
  double zp = .5*( 1.-x+Mui2+mui2-mui2 +
		   sqrt( sqr(1.-x+Mui2-mui2-mu2)-4.*mui2*mu2 ) ) /
    (1.-x+Mui2);

  return make_pair( zm>zm1 ? zm : zm1 , zp<zp1 ? zp : zp1 );
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMassiveInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FIMassiveInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FIMassiveInvertedTildeKinematics::Init() {

  static ClassDocumentation<FIMassiveInvertedTildeKinematics> documentation
    ("FIMassiveInvertedTildeKinematics inverts the final-initial tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMassiveInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFIMassiveInvertedTildeKinematics("Herwig::FIMassiveInvertedTildeKinematics", "HwMatchbox.so");
