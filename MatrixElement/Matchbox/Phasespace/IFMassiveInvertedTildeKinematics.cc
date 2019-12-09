// -*- C++ -*-
//
// IFMassiveInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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

  // Compute dipole scale
  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();
  Energy2 scale = 2.*(spectator*emitter);

  // Generate pt and z
  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ){
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  double z = ptz.second;
 
    // Compute x and u
    double ratio = sqr(pt)/scale;
    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);

    double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));
    double u = x*ratio / (1.-z);
      
    // Following Catani-Seymour paper
    double muk2CS = x*muk2;
    double up = (1.-x) /
      ( 1.-x + muk2CS );
      
    if ( x < emitterX() || x > 1. ||
	 u < 0. || u > up ) {
      jacobian(0.0);
      return false;
    }
   

    // Store x and u
    subtractionParameters().resize(2);
    subtractionParameters()[0] = x;
    subtractionParameters()[1] = u;
    
    // jac = sajk*(1./x^2)*dx*du
    // Note - lastScale() is not equal to scale!!!!!!!
    double jac = u/x/(u + x - 2.*u*x*(1.-muk2))*scale/sqr(pt);
    mapping *= jac;
    jacobian( mapping*(sqr(lastScale())/sHat()) / (16.*sqr(Constants::pi)) );

    // Compute the new momenta
    double phi = 2.*Constants::pi*r[2];
    Lorentz5Momentum kt = getKt(emitter,spectator,pt,phi,true);
    
    realEmitterMomentum() = (1./x)*emitter;
    realEmissionMomentum() = ((1.-x)*(1.-u)/x - 2.*u*muk2)*emitter + u*spectator + kt;
    realSpectatorMomentum() = ((1.-x)*u/x + 2.*u*muk2)*emitter + (1.-u)*spectator - kt;

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
    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    double x = subtractionParameters()[0];
    double u = subtractionParameters()[1];
    return sqrt(scale * ( u*(1.-u)*(1.-x)/x - u*u*muk2 ));
}

double IFMassiveInvertedTildeKinematics::lastZ() const {
    Energy2 scale = 2.*(bornEmitterMomentum()*bornSpectatorMomentum());
    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    double x = subtractionParameters()[0];
    double u = subtractionParameters()[1];  
    return u + x + u*x*(muk2-1.);
}

Energy IFMassiveInvertedTildeKinematics::ptMax() const {
    double xe = emitterX();
    Energy2 scale = 2.*(bornEmitterMomentum()*bornSpectatorMomentum());
    Energy2 A = scale*(1.-xe)/xe;
    Energy2 mk2 = sqr(bornSpectatorData()->hardProcessMass());
    return 0.5*A/sqrt(mk2+A);
}

pair<double,double> IFMassiveInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  double xe = emitterX();
  return make_pair(0.5*(1.+xe-(1.-xe)*s),0.5*(1.+xe+(1.-xe)*s));
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
     "kinematics involving a massive particle.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMassiveInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigIFMassiveInvertedTildeKinematics("Herwig::IFMassiveInvertedTildeKinematics", "Herwig.so");
