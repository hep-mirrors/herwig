// -*- C++ -*-
//
// IFMassiveKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMassiveKinematics class.
//

#include "IFMassiveKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/DipoleShower/Base/DipoleSplittingInfo.h"

using namespace Herwig;

IFMassiveKinematics::IFMassiveKinematics() 
  : DipoleSplittingKinematics(), theCollinearScheme(true) {}

IFMassiveKinematics::~IFMassiveKinematics() {}

IBPtr IFMassiveKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFMassiveKinematics::fullclone() const {
  return new_ptr(*this);
}

pair<double,double> IFMassiveKinematics::kappaSupport(const DipoleSplittingInfo&) const {
  return make_pair(0.0,1.0);
}

pair<double,double> IFMassiveKinematics::xiSupport(const DipoleSplittingInfo& split) const {

  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));

  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emitterData()->id() == ParticleID::g ) {
      double b = log((1.+c)/(1.-c));
      return make_pair(-b,b);
    } else {
      return make_pair(log(0.5*(1.-c)),log(0.5*(1.+c)));
    }
  }

  if ( split.index().emitterData()->id() != ParticleID::g &&
       split.emitterData()->id() != ParticleID::g ) {
    return make_pair(-log(0.5*(1.+c)),-log(0.5*(1.-c)));
  }

  return make_pair(0.5*(1.-c),0.5*(1.+c));

}

// sbar
Energy IFMassiveKinematics::dipoleScale(const Lorentz5Momentum& pEmitter,
				      const Lorentz5Momentum& pSpectator) const {
  return sqrt(2.*(pEmitter*pSpectator));
}

Energy IFMassiveKinematics::ptMax(Energy dScale, 
				double emX, double,
				const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  return dScale * sqrt(1.-emX) /2.;
}

Energy IFMassiveKinematics::QMax(Energy, 
				 double, double,
				 const DipoleIndex&,
				 const DipoleSplittingKernel&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy IFMassiveKinematics::PtFromQ(Energy, const DipoleSplittingInfo&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy IFMassiveKinematics::QFromPt(Energy, const DipoleSplittingInfo&) const {
  assert(false && "add this");
  return 0.0*GeV;
}


double IFMassiveKinematics::ptToRandom(Energy pt, Energy,
				       double,double,
				       const DipoleIndex&,
				       const DipoleSplittingKernel&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

bool IFMassiveKinematics::generateSplitting(double kappa, double xi, double rphi,
					    DipoleSplittingInfo& info,
					    const DipoleSplittingKernel&) {

  if ( info.emitterX() < xMin() ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = IRCutoff() * pow(0.5 * generator()->maximumCMEnergy()/IRCutoff(),kappa);

  if ( sqr(pt) > sqr(info.hardPt())/(1.+4.*sqr(info.hardPt()/info.scale())) ) {
    jacobian(0.0);
    return false;
  }

  double z = 0.;
  double mapZJacobian = 0.;

  if ( info.index().emitterData()->id() == ParticleID::g ) {
    if ( info.emitterData()->id() == ParticleID::g ) {
      z = exp(xi)/(1.+exp(xi));
      mapZJacobian = z*(1.-z);
    } else {
      z = exp(xi);
      mapZJacobian = z;
    }
  }

  if ( info.index().emitterData()->id() != ParticleID::g ) {
    if ( info.emitterData()->id() != ParticleID::g ) {
      z = 1.-exp(-xi);
      mapZJacobian = 1.-z;
    } else {
      z = xi;
      mapZJacobian = 1.;
    }
  }

  double ratio = sqr(pt/info.scale());

  double x = ( z*(1.-z) - ratio ) / ( 1. - z - ratio );
  double u = ratio/(1.-z);

  double up = (1.-x) /
    ( 1.-x + x*sqr(info.spectatorData()->mass()/info.scale()) );

  if ( x < 0. || x > 1. || u > up ) {
    jacobian(0.0);
    return false;
  }

  double xe = info.emitterX();
  double zmx = 0.5*(1.+xe-(1.-xe)*sqrt(1.-sqr(pt/info.hardPt())) );
  double zpx = 0.5*(1.+xe+(1.-xe)*sqrt(1.-sqr(pt/info.hardPt())) );

  double xq = sqr(pt/info.hardPt());
  double zpq = 0.5*( 1.+ xq + (1.-xq)*sqrt(1.-sqr(2.*pt/info.scale())/(1.-xq) ) );
  double zmq = 0.5*( 1.+ xq - (1.-xq)*sqrt(1.-sqr(2.*pt/info.scale())/(1.-xq) ) );

  double zp = min(zpx,zpq);
  double zm = max(zmx,zmq);

  if ( pt < IRCutoff() || 
       pt > info.hardPt() ||
       z > zp || z < zm ||
       x < xe ) {
    jacobian(0.0);
    return false;
  }

  double phi = 2.*Constants::pi*rphi;

  jacobian(2. * mapZJacobian * log(0.5 * generator()->maximumCMEnergy()/IRCutoff()) * (1.-z-ratio)/sqr(z*(1.-z)-ratio) );

  lastPt(pt);
  lastZ(z);
  lastPhi(phi);
  lastEmitterZ(x);

  if ( theMCCheck )
    theMCCheck->book(info.emitterX(),1.,info.scale(),info.hardPt(),pt,z,jacobian());

  return true;

}

void IFMassiveKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {
  Energy2 sbar = 2.*pEmitter*pSpectator;
  Energy pt = dInfo.lastPt();
  double ratio = pt*pt/sbar;
  double z = dInfo.lastZ();
  double x = (z*(1.-z)-ratio)/(1.-z-ratio);
  double u = ratio / (1.-z);

  pt = sqrt(sbar*u*(1.-u)*(1.-x));
  Energy magKt = 
    sqrt(sbar*u*(1.-u)*(1.-x)/x - sqr(u*dInfo.spectatorData()->mass()));
  Lorentz5Momentum kt =
    getKt (pEmitter, pSpectator, magKt, dInfo.lastPhi(),true);

  Lorentz5Momentum em;
  Lorentz5Momentum emm;
  Lorentz5Momentum spe;

  Energy2 mj2 = dInfo.spectatorData()->mass()*dInfo.spectatorData()->mass();
  double alpha = 1. - 2.*mj2/sbar;

  // TODO: adjust phasespace boundary condition
  if (!theCollinearScheme)
    assert(false);
  if ( !theCollinearScheme &&
       x > u && (1.-x)/(x-u) < 1. ) {

    double fkt = sqrt(sqr(x-u)+4.*x*u*mj2/sbar);

    //    em =
    //      ((1.-u)/(x-u))*pEmitter + ((u/x)*(1.-x)/(x-u))*pSpectator - kt/(x-u);
    Energy2 fa = (sbar*(x+u-2.*x*z)+2.*mj2*x*u) / sqrt(sqr(x-u)+4.*x*u*mj2/sbar);
    double a = (-sbar+fa) / (2.*x*(sbar-mj2));
    double ap = (sbar+alpha*fa) / (2.*x*(sbar-mj2));
    em = ap*pEmitter + a*pSpectator - fkt*kt;
    em.setMass(ZERO);
    em.rescaleEnergy();

    //    emm =
    //      ((1.-x)/(x-u))*pEmitter + ((u/x)*(1.-u)/(x-u))*pSpectator - kt/(x-u);
    Energy2 fb = abs(sbar*(u*(1.-u)-x*(1.-x))+2.*mj2*x*u) / sqrt(sqr(x-u)+4.*x*u*mj2/sbar);
    double b = (-sbar*(1.-x-u)+fb) / (2.*x*(sbar-mj2));
    double bp = (sbar*(1.-x-u)+alpha*fb) / (2.*x*(sbar-mj2));
    emm = bp*pEmitter + b*pSpectator + fkt*kt;
    emm.setMass(ZERO);
    emm.rescaleEnergy();

    //    spe =
    //      (1.-u/x)*pSpectator;
    Energy2 fc = sqrt(sqr(sbar*(x-u))+4.*sbar*mj2*x*u);
    double c = (sbar*(x-u)-2.*x*mj2+fc) / (2.*x*(sbar-mj2));
    double cp = (-sbar*(x-u)+2.*x*mj2+alpha*fc) / (2.*x*(sbar-mj2));
    spe = cp*pEmitter + c*pSpectator;
    spe.setMass(dInfo.spectatorData()->mass());
    spe.rescaleEnergy();

  } else {

    em = (1./x)*pEmitter;
    em.setMass(ZERO);
    em.rescaleEnergy();

    emm = (-kt*kt-u*u*mj2)/(u*sbar)*pEmitter +
      u*pSpectator - kt;
    emm.setMass(ZERO);
    emm.rescaleEnergy();

    spe = (-kt*kt + mj2*u*(2.-u))/((1.-u)*sbar)*pEmitter +
      (1.-u)*pSpectator + kt;
    spe.setMass(dInfo.spectatorData()->mass());
    spe.rescaleEnergy();

  }

  emitterMomentum(em);
  emissionMomentum(emm);
  spectatorMomentum(spe);

}

// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMassiveKinematics::persistentOutput(PersistentOStream &) const {
  //os << theCollinearScheme;
}

void IFMassiveKinematics::persistentInput(PersistentIStream &, int) {
  //is >> theCollinearScheme;
}

ClassDescription<IFMassiveKinematics> IFMassiveKinematics::initIFMassiveKinematics;
// Definition of the static class description member.

void IFMassiveKinematics::Init() {

  static ClassDocumentation<IFMassiveKinematics> documentation
    ("IFMassiveKinematics implements massless splittings "
     "off a initial-final dipole.");

  /*
  static Switch<IFMassiveKinematics,bool> interfaceCollinearScheme
    ("CollinearScheme",
     "[experimental] Switch on or off the collinear scheme",
     &IFMassiveKinematics::theCollinearScheme, false, false, false);
  static SwitchOption interfaceCollinearSchemeOn
    (interfaceCollinearScheme,
     "On",
     "Switch on the collinear scheme.",
     true);
  static SwitchOption interfaceCollinearSchemeOff
    (interfaceCollinearScheme,
     "Off",
     "Switch off the collinear scheme",
     false);

  interfaceCollinearScheme.rank(-1);
  */

}

