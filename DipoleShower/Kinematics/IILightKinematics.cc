// -*- C++ -*-
//
// IILightKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IILightKinematics class.
//

#include "IILightKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/DipoleShower/Base/DipoleSplittingInfo.h"

using namespace Herwig;

IILightKinematics::IILightKinematics() 
  : DipoleSplittingKinematics(), theCollinearScheme(false), didCollinear(false) {}

IILightKinematics::~IILightKinematics() {}

IBPtr IILightKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IILightKinematics::fullclone() const {
  return new_ptr(*this);
}

pair<double,double> IILightKinematics::kappaSupport(const DipoleSplittingInfo&) const {
  return make_pair(0.0,1.0);
}

pair<double,double> IILightKinematics::xiSupport(const DipoleSplittingInfo& split) const {

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

Energy IILightKinematics::dipoleScale(const Lorentz5Momentum& pEmitter,
				      const Lorentz5Momentum& pSpectator) const {
  return sqrt(2.*(pEmitter*pSpectator));
}

Energy IILightKinematics::ptMax(Energy dScale, 
				double emX, double specX,
				const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  double tau = emX*specX;
  return (1.-tau) * dScale / 2.;
}

Energy IILightKinematics::QMax(Energy, 
			       double, double,
			       const DipoleIndex&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy IILightKinematics::PtFromQ(Energy, const DipoleSplittingInfo&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy IILightKinematics::QFromPt(Energy, const DipoleSplittingInfo&) const {
  assert(false && "add this");
  return 0.0*GeV;
}


double IILightKinematics::ptToRandom(Energy pt, Energy,
				     const DipoleIndex&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

bool IILightKinematics::generateSplitting(double kappa, double xi, double rphi,
					  DipoleSplittingInfo& info) {

  if ( info.emitterX() < xMin() ||
       info.spectatorX() < xMin() ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = IRCutoff() * pow(0.5 * generator()->maximumCMEnergy()/IRCutoff(),kappa);

  if ( pt > info.hardPt() ) {
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
      mapZJacobian = (1.-z);
    } else {
      z = xi;
      mapZJacobian = 1.;
    }
  }

  double ratio = sqr(pt/info.scale());

  double x = ( z*(1.-z) - ratio ) / ( 1. - z );
  double v = ratio / (1.-z);

  if ( x < 0. || x > 1. || v > 1. || v > 1.-x ) {
    jacobian(0.0);
    return false;
  }

  double tau = info.emitterX()*info.spectatorX();

  double zp = 0.5*( 1.+ tau + 
		    (1.-tau)*sqrt(1.-sqr(pt/info.hardPt()) ) );
  double zm = 0.5*( 1.+ tau -
		    (1.-tau)*sqrt(1.-sqr(pt/info.hardPt()) ) );

  if ( pt < IRCutoff() ||
       pt > info.hardPt() || 
       z > zp || z < zm ) {
    jacobian(0.0);
    return false;
  }

  if ( !theCollinearScheme &&
       (1.-v-x)/(v+x) < 1. ) {
    if ( (x+v) < info.emitterX() ||
	 x/(x+v) < info.spectatorX() ) {
      jacobian(0.0);
      return false;
    }
  } else {
    if ( x < info.emitterX() ) {
      jacobian(0.0);
      return false;
    }
  }

  double phi = 2.*Constants::pi*rphi;

  jacobian(2. * mapZJacobian * (1.-z)/(z*(1.-z)-ratio) * log(0.5 * generator()->maximumCMEnergy()/IRCutoff()));

  lastPt(pt);
  lastZ(z);
  lastPhi(phi);

  if ( !theCollinearScheme &&
       (1.-v-x)/(v+x) < 1. ) {
    lastEmitterZ(x+v);
    lastSpectatorZ(x/(x+v));
  } else {
    lastEmitterZ(x);
    lastSpectatorZ(1.);
  }

  if ( theMCCheck )
    theMCCheck->book(info.emitterX(),info.spectatorX(),info.scale(),info.hardPt(),pt,z,jacobian());

  return true;

}

InvEnergy2 IILightKinematics::setKinematics(DipoleSplittingInfo& split) const {

  Lorentz5Momentum emitter = split.splitEmitter()->momentum();
  Lorentz5Momentum emission = split.emission()->momentum();
  Lorentz5Momentum spectator = split.splitSpectator()->momentum();

  split.splittingKinematics(const_cast<IILightKinematics*>(this));

  Energy2 scale = 2.*(-emission*emitter - emission*spectator + emitter*spectator);
  split.scale(sqrt(scale));

  double x = scale/(2.*(emitter*spectator));
  double v = (emitter*emission)/(emitter*spectator);

  split.lastZ(v+x);
  split.lastPt(split.scale() * sqrt(v*(1.-x-v)));

  split.hardPt(split.lastPt());

  if ( split.hardPt() > IRCutoff() ) {
    split.continuesEvolving();
  } else {
    split.didStopEvolving();
  }

  return 1./(2.*x*(emitter*emission));

}

double IILightKinematics::
jacobianTimesPropagator(const DipoleSplittingInfo&,
			Energy) const {
  assert(false && "implementation missing");
  return 0.;
}

void IILightKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  Energy pt = dInfo.lastPt();
  double z = dInfo.lastZ();

  double ratio = sqr(pt/(pEmitter+pSpectator).m());

  double x = ( z*(1.-z) - ratio ) / ( 1. - z );
  double v = ratio / (1.-z);

  pt = sqrt(v*(1.-x-v)/x) * (pEmitter+pSpectator).m();

  Lorentz5Momentum kt =
    getKt (pEmitter, pSpectator, pt, dInfo.lastPhi());

  if ( !theCollinearScheme &&
       (1.-v-x)/(v+x) < 1. ) {

    Lorentz5Momentum em =
      (1./(v+x))*pEmitter+(v*(1.-v-x)/(x*(x+v)))*pSpectator+kt/(x+v);
    em.setMass(0.*GeV);
    em.rescaleEnergy();

    Lorentz5Momentum emm =
      ((1.-v-x)/(v+x))*pEmitter+(v/(x*(x+v)))*pSpectator+kt/(x+v);
    emm.setMass(0.*GeV);
    emm.rescaleEnergy();

    Lorentz5Momentum spe =
      (1.+v/x)*pSpectator;
    spe.setMass(0.*GeV);
    spe.rescaleEnergy();

    emitterMomentum(em);
    emissionMomentum(emm);
    spectatorMomentum(spe);

    didCollinear = false;

  } else {

    Lorentz5Momentum em =
      (1./x)*pEmitter;
    em.setMass(0.*GeV);
    em.rescaleEnergy();

    Lorentz5Momentum emm =
      ((1.-x-v)/x)*pEmitter+v*pSpectator+kt;
    emm.setMass(0.*GeV);
    emm.rescaleEnergy();

    Lorentz5Momentum spe =
      pSpectator;

    emitterMomentum(em);
    emissionMomentum(emm);
    spectatorMomentum(spe);

    K = em + spe - emm;
    K2 = K.m2();
    
    Ktilde = pEmitter + pSpectator;
    KplusKtilde = K + Ktilde;
    
    KplusKtilde2 = KplusKtilde.m2();

    didCollinear = true;

  }

}

// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IILightKinematics::persistentOutput(PersistentOStream & os) const {
  os << theCollinearScheme;
}

void IILightKinematics::persistentInput(PersistentIStream & is, int) {
  is >> theCollinearScheme;
}

ClassDescription<IILightKinematics> IILightKinematics::initIILightKinematics;
// Definition of the static class description member.

void IILightKinematics::Init() {

  static ClassDocumentation<IILightKinematics> documentation
    ("IILightKinematics implements massless splittings "
     "off an initial-initial dipole.");


  static Switch<IILightKinematics,bool> interfaceCollinearScheme
    ("CollinearScheme",
     "[experimental] Switch on or off the collinear scheme",
     &IILightKinematics::theCollinearScheme, false, false, false);
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

}

