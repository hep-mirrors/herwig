// -*- C++ -*-
//
// IFLightKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFLightKinematics class.
//

#include "IFLightKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/DipoleShower/Base/DipoleSplittingInfo.h"

using namespace Herwig;

IFLightKinematics::IFLightKinematics() 
  : DipoleSplittingKinematics(), theCollinearScheme(false) {}

IFLightKinematics::~IFLightKinematics() {}

IBPtr IFLightKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFLightKinematics::fullclone() const {
  return new_ptr(*this);
}

pair<double,double> IFLightKinematics::kappaSupport(const DipoleSplittingInfo&) const {
  return make_pair(0.0,1.0);
}

pair<double,double> IFLightKinematics::xiSupport(const DipoleSplittingInfo& split) const {

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

Energy IFLightKinematics::dipoleScale(const Lorentz5Momentum& pEmitter,
				      const Lorentz5Momentum& pSpectator) const {
  return sqrt(2.*(pEmitter*pSpectator));
}

Energy IFLightKinematics::ptMax(Energy dScale, 
				double emX, double,
				const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  return dScale * sqrt((1.-emX)/emX) /2.;
}

Energy IFLightKinematics::QMax(Energy, 
			       double, double,
			       const DipoleIndex&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy IFLightKinematics::PtFromQ(Energy, const DipoleSplittingInfo&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy IFLightKinematics::QFromPt(Energy, const DipoleSplittingInfo&) const {
  assert(false && "add this");
  return 0.0*GeV;
}


double IFLightKinematics::ptToRandom(Energy pt, Energy,
				     const DipoleIndex&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

bool IFLightKinematics::generateSplitting(double kappar, double xi, double rphi,
					  DipoleSplittingInfo& info) {

  if ( info.emitterX() < xMin() ) {
    jacobian(0.0);
    return false;
  }

  Energy qt = IRCutoff() * pow(0.5 * generator()->maximumCMEnergy()/IRCutoff(),kappar);

  if ( qt < IRCutoff() || qt > info.hardPt() ) {
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

  double kappa = sqr(qt/info.scale());

  double rho = 1. - 4.*kappa*z*(1.-z) / sqr(1. - z + kappa);
  if ( rho < 0. ) {
    jacobian(0.0);
    return false;
  }

  Energy2 pt2 = (sqr(info.scale())/2.)*(1. - z - kappa)*(1. - sqrt(rho));
  if ( pt2 < ZERO ) {
    jacobian(0.0);
    return false;
  }

  double ratio = pt2/sqr(info.scale());

  double x = pt2/sqr(qt);
  double u = ratio/(1.-z);

  if ( x < 0. || x > 1. || u > 1. ) {
    jacobian(0.0);
    return false;
  }

  double xe = info.emitterX();
  double s = sqrt(1.-sqr(qt/info.hardPt()));

  double zp = 0.5*(1.+xe-(1.-xe)*s);
  double zm = 0.5*(1.+xe+(1.-xe)*s);

  if ( z > zp || z < zm ||
       x < xe ) {
    jacobian(0.0);
    return false;
  }

  double phi = 2.*Constants::pi*rphi;

  jacobian(2. * mapZJacobian * ((1.-z)*(1.-z-ratio)/(sqr(z*(1.-z)-ratio) + z*(1.-z)*sqr(1.-z))) * 
	   log(0.5 * generator()->maximumCMEnergy()/IRCutoff()));

  lastPt(qt);
  lastZ(z);
  lastPhi(phi);
  lastEmitterZ(x);

  if ( theMCCheck )
    theMCCheck->book(info.emitterX(),1.,info.scale(),info.hardPt(),qt,z,jacobian());

  return true;

}

InvEnergy2 IFLightKinematics::setKinematics(DipoleSplittingInfo&) const {
  // this is not used anymore
  return ZERO;
}

double IFLightKinematics::
jacobianTimesPropagator(const DipoleSplittingInfo&,
			Energy) const {
  assert(false && "implementation missing");
  return 0.;
}

void IFLightKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  Energy qt = dInfo.lastPt();
  double z = dInfo.lastZ();

  double kappa = sqr(qt)/(-(pEmitter-pSpectator).m2());

  double rho = 1. - 4.*kappa*z*(1.-z) / sqr(1. - z + kappa);
  Energy2 pt2 = (-(pEmitter-pSpectator).m2()/2.)*(1. - z - kappa)*(1. - sqrt(rho));

  double ratio = pt2/(-(pEmitter-pSpectator).m2());

  double x = pt2/sqr(qt);
  double u = ratio/(1.-z);


  Lorentz5Momentum kt =
    getKt (pEmitter, pSpectator, qt, dInfo.lastPhi(),true);

  Lorentz5Momentum em;
  Lorentz5Momentum emm;
  Lorentz5Momentum spe;

  if ( !theCollinearScheme &&
       x > u && (1.-x)/(x-u) < 1. ) {

    em =
      ((1.-u)/(x-u))*pEmitter + ((u/x)*(1.-x)/(x-u))*pSpectator - kt/(x-u);
    em.setMass(0.*GeV);
    em.rescaleEnergy();

    emm =
      ((1.-x)/(x-u))*pEmitter + ((u/x)*(1.-u)/(x-u))*pSpectator - kt/(x-u);
    emm.setMass(0.*GeV);
    emm.rescaleEnergy();

    spe =
      (1.-u/x)*pSpectator;
    spe.setMass(0.*GeV);
    spe.rescaleEnergy();

  } else {

    em = (1./x)*pEmitter;

    emm = ((1.-x)*(1.-u)/x)*pEmitter + u*pSpectator + kt;
    emm.setMass(0.*GeV);
    emm.rescaleEnergy();

    spe = ((1.-x)*u/x)*pEmitter + (1.-u)*pSpectator - kt;
    spe.setMass(0.*GeV);
    spe.rescaleEnergy();

  }
    
  emitterMomentum(em);
  emissionMomentum(emm);
  spectatorMomentum(spe);

}

// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFLightKinematics::persistentOutput(PersistentOStream & os) const {
  os << theCollinearScheme;
}

void IFLightKinematics::persistentInput(PersistentIStream & is, int) {
  is >> theCollinearScheme;
}

ClassDescription<IFLightKinematics> IFLightKinematics::initIFLightKinematics;
// Definition of the static class description member.

void IFLightKinematics::Init() {

  static ClassDocumentation<IFLightKinematics> documentation
    ("IFLightKinematics implements massless splittings "
     "off a initial-final dipole.");


  static Switch<IFLightKinematics,bool> interfaceCollinearScheme
    ("CollinearScheme",
     "[experimental] Switch on or off the collinear scheme",
     &IFLightKinematics::theCollinearScheme, false, false, false);
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

