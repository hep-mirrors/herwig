// -*- C++ -*-
//
// IFMassiveKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"
#include "Herwig/Shower/Dipole/Kernels/DipoleSplittingKernel.h"

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
  return {0.0,1.0};
}

pair<double,double> IFMassiveKinematics::xiSupport(const DipoleSplittingInfo& split) const {

  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));

  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emitterData()->id() == ParticleID::g ) {
      double b = log((1.+c)/(1.-c));
      return {-b,b};
    } else {
      return {log(0.5*(1.-c)),log(0.5*(1.+c))};
    }
  }

  if ( split.index().emitterData()->id() != ParticleID::g &&
       split.emitterData()->id() != ParticleID::g ) {
    return {-log(0.5*(1.+c)),-log(0.5*(1.-c))};
  }

  return {0.5*(1.-c),0.5*(1.+c)};

}

// sbar
Energy IFMassiveKinematics::dipoleScale(const Lorentz5Momentum& pEmitter,
					const Lorentz5Momentum& pSpectator) const {
  return sqrt(2.*(pEmitter*pSpectator));
}

Energy IFMassiveKinematics::ptMax(Energy dScale, 
				  double emX, double,
				  const DipoleSplittingInfo& dInfo,	
				  const DipoleSplittingKernel&) const {

  Energy2 A = sqr(dScale) * (1.-emX)/emX;
  Energy2 mk2 = sqr(dInfo.spectatorMass());
  Energy ptMax = 0.5*A/sqrt(mk2+A);
  return ptMax;
}

Energy IFMassiveKinematics::ptMax(Energy dScale, 
				  double emX, double,
				  const DipoleIndex&,
				  const DipoleSplittingKernel&,
				  tPPtr, tPPtr spectator) const {

  Energy2 A = sqr(dScale) * (1.-emX)/emX;
  Energy2 mk2 = sqr(spectator->mass());
  Energy ptMax = 0.5*A/sqrt(mk2+A);
  return ptMax;
}

Energy IFMassiveKinematics::QMax(Energy, 
				 double, double,
				 const DipoleSplittingInfo&,
				 const DipoleSplittingKernel&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy IFMassiveKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return scale*sqrt(1.-z);
}

Energy IFMassiveKinematics::QFromPt(Energy pt, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return pt/sqrt(1.-z);
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

  // Check emitter x against xmin
  if ( info.emitterX() < xMin() ) {
    jacobian(0.0);
    return false;
  }

  // Generate pt and check it against max allowed
  Energy pt = IRCutoff() * pow(0.5 * generator()->maximumCMEnergy()/IRCutoff(),kappa);
  if ( pt < IRCutoff() || pt > info.hardPt() ) {
    jacobian(0.0);
    return false;
  }

  // Compute scales required
  Energy2 pt2 = sqr(pt);
  Energy2 saj = sqr(info.scale());
  Energy2 mk2 = sqr(info.spectatorMass());

  // Generate z
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

  // Check limits on z
  double xe = info.emitterX();
  Energy hard = info.hardPt();

  if(openZBoundaries()==1){
        Energy2 A = saj*(1.-xe)/xe;
        hard = 0.5*A/sqrt(mk2+A);          
  }
  if(openZBoundaries()==2){
        Energy2 A = saj*min(1.,(1.-xe)/xe);
        hard= 0.5*A/sqrt(mk2+A);
	assert(pt2<=sqr(hard));
  }

  double ptRatio = sqrt(1. - pt2/sqr(hard) );
  double zp = 0.5*(1.+xe + (1.-xe)*ptRatio);
  double zm = 0.5*(1.+xe - (1.-xe)*ptRatio);
  
  if ( z < zm || z > zp ) {
    jacobian(0.0);
    return false;
  }
    
  // Calculate x and u in terms of z and pt
  double r = pt2/saj;
  double muk2 = mk2/saj;
  double rho = 1. - 4.*r*(1.-muk2)*z*(1.-z)/sqr(1.-z+r);
  if ( rho < 0.0 ) {
    // This has never happened
    jacobian(0.0);
    return false;
  }

  double x = 0.5*((1.-z+r)/(r*(1.-muk2))) * (1. - sqrt(rho));
  double u = x*r / (1.-z);

  // Check limits on x and u      
  // Following Catani-Seymour paper
  double muk2CS = x*muk2;
  double up = (1.-x) / ( 1.-x + muk2CS );
  if ( x < xe || x > 1. ||
       u < 0. || u > up ) {
    jacobian(0.0);
    return false;
  }


  // Compute the Jacobian
  double jac = 1./(u + x - 2.*u*x*(1.-muk2));

  jacobian( jac * mapZJacobian * 2. * log(0.5 * generator()->maximumCMEnergy()/IRCutoff()));
    
  // Log results
  double phi = 2.*Constants::pi*rphi;
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

  // Initialise the momenta
  Lorentz5Momentum em;
  Lorentz5Momentum emm;
  Lorentz5Momentum spe;

  if (!theCollinearScheme) {
    assert(false);

    Energy2 sbar = 2.*pEmitter*pSpectator;
    Energy pt = dInfo.lastPt();
    double ratio = pt*pt/sbar;
    double z = dInfo.lastZ();
    double x = (z*(1.-z)-ratio)/(1.-z-ratio);
    double u = ratio / (1.-z);

    pt = sqrt(sbar*u*(1.-u)*(1.-x));
    Energy magKt = 
      sqrt(sbar*u*(1.-u)*(1.-x)/x - sqr(u*dInfo.spectatorMass()));
    Lorentz5Momentum kt =
      getKt (pSpectator, pEmitter, magKt, dInfo.lastPhi(),true);

    Energy2 mj2 = sqr(dInfo.spectatorMass());
    double alpha = 1. - 2.*mj2/sbar;

    if ( x > u && (1.-x)/(x-u) < 1. ) {

      double fkt = sqrt(sqr(x-u)+4.*x*u*mj2/sbar);

      //    em =
      //      ((1.-u)/(x-u))*pEmitter + ((u/x)*(1.-x)/(x-u))*pSpectator - kt/(x-u);
      Energy2 fa = (sbar*(x+u-2.*x*z)+2.*mj2*x*u) / sqrt(sqr(x-u)+4.*x*u*mj2/sbar);
      double a = (-sbar+fa) / (2.*x*(sbar-mj2));
      double ap = (sbar+alpha*fa) / (2.*x*(sbar-mj2));
      em = ap*pEmitter + a*pSpectator - fkt*kt;

      //    emm =
      //      ((1.-x)/(x-u))*pEmitter + ((u/x)*(1.-u)/(x-u))*pSpectator - kt/(x-u);
      Energy2 fb = abs(sbar*(u*(1.-u)-x*(1.-x))+2.*mj2*x*u) / sqrt(sqr(x-u)+4.*x*u*mj2/sbar);
      double b = (-sbar*(1.-x-u)+fb) / (2.*x*(sbar-mj2));
      double bp = (sbar*(1.-x-u)+alpha*fb) / (2.*x*(sbar-mj2));
      emm = bp*pEmitter + b*pSpectator + fkt*kt;

      //    spe =
      //      (1.-u/x)*pSpectator;
      Energy2 fc = sqrt(sqr(sbar*(x-u))+4.*sbar*mj2*x*u);
      double c = (sbar*(x-u)-2.*x*mj2+fc) / (2.*x*(sbar-mj2));
      double cp = (-sbar*(x-u)+2.*x*mj2+alpha*fc) / (2.*x*(sbar-mj2));
      spe = cp*pEmitter + c*pSpectator;

    }
  }

  else {
    
    // Get z, pt and the relevant scales
    double z = dInfo.lastZ();
    Energy pt = dInfo.lastPt();
    
    Energy2 pt2 = sqr(pt);
    Energy2 saj = 2.*pEmitter*pSpectator;

    double muk2 = sqr(dInfo.spectatorMass())/saj;
    double r = pt2/saj;
    
    // Calculate x and u
    double rho = 1. - 4.*r*(1.-muk2)*z*(1.-z)/sqr(1.-z+r);
    double x = 0.5*((1.-z+r)/(r*(1.-muk2))) * (1. - sqrt(rho));
    double u = x*r / (1.-z);
    
    // Generate kt
    Lorentz5Momentum kt = getKt(pEmitter, pSpectator, pt, dInfo.lastPhi(), true);
 
    // Set the momenta  
    em = (1./x)*pEmitter;
    emm = ((1.-x)*(1.-u)/x - 2.*u*muk2)*pEmitter + u*pSpectator + kt;
    spe = ((1.-x)*u/x + 2.*u*muk2)*pEmitter + (1.-u)*pSpectator - kt;
    
  }
  
  em.setMass(ZERO);
  em.rescaleEnergy();

  emm.setMass(ZERO);
  emm.rescaleEnergy();

  spe.setMass(dInfo.spectatorMass());
  spe.rescaleEnergy();

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
    static SwitchOption interfaceCollinearSchemeYes
    (interfaceCollinearScheme,
    "Yes",
    "Switch on the collinear scheme.",
    true);
    static SwitchOption interfaceCollinearSchemeNo
    (interfaceCollinearScheme,
    "No",
    "Switch off the collinear scheme",
    false);

    interfaceCollinearScheme.rank(-1);
  */

}

