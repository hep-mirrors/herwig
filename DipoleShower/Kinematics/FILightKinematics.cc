// -*- C++ -*-
//
// FILightKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FILightKinematics class.
//

#include "FILightKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/DipoleShower/Base/DipoleSplittingInfo.h"

using namespace Herwig;

FILightKinematics::FILightKinematics() 
  : DipoleSplittingKinematics() {}

FILightKinematics::~FILightKinematics() {}

IBPtr FILightKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FILightKinematics::fullclone() const {
  return new_ptr(*this);
}

pair<double,double> FILightKinematics::kappaSupport(const DipoleSplittingInfo&) const {
  return make_pair(0.0,1.0);
}

pair<double,double> FILightKinematics::xiSupport(const DipoleSplittingInfo& split) const {
  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));
  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emissionData()->id() != ParticleID::g )
      return make_pair(0.5*(1.-c),0.5*(1.+c));
    double b = log((1.+c)/(1.-c));
    return make_pair(-b,b);
  }
  return make_pair(-log(0.5*(1.+c)),-log(0.5*(1.-c)));
}

Energy FILightKinematics::dipoleScale(const Lorentz5Momentum& pEmitter,
				      const Lorentz5Momentum& pSpectator) const {
  return sqrt(2.*(pEmitter*pSpectator));
}

Energy FILightKinematics::ptMax(Energy dScale, 
				double, double specX,
				const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  return dScale * sqrt((1.-specX)/specX) /2.;
}

Energy FILightKinematics::QMax(Energy dScale, 
			       double, double specX,
			       const DipoleIndex&) const {
  return dScale * sqrt((1.-specX)/specX);
}

Energy FILightKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return scale*sqrt(z*(1.-z));
}

Energy FILightKinematics::QFromPt(Energy scale, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return scale/sqrt(z*(1.-z));
}


double FILightKinematics::ptToRandom(Energy pt, Energy,
				     const DipoleIndex&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

bool FILightKinematics::generateSplitting(double kappa, double xi, double rphi,
					  DipoleSplittingInfo& info) {

  if ( info.spectatorX() < xMin() ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = IRCutoff() * pow(0.5 * generator()->maximumCMEnergy()/IRCutoff(),kappa);

  if ( pt > info.hardPt() ) {
    jacobian(0.0);
    return false;
  }

  double z;
  double mapZJacobian;

  if ( info.index().emitterData()->id() == ParticleID::g ) {
    if ( info.emissionData()->id() != ParticleID::g ) {
      z = xi;
      mapZJacobian = 1.;
    } else {
      z = exp(xi)/(1.+exp(xi));
      mapZJacobian = z*(1.-z);
    }
  } else {
    z = 1.-exp(-xi);
    mapZJacobian = 1.-z;
  }

  double s = z*(1.-z);

  double xs = info.spectatorX();

  double x = 1. / ( 1. + sqr(pt/info.scale()) / s );

  double zp = 0.5*(1.+sqrt(1.-sqr(pt/info.hardPt())));
  double zm = 0.5*(1.-sqrt(1.-sqr(pt/info.hardPt())));

  if ( pt < IRCutoff() || 
       pt > info.hardPt() ||
       z > zp || z < zm ||
       x < xs ) {
    jacobian(0.0);
    return false;
  }

  double phi = 2.*Constants::pi*rphi;

  jacobian( 2. * mapZJacobian * log(0.5 * generator()->maximumCMEnergy()/IRCutoff()));

  lastPt(pt);
  lastZ(z);
  lastPhi(phi);
  lastSpectatorZ(x);

  if ( theMCCheck )
    theMCCheck->book(1.,info.spectatorX(),info.scale(),info.hardPt(),pt,z,jacobian());

  return true;

}

InvEnergy2 FILightKinematics::setKinematics(DipoleSplittingInfo& split) const {

  Lorentz5Momentum emitter = split.splitEmitter()->momentum();
  Lorentz5Momentum emission = split.emission()->momentum();
  Lorentz5Momentum spectator = split.splitSpectator()->momentum();

  split.splittingKinematics(const_cast<FILightKinematics*>(this));

  Energy2 scale = 2.*(- emission*emitter + emission*spectator + emitter*spectator);
  split.scale(sqrt(scale));

  double x = 
    scale / (2.*(emitter*spectator + emission*spectator));
  double z = emitter*spectator / (emitter*spectator + emission*spectator);

  split.lastPt(split.scale() * sqrt(z*(1.-z)*(1.-x)/x));
  split.lastZ(z);

  split.hardPt(split.lastPt());

  if ( split.hardPt() > IRCutoff() ) {
    split.continuesEvolving();
  } else {
    split.didStopEvolving();
  }

  return 1./(2.*x*(emitter*emission));

}

double FILightKinematics::
jacobianTimesPropagator(const DipoleSplittingInfo&,
			Energy) const {
  assert(false && "implementation missing");
  return 0.;
}


void FILightKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  Energy pt = dInfo.lastPt();
  double z = dInfo.lastZ();

  Lorentz5Momentum kt =
    getKt (pEmitter, pSpectator, pt, dInfo.lastPhi());

  double ratio = sqr(pt/(-(pEmitter-pSpectator).m()));
  double xInv = (1.+ratio/(z*(1.-z)));

  Lorentz5Momentum em = z*pEmitter + (ratio/z)*pSpectator + kt;
  em.setMass(0.*GeV);
  em.rescaleEnergy();

  Lorentz5Momentum emm = (1.-z)*pEmitter + (ratio/(1.-z))*pSpectator - kt;
  emm.setMass(0.*GeV);
  emm.rescaleEnergy();

  Lorentz5Momentum spe = xInv*pSpectator;
  spe.setMass(0.*GeV);
  spe.rescaleEnergy();

  emitterMomentum(em);
  emissionMomentum(emm);
  spectatorMomentum(spe);

}

// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FILightKinematics::persistentOutput(PersistentOStream & ) const {
}

void FILightKinematics::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FILightKinematics> FILightKinematics::initFILightKinematics;
// Definition of the static class description member.

void FILightKinematics::Init() {

  static ClassDocumentation<FILightKinematics> documentation
    ("FILightKinematics implements massless splittings "
     "off a final-initial dipole.");

}

