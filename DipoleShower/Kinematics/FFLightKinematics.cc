// -*- C++ -*-
//
// FFLightKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFLightKinematics class.
//

#include "FFLightKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/DipoleShower/Base/DipoleSplittingInfo.h"

using namespace Herwig;

FFLightKinematics::FFLightKinematics() 
  : DipoleSplittingKinematics() {}

FFLightKinematics::~FFLightKinematics() {}

IBPtr FFLightKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFLightKinematics::fullclone() const {
  return new_ptr(*this);
}

pair<double,double> FFLightKinematics::kappaSupport(const DipoleSplittingInfo&) const {
  return make_pair(0.0,1.0);
}

pair<double,double> FFLightKinematics::xiSupport(const DipoleSplittingInfo& split) const {
  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));
  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emissionData()->id() != ParticleID::g )
      return make_pair(0.5*(1.-c),0.5*(1.+c));
    double b = log((1.+c)/(1.-c));
    return make_pair(-b,b);
  }
  return make_pair(-log(0.5*(1.+c)),-log(0.5*(1.-c)));
}

Energy FFLightKinematics::dipoleScale(const Lorentz5Momentum& pEmitter,
				      const Lorentz5Momentum& pSpectator) const {
  return sqrt(2.*(pEmitter*pSpectator));
}

Energy FFLightKinematics::ptMax(Energy dScale, 
				double, double,
				const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  return dScale/2.;
}

Energy FFLightKinematics::QMax(Energy dScale, 
			       double, double,
			       const DipoleIndex&) const {
  return dScale;
}

Energy FFLightKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return scale*sqrt(z*(1.-z));
}

Energy FFLightKinematics::QFromPt(Energy scale, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return scale/sqrt(z*(1.-z));
}

double FFLightKinematics::ptToRandom(Energy pt, Energy,
				     const DipoleIndex&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

bool FFLightKinematics::generateSplitting(double kappa, double xi, double rphi,
					  DipoleSplittingInfo& info) {

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
  double zp = 0.5*(1.+sqrt(1.-sqr(pt/info.hardPt())));
  double zm = 0.5*(1.-sqrt(1.-sqr(pt/info.hardPt())));

  if ( pt < IRCutoff() || 
       pt > info.hardPt() ||
       z > zp || z < zm ) {
    jacobian(0.0);
    return false;
  }

  double phi = 2.*Constants::pi*rphi;

  jacobian( 2. * mapZJacobian * (1. - sqr(pt) / (s * sqr(info.scale())) ) * 
	    log(0.5 * generator()->maximumCMEnergy()/IRCutoff()) );

  lastPt(pt);
  lastZ(z);
  lastPhi(phi);

  if ( theMCCheck )
    theMCCheck->book(1.,1.,info.scale(),info.hardPt(),pt,z,jacobian());

  return true;

}

InvEnergy2 FFLightKinematics::setKinematics(DipoleSplittingInfo&) const {
  // this is not used anymore
  return ZERO;
}

double FFLightKinematics::
jacobianTimesPropagator(const DipoleSplittingInfo& split,
			Energy scale) const {

  Energy pt = split.lastPt();
  double z = split.lastZ();
  double s = z*(1.-z);
  double zp = 0.5*(1.+sqrt(1.-sqr(pt/split.hardPt())));
  double zm = 0.5*(1.-sqrt(1.-sqr(pt/split.hardPt())));

  if ( pt < IRCutoff() || 
       pt > split.hardPt() ||
       z > zp || z < zm ) {
    return 0.;
  }

  return (2.*scale/pt)*(1.-sqr(pt)/(s*sqr(scale)));
  
}


void FFLightKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  double y = sqr(pt / (pEmitter+pSpectator).m()) / (z*(1.-z));

  Lorentz5Momentum kt =
    getKt(pEmitter, pSpectator, pt, dInfo.lastPhi());

  Lorentz5Momentum em = z*pEmitter + y*(1.-z)*pSpectator + kt;
  em.setMass(0.*GeV);
  em.rescaleEnergy();

  Lorentz5Momentum emm = (1.-z)*pEmitter + z*y*pSpectator - kt;
  emm.setMass(0.*GeV);
  emm.rescaleEnergy();

  Lorentz5Momentum spe = (1.-y)*pSpectator;
  spe.setMass(0.*GeV);
  spe.rescaleEnergy();

  emitterMomentum(em);
  emissionMomentum(emm);
  spectatorMomentum(spe);

}

// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFLightKinematics::persistentOutput(PersistentOStream & ) const {
}

void FFLightKinematics::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFLightKinematics> FFLightKinematics::initFFLightKinematics;
// Definition of the static class description member.

void FFLightKinematics::Init() {

  static ClassDocumentation<FFLightKinematics> documentation
    ("FFLightKinematics implements massless splittings "
     "off a final-final dipole.");

}

