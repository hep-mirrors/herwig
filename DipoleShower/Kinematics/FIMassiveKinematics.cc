// -*- C++ -*-
//
// FIMassiveKinematics.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMassiveKinematics class.
//

#include "FIMassiveKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/DipoleShower/Base/DipoleSplittingInfo.h"
#include "Herwig++/DipoleShower/Kernels/DipoleSplittingKernel.h"

using namespace Herwig;

FIMassiveKinematics::FIMassiveKinematics() 
  : DipoleSplittingKinematics() {}

FIMassiveKinematics::~FIMassiveKinematics() {}

IBPtr FIMassiveKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FIMassiveKinematics::fullclone() const {
  return new_ptr(*this);
}

pair<double,double> FIMassiveKinematics::kappaSupport(const DipoleSplittingInfo&) const {
  return make_pair(0.0,1.0);
}

pair<double,double> FIMassiveKinematics::xiSupport(const DipoleSplittingInfo& split) const {
  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));
  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emissionData()->id() != ParticleID::g )
      return make_pair(0.5*(1.-c),0.5*(1.+c));
    double b = log((1.+c)/(1.-c));
    return make_pair(-b,b);
  }
  return make_pair(-log(0.5*(1.+c)),-log(0.5*(1.-c)));
}

// sbar
Energy FIMassiveKinematics::dipoleScale(const Lorentz5Momentum& pEmitter,
				      const Lorentz5Momentum& pSpectator) const {
  return sqrt(2.*(pEmitter*pSpectator));
}

Energy FIMassiveKinematics::ptMax(Energy dScale, 
				double, double specX,
				const DipoleIndex& ind,
				const DipoleSplittingKernel& split) const {
  Energy2 mi2 = sqr(split.emitter(ind)->mass());
  Energy2 m2  = sqr(split.emission(ind)->mass());
  // Energy2 Mi2 = split.emitter(int)->id() + split.emission(int)->id() == 0 ?
  //   0.*GeV2 : mi2;
  Energy2 Mi2 = mi2 == m2 ? 0.*GeV2 : mi2;

  // s^star/x
  Energy2 s = sqr(dScale) * (1.-specX)/specX + Mi2;

  Energy ptmax = .5 * sqrt(s) * rootOfKallen( s/s, mi2/s, m2/s );
  return ptmax > 0.*GeV ? ptmax : 0.*GeV;
}

// what is this? in FF it is given by y+*dScale = sqrt( 2qi*q / bar )->max
Energy FIMassiveKinematics::QMax(Energy dScale, 
			       double, double specX,
			       const DipoleIndex&) const {
  assert(false && "implementation missing");
  cout << "FIMassiveKinematics::QMax called.\n" << flush;
  // this is sqrt( 2qi*q ) -> max;
  return dScale * sqrt((1.-specX)/specX);
}

Energy FIMassiveKinematics::PtFromQ(Energy, const DipoleSplittingInfo&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy FIMassiveKinematics::QFromPt(Energy, const DipoleSplittingInfo&) const {
  assert(false && "add this");
  return 0.0*GeV;
}


double FIMassiveKinematics::ptToRandom(Energy pt, Energy,
				     const DipoleIndex&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

bool FIMassiveKinematics::generateSplitting(double kappa, double xi, double rphi,
					  DipoleSplittingInfo& info) {

  if ( info.spectatorX() < xMin() ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = IRCutoff() * pow(0.5 * generator()->maximumCMEnergy()/IRCutoff(),kappa);

  if ( pt > info.hardPt() || pt < IRCutoff() ) {
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

  // double s = z*(1.-z);

  // double xs = info.spectatorX();

  // double x = 1. / ( 1. + sqr(pt/info.scale()) / s );

  // double zp = 0.5*(1.+sqrt(1.-sqr(pt/info.hardPt())));
  // double zm = 0.5*(1.-sqrt(1.-sqr(pt/info.hardPt())));

  Energy2 mi2 = sqr(info.emitterData()->mass());
  Energy2 m2  = sqr(info.emissionData()->mass());
  Energy2 Mi2 = info.emitterData()->id()+info.emissionData()->id() == 0 ?
    0.*GeV2 : mi2;

  // s^star/x
  Energy2 s = sqr(info.scale()) * (1.-info.spectatorX())/info.spectatorX() + Mi2;

  double xs = info.spectatorX();
  double x = 1. / ( 1. +
		    ( sqr(pt) + (1.-z)*mi2 + z*m2 - z*(1.-z)*Mi2 ) /
		    ( z*(1.-z)*s ) );

  double zm1 = .5*( 1.+(mi2-m2)/s - rootOfKallen(s/s,mi2/s,m2/s) *
		   sqrt( 1.-sqr(pt/info.hardPt()) ) );
  double zp1 = .5*( 1.+(mi2-m2)/s + rootOfKallen(s/s,mi2/s,m2/s) *
		   sqrt( 1.-sqr(pt/info.hardPt()) ) );

  if ( // pt < IRCutoff() || 
       // pt > info.hardPt() ||
       z > zp1 || z < zm1 ||
       x < xs ) {
    jacobian(0.0);
    return false;
  }

  // additional purely kinematic constraints
  double mui2 = x*mi2/sqr(info.scale());
  double mu2  = x*m2/sqr(info.scale());
  double Mui2 = x*Mi2/sqr(info.scale());
  double xp = 1. + Mui2 - sqr(sqrt(mui2)+sqrt(mu2));
  double zm = .5*( 1.-x+Mui2+mui2-mui2 -
		   sqrt( sqr(1.-x+Mui2-mui2-mu2)-4.*mui2*mu2 ) ) /
    (1.-x+Mui2);
  double zp = .5*( 1.-x+Mui2+mui2-mui2 +
		   sqrt( sqr(1.-x+Mui2-mui2-mu2)-4.*mui2*mu2 ) ) /
    (1.-x+Mui2);

  if (x > xp ||
      z > zp || z < zm ) {
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

InvEnergy2 FIMassiveKinematics::setKinematics(DipoleSplittingInfo& split) const {

  Energy2 mi2 = sqr(split.emitterData()->mass());
  Energy2 m2  = sqr(split.emissionData()->mass());
  Energy2 Mi2 = split.emitterData()->id() + split.emissionData()->id() == 0 ?
    0.*GeV2 : mi2;

  Lorentz5Momentum emitter = split.splitEmitter()->momentum();
  Lorentz5Momentum emission = split.emission()->momentum();
  Lorentz5Momentum spectator = split.splitSpectator()->momentum();

  split.splittingKinematics(const_cast<FIMassiveKinematics*>(this));

  // sbar is the relevant scale
  // Energy2 scale = 2.*(- emission*emitter + emission*spectator + emitter*spectator);
  Energy2 scale = Mi2 - (emitter+emission-spectator).m2();
  split.scale(sqrt(scale));

  double x = 
    scale / (2.*(emitter*spectator + emission*spectator));
  double z = emitter*spectator / (emitter*spectator + emission*spectator);

  split.lastPt( sqrt( z*(1.-z)*(1.-x)/x*scale -
				     ((1.-z)*mi2+z*m2-z*(1.-z)*Mi2) ) );
  split.lastZ(z);

  split.hardPt(split.lastPt());

  if ( split.hardPt() > IRCutoff() ) {
    split.continuesEvolving();
  } else {
    split.didStopEvolving();
  }

  return 1./(2.*x*((emitter+emission).m2()-Mi2));

}

double FIMassiveKinematics::
jacobianTimesPropagator(const DipoleSplittingInfo&,
			Energy) const {
  assert(false && "implementation missing");
  return 0.;
}

void FIMassiveKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  Energy pt = dInfo.lastPt();
  double z = dInfo.lastZ();

  Lorentz5Momentum kt =
    getKt (pSpectator, pEmitter, pt, dInfo.lastPhi(),true);

  Energy2 mi2 = sqr(dInfo.emitterData()->mass());
  Energy2 m2  = sqr(dInfo.emissionData()->mass());
  Energy2 Mi2 = dInfo.emitterData()->id() + dInfo.emissionData()->id() == 0 ?
    0.*GeV2 : mi2;

  double xInv = ( 1. +
		  (pt*pt+(1.-z)*mi2+z*m2-z*(1.-z)*Mi2) /
		  (z*(1.-z)*sqr(dInfo.scale())) );

  Lorentz5Momentum em = z*pEmitter +
    (sqr(pt)+mi2-z*z*Mi2)/(z*sqr(dInfo.scale()))*pSpectator + kt;
  em.setMass(sqrt(mi2));
  em.rescaleEnergy();

  Lorentz5Momentum emm = (1.-z)*pEmitter +
    (pt*pt+m2-sqr(1.-z)*Mi2)/((1.-z)*sqr(dInfo.scale()))*pSpectator - kt;
  emm.setMass(sqrt(m2));
  emm.rescaleEnergy();

  Lorentz5Momentum spe = xInv*pSpectator;
  spe.setMass(ZERO);
  spe.rescaleEnergy();

  emitterMomentum(em);
  emissionMomentum(emm);
  spectatorMomentum(spe);

}

// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMassiveKinematics::persistentOutput(PersistentOStream & ) const {
}

void FIMassiveKinematics::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMassiveKinematics> FIMassiveKinematics::initFIMassiveKinematics;
// Definition of the static class description member.

void FIMassiveKinematics::Init() {

  static ClassDocumentation<FIMassiveKinematics> documentation
    ("FIMassiveKinematics implements massless splittings "
     "off a final-initial dipole.");

}

