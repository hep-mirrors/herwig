// -*- C++ -*-
//
// FIMassiveKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"
#include "Herwig/Shower/Dipole/Kernels/DipoleSplittingKernel.h"

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
  return {0.0,1.0};
}

pair<double,double> FIMassiveKinematics::xiSupport(const DipoleSplittingInfo& split) const {
  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));
  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emissionData()->id() != ParticleID::g )
    return {0.5*(1.-c),0.5*(1.+c)};
    double b = log((1.+c)/(1.-c));
        return {-b,b};
  }
    return {-log(0.5*(1.+c)),-log(0.5*(1.-c))};
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
  Energy mi = split.emitter(ind)->mass(), mj = split.emission(ind)->mass();
  Energy2 mi2 = sqr(mi), mj2 = sqr(mj);
  Energy2 mij2 = split.emitter(ind)->id() + split.emission(ind)->id() == 0 ?
    0.*GeV2 : mi2;

  Energy2 sPrime = sqr(dScale) * (1.-specX)/specX + mij2;
  return .5 * sqrt(sPrime) * rootOfKallen( sPrime/sPrime, mi2/sPrime, mj2/sPrime );
}

// what is this? in FF it is given by y+*dScale = sqrt( 2qi*q / bar )->max
Energy FIMassiveKinematics::QMax(Energy dScale, 
			       double, double specX,
			       const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  generator()->log() << "FIMassiveKinematics::QMax called.\n" << flush;
  assert(false && "implementation missing");
  // this is sqrt( 2qi*q ) -> max;
  return dScale * sqrt((1.-specX)/specX);
}

Energy FIMassiveKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  // from Martin's thesis
  double z = split.lastZ();
  Energy mi = split.emitterData()->mass();
  Energy mj = split.emissionData()->mass();
  Energy2 pt2 = z*(1.-z)*sqr(scale) - (1-z)*sqr(mi) - z*sqr(mj);
  assert(pt2 >= ZERO);
  return sqrt(pt2);
}

Energy FIMassiveKinematics::QFromPt(Energy pt, const DipoleSplittingInfo& split) const {
  // from Martin's thesis
  double z = split.lastZ();
  Energy mi = split.emitterData()->mass();
  Energy mj = split.emissionData()->mass();
  Energy2 Q2 = (sqr(pt) + (1-z)*sqr(mi) + z*sqr(mj))/(z*(1.-z));
  return sqrt(Q2);
}

double FIMassiveKinematics::ptToRandom(Energy pt, Energy,
				       double,double,
				       const DipoleIndex&,
				       const DipoleSplittingKernel&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

bool FIMassiveKinematics::generateSplitting(double kappa, double xi, double rphi,
					  DipoleSplittingInfo& info,
					    const DipoleSplittingKernel&) {

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

  // Construct mass squared variables
  Energy2 mi2 = sqr(info.emitterData()->mass());
  Energy2 mj2  = sqr(info.emissionData()->mass());
  Energy2 mij2 = info.emitterData()->id()+info.emissionData()->id() == 0 ?
    0.*GeV2 : mi2;
  Energy2 pt2 = sqr(pt);

  
  // 2 pij.pb
  Energy2 sbar = sqr(info.scale());

  // Compute x
  double x = 1. / ( 1. +
		    ( pt2 + (1.-z)*mi2 + z*mj2 - z*(1.-z)*mij2 ) /
		    ( z*(1.-z)*sbar ) );

  // Check the limit on x
  double xs = info.spectatorX();
  
  if ( x < xs ) {
    jacobian(0.0);
    return false;
  }


  // Compute and check the z limits
  Energy2 sPrime = sbar * (1.-xs)/xs + mij2;
  Energy hard=info.hardPt();

  if(openZBoundaries()==1){
    hard=.5 * sqrt(sPrime) * rootOfKallen( sPrime/sPrime, mi2/sPrime, mj2/sPrime );
  }

  if(openZBoundaries()==2){
    Energy2 s = mij2 - sbar;
	hard=min(0.5*sqrt(sPrime) * 
		rootOfKallen( sPrime/sPrime, mi2/sPrime, mj2/sPrime ) , 
		 0.5*sqrt(s)  * 
		 rootOfKallen( s/s, mi2/s, mj2/s ));
 }

  double ptRatio = sqrt(1.-sqr(pt/hard));

  double zm1 = .5*( 1.+(mi2-mj2)/sPrime - rootOfKallen(sPrime/sPrime,mi2/sPrime,mj2/sPrime) * ptRatio);
  double zp1 = .5*( 1.+(mi2-mj2)/sPrime + rootOfKallen(sPrime/sPrime,mi2/sPrime,mj2/sPrime) * ptRatio);

  if ( z > zp1 || z < zm1 ) {
    jacobian(0.0);
    return false;
  }

  // additional purely kinematic constraints from
  // the integration limits in Catani-Seymour
  double mui2CS = x*mi2/sbar;
  double muj2CS  = x*mj2/sbar;
  double muij2CS = x*mij2/sbar;

  // Limit on x
  double xp = 1. + muij2CS - sqr(sqrt(mui2CS)+sqrt(muj2CS));  
  if (x > xp ) {
    jacobian(0.0);
    return false;
  }

  // Limit on z
  double root = sqr(1.-x+muij2CS-mui2CS-muj2CS)-4.*mui2CS*muj2CS;

  if( root < 0. && root>-1e-10 ) {
    //    assert(false);
    root = 0.;
  }
  else if (root <0. ) {
    jacobian(0.0);
    return false;
  }

  root = sqrt(root);
  double zm2 = .5*( 1.-x+muij2CS+mui2CS-muj2CS - root ) / (1.-x+muij2CS);
  double zp2 = .5*( 1.-x+muij2CS+mui2CS-muj2CS + root ) / (1.-x+muij2CS);

  if ( z > zp2 || z < zm2 ) {
    jacobian(0.0);
    return false;
  }

  // Store the splitting variables
  double phi = 2.*Constants::pi*rphi;
  
  // Compute and store the jacobian
  double jacPt2 = 1. / ( 1. + (1.-z)*mi2/pt2 + z*mj2/pt2 - z*(1.-z)*mij2/pt2 );
  jacobian( jacPt2 * mapZJacobian * 2.*log(0.5 * generator()->maximumCMEnergy()/IRCutoff()));

  lastPt(pt);
  lastZ(z);
  lastPhi(phi);
  lastSpectatorZ(x);

  if ( theMCCheck )
    theMCCheck->book(1.,info.spectatorX(),info.scale(),info.hardPt(),pt,z,jacobian());

  return true;

}

void FIMassiveKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  // Get splitting variables
  Energy pt = dInfo.lastPt();
  double z = dInfo.lastZ();

  // Compute sqr scales
  Energy2 pt2 = sqr(pt);
  Energy2 sbar = sqr(dInfo.scale());
  
  Lorentz5Momentum kt =
    getKt (pSpectator, pEmitter, pt, dInfo.lastPhi(),true);

  Energy2 mi2 = sqr(dInfo.emitterData()->mass());
  Energy2 mj2  = sqr(dInfo.emissionData()->mass());
  Energy2 mij2 = dInfo.emitterData()->id() + dInfo.emissionData()->id() == 0 ?
    0.*GeV2 : mi2;

  double xInv = ( 1. +
		  (pt2+(1.-z)*mi2+z*mj2-z*(1.-z)*mij2) /
		  (z*(1.-z)*sbar) );

  Lorentz5Momentum em = z*pEmitter +
    (pt2+mi2-z*z*mij2)/(z*sbar)*pSpectator + kt;
  em.setMass(sqrt(mi2));
  em.rescaleEnergy();

  Lorentz5Momentum emm = (1.-z)*pEmitter +
    (pt2+mj2-sqr(1.-z)*mij2)/((1.-z)*sbar)*pSpectator - kt;
  emm.setMass(sqrt(mj2));
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

