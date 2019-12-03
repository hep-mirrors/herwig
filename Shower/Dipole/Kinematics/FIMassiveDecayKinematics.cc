// -*- C++ -*-
//
// FIMassiveDecayKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMassiveDecayKinematics class.
//

#include "FIMassiveDecayKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"
#include "Herwig/Shower/Dipole/Kernels/DipoleSplittingKernel.h"

#include "ThePEG/Interface/Switch.h"


using namespace Herwig;

FIMassiveDecayKinematics::FIMassiveDecayKinematics() 
  : DipoleSplittingKinematics() {}

FIMassiveDecayKinematics::~FIMassiveDecayKinematics() {}

IBPtr FIMassiveDecayKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FIMassiveDecayKinematics::fullclone() const {
  return new_ptr(*this);
}

pair<double,double> FIMassiveDecayKinematics::kappaSupport(const DipoleSplittingInfo&) const {
  return {0.0,1.0};
}

pair<double,double> FIMassiveDecayKinematics::xiSupport(const DipoleSplittingInfo& split) const {

  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));

  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emissionData()->id() != ParticleID::g ){
      return {0.5*(1.-c),0.5*(1.+c)};
    }
    double b = log((1.+c)/(1.-c));
    return {-b,b};
  }
  return {-log(0.5*(1.+c)),-log(0.5*(1.-c))};
}

Energy FIMassiveDecayKinematics::dipoleScale(const Lorentz5Momentum&,
                                             const Lorentz5Momentum& pSpectator) const {
  return pSpectator.m();
}

Energy FIMassiveDecayKinematics::recoilMassKin(const Lorentz5Momentum& pEmitter,
                                               const Lorentz5Momentum& pSpectator) const {
  Lorentz5Momentum pk = pSpectator - pEmitter;
  Energy pkmass = pk.m();
  return pkmass;
}

Energy FIMassiveDecayKinematics::ptMax(Energy dScale, 
                                       double, double,
                                       const DipoleSplittingInfo& dInfo,
                                       const DipoleSplittingKernel& split) const {


  DipoleIndex ind = dInfo.index();
  double mui2 = 0.;

  // g->gg and g->qqbar
  if ( abs(split.emitter(ind)->id()) == abs(split.emission(ind)->id()) ) {
    mui2 = sqr(split.emitter(ind)->mass() / dScale);
  }
  // Otherwise have X->Xg (should work for SUSY)
  else {
    mui2 = sqr(dInfo.emitterMass()/dScale);
  }
  double muj2  = sqr( split.emission(ind)->mass() / dScale );
  
  // Mass of recoil system
  double muk = dInfo.recoilMass() / dScale;

  return rootOfKallen( mui2, muj2, sqr(1.-muk) ) / ( 2.-2.*muk ) * dScale;
}


Energy FIMassiveDecayKinematics::ptMax(Energy dScale, 
                                       double, double,
                                       const DipoleIndex& ind,
                                       const DipoleSplittingKernel& split,
                                       tPPtr emitter, tPPtr spectator) const {
  double mui2 = 0.;
  // g->gg and g->qqbar
  if ( abs(split.emitter(ind)->id()) == abs(split.emission(ind)->id()) ) {
    mui2 = sqr(split.emitter(ind)->mass() / dScale);
  }
  // Otherwise have X->Xg (should work for SUSY)
  else {
    mui2 = sqr(emitter->mass()/dScale);
  }
  double muj2  = sqr(split.emission(ind)->mass() / dScale);
  
  // Mass of recoil system
  double muk = recoilMassKin(emitter->momentum(),spectator->momentum()) / dScale;

  return rootOfKallen( mui2, muj2, sqr(1.-muk) ) / ( 2.-2.*muk ) * dScale;
}


Energy FIMassiveDecayKinematics::QMax(Energy dScale, 
                                      double, double,
                                      const DipoleSplittingInfo& dInfo,
                                      const DipoleSplittingKernel&) const {
  assert(false && "implementation missing");
  // Mass of recoil system
  double Muk2 = sqr( dInfo.recoilMass() / dScale );
  double Muk = sqrt( Muk2 );

  return dScale * ( 1.-2.*Muk+Muk2 );
}

// The name of this function is misleading
// scale here is defined as sqr(scale) = sqr(qi+qj)
// Here, scale is Q
Energy FIMassiveDecayKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  double zPrime = split.lastSplittingParameters()[0];
  
  // masses
  Energy2 mi2 = ZERO;
  // g->gg and g->qqbar
  if ( abs(split.emitterData()->id()) == abs(split.emissionData()->id()) ) {
    mi2 = sqr( split.emitterData()->mass());
  }
  // Otherwise have X->Xg (should work for SUSY)
  else {
    mi2 = sqr(split.emitterMass());
  }
  Energy2 m2  = sqr(split.emissionData()->mass());

  Energy2 pt2 = zPrime*(1.-zPrime)*sqr(scale) - (1-zPrime)*mi2 - zPrime*m2;
  assert(pt2 >= ZERO);
  return sqrt(pt2);
}


// This is simply the inverse of PtFromQ
Energy FIMassiveDecayKinematics::QFromPt(Energy pt, const DipoleSplittingInfo& split) const {
  double zPrime = split.lastSplittingParameters()[0];
  
  // masses
  Energy2 mi2 = ZERO;
  // g->gg and g->qqbar
  if ( abs(split.emitterData()->id()) == abs(split.emissionData()->id()) ) {
    mi2 = sqr( split.emitterData()->mass());
  }
  // Otherwise have X->Xg (should work for SUSY)
  else {
    mi2 = sqr(split.emitterMass());
  }
  Energy2 m2  = sqr(split.emissionData()->mass());

  Energy2 Q2 = (sqr(pt) + (1-zPrime)*mi2 + zPrime*m2)/(zPrime*(1.-zPrime));
  return sqrt(Q2);
}

double FIMassiveDecayKinematics::ptToRandom(Energy pt, Energy,
                                            double,double,
                                            const DipoleIndex&,
                                            const DipoleSplittingKernel&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

// SW, 14/02/2019: Tidied to match thesis
bool FIMassiveDecayKinematics::generateSplitting(double kappa, double xi, double rphi,
                                                 DipoleSplittingInfo& info,
                                                 const DipoleSplittingKernel&) {
  
  // Scale 's' and masses
  Energy2 Qijk = sqr(info.scale());
  Energy2 mij2 = sqr(info.emitterMass());
  Energy2 Mk2 = sqr(info.recoilMass());

  // To solve issue with scale during presampling
  // need to enforce that Qijk-mij2-mk2 = 2*pij.pk > 0.
  // Combine checks by comparing against square root
  if ( Qijk-mij2-Mk2 < sqrt(4.*mij2*Mk2) ) {
    jacobian(0.0);
    return false;
  }

  Energy2 mk2 = Mk2;  
  Energy2 mi2 = ZERO;
  // g->gg and g->qqbar
  if ( abs(info.emitterData()->id()) == abs(info.emissionData()->id()) ) {
    mi2 = sqr(info.emitterData()->mass());
  }
  // Otherwise have X->Xg (should work for SUSY)
  else {
    mi2 = mij2;
  }
  Energy2 mj2 = sqr(info.emissionData()->mass());

  // Calculate pt 
  Energy pt = IRCutoff() * pow(0.5 * generator()->maximumCMEnergy()/IRCutoff(),kappa);
  Energy2 pt2 = sqr(pt);
  
  if ( pt > info.hardPt() || pt < IRCutoff() ) {
    jacobian(0.0);
    return false;
  }

  // Generate z
  double z;
  if ( info.index().emitterData()->id() == ParticleID::g ) {
    if ( info.emissionData()->id() != ParticleID::g ) {
      z = xi;
    } 
    else {
      z = exp(xi)/(1.+exp(xi));
    }
  } 
  else {
    z = 1.-exp(-xi);
  }

  // new: 2011-08-31
  // 2011-11-08: this does happen
  if( (sqrt(mi2)+sqrt(mj2)+sqrt(mk2))/ sqrt(Qijk) > 1. ){
    jacobian(0.0);
    return false;
  }

  // Limits on z.
  // Phasespace constraint to incorporate ptMax.
  Energy hard = info.hardPt();
  Energy2 sqrRootQijkMk = sqr(sqrt(Qijk)-sqrt(mk2));

  if(openZBoundaries()>0){
	// From ptMax(..)
	hard = rootOfKallen(mi2, mj2, sqrRootQijkMk) / ( 2.*sqrt(sqrRootQijkMk) );
	assert(pt<=hard);
  }

  double ptRatio = sqrt(1.-sqr(pt/hard));  
  double zp1 = ( mi2 - mj2 + sqrRootQijkMk +
                 rootOfKallen(mi2,mj2,sqrRootQijkMk) * ptRatio )
    / 2. / sqrRootQijkMk ;
  double zm1 = ( mi2 - mj2 + sqrRootQijkMk -
                 rootOfKallen(mi2,mj2,sqrRootQijkMk) * ptRatio )
    / 2. / sqrRootQijkMk ;
  
  if ( z > zp1 || z < zm1 ) {
    jacobian(0.0);
    return false;
  }

  // Calculate y
  Energy2 sbar = Qijk - mi2 - mj2 - mk2;
  double y = (pt2 + sqr(1.-z)*mi2 + sqr(z)*mj2)
    / sbar / z / (1.-z);

  // Kinematic phasespace boundaries for y.
  // Same as in Dittmaier hep-ph/9904440v2 (equivalent to CS).
  double ym = 2.*sqrt(mi2)*sqrt(mj2)/sbar;
  double yp = 1. - 2.*sqrt(mk2)*sqrt(sqrRootQijkMk) / sbar;
  if ( y < ym || y > yp ) {
    jacobian(0.0);
    return false;
  }

  // Virtuality of emitted pair and other invariant scale
  Energy2 Qij2 = (pt2 + (1.-z)*mi2 + z*mj2) / z / (1.-z);
  Energy2 sijk = 0.5*( Qijk - mij2 - Mk2 + rootOfKallen(Qijk,mij2,mk2) );
  
  // Calculate xk and xij
  double lambdaIJ = 1. + (mij2/sijk);
  double lambdaK = 1. + (mk2/sijk);
  double fac1 = lambdaIJ*lambdaK + (mk2 - Qij2)/sijk;
  double xk =
    ( fac1 + sqrt( sqr(fac1) - 4.*lambdaIJ*lambdaK*mk2/sijk ) )
    / 2. / lambdaK ;
  double xij = 1. - mk2*(1.-xk) / xk / sijk;

  // Calculate zi
  double zi =
    ( z*xij*xk*sijk + mk2*(pt2+mi2) / (z*xij*xk*sijk) )
    / (1.-y) / sbar;

  // Limits on zi
  double facA = (2.*mi2 + sbar*y) / 2. / (mi2 + mj2 + sbar*y);
  // viji*vijk
  double facB =
    sqrt( (sqr(2.*mk2 + sbar*(1.-y)) - 4.*mk2*Qijk) *
          (sqr(sbar)*sqr(y) - 4.*mi2*mj2))
    / sbar / (1.-y) / (sbar*y + 2.*mi2);
  double zim = facA * (1. - facB);
  double zip = facA * (1. + facB);

  if ( zi < zim || zi > zip ) {
    jacobian(0.0);
    return false;
  }
  
  double mapZJacobian;
  if ( info.index().emitterData()->id() == ParticleID::g ) {
    if ( info.emissionData()->id() != ParticleID::g ) {
      mapZJacobian = 1.;
    } 
    else {
      mapZJacobian = z*(1.-z);
    }
  } 
  else {
    mapZJacobian = 1.-z;
  }

  // Compute and store the jacobian
  double jac = 0.0;
  jac = sbar / rootOfKallen(Qijk,mij2,mk2) * (1.-y) / ( 1. + (mi2 + mj2 - mij2)/sbar/y )
    * (pt2 / (pt2 + sqr(1.-z)*mi2+sqr(z)*mj2))
    * abs(1. - 2.*mk2*Qij2 / (sbar*(1.-y)*xij*xk*sijk));

  jacobian(jac * mapZJacobian * 2. * log(0.5 * generator()->maximumCMEnergy()/IRCutoff()) );
  
  // Record the physical variables, as used by the CS kernel definitions
  double phi = 2.*Constants::pi*rphi;
   
  lastPt(pt);
  lastZ(z);
  lastPhi(phi);

  // Record zi for use in kinematics generation and kernel evaluation
  splittingParameters().clear();
  splittingParameters().push_back(zi);
  
  if ( theMCCheck ) {
    theMCCheck->book(1.,1.,info.scale(),info.hardPt(),pt,z,jacobian());
  }
  return true;
  
}

// SW, 14/02/2019: Tidied to match thesis
void FIMassiveDecayKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
                                                  const Lorentz5Momentum& pSpectator,
                                                  const DipoleSplittingInfo& dInfo) {

  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();
  Energy2 pt2 = sqr(pt);

  // Momentum of the recoil system
  Lorentz5Momentum pk = pSpectator - pEmitter;
  Lorentz5Momentum pij = pEmitter;

  // scaled masses
  Energy2 mij2 = sqr(dInfo.emitterMass());
  Energy2 Mk2 = sqr(dInfo.recoilMass());
  Energy2 mk2 = Mk2;

  Energy2 mi2 = ZERO;
  // g->gg and g->qqbar
  if ( abs(dInfo.emitterData()->id()) == abs(dInfo.emissionData()->id()) ) {
    mi2 = sqr(dInfo.emitterData()->mass());
  }
  // Otherwise have X->Xg (should work for SUSY)
  else {
    mi2 = mij2;
  }
  Energy2 mj2 = sqr(dInfo.emissionData()->mass());

  // Scales
  Energy2 Qijk = sqr(dInfo.scale());
  Energy2 Qij2 = (pt2 + (1.-z)*mi2 + z*mj2) / z / (1.-z);
  Energy2 sijk = 0.5*( Qijk - mij2 - mk2 + rootOfKallen(Qijk,mij2,mk2) );
  Energy4 sijk2 = sqr(sijk);

  // Calculate xk and xij
  double lambdaIJ = 1. + (mij2/sijk);
  double lambdaK = 1. + (mk2/sijk);
  double fac1 = lambdaIJ*lambdaK + (mk2 - Qij2)/sijk;
  double xk =
    ( fac1 + sqrt( sqr(fac1) - 4.*lambdaIJ*lambdaK*mk2/sijk ) )
    / 2. / lambdaK ;
  double xij = 1. - mk2*(1.-xk) / xk / sijk;
  
  // Construct reference momenta nk and nij
  Lorentz5Momentum nij = ( sijk2 / (sijk2-mij2*Mk2) ) * (pij - (mij2/sijk)*pk);
  Lorentz5Momentum nk = ( sijk2 / (sijk2-mij2*Mk2) ) * (pk - (Mk2/sijk)*pij);

  // Construct qij, qk, qi and qj
  Lorentz5Momentum qij = xij*nij + (mij2/(xij*sijk))*nk;
  Lorentz5Momentum qk = xk*nk + (Mk2/(xk*sijk))*nij;
  
  Lorentz5Momentum qt = getKt(pij, pk, pt, dInfo.lastPhi());

  // No need to actually calculate nt and wt:
  Lorentz5Momentum qi = z*qij + ((pt2 + mi2 - z*z*mij2)/(xij*sijk*z))*nk + qt;
  Lorentz5Momentum qj = (1.-z)*qij + ((pt2 + mj2 - sqr(1.-z)*mij2)/(xij*sijk*(1.-z)))*nk - qt;

  qi.setMass(sqrt(mi2));
  qi.rescaleEnergy();
  
  qj.setMass(sqrt(mj2));
  qj.rescaleEnergy();

  emitterMomentum(qi);
  emissionMomentum(qj);
  spectatorMomentum(pSpectator);

  // Required for absorbing recoil in DipoleEventRecord::update
  splitRecoilMomentum(qk);

}

Lorentz5Momentum FIMassiveDecayKinematics::nVector(const Lorentz5Momentum& pEmitter, const Lorentz5Momentum& pSpectator, const DipoleSplittingInfo&) const {
  return (pSpectator-pEmitter);
}
  
// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void FIMassiveDecayKinematics::persistentOutput(PersistentOStream &) const {
  //os << ;
}

void FIMassiveDecayKinematics::persistentInput(PersistentIStream &, int) {
  //is >> ;
}

ClassDescription<FIMassiveDecayKinematics> FIMassiveDecayKinematics::initFIMassiveDecayKinematics;
// Definition of the static class description member.

void FIMassiveDecayKinematics::Init() {

  static ClassDocumentation<FIMassiveDecayKinematics> documentation
    ("FIMassiveDecayKinematics implements implements massive splittings "
     "off a final-initial decay dipole.");

}
