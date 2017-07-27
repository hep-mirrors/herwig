// -*- C++ -*-
//
// FFMassiveKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMassiveKinematics class.
//

#include "FFMassiveKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"
#include "Herwig/Shower/Dipole/Kernels/DipoleSplittingKernel.h"

#include "ThePEG/Interface/Switch.h"

// TODO: remove after verification
// only for checking for NaN or inf
#include <gsl/gsl_math.h>

using namespace Herwig;

FFMassiveKinematics::FFMassiveKinematics() 
  : DipoleSplittingKinematics(),
    theFullJacobian(true) {}

FFMassiveKinematics::~FFMassiveKinematics() {}

IBPtr FFMassiveKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFMassiveKinematics::fullclone() const {
  return new_ptr(*this);
}

pair<double,double> FFMassiveKinematics::kappaSupport(const DipoleSplittingInfo&) const {
  return {0.0,1.0};
}

pair<double,double> FFMassiveKinematics::xiSupport(const DipoleSplittingInfo& split) const {

  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));
  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emissionData()->id() != ParticleID::g )
      return {0.5*(1.-c),0.5*(1.+c)};
    double b = log((1.+c)/(1.-c));
    return {-b,b};
  }
  return {-log(0.5*(1.+c)),-log(0.5*(1.-c))};
}

Energy FFMassiveKinematics::dipoleScale(const Lorentz5Momentum& pEmitter,
					const Lorentz5Momentum& pSpectator) const {
  return (pEmitter+pSpectator).m();
}

Energy FFMassiveKinematics::ptMax(Energy dScale, 
				  double, double,
				  const DipoleIndex& ind,
				  const DipoleSplittingKernel& split) const {

  double mui = split.emitter(ind)->mass() / dScale;
  double mu  = split.emission(ind)->mass() / dScale;
  double muj = split.spectator(ind)->mass() / dScale;
  double mui2 = sqr( mui ), mu2  = sqr( mu );
  return rootOfKallen( mui2, mu2, sqr(1.-muj) ) / ( 2.-2.*muj ) * dScale;
}

Energy FFMassiveKinematics::QMax(Energy dScale, 
				 double, double,
				 const DipoleIndex& ind,
				 const DipoleSplittingKernel&) const {
  assert(false && "implementation missing");
  double Muj = ind.spectatorData()->mass() / dScale;
  return dScale * ( 1.-2.*Muj+sqr(Muj) );
}

// The name of this function is misleading
// scale here is defined as sqr(scale) = sqr(qi+qj)
// Here, scale is Q
Energy FFMassiveKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  double zPrime=split.lastSplittingParameters()[0];
  Energy mi = split.emitterData()->mass();
  Energy m = split.emissionData()->mass();
  Energy2 pt2 = zPrime*(1.-zPrime)*sqr(scale) - (1-zPrime)*sqr(mi) - zPrime*sqr(m);  
  assert(pt2 >= ZERO);
  return sqrt(pt2);
}

// This is simply the inverse of PtFromQ
Energy FFMassiveKinematics::QFromPt(Energy pt, const DipoleSplittingInfo& split) const {
  double zPrime=split.lastSplittingParameters()[0];
  Energy mi = split.emitterData()->mass();
  Energy m = split.emissionData()->mass();
  Energy2 Q2 = (sqr(pt) + (1-zPrime)*sqr(mi) + zPrime*sqr(m))/(zPrime*(1.-zPrime));
  return sqrt(Q2);
}

double FFMassiveKinematics::ptToRandom(Energy pt, Energy,
				       double,double,
				       const DipoleIndex&,
				       const DipoleSplittingKernel&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

// TODO - SW: Implement if statement to check 
// for spectator mass when tidying up the notation.
bool FFMassiveKinematics::generateSplitting(double kappa, double xi, double rphi,
					    DipoleSplittingInfo& info,
					    const DipoleSplittingKernel&) {

  // scaled masses
  double mui2 = sqr( info.emitterData()->mass() / info.scale() );
  double mu2 = sqr( info.emissionData()->mass() / info.scale() );
  double muj2 = sqr( info.spectatorData()->mass() / info.scale() );

  double Mui2 = 0.;
  if ( info.emitterData()->id() + info.emissionData()->id() == 0 ) Mui2 = 0.; // gluon
  else Mui2   = mui2; // (anti)quark 
  double Muj2 = muj2;

  // To solve issue with scale during presampling
  // need to enforce that Qijk-mij2-mk2 = 2*pij.pk > 0,
  // so combine checks by comparing against square root.
  if ( 1.-Mui2-Muj2 < sqrt(4.*Mui2*Muj2) ) {
    jacobian(0.0);
    return false;
  }

  Energy2 Qijk = sqr(info.scale());
  double suijk = 0.5*( 1. - Mui2 - Muj2 + sqrt( sqr(1.-Mui2-Muj2) - 4.*Mui2*Muj2 ) );

  // Calculate pt 
  Energy pt = IRCutoff() * pow(0.5 * generator()->maximumCMEnergy()/IRCutoff(),kappa);
  Energy2 pt2 = sqr(pt);

  if ( pt > info.hardPt() || pt < IRCutoff() ) {
    jacobian(0.0);
    return false;
  }

  // Generate zPrime (i.e. the new definition of z specific to massive FF)
  double zPrime;
  if ( info.index().emitterData()->id() == ParticleID::g ) {
    if ( info.emissionData()->id() != ParticleID::g ) {
      zPrime = xi;
    } 
    else {
      zPrime = exp(xi)/(1.+exp(xi));
    }
  } 
  else {
    zPrime = 1.-exp(-xi);
  }

  // new: 2011-08-31
  // 2011-11-08: this does happen
  if( sqrt(mui2)+sqrt(mu2)+sqrt(muj2) > 1. ){
    jacobian(0.0);
    return false;
  }

  // These apply to zPrime
  // phasespace constraint to incorporate ptMax
  Energy hard = info.hardPt();

  if(openZBoundaries()>0){
	// From ptMax(..)
	hard = rootOfKallen( mui2, mu2, sqr(1.-sqrt(muj2)) ) /
	       ( 2.-2.*sqrt(muj2) ) * info.scale();
	assert(pt<=hard);
  }

  double ptRatio = sqrt(1.-sqr(pt/hard));  
  double zp1 = ( 1.+mui2-mu2+muj2-2.*sqrt(muj2) +
		 rootOfKallen(mui2,mu2,sqr(1-sqrt(muj2))) * ptRatio) /
    ( 2.*sqr(1.-sqrt(muj2)) );
  double zm1 = ( 1.+mui2-mu2+muj2-2.*sqrt(muj2) -
		 rootOfKallen(mui2,mu2,sqr(1-sqrt(muj2))) * ptRatio) /
    ( 2.*sqr(1.-sqrt(muj2)) );

  if ( zPrime > zp1 || zPrime < zm1 ) {
    jacobian(0.0);
    return false;
  }
    
  // Calculate A:=xij*w
  double A = (1./(suijk*zPrime*(1.-zPrime))) * ( pt2/Qijk + zPrime*mu2 + (1.-zPrime)*mui2 - zPrime*(1.-zPrime)*Mui2 );

  // Calculate y from A (can also write explicitly in terms of qt, zPrime and masses however we need A anyway)
  double bar = 1.-mui2-mu2-muj2;
  double y = (1./bar) * (A*suijk + Mui2 - mui2 - mu2 );

  // kinematic phasespace boundaries for y
  // same as in Dittmaier hep-ph/9904440v2 (equivalent to CS)
  double ym = 2.*sqrt(mui2)*sqrt(mu2)/bar;
  double yp = 1. - 2.*sqrt(muj2)*(1.-sqrt(muj2))/bar;
  if ( y < ym || y > yp ) {
    jacobian(0.0);
    return false;
  }

  // Calculate xk and xij
  double lambdaK = 1. + (Muj2/suijk);
  double lambdaIJ = 1. + (Mui2/suijk);
  double xk = (1./(2.*lambdaK)) * ( (lambdaK + (Muj2/suijk)*lambdaIJ - A) + sqrt( sqr(lambdaK + (Muj2/suijk)*lambdaIJ - A) - 4.*lambdaK*lambdaIJ*Muj2/suijk) );
  double xij = 1. - ( (Muj2/suijk) * (1.-xk) / xk );

  // Transform to standard z definition as used in the kernels (i.e. that used in CS and standard sudakov parametrisations)
  double z = ( (zPrime*xij*xk*suijk/2.) + (Muj2/ ( 2.*xk*xij*suijk*zPrime))*(pt2/Qijk + mui2) ) /
    ( (xij*xk*suijk/2.) + (Muj2/(2.*xk*xij))*(Mui2/suijk + A) );

  // These apply to z, not zPrime
  double zm = ( (2.*mui2+bar*y)*(1.-y) - sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );
  double zp = ( (2.*mui2+bar*y)*(1.-y) + sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );

  if ( z < zm || z > zp ) {
    jacobian(0.0);
    return false;
  }

  double phi = 2.*Constants::pi*rphi;

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
  double propCntrb = 1./ ( 1. + (mui2+mu2-Mui2)/(bar*y) );  

  // TODO - SW - Tidy up notation everywhere alongside writing for the manual
  // The full jacobian including the z->zprime jacobian  
  if ( theFullJacobian && info.spectatorData()->mass() != ZERO ) {

    // Sort out variables
    Energy2 sbar = Qijk*bar;
    Energy2 sijk = Qijk*suijk;
    Energy2 mi2 = Qijk*mui2;
    Energy2 mj2 = Qijk*mu2;
    Energy2 mk2 = Qijk*muj2;    
    
    // Compute dy/dzPrime and pt2* dy/dpt2
    double dyBydzPrime = (1./sbar) * ( -pt2*(1.-2.*zPrime)/sqr(zPrime*(1.-zPrime)) - mi2/sqr(zPrime) + mj2/sqr(1.-zPrime) );
    InvEnergy2 dyBydpt2 = 1./(sbar*zPrime*(1.-zPrime));

    // Compute dA/dzPrime and dA/dpt2
    double dABydzPrime = (sbar/sijk) * dyBydzPrime;
    InvEnergy2 dABydpt2 = (sbar/sijk) * dyBydpt2;

    // Compute dxk/dzPrime, dxk/dpt2, dxij/dzPrime and dxij/dpt2
    double factor = (0.5/lambdaK) * (-1. - (1./sqrt( sqr(lambdaK + (mk2/sijk)*lambdaIJ - A) - 4.*lambdaK*lambdaIJ*mk2/sijk)) * (lambdaK + (mk2/sijk)*lambdaIJ - A));
    
    double dxkBydzPrime = factor * dABydzPrime;
    InvEnergy2 dxkBydpt2 = factor * dABydpt2;

    double dxijBydzPrime = (mk2/sijk) * (1./sqr(xk)) * dxkBydzPrime;
    InvEnergy2 dxijBydpt2 = (mk2/sijk) * (1./sqr(xk)) * dxkBydpt2;

    Energy2 dqiDotqkBydzPrime = xij*xk*0.5*sijk + zPrime*dxijBydzPrime*xk*0.5*sijk + zPrime*xij*dxkBydzPrime*0.5*sijk
      + 0.5*(mk2/sijk)*(pt2 + mi2) * (-1./(xk*xij*sqr(zPrime)) - dxkBydzPrime/(zPrime*xij*sqr(xk)) - dxijBydzPrime/(zPrime*xk*sqr(xij)));

    double dqiDotqkBydpt2 =  dxijBydpt2*zPrime*xk*0.5*sijk + zPrime*xij*dxkBydpt2*0.5*sijk
      + (0.5*mk2/sijk) * (1./(zPrime*xk*xij)) * (1. + (pt2+mi2)*(-dxkBydpt2/xk - dxijBydpt2/xij) );

    
    // Compute dzBydzPrime and dzBydpt2
    Energy2 qiDotqk = (zPrime*xij*xk*sijk*0.5) + (mk2/ ( 2.*xk*xij*sijk*zPrime))*(pt2 + mi2);

    double dzBydzPrime = (1./sbar) * ( 2.*qiDotqk*dyBydzPrime/sqr(1.-y) + (1./(1.-y)) * 2.*dqiDotqkBydzPrime );
    InvEnergy2 dzBydpt2    = (1./sbar) * ( 2.*qiDotqk*dyBydpt2/sqr(1.-y) + (1./(1.-y)) * 2.*dqiDotqkBydpt2 );

    double pt2Jac = pt2*abs(dzBydpt2*dyBydzPrime - dzBydzPrime*dyBydpt2);

    // Include the other terms and calculate the jacobian
    jac = propCntrb * bar / rootOfKallen(1.,Mui2,Muj2) * (1.-y)/y * pt2Jac;
  }


  // the exact result when the spectator is massless.
  else {
    double jacPt2 = 1. / ( 1. + sqr(1.-zPrime)*Qijk*mui2/pt2 + zPrime*zPrime*Qijk*mu2/pt2 );
    jac = propCntrb * bar / rootOfKallen(1.,Mui2,Muj2) * (1.-y) * jacPt2;    
  }

  jacobian(jac * mapZJacobian * 2. * log(0.5 * generator()->maximumCMEnergy()/IRCutoff()) );

  // Record the physical variables, as used by the CS kernel definitions
  lastPt(pt);
  lastZ(z);
  lastPhi(phi);

  // Record zPrime for use in kinematics generation and kernel evaluation
  splittingParameters().clear();
  splittingParameters().push_back(zPrime);

  if ( theMCCheck ) {
    theMCCheck->book(1.,1.,info.scale(),info.hardPt(),pt,z,jacobian());
  }
  return true;

}

// revised 2011-08-22
// revised 2011-11-06
// revised with new FF kinematics in 2016 - SW
void FFMassiveKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					     const Lorentz5Momentum& pSpectator,
					     const DipoleSplittingInfo& dInfo) {

  // The only value stored in dInfo.lastSplittingParameters() should be zPrime
  assert(dInfo.lastSplittingParameters().size() == 1 );
  double zPrime = dInfo.lastSplittingParameters()[0];
  Energy pt = dInfo.lastPt();
  Energy2 pt2 = sqr(pt);

  // scaled masses
  double mui2 = sqr( dInfo.emitterData()->mass() / dInfo.scale() );
  double mu2 = sqr( dInfo.emissionData()->mass() / dInfo.scale() );
  double muj2 = sqr( dInfo.spectatorData()->mass() / dInfo.scale() );

  double Mui2 = 0.;
  if ( dInfo.emitterData()->id() + dInfo.emissionData()->id() == 0 ) Mui2 = 0.; // gluon
  else Mui2   = mui2; // (anti)quark 
  double Muj2 = muj2;

  Energy2 Qijk = sqr(dInfo.scale());
  double suijk = 0.5*( 1. - Mui2 - Muj2 + sqrt( sqr(1.-Mui2-Muj2) - 4.*Mui2*Muj2 ) );
  double suijk2 = sqr(suijk);

  // Calculate A:=xij*w
  double A = (1./(suijk*zPrime*(1.-zPrime))) * ( pt2/Qijk + zPrime*mu2 + (1.-zPrime)*mui2 - zPrime*(1.-zPrime)*Mui2 );

  // Calculate the scaling factors, xk and xij
  double lambdaK = 1. + (Muj2/suijk);
  double lambdaIJ = 1. + (Mui2/suijk);
  double xk = (1./(2.*lambdaK)) * ( (lambdaK + (Muj2/suijk)*lambdaIJ - A) + sqrt( sqr(lambdaK + (Muj2/suijk)*lambdaIJ - A) - 4.*lambdaK*lambdaIJ*Muj2/suijk) );
  double xij = 1. - ( (Muj2/suijk) * (1.-xk) / xk );

  // Construct reference momenta nk, nij, nt
  Lorentz5Momentum nij = ( suijk2 / (suijk2-Mui2*Muj2) ) * (pEmitter - (Mui2/suijk)*pSpectator);
  Lorentz5Momentum nk = ( suijk2 / (suijk2-Mui2*Muj2) ) * (pSpectator - (Muj2/suijk)*pEmitter);

  // Following notation in notes, qt = sqrt(wt)*nt
  Lorentz5Momentum qt = getKt(nij,nk,pt,dInfo.lastPhi());

  // Construct qij, qk, qi and qj
  Lorentz5Momentum qij = xij*nij + (Mui2/(xij*suijk))*nk;
  Lorentz5Momentum qk = xk*nk + (Muj2/(xk*suijk))*nij;

  // For clarity, following notation in notes:
  //Lorentz5Momentum qi = zPrime*qij + ((pt2 + mi2 - zPrime*zPrime*mij2)/(xij*sijk*zPrime))*nk + sqrt(wt)*nt;
  //Lorentz5Momentum qj = (1.-zPrime)*qij + ((pt2 + mj2 - sqr(1.-zPrime)*mij2)/(xij*sijk*(1.-zPrime)))*nk - sqrt(wt)*nt;

  // No need to actually calculate nt and wt:
  Lorentz5Momentum qi = zPrime*qij + ((pt2/Qijk + mui2 - zPrime*zPrime*Mui2)/(xij*suijk*zPrime))*nk + qt;
  Lorentz5Momentum qj = (1.-zPrime)*qij + ((pt2/Qijk + mu2 - sqr(1.-zPrime)*Mui2)/(xij*suijk*(1.-zPrime)))*nk - qt;

  emitterMomentum(qi);
  emissionMomentum(qj);
  spectatorMomentum(qk);
  
}


// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void FFMassiveKinematics::persistentOutput(PersistentOStream & os) const {
  os << theFullJacobian;
}

void FFMassiveKinematics::persistentInput(PersistentIStream & is, int) {
  is >> theFullJacobian;
}

ClassDescription<FFMassiveKinematics> FFMassiveKinematics::initFFMassiveKinematics;
// Definition of the static class description member.

void FFMassiveKinematics::Init() {

  static ClassDocumentation<FFMassiveKinematics> documentation
    ("FFMassiveKinematics implements massive splittings "
     "off a final-final dipole.");

  
  static Switch<FFMassiveKinematics,bool> interfaceFullJacobian
    ("FullJacobian",
     "Use the full jacobian expression for the FF kinematics.",
     &FFMassiveKinematics::theFullJacobian, true, false, false);
  static SwitchOption interfaceFullJacobianYes
    (interfaceFullJacobian,
     "Yes",
     "Use the full jacobian.",
     true);
  static SwitchOption interfaceFullJacobianNo
    (interfaceFullJacobian,
     "No",
     "Do not use the full jacobian.",
     false);
  interfaceFullJacobian.rank(-1);
}

