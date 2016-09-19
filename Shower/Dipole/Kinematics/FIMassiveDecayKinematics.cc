// -*- C++ -*-
//
// FIMassiveDecayKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
  return make_pair(0.0,1.0);
}

pair<double,double> FIMassiveDecayKinematics::xiSupport(const DipoleSplittingInfo& split) const {

  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));

  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emissionData()->id() != ParticleID::g ){
      return make_pair(0.5*(1.-c),0.5*(1.+c));}
    double b = log((1.+c)/(1.-c));
    return make_pair(-b,b);
  }
  return make_pair(-log(0.5*(1.+c)),-log(0.5*(1.-c)));
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

// Look at derivation of this and make sure it still applies
Energy FIMassiveDecayKinematics::ptMax(Energy dScale, 
				 double, double,
				 const DipoleSplittingInfo& dInfo,
				 const DipoleSplittingKernel& split) const {

  DipoleIndex ind = dInfo.index();
  double mui = split.emitter(ind)->mass() / dScale;
  double mu  = split.emission(ind)->mass() / dScale;

  // Mass of recoil system
  // Use abs() due to generation of negative
  // recoilMass during sampling.

  if (dInfo.recoilMass() <= ZERO) {
    std::cerr << "Test triggered in FIMassiveDecayKinematics\n";
  }

  double muj = abs(dInfo.recoilMass() / dScale); 

  double mui2 = sqr( mui ), mu2  = sqr( mu );

  return rootOfKallen( mui2, mu2, sqr(1.-muj) ) / ( 2.-2.*muj ) * dScale;
}

// Where does this come from? - Ask Simon
Energy FIMassiveDecayKinematics::QMax(Energy dScale, 
				double, double,
				const DipoleSplittingInfo& dInfo,
				const DipoleSplittingKernel&) const {
  assert(false && "implementation missing");
  // Mass of recoil system
  double Muj = abs(dInfo.recoilMass() / dScale ); 

  return dScale * ( 1.-2.*Muj+sqr(Muj) );
}

Energy FIMassiveDecayKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  // from Martin's thesis
  double zPrime = split.lastSplittingParameters()[0];
  Energy mi = split.emitterData()->mass();
  Energy m = split.emissionData()->mass();
  Energy2 pt2 = zPrime*(1.-zPrime)*sqr(scale) - (1-zPrime)*sqr(mi) - zPrime*sqr(m);
  assert(pt2 >= ZERO);
  return sqrt(pt2);
}

Energy FIMassiveDecayKinematics::QFromPt(Energy scale, const DipoleSplittingInfo& split) const {
  // from Martin's thesis
  double zPrime = split.lastSplittingParameters()[0];
  Energy mi = split.emitterData()->mass();
  Energy m = split.emissionData()->mass();
  Energy2 Q2 = (sqr(scale) + (1-zPrime)*sqr(mi) + zPrime*sqr(m))/(zPrime*(1.-zPrime));
  return sqrt(Q2);
}

double FIMassiveDecayKinematics::ptToRandom(Energy pt, Energy,
				       double,double,
				       const DipoleIndex&,
				       const DipoleSplittingKernel&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}


bool FIMassiveDecayKinematics::generateSplitting(double kappa, double xi, double rphi,
					    DipoleSplittingInfo& info,
					    const DipoleSplittingKernel&) {

  // Construct mass squared variables
  Energy2 mi2 = sqr( info.emitterData()->mass() );
  Energy2 mj2 = sqr( info.emissionData()->mass() );

  Energy2 mij2 = ZERO;
  if ( info.emitterData()->id() + info.emissionData()->id() == 0 ) mij2 = ZERO; // g -> qq
  else mij2 = mi2; // g -> gg or q -> qg

  Energy2 mk2 = sqr( info.recoilMass() );  
  Energy2 Qijk = sqr(info.scale());

  // To solve issue with scale during presampling
  // need to enforce that Qijk-mij2-mk2 = 2*pij.pk > 0,
  // so combine checks by comparing against square root.
  if ( Qijk-mij2-mk2 < sqrt(4.*mij2*mk2) ) {
    jacobian(0.0);
    return false;
  }

  Energy2 sijk = 0.5*( Qijk - mij2 - mk2 + sqrt( sqr(Qijk-mij2-mk2) - 4.*mij2*mk2 ) );

  // Calculate pt
  Energy pt = IRCutoff() * pow(0.5 * generator()->maximumCMEnergy()/IRCutoff(),kappa);
  Energy2 pt2 = sqr(pt);

  if ( pt > info.hardPt() || pt < IRCutoff() ) {
    jacobian(0.0);
    return false;
  }

  // Generate zPrime (i.e. the new definition of z specific to massive FF and decays)
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

  // scaled masses *** TODO: rewrite above in terms of mu to avoid this calculation ***, and mix up of notation
  double mui2 = mi2 / Qijk;
  double mu2  = mj2 / Qijk;
  double muj2 = mk2 / Qijk;

  double Mui2 = mij2 / Qijk;
  double Muj2 = muj2;

  // Check limit on pt
  Energy ptmax1 = rootOfKallen( mui2, mu2, sqr(1.-sqrt(muj2)) ) /
    ( 2.-2.*sqrt(muj2) ) * info.scale();
  Energy auxHardPt = ptmax1 > info.hardPt() ? info.hardPt() : ptmax1;

  // 24/05/2015: Moved this check from the zPrime limit checks
  if ( pt > auxHardPt ){
    jacobian(0.0);
    return false;
  }

  // 2011-11-09
  //assert(ptmax1>info.hardPt());
  // 24/05/2015:
  // The simple >= assert above is triggered 
  // during sampling due to precision.
  // Have added a tolerance to deal with this.
  assert( abs(ptmax1 - info.hardPt()) <= 1e-8*GeV || ptmax1>=info.hardPt() );
  
  // new: 2011-08-31
  // 2011-11-08: this does happen
  if( sqrt(mui2)+sqrt(mu2)+sqrt(muj2) > 1. ){
    jacobian(0.0);
    return false;
  }
  
  // Have derived and checked the equations for zp1 and zm1, these apply to zPrime
  // phasespace constraint to incorporate ptMax
  double zp1 = ( 1.+mui2-mu2+muj2-2.*sqrt(muj2) +
		 rootOfKallen(mui2,mu2,sqr(1-sqrt(muj2))) *
		 sqrt( 1.-sqr(pt/auxHardPt) ) ) /
    ( 2.*sqr(1.-sqrt(muj2)) );
  double zm1 = ( 1.+mui2-mu2+muj2-2.*sqrt(muj2) -
		 rootOfKallen(mui2,mu2,sqr(1-sqrt(muj2))) *
		 sqrt( 1.-sqr(pt/auxHardPt) ) ) /
    ( 2.*sqr(1.-sqrt(muj2)) );
  
  if ( zPrime > zp1 || zPrime < zm1 ) {
    jacobian(0.0);
    return false;
  }

  // Calculate A:=xij*w
  double A = (1./(sijk*zPrime*(1.-zPrime))) * ( pt2 + zPrime*mj2 + (1.-zPrime)*mi2 - zPrime*(1.-zPrime)*mij2 );

  // Calculate y from A (can also write explicitly in terms of qt, zPrime and masses however we need A anyway)
  Energy2 sbar = Qijk - mi2 - mj2 - mk2;
  double y = (1./sbar) * (A*sijk + mij2 - mi2 - mj2);

  // kinematic phasespace boundaries for y
  // same as in Dittmaier hep-ph/9904440v2 (equivalent to CS)
  double bar = 1.-mui2-mu2-muj2;
  double ym = 2.*sqrt(mui2)*sqrt(mu2)/bar;
  double yp = 1. - 2.*sqrt(muj2)*(1.-sqrt(muj2))/bar;
  if ( y < ym || y > yp ) {
    jacobian(0.0);
    return false;
  }

  // Calculate xk and xij
  double lambdaK = 1. + (mk2/sijk);
  double lambdaIJ = 1. + (mij2/sijk);

  double xk = (1./(2.*lambdaK)) * ( (lambdaK + (mk2/sijk)*lambdaIJ - A) + sqrt( sqr(lambdaK + (mk2/sijk)*lambdaIJ - A) - 4.*lambdaK*lambdaIJ*mk2/sijk) );
  double xij = 1. - ( (mk2/sijk) * (1.-xk) / xk );

  // Transform to standard z definition as used in the kernels (i.e. that used in CS and standard sudakov parametrisations)
  double z = 
    ( (zPrime*xij*xk*sijk/2.) + (mk2/ ( 2.*xk*xij*sijk*zPrime))*(pt2 + mi2) ) /
    ( (xij*xk*sijk/2.) + (mk2*mij2/(2.*xk*xij*sijk)) + (mk2/(2.*xk*xij))*A );

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

  // TODO: This may need changing due to different definitions of z
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

  // TODO: May need a redefinition due to different definitions of z
  jacobian( 2. * mapZJacobian * (1.-y) * 
	    log(0.5 * generator()->maximumCMEnergy()/IRCutoff()) *
	    bar / rootOfKallen(1.,Mui2,Muj2) );

  // Record the physical variables, as used by the CS kernel definitions
  lastPt(pt);
  lastZ(z);
  lastPhi(phi);

  // Record zPrime for use in kinematics generation
  splittingParameters().clear();
  splittingParameters().push_back(zPrime);

  if ( theMCCheck ) {
    theMCCheck->book(1.,1.,info.scale(),info.hardPt(),pt,z,jacobian());
  }
  return true;

}

void FIMassiveDecayKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					    const Lorentz5Momentum& pSpectator,
					    const DipoleSplittingInfo& dInfo) {


  // The only value stored in dInfo.lastSplittingParameters() should be zPrime
  assert(dInfo.lastSplittingParameters().size() == 1 );
  double zPrime = dInfo.lastSplittingParameters()[0];
  Energy pt = dInfo.lastPt();
  Energy2 pt2 = sqr(pt);

  // Momentum of the recoil system
  Lorentz5Momentum pk = pSpectator - pEmitter;
  Lorentz5Momentum pij = pEmitter;

  // Masses - Currently not using the mu-ratio formalism and using a different notation to Simon
  Energy2 mi2 = sqr( dInfo.emitterData()->mass() ); 
  Energy2 mj2 = sqr( dInfo.emissionData()->mass() );

  Energy2 mij2 = ZERO;
  if ( dInfo.emitterData()->id() + dInfo.emissionData()->id() == 0 ) mij2 = ZERO; // g -> qqbar
  else mij2 = mi2; // g -> gg or q -> qg
  Energy2 mk2 = sqr(dInfo.recoilMass());
    
  Energy2 Qijk = sqr(dInfo.scale());

  Energy2 sijk = 0.5*( Qijk - mij2 - mk2 + sqrt( sqr(Qijk-mij2-mk2) - 4.*mij2*mk2 ) );
  Energy4 sijk2 = sqr(sijk);

  // Calculate A:=xij*w
  double A = (1./(sijk*zPrime*(1.-zPrime))) * ( pt2 + zPrime*mj2 + (1.-zPrime)*mi2 - zPrime*(1.-zPrime)*mij2 );

  // Calculate xk and xij
  double lambdaK = 1. + (mk2/sijk);
  double lambdaIJ = 1. + (mij2/sijk);

  double xk = (1./(2.*lambdaK)) * ( (lambdaK + (mk2/sijk)*lambdaIJ - A) + sqrt( sqr(lambdaK + (mk2/sijk)*lambdaIJ - A) - 4.*lambdaK*lambdaIJ*mk2/sijk) );
  double xij = 1. - ( (mk2/sijk) * (1.-xk) / xk );

  // Construct reference momenta nk, nij, nt
  Lorentz5Momentum nij = ( sijk2 / (sijk2-mij2*mk2) ) * (pij - (mij2/sijk)*pk);
  Lorentz5Momentum nk = ( sijk2 / (sijk2-mij2*mk2) ) * (pk - (mk2/sijk)*pij);

  // Following notation in notes, qt = sqrt(wt)*nt
  Lorentz5Momentum qt = getKt(nij,nk,pt,dInfo.lastPhi());

  // Construct qij, qk, qi and qj
  Lorentz5Momentum qij = xij*nij + (mij2/(xij*sijk))*nk;
  Lorentz5Momentum qk = xk*nk + (mk2/(xk*sijk))*nij;

  // No need to actually calculate nt and wt:
  Lorentz5Momentum qi = zPrime*qij + ((pt2 + mi2 - zPrime*zPrime*mij2)/(xij*sijk*zPrime))*nk + qt;
  Lorentz5Momentum qj = (1.-zPrime)*qij + ((pt2 + mj2 - sqr(1.-zPrime)*mij2)/(xij*sijk*(1.-zPrime)))*nk - qt;

  emitterMomentum(qi);
  emissionMomentum(qj);
  spectatorMomentum(pSpectator);

  //recoilMomentum is not currently used
  //recoilMomentum(pk); 
  splitRecoilMomentum(qk);

}


  // If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMassiveDecayKinematics::persistentOutput(PersistentOStream & ) const {
}

void FIMassiveDecayKinematics::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FIMassiveDecayKinematics> FIMassiveDecayKinematics::initFIMassiveDecayKinematics;
// Definition of the static class description member.

void FIMassiveDecayKinematics::Init() {

  static ClassDocumentation<FIMassiveDecayKinematics> documentation
    ("FIMassiveDecayKinematics implements implements massive splittings "
     "off a final-initial decay dipole.");
}
