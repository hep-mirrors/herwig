// -*- C++ -*-
//
// FFMassiveKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig/DipoleShower/Base/DipoleSplittingInfo.h"
#include "Herwig/DipoleShower/Kernels/DipoleSplittingKernel.h"

// TODO: remove after verification
// only for checking for NaN or inf
#include <gsl/gsl_math.h>

using namespace Herwig;

FFMassiveKinematics::FFMassiveKinematics() 
  : DipoleSplittingKinematics() {}

FFMassiveKinematics::~FFMassiveKinematics() {}

IBPtr FFMassiveKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFMassiveKinematics::fullclone() const {
  return new_ptr(*this);
}

pair<double,double> FFMassiveKinematics::kappaSupport(const DipoleSplittingInfo&) const {
  return make_pair(0.0,1.0);
}

pair<double,double> FFMassiveKinematics::xiSupport(const DipoleSplittingInfo& split) const {
  double c = sqrt(1.-4.*sqr(IRCutoff()/generator()->maximumCMEnergy()));
  if ( split.index().emitterData()->id() == ParticleID::g ) {
    if ( split.emissionData()->id() != ParticleID::g )
      return make_pair(0.5*(1.-c),0.5*(1.+c));
    double b = log((1.+c)/(1.-c));
    return make_pair(-b,b);
  }
  return make_pair(-log(0.5*(1.+c)),-log(0.5*(1.-c)));
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
  double mui2 = sqr( mui ), mu2  = sqr( mu ), muj2 = sqr( muj );
  return rootOfKallen( mui2, mu2, sqr(1.-muj) ) / ( 2.-2.*sqrt(muj2) ) * dScale;
}

Energy FFMassiveKinematics::QMax(Energy dScale, 
			       double, double,
			       const DipoleIndex& ind,
				const DipoleSplittingKernel&) const {
  assert(false && "implementation missing");
  double Muj = ind.spectatorData()->mass() / dScale;
  return dScale * ( 1.-2.*Muj+sqr(Muj) );
}

Energy FFMassiveKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  // from Martin's thesis
  double z = split.lastZ();
  Energy mi = split.emitterData()->mass();
  Energy m = split.emissionData()->mass();
  Energy2 pt2 = z*(1.-z)*sqr(scale) - (1-z)*sqr(mi) - z*sqr(m);
  assert(pt2 >= ZERO);
  return sqrt(pt2);
}

Energy FFMassiveKinematics::QFromPt(Energy scale, const DipoleSplittingInfo& split) const {
  // from Martin's thesis
  double z = split.lastZ();
  Energy mi = split.emitterData()->mass();
  Energy m = split.emissionData()->mass();
  Energy2 Q2 = (sqr(scale) + (1-z)*sqr(mi) + z*sqr(m))/(z*(1.-z));
  return sqrt(Q2);
}

double FFMassiveKinematics::ptToRandom(Energy pt, Energy,
				       double,double,
				       const DipoleIndex&,
				       const DipoleSplittingKernel&) const {
  return log(pt/IRCutoff()) / log(0.5 * generator()->maximumCMEnergy()/IRCutoff());
}

// own kinematics close to Dinsdale,Ternick,Weinzierl
bool FFMassiveKinematics::generateSplitting(double kappa, double xi, double rphi,
					    DipoleSplittingInfo& info,
					    const DipoleSplittingKernel&) {
  
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

  // masses
  double mui2 = sqr( info.emitterData()->mass() / info.scale() );
  double mu2  = sqr( info.emissionData()->mass() / info.scale() );
  double muj2 = sqr( info.spectatorData()->mass() / info.scale() );
  double Mui2 = 0.;
  if ( info.emitterData()->id() + info.emissionData()->id() == 0 ) Mui2 = 0.; // gluon
  else Mui2   = mui2; // (anti)quark 
  double Muj2 = muj2;
  
  // new: 2011-08-31
  // 2011-11-08: this does happen
  if( sqrt(mui2)+sqrt(mu2)+sqrt(muj2) > 1. ){
    jacobian(0.0);
    return false;
  }

  double bar = 1.-mui2-mu2-muj2;
  double y = ( sqr( pt / info.scale() ) + sqr(1.-z)*mui2 + z*z*mu2 ) /
    (z*(1.-z)*bar);
  
  // do this here for simplicity
  Energy ptmax1 = rootOfKallen( mui2, mu2, sqr(1.-sqrt(muj2)) ) /
    ( 2.-2.*sqrt(muj2) ) * info.scale();
  Energy auxHardPt = ptmax1 > info.hardPt() ? info.hardPt() : ptmax1;
  // 2011-11-09
  assert(ptmax1>=info.hardPt());

  // phasespace constraint to incorporate ptMax
  double zp1 = ( 1.+mui2-mu2+muj2-2.*sqrt(muj2) +
    rootOfKallen(mui2,mu2,sqr(1-sqrt(muj2))) *
    sqrt( 1.-sqr(pt/auxHardPt) ) ) /
    ( 2.*sqr(1.-sqrt(muj2)) );
  double zm1 = ( 1.+mui2-mu2+muj2-2.*sqrt(muj2) -
    rootOfKallen(mui2,mu2,sqr(1-sqrt(muj2))) *
    sqrt( 1.-sqr(pt/auxHardPt) ) ) /
    ( 2.*sqr(1.-sqrt(muj2)) );
  if ( z > zp1 || z < zm1 ||
       pt > auxHardPt) {
    jacobian(0.0);
    return false;
  }
  
  // kinematic phasespace boundaries for (y,z)
  // same as in Dittmaier hep-ph/9904440v2 (equivalent to CS)
  double ym = 2.*sqrt(mui2)*sqrt(mu2)/bar;
  double yp = 1. - 2.*sqrt(muj2)*(1.-sqrt(muj2))/bar;
  if ( y < ym || y > yp ) {
    jacobian(0.0);
    return false;
  }
  
  double zm = ( (2.*mui2+bar*y)*(1.-y) - sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );
  double zp = ( (2.*mui2+bar*y)*(1.-y) + sqrt(y*y-ym*ym)*sqrt(sqr(2.*muj2+bar-bar*y)-4.*muj2) ) /
    ( 2.*(1.-y)*(mui2+mu2+bar*y) );
  
  if ( z < zm || z > zp ) {
    jacobian(0.0);
    return false;
  }

  double phi = 2.*Constants::pi*rphi;

  jacobian( 2. * mapZJacobian * (1.-y) * 
	    log(0.5 * generator()->maximumCMEnergy()/IRCutoff()) *
	    bar / rootOfKallen(1.,Mui2,Muj2) );

  lastPt(pt);
  lastZ(z);
  lastPhi(phi);

  if ( theMCCheck )
    theMCCheck->book(1.,1.,info.scale(),info.hardPt(),pt,z,jacobian());
  
  return true;

}

// kinematics close to Dinsdale,Ternick,Weinzierl
// revised 2011-08-22
// revised 2011-11-06
void FFMassiveKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  double z = dInfo.lastZ();
  Energy pt = dInfo.lastPt();

  // masses
  double mui2 = sqr( dInfo.emitterData()->mass() / dInfo.scale() );
  double mu2  = sqr( dInfo.emissionData()->mass() / dInfo.scale() );
  double muj2 = sqr( dInfo.spectatorData()->mass() / dInfo.scale() );

  double y = ( sqr( pt / dInfo.scale() ) + sqr(1.-z)*mui2 + z*z*mu2 ) /
      (z*(1.-z)*(1.-mui2-mu2-muj2));

  Energy2 sbar = sqr(dInfo.scale()) *(1.-mui2-mu2-muj2);

  // CMF: particle energies
  Energy Ei = ( sbar*(1.-(1.-z)*(1.-y)) + 2.*sqr(dInfo.scale())*mui2 ) / (2.*dInfo.scale());
  Energy E  = ( sbar*(1.-    z *(1.-y)) + 2.*sqr(dInfo.scale())*mu2  ) / (2.*dInfo.scale());
  Energy Ej = ( sbar*(1.-           y ) + 2.*sqr(dInfo.scale())*muj2 ) / (2.*dInfo.scale());
  // CMF: momenta in z-direction (axis of pEmitter & pSpectator)  
  Energy qi3 = (2.*Ei*Ej-z*(1.-y)*sbar     ) / 2./sqrt(Ej*Ej-sqr(dInfo.scale())*muj2);
  Energy q3  = (2.*E *Ej-(1.-z)*(1.-y)*sbar) / 2./sqrt(Ej*Ej-sqr(dInfo.scale())*muj2);
  Energy qj3 = sqrt(   sqr(Ej) - sqr(dInfo.scale())*muj2 );

  // get z axis in the dipole's CMF which is parallel to pSpectator
  Boost toCMF = (pEmitter+pSpectator).findBoostToCM();
  Lorentz5Momentum pjAux = pSpectator; pjAux.boost(toCMF);
  ThreeVector<double> pjAxis = pjAux.vect().unit();
  
  // set the momenta in this special reference frame
  // note that pt might in some cases differ from the physical pt!
  // phi is defined exactly as in getKt
  Energy ptResc = sqrt( sqr(Ei)-sqr(dInfo.scale())*mui2-sqr(qi3) );
  Lorentz5Momentum em  ( ptResc*cos(dInfo.lastPhi()), -ptResc*sin(dInfo.lastPhi()), qi3, Ei );
  Lorentz5Momentum emm ( -ptResc*cos(dInfo.lastPhi()), ptResc*sin(dInfo.lastPhi()), q3, E );
  Lorentz5Momentum spe ( 0.*GeV, 0.*GeV, qj3, Ej );
  
  // output the mismatch between pt and physical pt
//   ofstream output1("ptDiffOnPtAxis-uub-m.dat",ofstream::app);
//   ofstream output2("ptDiffOnCosAxis-uub-m.dat",ofstream::app);
//   if( abs(dInfo.spectatorData()->id())==5 && dInfo.emitterData()->id()+dInfo.emissionData()->id()==0 &&
//       abs(dInfo.emitterData()->id())==1 ) {
//     output1 << pt/dInfo.scale() << " " << abs(ptResc-pt)/(ptResc+pt) << " " << endl;
//     output2 << em.vect().unit()*emm.vect().unit() << " " << abs(ptResc-pt)/(ptResc+pt) << " " << endl;
//   }
//   output1.close(); output2.close();
  
  // rotate back
  em.rotateUz (pjAxis);
  emm.rotateUz(pjAxis);
  spe.rotateUz(pjAxis);

  // boost back
  em.boost (-toCMF);
  emm.boost(-toCMF);
  spe.boost(-toCMF);

  // mass shells, rescale energy
  em.setMass(dInfo.scale()*sqrt(mui2));
  em.rescaleEnergy();
  emm.setMass(dInfo.scale()*sqrt(mu2));
  emm.rescaleEnergy();
  spe.setMass(dInfo.scale()*sqrt(muj2));
  spe.rescaleEnergy();
  
  // book
  emitterMomentum(em);
  emissionMomentum(emm);
  spectatorMomentum(spe);
  
  // TODO: remove
  // 2011-11-09: never occurred
  if(em.t()/GeV>=0. && emm.t()/GeV>=0. && spe.t()/GeV>=0.);
  else cout << "FFMassiveKinematics::generateKinematics momenta corrupt" << endl;
  
  // 2011-11-03 LEP run with full masses:
  // x,y,t no problem
  // z order > 5.e-7 happend 41 times in 10000 runs
  // maximum was 2.5e-6
  
//   double order=2.4e-6;
//   Lorentz5Momentum pDif=em+emm+spe-pEmitter-pSpectator, pSum=em+emm+spe+pEmitter+pSpectator;
//   if(abs(pDif.x()/(pSum.x()==ZERO ? 1.*GeV : pSum.x())) <order || (pDif-pSum).x()==ZERO);
//   else cout << "FFMassiveKinematics::generateKinematics momenta corrupt: x  " <<
//     abs(pDif.x()/(pSum.x()==ZERO ? 1.*GeV : pSum.x())) << endl <<
//     "  " << em/GeV << "  " << emm/GeV << "  " << spe/GeV << endl <<
//     "  " << pEmitter/GeV << "  " << pSpectator/GeV << endl;
//   if(abs(pDif.y()/(pSum.y()==ZERO ? 1.*GeV : pSum.y())) <order || (pDif-pSum).y()==ZERO);
//   else cout << "FFMassiveKinematics::generateKinematics momenta corrupt: y  " <<
//     abs(pDif.y()/(pSum.y()==ZERO ? 1.*GeV : pSum.y())) << endl;
//   if(abs(pDif.z()/(pSum.z()==ZERO ? 1.*GeV : pSum.z())) <order | (pDif-pSum).z()==ZERO);
//   else cout << "FFMassiveKinematics::generateKinematics momenta corrupt: z  " <<
//     abs(pDif.z()/(pSum.z()==ZERO ? 1.*GeV : pSum.z())) << endl <<
//     "  " << em/GeV << "  " << emm/GeV << "  " << spe/GeV << endl <<
//     "  " << pEmitter/GeV << "  " << pSpectator/GeV << endl;
//   if(abs(pDif.t()/(pSum.t()==ZERO ? 1.*GeV : pSum.t())) <order || (pDif-pSum).t()==ZERO);
//   else cout << "FFMassiveKinematics::generateKinematics momenta corrupt: t  " <<
// //     abs(pDif.t()/(pSum.t()==ZERO ? 1.*GeV : pSum.t())) << endl;
//   cout << endl;

}

// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMassiveKinematics::persistentOutput(PersistentOStream & ) const {
}

void FFMassiveKinematics::persistentInput(PersistentIStream & , int) {
}

ClassDescription<FFMassiveKinematics> FFMassiveKinematics::initFFMassiveKinematics;
// Definition of the static class description member.

void FFMassiveKinematics::Init() {

  static ClassDocumentation<FFMassiveKinematics> documentation
    ("FFMassiveKinematics implements massive splittings "
     "off a final-final dipole.");

}

