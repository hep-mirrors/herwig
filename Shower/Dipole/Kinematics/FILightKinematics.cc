// -*- C++ -*-
//
// FILightKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"

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

Energy FILightKinematics::ptMax(Energy dScale, 
				double, double specX,
				const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  return dScale * sqrt((1.-specX)/specX) /2.;
}

Energy FILightKinematics::QMax(Energy dScale, 
			       double, double specX,
			       const DipoleIndex&,
			       const DipoleSplittingKernel&) const {
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

pair<double,double> FILightKinematics::zBoundaries(Energy pt,
						   const DipoleSplittingInfo& dInfo,
						   const DipoleSplittingKernel&) const {
  Energy hard=dInfo.hardPt();
  if(openZBoundaries()==1)
	hard=dInfo.scale()*sqrt((1.-dInfo.spectatorX())/dInfo.spectatorX()/2.);
  if(openZBoundaries()==2)
	hard=dInfo.scale()*min(1.,sqrt((1.-dInfo.spectatorX())/dInfo.spectatorX())/2.);
  if(hard<pt)return {0.5,0.5};
  double s = sqrt(1.-sqr(pt/hard));
  return {0.5*(1.-s),0.5*(1.+s)};
}

bool FILightKinematics::generateSplitting(double kappa, double xi, double rphi,
					  DipoleSplittingInfo& info,
					  const DipoleSplittingKernel& split) {

  if ( info.spectatorX() < xMin() ) {
    jacobian(0.0);
    return false;
  }

  double weight = 1.0;

  Energy pt = generatePt(kappa,info.scale(),
			 info.emitterX(),info.spectatorX(),
			 info.index(),split,
			 weight);

  if ( pt < IRCutoff() || pt > info.hardPt() ) {
    jacobian(0.0);
    return false;
  }

  double z = 0.0;

  if ( info.index().emitterData()->id() == ParticleID::g ) {
    if ( info.emissionData()->id() != ParticleID::g ) {
      z = generateZ(xi,pt,FlatZ,
		    info,split,weight);
    } else {
      z = generateZ(xi,pt,OneOverZOneMinusZ,
		    info,split,weight);
    }
  } else {
    z = generateZ(xi,pt,OneOverOneMinusZ,
		  info,split,weight);
  }

  double x = 1./(1.+sqr(pt/info.scale())/(z*(1.-z)));

  if ( z < 0.0 || z > 1.0 ||
       x < info.spectatorX() || x > 1.0 ) {
    jacobian(0.0);
    return false;
  }

  double phi = 2.*Constants::pi*rphi;

  jacobian(weight);

  lastPt(pt);
  lastZ(z);
  lastPhi(phi);
  lastSpectatorZ(x);

  if ( theMCCheck )
    theMCCheck->book(1.,info.spectatorX(),info.scale(),info.hardPt(),pt,z,jacobian());

  return true;

}

void FILightKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  Energy pt = dInfo.lastPt();
  double z = dInfo.lastZ();
  double x = 1./(1.+sqr(pt/dInfo.scale())/(z*(1.-z)));

  Lorentz5Momentum kt =
    getKt (pSpectator, pEmitter, pt, dInfo.lastPhi(),true);

  Lorentz5Momentum em = z*pEmitter + (1.-z)*((1.-x)/x)*pSpectator + kt;
  em.setMass(0.*GeV);
  em.rescaleEnergy();

  Lorentz5Momentum emm = (1.-z)*pEmitter + z*((1.-x)/x)*pSpectator - kt;
  emm.setMass(0.*GeV);
  emm.rescaleEnergy();

  Lorentz5Momentum spe = (1./x)*pSpectator;
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

