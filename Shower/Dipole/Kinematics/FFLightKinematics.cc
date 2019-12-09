// -*- C++ -*-
//
// FFLightKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"

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

Energy FFLightKinematics::ptMax(Energy dScale, 
				double, double,
				const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  return dScale/2.;
}

Energy FFLightKinematics::QMax(Energy dScale, 
			       double, double,
			       const DipoleIndex&,
			       const DipoleSplittingKernel&) const {
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

pair<double,double> FFLightKinematics::zBoundaries(Energy pt,
						   const DipoleSplittingInfo& dInfo,
						   const DipoleSplittingKernel&) const {
  Energy hard=dInfo.hardPt();
  if(openZBoundaries()>0)hard=dInfo.scale()/2.;
  const double s = sqrt(1.-sqr(pt/hard));
  return {0.5*(1.-s),0.5*(1.+s)};
}

bool FFLightKinematics::generateSplitting(double kappa, double xi, double rphi,
					  DipoleSplittingInfo& info,
					  const DipoleSplittingKernel& split) {

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

  double y = sqr(pt/info.scale())/(z*(1.-z));

  if ( z < 0.0 || z > 1.0 ||
       y < 0.0 || y > 1.0 ) {
    jacobian(0.0);
    return false;
  }

  double phi = 2.*Constants::pi*rphi;

  jacobian(weight*(1.-y));

  lastPt(pt);
  lastZ(z);
  lastPhi(phi);

  if ( theMCCheck )
    theMCCheck->book(1.,1.,info.scale(),info.hardPt(),pt,z,jacobian());

  return true;
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
  Lorentz5Momentum emm = (1.-z)*pEmitter + z*y*pSpectator - kt;
  Lorentz5Momentum spe = (1.-y)*pSpectator;
  
  em.setMass(ZERO);
  em.rescaleEnergy();

  emm.setMass(ZERO);
  emm.rescaleEnergy();
  
  spe.setMass(ZERO);
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

