// -*- C++ -*-
//
// IFLightKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFLightKinematics class.
//

#include "IFLightKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"

using namespace Herwig;

IFLightKinematics::IFLightKinematics() 
  : DipoleSplittingKinematics(), theCollinearScheme(true) {}

IFLightKinematics::~IFLightKinematics() {}

IBPtr IFLightKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFLightKinematics::fullclone() const {
  return new_ptr(*this);
}

Energy IFLightKinematics::ptMax(Energy dScale, 
				double emX, double,
				const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  return dScale * sqrt((1.-emX)/emX) /2.;
}

Energy IFLightKinematics::QMax(Energy, 
			       double, double,
			       const DipoleIndex&,
			       const DipoleSplittingKernel&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy IFLightKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return scale*sqrt(1.-z);
}

Energy IFLightKinematics::QFromPt(Energy scale, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return scale/sqrt(1.-z);
}

pair<double,double> IFLightKinematics::zBoundaries(Energy pt,
						   const DipoleSplittingInfo& dInfo,
						   const DipoleSplittingKernel&) const {
  double x = dInfo.emitterX();

  Energy hard=dInfo.hardPt();
  if(openZBoundaries()==1)hard=dInfo.scale() * sqrt((1.-x)/x) /2.;
  if(openZBoundaries()==2)hard=dInfo.scale() * min(1.,sqrt((1.-x)/x) /2.);
  if(hard<pt)return {0.5*(1.+x),0.5*(1.+x)};

  double s = sqrt(1.-sqr(pt/hard));

  return {0.5*(1.+x-(1.-x)*s),0.5*(1.+x+(1.-x)*s)};
}


bool IFLightKinematics::generateSplitting(double kappa, double xi, double rphi,
					  DipoleSplittingInfo& info,
					  const DipoleSplittingKernel& split) {

  if ( info.emitterX() < xMin() ) {
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
    if ( info.emitterData()->id() == ParticleID::g ) {
      z = generateZ(xi,pt,OneOverZOneMinusZ,
		    info,split,weight);
    } else {
      z = generateZ(xi,pt,OneOverZ,
		    info,split,weight);
    }
  }

  if ( info.index().emitterData()->id() != ParticleID::g ) {
    if ( info.emitterData()->id() != ParticleID::g ) {
      z = generateZ(xi,pt,OneOverOneMinusZ,
		    info,split,weight);
    } else {
      z = generateZ(xi,pt,FlatZ,
		    info,split,weight);
    }
  }

  if ( weight == 0. && z == -1. ) {
    jacobian(0.0);
    return false;
  }

  double ratio = sqr(pt/info.scale());
  
  double rho = 1. - 4.*ratio*z*(1.-z)/sqr(1.-z+ratio);
  if ( rho < 0.0 ) {
    jacobian(0.0);
    return false;
  }

  double x = 0.5*((1.-z+ratio)/ratio)*(1.-sqrt(rho));
  double u = 0.5*((1.-z+ratio)/(1.-z))*(1.-sqrt(rho));

  if ( x < info.emitterX() || x > 1. ||
       u < 0. || u > 1. ) {
    jacobian(0.0);
    return false;
  }

  double phi = 2.*Constants::pi*rphi;

  jacobian(weight*(1./(u+x-2.*u*x)));
  
  lastPt(pt);
  lastZ(z);
  lastPhi(phi);
  lastEmitterZ(x);

  if ( theMCCheck )
    theMCCheck->book(info.emitterX(),1.,info.scale(),info.hardPt(),pt,z,jacobian());

  return true;

}

void IFLightKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  Energy pt = dInfo.lastPt();
  double z = dInfo.lastZ();

  double ratio = sqr(pt)/(2.*pEmitter*pSpectator);
  double rho = 1. - 4.*ratio*z*(1.-z)/sqr(1.-z+ratio);
  
  double x = 0.5*((1.-z+ratio)/ratio)*(1.-sqrt(rho));
  double u = 0.5*((1.-z+ratio)/(1.-z))*(1.-sqrt(rho));

  Lorentz5Momentum kt =
    getKt(pEmitter, pSpectator, pt, dInfo.lastPhi(), true);

  // Initialise the momenta
  Lorentz5Momentum em;
  Lorentz5Momentum emm;
  Lorentz5Momentum spe;

  if ( !theCollinearScheme &&
       x > u && (1.-x)/(x-u) < 1. ) {
    
    assert(false);

    em = ((1.-u)/(x-u))*pEmitter + ((u/x)*(1.-x)/(x-u))*pSpectator - kt/(x-u);
    emm = ((1.-x)/(x-u))*pEmitter + ((u/x)*(1.-u)/(x-u))*pSpectator - kt/(x-u);
    spe = (1.-u/x)*pSpectator;

  } else {

    em = (1./x)*pEmitter;
    emm = ((1.-x)*(1.-u)/x)*pEmitter + u*pSpectator + kt;
    spe = ((1.-x)*u/x)*pEmitter + (1.-u)*pSpectator - kt;
  }

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


void IFLightKinematics::persistentOutput(PersistentOStream &) const {
  //os << theCollinearScheme;
}

void IFLightKinematics::persistentInput(PersistentIStream &, int) {
  //is >> theCollinearScheme;
}

ClassDescription<IFLightKinematics> IFLightKinematics::initIFLightKinematics;
// Definition of the static class description member.

void IFLightKinematics::Init() {

  static ClassDocumentation<IFLightKinematics> documentation
    ("IFLightKinematics implements massless splittings "
     "off a initial-final dipole.");

  /*
  static Switch<IFLightKinematics,bool> interfaceCollinearScheme
    ("CollinearScheme",
     "[experimental] Switch on or off the collinear scheme",
     &IFLightKinematics::theCollinearScheme, false, false, false);
  static SwitchOption interfaceCollinearSchemeYes
    (interfaceCollinearScheme,
     "Yes",
     "Switch on the collinear scheme.",
     true);
  static SwitchOption interfaceCollinearSchemeNo
    (interfaceCollinearScheme,
     "No",
     "Switch off the collinear scheme",
     false);
  
  interfaceCollinearScheme.rank(-1);
  */

}

