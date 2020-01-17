// -*- C++ -*-
//
// IILightKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IILightKinematics class.
//

#include "IILightKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"

using namespace Herwig;

IILightKinematics::IILightKinematics() 
  : DipoleSplittingKinematics(), theCollinearScheme(true), didCollinear(false),
    theTransformationCalculated(false) {}
//theTransformHardOnly(false) {}

IILightKinematics::~IILightKinematics() {}

IBPtr IILightKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IILightKinematics::fullclone() const {
  return new_ptr(*this);
}

Energy IILightKinematics::ptMax(Energy dScale, 
				double emX, double specX,
				const DipoleIndex&,
				const DipoleSplittingKernel&) const {
  double tau = 
    !theCollinearScheme ? emX*specX : emX;
  return (1.-tau) * dScale / (2.*sqrt(tau));
}

Energy IILightKinematics::QMax(Energy, 
			       double, double,
			       const DipoleIndex&,
			       const DipoleSplittingKernel&) const {
  assert(false && "add this");
  return 0.0*GeV;
}

Energy IILightKinematics::PtFromQ(Energy scale, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return scale*sqrt(1.-z);
}

Energy IILightKinematics::QFromPt(Energy scale, const DipoleSplittingInfo& split) const {
  double z = split.lastZ();
  return scale/sqrt(1.-z);
}

pair<double,double> IILightKinematics::zBoundaries(Energy pt,
						   const DipoleSplittingInfo& dInfo,
						   const DipoleSplittingKernel&) const {
  double x = 
    !theCollinearScheme ?
    dInfo.emitterX()*dInfo.spectatorX() :
    dInfo.emitterX();

  Energy hard=dInfo.hardPt();
  if(openZBoundaries()==1)hard=(1.-x) *dInfo.scale()/(2.*sqrt(x));
  if(openZBoundaries()==2)hard=min(dInfo.scale(),(1.-x) *dInfo.scale()/(2.*sqrt(x)));
  if(hard<pt)return {0.5*(1.+x),0.5*(1.+x)};

  double s = sqrt(1.-sqr(pt/hard));

  return {0.5*(1.+x-(1.-x)*s),0.5*(1.+x+(1.-x)*s)};
}

bool IILightKinematics::generateSplitting(double kappa, double xi, double rphi,
					  DipoleSplittingInfo& info,
					  const DipoleSplittingKernel& split) {

  if ( info.emitterX() < xMin() ||
       info.spectatorX() < xMin() ) {
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
  double x = z*(1.-z)/(1.-z+ratio);
  double v = ratio*z /(1.-z+ratio);

  if ( x < 0. || x > 1. || v < 0. || v > 1.-x ) {
    jacobian(0.0);
    return false;
  }
  
  if ( !theCollinearScheme &&
       (1.-v-x)/(v+x) < 1. ) {
    if ( (x+v) < info.emitterX() ||
         x/(x+v) < info.spectatorX() ) {
      jacobian(0.0);
      return false;
    }
  } else {
    if ( x < info.emitterX() ) {
      jacobian(0.0);
      return false;
    }
  }
  
  double phi = 2.*Constants::pi*rphi;

  jacobian(weight*(1./z));

  lastPt(pt);
  lastZ(z);
  lastPhi(phi);

  if ( !theCollinearScheme &&
       (1.-v-x)/(v+x) < 1. ) {
    lastEmitterZ(x+v);
    lastSpectatorZ(x/(x+v));
  } else {
    lastEmitterZ(x);
    lastSpectatorZ(1.);
  }

  if ( theMCCheck )
    theMCCheck->book(info.emitterX(),info.spectatorX(),info.scale(),info.hardPt(),pt,z,jacobian());

  return true;

}

void IILightKinematics::generateKinematics(const Lorentz5Momentum& pEmitter,
					   const Lorentz5Momentum& pSpectator,
					   const DipoleSplittingInfo& dInfo) {

  Energy pt = dInfo.lastPt(); 
  double z = dInfo.lastZ();

  double ratio = sqr(pt)/(2.*pEmitter*pSpectator);

  double x = z*(1.-z)/(1.-z+ratio);
  double v = ratio*z /(1.-z+ratio);
    
  Lorentz5Momentum kt =
    getKt(pEmitter, pSpectator, pt, dInfo.lastPhi());

  // Initialise the momenta
  Lorentz5Momentum em;
  Lorentz5Momentum emm;
  Lorentz5Momentum spe;
  
  if ( !theCollinearScheme &&
       (1.-v-x)/(v+x) < 1. ) {

    assert(false);

    em = (1./(v+x))*pEmitter+(v*(1.-v-x)/(x*(x+v)))*pSpectator+kt/(x+v);
    emm = ((1.-v-x)/(v+x))*pEmitter+(v/(x*(x+v)))*pSpectator+kt/(x+v);
    spe = (1.+v/x)*pSpectator;
    
    didCollinear = false;

  } else {

    em = (1./x)*pEmitter;
    emm = ((1.-x-v)/x)*pEmitter+v*pSpectator+kt;
    spe = pSpectator;

    K = em + spe - emm;
    K2 = K.m2();
    
    Ktilde = pEmitter + pSpectator;
    KplusKtilde = K + Ktilde;
    KplusKtilde2 = KplusKtilde.m2();

    didCollinear = true;

    // Set indicator that the transformation,
    // required for spin correlations, hasn't
    // been calculated.
    theTransformationCalculated = false;
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

void IILightKinematics::setTransformation () {
  
    // Construct transformation for spin correlations
    // Clear the rotation
    theRecoilTransformation = LorentzRotation();
  
    // Construct boost part
    Energy KplusKtildeT = KplusKtilde.t();
    Energy KT = K.t();
  
    double tt = 1. - 2.*sqr(KplusKtildeT)/KplusKtilde2            + 2.*KT*Ktilde.t()/K2;
    double tx =      2.*KplusKtildeT*KplusKtilde.x()/KplusKtilde2 - 2.*KT*Ktilde.x()/K2;
    double ty =      2.*KplusKtildeT*KplusKtilde.y()/KplusKtilde2 - 2.*KT*Ktilde.y()/K2;  
    double tz =      2.*KplusKtildeT*KplusKtilde.z()/KplusKtilde2 - 2.*KT*Ktilde.z()/K2;
    
    theRecoilTransformation.boost(tx/tt, ty/tt, tz/tt, tt);
    
    // Rotate KtildeSpinCorrTrans to z-axis
    Lorentz5Momentum KtildeSpinCorrTrans = theRecoilTransformation*Lorentz5Momentum(Ktilde);
    Axis axis(KtildeSpinCorrTrans.vect().unit());
    if( axis.perp2() > 1e-12 ) {
      double sinth(sqrt(1.-sqr(axis.z())));
      theRecoilTransformation.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    }
    else if( axis.z() < 0. ) {
      theRecoilTransformation.rotate(Constants::pi,Axis(1.,0.,0.));
    }
    
    // Rotate from z-axis to K
    Axis axis2(K.vect().unit());
    if( axis2.perp2() > 1e-12 ) {
      double sinth(sqrt(1.-sqr(axis2.z())));
      theRecoilTransformation.rotate(acos(axis2.z()),Axis(-axis2.y()/sinth,axis2.x()/sinth,0.));
    }
    else if( axis2.z() < 0. ) {
      theRecoilTransformation.rotate(Constants::pi,Axis(1.,0.,0.));
    }
  
    // SW 30/01/2019: Test feature only, not for release.
    // Required for absorbing recoil in hard particles only
    //splitRecoilMomentum(K);

    // Set indicator that the transformation,
    // has been calculated
    theTransformationCalculated = true;
  }


// If needed, insert default implementations of function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IILightKinematics::persistentOutput(PersistentOStream &) const {
  //os << theTransformHardOnly;
  //os << theCollinearScheme;
}

void IILightKinematics::persistentInput(PersistentIStream &, int) {
  //is >> theTransformHardOnly;
  //is >> theCollinearScheme;
}

ClassDescription<IILightKinematics> IILightKinematics::initIILightKinematics;
// Definition of the static class description member.

void IILightKinematics::Init() {

  static ClassDocumentation<IILightKinematics> documentation
    ("IILightKinematics implements massless splittings "
     "off an initial-initial dipole.");

  /*
  static Switch<IILightKinematics,bool> interfaceTransformHardOnly
    ("TransformHardOnly",
     "[Dev only] Apply recoil transformation to colourless particles only.",
     &IILightKinematics::theTransformHardOnly, false, false, false);
  static SwitchOption interfaceTransformHardOnlyYes
    (interfaceTransformHardOnly,
     "Yes",
     "Transform only colourless particles.",
     true);
  static SwitchOption interfaceTransformHardOnlyNo
    (interfaceTransformHardOnly,
     "No",
     "Transform all particles.",
     false);

  interfaceTransformHardOnly.rank(-1);
  */
  /*
  static Switch<IILightKinematics,bool> interfaceCollinearScheme
    ("CollinearScheme",
     "[experimental] Switch on or off the collinear scheme",
     &IILightKinematics::theCollinearScheme, false, false, false);
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

