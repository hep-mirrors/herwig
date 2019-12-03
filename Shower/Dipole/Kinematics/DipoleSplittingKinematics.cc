// -*- C++ -*-
//
// DipoleSplittingKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleSplittingKinematics class.
//

#include "DipoleSplittingKinematics.h"
#include "Herwig/Shower/Dipole/Base/DipoleSplittingInfo.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

#include <limits>

using namespace Herwig;

DipoleSplittingKinematics::DipoleSplittingKinematics()
  : HandlerBase(), theIRCutoff(1.0*GeV), 
    theXMin(1.e-5), theJacobian(0.0),
    theLastPt(0.0*GeV), theLastZ(0.0), theLastPhi(0.0),
    theLastEmitterZ(1.0), theLastSpectatorZ(1.0),
    theLastSplittingParameters(),theOpenZBoundaries(1) {}

DipoleSplittingKinematics::~DipoleSplittingKinematics() {}

void DipoleSplittingKinematics::persistentOutput(PersistentOStream & os) const {
  os << ounit(theIRCutoff,GeV) << theXMin << theMCCheck<<theOpenZBoundaries;
}

void DipoleSplittingKinematics::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theIRCutoff,GeV) >> theXMin >> theMCCheck>>theOpenZBoundaries;
}

void DipoleSplittingKinematics::prepareSplitting(DipoleSplittingInfo& dInfo) {

  dInfo.splittingKinematics(this);

  if ( lastPt() > IRCutoff() )
    dInfo.lastPt(lastPt());
  else {
    dInfo.lastPt(0.0*GeV);
    dInfo.didStopEvolving();
  }

  dInfo.lastZ(lastZ());
  dInfo.lastPhi(lastPhi());
  dInfo.lastEmitterZ(lastEmitterZ());
  dInfo.lastSpectatorZ(lastSpectatorZ());
  dInfo.splittingParameters().resize(lastSplittingParameters().size());
  copy(lastSplittingParameters().begin(),lastSplittingParameters().end(),
       dInfo.splittingParameters().begin());  
}


Energy DipoleSplittingKinematics::ptMax(Energy dScale, 
                                        double emX, double specX,
                                        const DipoleSplittingInfo& dInfo,
                                        const DipoleSplittingKernel& split) const {
  return ptMax(dScale, emX, specX, dInfo.index(), split);
} 


Energy DipoleSplittingKinematics::ptMax(Energy dScale, 
                                        double emX, double specX,
                                        const DipoleIndex& dIndex,
                                        const DipoleSplittingKernel& split,
                                        tPPtr, tPPtr) const {
  return ptMax(dScale, emX, specX, dIndex, split);
} 

Energy DipoleSplittingKinematics::QMax(Energy dScale, 
                                       double emX, double specX,
                                       const DipoleSplittingInfo& dInfo,
                                       const DipoleSplittingKernel& split) const {
  return QMax(dScale, emX, specX, dInfo.index(), split);
}

Energy DipoleSplittingKinematics::QMax(Energy dScale, 
                                       double emX, double specX,
                                       const DipoleIndex& dIndex,
                                       const DipoleSplittingKernel& split,
                                       tPPtr, tPPtr) const {
  return QMax(dScale, emX, specX, dIndex, split);
}

Energy DipoleSplittingKinematics::generatePt(double r, Energy dScale,
                                             double emX, double specX,
                                             const DipoleIndex& dIndex,
                                             const DipoleSplittingKernel& split,
                                             double& weight) const {

  Energy maxPt = ptMax(dScale,emX,specX,dIndex,split);
  if ( maxPt <= IRCutoff() ) {
    weight = 0.0;
    return ZERO;
  }

  weight *= log(sqr(maxPt/IRCutoff()));

  return IRCutoff()*pow(maxPt/IRCutoff(),r);

}

double DipoleSplittingKinematics::ptToRandom(Energy pt, Energy dScale,
                                             double emX, double specX,
                                             const DipoleIndex& dIndex,
                                             const DipoleSplittingKernel& split) const {

  Energy maxPt = ptMax(dScale,emX,specX,dIndex,split);
  assert(pt >= IRCutoff() && pt <= maxPt);

  return log(pt/IRCutoff())/log(maxPt/IRCutoff());

}

double DipoleSplittingKinematics::generateZ(double r, Energy pt, int sampling,
                                            const DipoleSplittingInfo& dInfo,
                                            const DipoleSplittingKernel& split,
                                            double& weight) const {

  pair<double,double> zLims = zBoundaries(pt,dInfo,split);

  if(zLims.first==zLims.second){  
	weight = 0.0;
    return 0.0;
  }

  using namespace RandomHelpers;

  if ( sampling == FlatZ ) {
    pair<double,double> kw = generate(flat(zLims.first,zLims.second),r);

    if ( kw.second != 0. ) {
      weight *= kw.second;
      return kw.first;
    }
    else {
      assert( kw.first < zLims.first || kw.first > zLims.second );
      weight *= kw.second;
      return -1.;
    }

  }

  if ( sampling == OneOverZ ) {
    pair<double,double> kw = generate(inverse(0.0,zLims.first,zLims.second),r);

    if ( kw.second != 0. ) {
      weight *= kw.second;
      return kw.first;
    }
    else {
      assert( kw.first < zLims.first || kw.first > zLims.second );
      weight *= kw.second;
      return -1.;
    }

  }

  if ( sampling == OneOverOneMinusZ ) {
    pair<double,double> kw = generate(inverse(1.0,zLims.first,zLims.second),r);

    if ( kw.second != 0. ) {
      weight *= kw.second;
      return kw.first;
    }
    else {
      assert( kw.first < zLims.first || kw.first > zLims.second );
      weight *= kw.second;
      return -1.;
    }

  }

  if ( sampling == OneOverZOneMinusZ ) {
    pair<double,double> kw = generate(inverse(0.0,zLims.first,zLims.second) + 
                                      inverse(1.0,zLims.first,zLims.second),r);

    if ( kw.second != 0. ) {
      weight *= kw.second;
      return kw.first;
    }
    else {
      assert( kw.first < zLims.first || kw.first > zLims.second );
      weight *= kw.second;
      return -1.;
    }
  }

  weight = 0.0;
  return 0.0;

}

Lorentz5Momentum DipoleSplittingKinematics::getKt(const Lorentz5Momentum& p1,
                                                  const Lorentz5Momentum& p2,
                                                  Energy pt,
                                                  double phi,
                                                  bool spacelike) const {

  Lorentz5Momentum P;
  // CoM frame
  if ( !spacelike )
    P = p1 + p2;
  // Breit frame
  else
    P = p1 - p2;
  
  Energy mag = sqrt(abs(P.m2()));

  // Define Q. The 'boost' part of this transforms from the current
  // frame into a frame (') in which P' = Q
  Lorentz5Momentum Q = 
    !spacelike ? 
    Lorentz5Momentum(ZERO,ZERO,ZERO,mag,mag) :
    Lorentz5Momentum(ZERO,ZERO,mag,ZERO,-mag);

  // This is required to make the boost parameter
  // gamma positive (construct the 00 term of the
  // transformtion below to see this)
  //if ( spacelike && P.z() < -mag )
  //Q.setZ(-Q.z());
  // Below is safer than above as it avoids
  // cases of very small positive (P*Q + Q2)
  if ( spacelike && P.z() < ZERO )
    Q.setZ(-Q.z());

  Energy2 Q2 = Q.m2();

  // Establish if we need to boost
  bool boost =
    abs((P-Q).vect().mag2()/GeV2) > 1e-10 ||
    abs((P-Q).t()/GeV) > 1e-5;
  
  // Initialise copy of p1 to transform in the following
  Lorentz5Momentum inFrame1(p1);
  if ( boost )
    inFrame1 = inFrame1 - ((P*inFrame1+Q*inFrame1)/(Q2+P*Q))*(P+Q) + 2.*((P*inFrame1)/Q2)*Q; 

  // Compute components of kt
  double cPhi = cos(phi);
  double sPhi = sqrt(1.-sqr(cPhi));
  if ( phi > Constants::pi )
    sPhi = -sPhi;
  
  // Initialise kt
  Lorentz5Momentum kt;
  
  // By 'timelike' case we mean we work in the centre-of-momentum frame
  // The boost to the com frame is defined upto some rotation,
  // here we do the rotation to/from the frame with boosted p1 along the +ve z-axis
  if ( !spacelike ) {

    Axis inFrame1Unit = inFrame1.vect().unit();
    if ( inFrame1Unit.perp2() > 1e-12 ) {
      // 'n' indicates normalised momenta components
      double pxn = inFrame1Unit.x();
      double pyn = inFrame1Unit.y();
      double pzn = inFrame1Unit.z();
      double den = 1./(1.+pzn);
      
      kt.setT(ZERO);
      kt.setX( pt * ( (sqr(pyn)*den + pzn)*cPhi - pxn*pyn*den*sPhi ) );
      kt.setY( pt * ( -pxn*pyn*den*cPhi + (sqr(pxn)*den + pzn)*sPhi) );
      kt.setZ( -pt * ( pxn*cPhi + pyn*sPhi ) );
    }

    // If boosted p1 already lies along the z-axis, construct the pt
    // in this frame, rotating to put boosted p1 along the *+ve* z-axis
    // if required
    else {
      
      // Note pzn will simply be +1 or -1 in this case
      double pzn = inFrame1Unit.z();
      
      // Multiply y component by pzn:
      // In the case of pzn = -1, this corresponds to a rotation
      // about the x-axis as done in boostToSplitting
      kt.setT(ZERO);
      kt.setX( pt * cPhi );
      kt.setY( pt * pzn*sPhi );
      kt.setZ(ZERO);
    }
  }

  // By 'spacelike' we mean we work in the Breit frame.
  // The transformation to the breit frame above is
  // defined up to boosts in the x- and y-directions,
  // here we do the boosts to put the momenta along the z-axis
  else {
    Energy ptx = inFrame1.x();
    Energy pty = inFrame1.y();

    // q/2 = energy component of inFrame1 AFTER applying
    // boosts to eliminate the x and y components. Therefore
    // we calculate q from the mass and the z-component as
    // these will not change due to these boosts.
    Energy q = 2.*sqrt(sqr(inFrame1) + sqr(inFrame1.z()));
  
    Energy Qp = sqrt(4.*(sqr(ptx)+sqr(pty))+sqr(q));
    Energy Qy = sqrt(4.*sqr(pty)+sqr(q));
    
    // Most straightforward way to construct kt in frame
    // where p1 lies along the positive z-axis
    double pzn = inFrame1.z()/abs(inFrame1.z());

    kt.setT(2.*pt*(ptx*q*cPhi+pty*Qp*pzn*sPhi)/(q*Qy));
    kt.setX(pt*(Qp*q*cPhi+4.*ptx*pty*pzn*sPhi)/(q*Qy));
    kt.setY(pt*Qy*pzn*sPhi/q);      
    kt.setZ(ZERO);
  }

  // Transform back to the lab frame
  // Note Q*kt = 0
  if ( boost )
    kt = kt - ((P*kt+Q*kt)/(Q2+P*Q))*(P+Q);// + 2.*((Q*kt)/Q2)*P;
  
  kt.setMass(-pt);
  kt.rescaleRho();

  return kt;  
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

AbstractClassDescription<DipoleSplittingKinematics> DipoleSplittingKinematics::initDipoleSplittingKinematics;
// Definition of the static class description member.

void DipoleSplittingKinematics::Init() {

  static ClassDocumentation<DipoleSplittingKinematics> documentation
    ("DipoleSplittingKinematics is the base class for dipole splittings "
     "as performed in the dipole shower.");


  static Parameter<DipoleSplittingKinematics,Energy> interfaceIRCutoff
    ("IRCutoff",
     "The IR cutoff to be used by this splitting kinematics.",
     &DipoleSplittingKinematics::theIRCutoff, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


  static Parameter<DipoleSplittingKinematics,double> interfaceXMin
    ("XMin",
     "The minimum momentum fraction for incoming partons",
     &DipoleSplittingKinematics::theXMin, 1.0e-5, 0.0, 1.0,
     false, false, Interface::limited);


  static Reference<DipoleSplittingKinematics,DipoleMCCheck> interfaceMCCheck
    ("MCCheck",
     "[debug option] MCCheck",
     &DipoleSplittingKinematics::theMCCheck, false, false, true, true, false);

  interfaceMCCheck.rank(-1);
  
  
  static Switch<DipoleSplittingKinematics,int> interfaceOpenZBoundaries
    ("OpenZBoundaries", "",
     &DipoleSplittingKinematics::theOpenZBoundaries, 0, false, false);
  static SwitchOption interfaceOpenZBoundarieshardScale
    (interfaceOpenZBoundaries,   "Hard",   "",   0);
  static SwitchOption interfaceOpenZBoundariesfull
    (interfaceOpenZBoundaries,   "Full",   "",   1);
  static SwitchOption interfaceOpenZBoundariesDipoleScale
    (interfaceOpenZBoundaries,   "DipoleScale",   "",   2);
  
  
  
  
  

}

