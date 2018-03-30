// -*- C++ -*-
//
// DipoleSplittingKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
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


Energy DipoleSplittingKinematics::QMax(Energy dScale, 
				       double emX, double specX,
				       const DipoleSplittingInfo& dInfo,
				       const DipoleSplittingKernel& split) const {
  return QMax(dScale, emX, specX, dInfo.index(), split);
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
  if ( !spacelike )
    P = p1 + p2;
  else
    P = p1 - p2;

  Energy2 Q2 = abs(P.m2());

  Lorentz5Momentum Q = 
    !spacelike ? 
    Lorentz5Momentum(ZERO,ZERO,ZERO,sqrt(Q2),sqrt(Q2)) :
    Lorentz5Momentum(ZERO,ZERO,sqrt(Q2),ZERO,-sqrt(Q2));

  if ( spacelike && Q.z() < P.z() )
    Q.setZ(-Q.z());

  bool boost =
    abs((P-Q).vect().mag2()/GeV2) > 1e-10 ||
    abs((P-Q).t()/GeV) > 1e-5;
  boost &= (P*Q-Q.mass2())/GeV2 > 1e-8;

  Lorentz5Momentum inFrame1;
  if ( boost )
    inFrame1 = p1 + ((P*p1-Q*p1)/(P*Q-Q.mass2()))*(P-Q);
  else
    inFrame1 = p1;

  Energy ptx = inFrame1.x();
  Energy pty = inFrame1.y();
  Energy q = 2.*inFrame1.z();

  Energy Qp = sqrt(4.*(sqr(ptx)+sqr(pty))+sqr(q));
  Energy Qy = sqrt(4.*sqr(pty)+sqr(q));

  double cPhi = cos(phi);
  double sPhi = sqrt(1.-sqr(cPhi));
  if ( phi > Constants::pi )
    sPhi = -sPhi;

  Lorentz5Momentum kt;

  if ( !spacelike ) {
    kt.setT(ZERO);
    kt.setX(pt*Qy*cPhi/Qp);
    kt.setY(-pt*(4*ptx*pty*cPhi/Qp+q*sPhi)/Qy);
    kt.setZ(2.*pt*(-ptx*q*cPhi/Qp + pty*sPhi)/Qy);
  } else {
    kt.setT(2.*pt*(ptx*q*cPhi+pty*Qp*sPhi)/(q*Qy));
    kt.setX(pt*(Qp*q*cPhi+4.*ptx*pty*sPhi)/(q*Qy));
    kt.setY(pt*Qy*sPhi/q);
    kt.setZ(ZERO);
  }

  if ( boost )
    kt = kt + ((P*kt-Q*kt)/(P*Q-Q.mass2()))*(P-Q);
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

