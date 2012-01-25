// -*- C++ -*-
//
// DipoleSplittingKinematics.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleSplittingKinematics class.
//

#include "DipoleSplittingKinematics.h"
#include "Herwig++/DipoleShower/Base/DipoleSplittingInfo.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <limits>

using namespace Herwig;

DipoleSplittingKinematics::DipoleSplittingKinematics()
  : HandlerBase(), theIRCutoff(1.0*GeV), 
    theXMin(1.e-5), theJacobian(0.0),
    theLastPt(0.0*GeV), theLastZ(0.0), theLastPhi(0.0),
    theLastEmitterZ(1.0), theLastSpectatorZ(1.0),
    theLastSplittingParameters() {}

DipoleSplittingKinematics::~DipoleSplittingKinematics() {}

void DipoleSplittingKinematics::persistentOutput(PersistentOStream & os) const {
  os << ounit(theIRCutoff,GeV) << theXMin << theMCCheck;
}

void DipoleSplittingKinematics::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theIRCutoff,GeV) >> theXMin >> theMCCheck;
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

Lorentz5Momentum DipoleSplittingKinematics::getKt (const Lorentz5Momentum& p1,
						   const Lorentz5Momentum& p2,
						   Energy pt,
						   double phi) const {

  Boost beta = (p1+p2).findBoostToCM();
      
  Lorentz5Momentum p1c = p1;

  if (beta.mag2() > Constants::epsilon) {
    p1c.boost(beta);
  }
      
  Lorentz5Momentum k (0.*GeV,0.*GeV,0.*GeV,0.*GeV);
      
  double ct = p1c.vect().unit().z();
  double st = sqrt(1.-ct*ct);
      
  double cphi = cos(phi);
  double sphi = sqrt(1.-cphi*cphi);
  if (phi  > Constants::pi) sphi = -sphi;
      
  if (st > Constants::epsilon) {
    double cchi = p1c.vect().unit().x()/st;
    double schi = p1c.vect().unit().y()/st;
    k.setX((cphi*cchi*ct-sphi*schi)*pt);
    k.setY((cphi*schi*ct+sphi*cchi)*pt);
    k.setZ(-cphi*st*pt);
  } else {
    k.setX(pt*cphi);
    k.setY(pt*sphi);
    k.setZ(0.*GeV);
  }
      
  if (beta.mag2() > Constants::epsilon)
    k.boost(-beta);
      
  return k;

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

}

