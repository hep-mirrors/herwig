// -*- C++ -*-
//
// SSWWHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWWHVertex class.
//

#include "SSWWHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWWHVertex::SSWWHVertex() : theh0Wfact(ZERO), theH0Wfact(ZERO), 
			     theh0Zfact(ZERO), theH0Zfact(ZERO),
			     theCoupLast(ZERO), theElast(0.0),
			     theq2last(ZERO), theHlast(0), 
			     theGBlast(0) {
  orderInGem(1);
  orderInGs(0);
}

void SSWWHVertex::doinit() {
  addToList(23,23,25);
  addToList(-24,24,25);
  addToList(23,23,35);
  addToList(-24,24,35);
  VVSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw InitException() 
      << "SSWWHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;

  Energy mw = getParticleData(ParticleID::Wplus)->mass();
  Energy mz = getParticleData(ParticleID::Z0)->mass();
  double sw = sqrt(sin2ThetaW());
  double sinalp = sin(model->higgsMixingAngle());
  double cosalp = sqrt(1. - sqr(sinalp));
  double tanbeta = model->tanBeta();
  double sinbeta = tanbeta/sqrt(1. + sqr(tanbeta));
  double cosbeta = sqrt( 1. - sqr(sinbeta) );
  double sinbma = sinbeta*cosalp - cosbeta*sinalp;
  double cosbma = cosbeta*cosalp + sinbeta*sinalp;
  
  theh0Wfact = mw*sinbma/sw;
  theH0Wfact = mw*cosbma/sw;
  theh0Zfact = mz*sinbma/sw/sqrt(1. - sw*sw);
  theH0Zfact = mz*cosbma/sw/sqrt(1. - sw*sw);
}

void SSWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theh0Wfact,GeV) << ounit(theH0Wfact,GeV) 
     << ounit(theh0Zfact,GeV) << ounit(theH0Zfact,GeV);
}

void SSWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theh0Wfact,GeV) >> iunit(theH0Wfact,GeV) 
     >> iunit(theh0Zfact,GeV) >> iunit(theH0Zfact,GeV);
}

ClassDescription<SSWWHVertex> SSWWHVertex::initSSWWHVertex;
// Definition of the static class description member.

void SSWWHVertex::Init() {

  static ClassDocumentation<SSWWHVertex> documentation
    ("This is the coupling of a pair of SM gauge bosons"
     "to the higgs particles in the MSSM");

}

void SSWWHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr,
			      tcPDPtr particle3) {
  long bosonID = abs(particle1->id());
  long higgsID =     particle3->id();
  assert( higgsID == ParticleID::h0    || higgsID == ParticleID::H0 );
  assert( bosonID == ParticleID::Wplus || bosonID == ParticleID::Z0 );
  if( higgsID != theHlast || bosonID != theGBlast) {
    if( higgsID == ParticleID::h0 )
      theCoupLast = (bosonID == ParticleID::Z0) ? theh0Zfact : theh0Wfact;
    else
      theCoupLast = (bosonID == ParticleID::Z0) ? theH0Zfact : theH0Wfact;
  }
  if( q2 != theq2last ) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);
  }
  norm(theElast*theCoupLast*UnitRemoval::InvE);
}
