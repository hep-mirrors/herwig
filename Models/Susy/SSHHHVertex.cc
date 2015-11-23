// -*- C++ -*-
//
// SSHHHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHHHVertex class.
//

#include "SSHHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHHHVertex::SSHHHVertex() : theMw(ZERO), theZfact(ZERO), theSw(0.),
			     theSbpa(0.), theCbpa(0.), theSbma(0.),
			     theCbma(0.), theS2a(0.), theC2a(0.),
			     theS2b(0.), theC2b(0.), theElast(0.),
			     theq2last(ZERO) {
  orderInGem(1);
  orderInGs(0);
}

void SSHHHVertex::doinit() {
  long sec = 35;
  for(long h = 25; h < 36; h += 10) {
    //self-coupling
    addToList(h, h, h);
    //first-second
    addToList(h,sec,sec);
    //pseudo-scalar
    addToList(h, 36, 36);
    //charged higgs
    addToList(h, 37,-37);
    
    sec = 25;
  }
  SSSVertex::doinit();
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSHHHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  
  theMw = getParticleData(ParticleID::Wplus)->mass();
  theSw = sqrt(sin2ThetaW());
  theZfact = getParticleData(ParticleID::Z0)->mass()/2./
    theSw/sqrt(1. - sqr(theSw));
  
  double tanbeta = theMSSM->tanBeta();
  double sinbeta = tanbeta/sqrt(1. + sqr(tanbeta));
  double cosbeta = sqrt(1. - sqr(sinbeta));
  double sinalpha = sin(theMSSM->higgsMixingAngle());
  double cosalpha = sqrt( 1. - sqr(sinalpha) );
  
  theS2a = 2.*sinalpha*cosalpha;
  theS2b = 2.*sinbeta*cosbeta;
  theC2a = cosalpha*cosalpha - sinalpha*sinalpha;
  theC2b = cosbeta*cosbeta - sinbeta*sinbeta;
  theSbpa = sinbeta*cosalpha + sinalpha*cosbeta;
  theCbpa = cosbeta*cosalpha - sinbeta*sinalpha;
  theSbma = sinbeta*cosalpha - sinalpha*cosbeta;
  theCbma = cosbeta*cosalpha + sinbeta*sinalpha;

}

void SSHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMw,GeV) << ounit(theZfact,GeV) << theSw 
     << theSbpa << theCbpa << theSbma << theCbma << theS2a << theC2a 
     << theS2b << theC2b; 
}

void SSHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMw,GeV) >> iunit(theZfact,GeV) >> theSw 
     >> theSbpa >> theCbpa >> theSbma >> theCbma >> theS2a >> theC2a 
     >> theS2b >> theC2b;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSHHHVertex,SSSVertex>
describeHerwigSSHHHVertex("Herwig::SSHHHVertex", "HwSusy.so");

void SSHHHVertex::Init() {

  static ClassDocumentation<SSHHHVertex> documentation
    ("This is the coupling of a higgs to a pair of higgs bosons "
     "in the MSSM.");

}

void SSHHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
			      tcPDPtr particle2,tcPDPtr particle3) {
  long ids[3] = { abs(particle1->id()), abs(particle2->id()),
		  abs(particle3->id()) };
  long h1(0), h2(0), h3(0), hc(0);
  for(unsigned int i = 0; i < 3; ++i) {
    if( ids[i] == ParticleID::h0) ++h1;
    else if( ids[i] == ParticleID::H0) ++h2;
    else if( ids[i] == ParticleID::A0) ++h3;
    else if( ids[i] == ParticleID::Hplus) ++hc;
    else assert(false);
  }
  assert(h1 + h2 + h3 + hc == 3);
  
  complex<Energy> coupling;
  if( h1 == 3 || h2 == 3 ) {
    coupling = -3.*theZfact*theC2a;
    if( h1 == 3 )
      coupling *= theSbpa;
    else
      coupling *= theCbpa;
  }
  else if( h1 == 1 ) {
    if( h2 == 2 )
      coupling = theZfact*( 2.*theS2a*theCbpa + theSbpa*theC2a );
    else if( h3 == 2 )
      coupling = -theZfact*theC2b*theSbpa;
    else if( hc == 2 )
      coupling = -theMw*theSbma/theSw - theZfact*theC2b*theSbpa;
    else assert(false);
  }
  else if( h2 == 1 ) {
    if( h1 == 2 )
      coupling = -theZfact*( 2.*theS2a*theSbpa - theCbpa*theC2a );
    else if( h3 == 2 )
      coupling = theZfact*theC2b*theCbpa;
    else if( hc == 2 )
      coupling = -theMw*theCbma/theSw + theZfact*theC2b*theCbpa;
    else assert(false);
  }
  
  if( q2 != theq2last || theElast==0. ) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);
  }
  norm(theElast*coupling*UnitRemoval::InvE);
}
			      
