// -*- C++ -*-
//
// SSHHHVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHHHVertex class.
//

#include "SSHHHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHHHVertex::SSHHHVertex() : theMw(0.*MeV), theZfact(0.*MeV), theSw(0.),
			     theSbpa(0.), theCbpa(0.), theSbma(0.),
			     theCbma(0.), theS2a(0.), theC2a(0.),
			     theS2b(0.), theC2b(0.), theElast(0.),
			     theq2last(0.*MeV2) 
{
  vector<long> first, second, third;
  int sec = 35;
  for(long h = 25; h < 36; h += 10) {
    //self-coupling
    first.push_back(h);
    second.push_back(h);
    third.push_back(h);
    //first-second
    first.push_back(h);
    second.push_back(sec);
    third.push_back(sec);
    //pseudo-scalar
    first.push_back(h);
    second.push_back(36);
    third.push_back(36);
    //charged higgs
    first.push_back(h);
    second.push_back(37);
    third.push_back(37);
    
    sec = 25;
  }
  
  setList(first, second, third);

}

void SSHHHVertex::doinit() throw(InitException) {
  SSSVertex::doinit();
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSHHHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  
  theMw = getParticleData(ParticleID::Wplus)->mass();
  theSw = sqrt(theMSSM->sin2ThetaW());
  theZfact = getParticleData(ParticleID::Z0)->mass()/2./theSw/sqrt(1. - sqr(theSw));
  
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

  orderInGem(1);
  orderInGs(0);
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

ClassDescription<SSHHHVertex> SSHHHVertex::initSSHHHVertex;
// Definition of the static class description member.

void SSHHHVertex::Init() {

  static ClassDocumentation<SSHHHVertex> documentation
    ("This is the coupling of a higgs to a pair of higgs bosons "
     "in the MSSM.");

}

void SSHHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			      tcPDPtr particle3) {
  long ids[3] = { abs(particle1->id()), abs(particle2->id()),
		  abs(particle3->id()) };
  long h1(0), h2(0), h3(0), hc(0);
  for(unsigned int i = 0; i < 3; ++i) {
    if( ids[i] == ParticleID::h0) ++h1;
    else if( ids[i] == ParticleID::H0) ++h2;
    else if( ids[i] == ParticleID::A0) ++h3;
    else if( ids[i] == ParticleID::Hplus) ++hc;
    else {
      throw HelicityConsistencyError()
	<< "SSHHHVertex::setCoupling - There is an unrecognised particle in "
	<< "this vertex! " << ids[i] << Exception::runerror; 
	break;
    }
  }
  assert(h1 + h2 + h3 + hc == 3);
  
  complex<Energy> coupling;
  bool unrec(false);
  if( h1 == 3 || h2 == 3 ) {
    coupling = -3.*theZfact*theC2a;
    if( h1 == 3 )
      coupling *= theCbpa;
    else
      coupling *= theSbpa;
  }
  else if( h1 == 1 ) {
    if( h2 == 2 )
      coupling = theZfact*( 2.*theS2a*theCbpa + theSbpa*theC2a );
    else if( h3 == 2 )
      coupling = -theZfact*theC2b*theSbpa;
    else if( hc == 2 )
      coupling = -theMw*theSbma/theSw - theZfact*theC2b*theSbpa;
    else unrec = true;
  }
  else if( h2 == 1 ){
    if( h1 == 2 )
      coupling = -theZfact*( 2.*theS2a*theSbpa - theCbpa*theC2a );
    else if( h3 == 2 )
      coupling = theZfact*theC2b*theCbpa;
    else if( hc == 2 )
      coupling = -theMw*theCbma/theSw + theZfact*theC2b*theCbpa;
    else unrec = true;
  }
  
  if( unrec ) {
    throw HelicityConsistencyError() 
      << "SSHHHVertex::setCoupling - Trying to calculate the coupling "
      << "for an unrecognised combination of particles "
      << ids[0] << " " << ids[1] << " " << ids[2] << Exception::runerror;
    setNorm(0.);
    return;
  }
  
  if( q2 != theq2last ) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);;
  }

  setNorm(theElast*coupling*UnitRemoval::InvE);
  
}
			      
