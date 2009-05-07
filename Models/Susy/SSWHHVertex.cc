// -*- C++ -*-
//
// SSWHHVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWHHVertex class.
//

#include "SSWHHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWHHVertex::SSWHHVertex() : 
  theSw(0.), theS2w(0.), theC2w(0.), thesbma(0.), thecbma(0.), 
  theGBlast(0), theHlast(0), theCouplast(0.), theq2last(ZERO), 
  theElast(0.) {
  vector<long> first, second, third;
  //photon
  first.push_back(22);
  second.push_back(37);
  third.push_back(-37);
  //Z
  first.push_back(23);
  second.push_back(36);
  third.push_back(25);
  first.push_back(23);
  second.push_back(36);
  third.push_back(35);
  first.push_back(23);
  second.push_back(37);
  third.push_back(-37);
  //outgoing W+
  first.push_back(24);
  second.push_back(-37);
  third.push_back(25);
  first.push_back(24);
  second.push_back(-37);
  third.push_back(35);
  first.push_back(24);
  second.push_back(-37);
  third.push_back(36);
  //outgoing W-
  first.push_back(-24);
  second.push_back(37);
  third.push_back(25);
  first.push_back(-24);
  second.push_back(37);
  third.push_back(35);
  first.push_back(-24);
  second.push_back(37);
  third.push_back(36);
  
  setList(first, second, third);
}

void SSWHHVertex::doinit() {
  VSSVertex::doinit();
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSWHHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  theSw = sqrt(theMSSM->sin2ThetaW());
  double cw = sqrt(1. - sqr(theSw));
  theS2w = 2.*theSw*cw;
  theC2w = cw*cw - theSw*theSw;

  double sina = sin(theMSSM->higgsMixingAngle());
  double cosa =  sqrt(1. - sqr(sina));
  double tanb = theMSSM->tanBeta();
  double sinb = tanb/sqrt(1. + sqr(tanb));
  double cosb = sqrt( 1. - sqr(sinb) );
  thesbma = sinb*cosa - sina*cosb;
  thecbma = cosa*cosb + sina*sinb;

  orderInGs(0);
  orderInGem(1);
}

void SSWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << theSw << theS2w << theC2w << thesbma << thecbma;
}

void SSWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> theSw >> theS2w >> theC2w >> thesbma >> thecbma;
}

ClassDescription<SSWHHVertex> SSWHHVertex::initSSWHHVertex;
// Definition of the static class description member.

void SSWHHVertex::Init() {

  static ClassDocumentation<SSWHHVertex> documentation
    ("The coupling of a pair of higgs bosons and a SM gauge boson");

}

void SSWHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			      tcPDPtr particle3) {
  long id1(abs(particle1->id())), id2(abs(particle2->id())),
    id3(abs(particle3->id())), gboson(0), h1ID(0), h2ID(0);
  if( id1 == ParticleID::Z0 || id1 == ParticleID::gamma || 
      id1 == ParticleID::Wplus ) {
    gboson = id1;
    h1ID = id2;
    h2ID = id3;
  }
  else if( id2 == ParticleID::Z0 || id2 == ParticleID::gamma || 
	   id2 == ParticleID::Wplus ) {
    gboson = id2;
    h1ID = id1;
    h2ID = id3;
  }
  else if( id3 == ParticleID::Z0 || id3 == ParticleID::gamma || 
	   id3 == ParticleID::Wplus ) {
    gboson = id3;
    h1ID = id1;
    h2ID = id2;
  }
  else {
    throw HelicityConsistencyError() 
      << "SSWHHVertex::setCoupling - There is no gauge boson particle in "
      << "this vertex. Particles: " << id1 << " " << id2 << " " << id3
      << Exception::warning;
    return;
  }
  long higgs(0);
  if( gboson == ParticleID::Wplus )
    higgs = (h1ID == ParticleID::Hplus) ? h2ID : h1ID;
  else
    higgs = (h1ID >= ParticleID::A0) ? h2ID : h1ID;

  if( higgs != ParticleID::h0 && higgs != ParticleID::H0 &&
      higgs != ParticleID::A0 && higgs != ParticleID::Hplus) {
    throw HelicityConsistencyError() 
      << "SSWHHVertex::setCoupling - There is no higgs in this "
      << "this vertex. " << higgs << Exception::warning;
    setNorm(0.0);
    return;
  }

  if( gboson != theGBlast || higgs != theHlast ) {
    theGBlast = gboson;
    theHlast = higgs;
    //photon
    if( gboson == ParticleID::gamma )
      theCouplast = 1.;
    else if( gboson == ParticleID::Z0 ) {
      if( higgs == ParticleID::Hplus ) 
	theCouplast = theC2w/theS2w;
      else if( higgs == ParticleID::h0 ) 
	theCouplast = -Complex(0., 1.)*thecbma/theS2w;
      else 
	theCouplast = Complex(0., 1.)*thesbma/theS2w;
    }
    else {
      if( higgs == ParticleID::h0 ) {
	theCouplast = -0.5*thecbma/theSw;
      }
      else if( higgs == ParticleID::H0) 
	theCouplast = 0.5*thesbma/theSw;
      else 
	theCouplast = Complex(0., 0.5)/theSw;
    }
  }
  if( q2 != theq2last || theElast==0.) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);
  }

  setNorm(theElast*theCouplast);
}
