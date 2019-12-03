// -*- C++ -*-
//
// SSWHHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWHHVertex class.
//

#include "SSWHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWHHVertex::SSWHHVertex() : 
  theSw(0.), theS2w(0.), theC2w(0.), thesbma(0.), thecbma(0.), 
  theq2last(ZERO), theElast(0.) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void SSWHHVertex::doinit() {
  addToList(22,37,-37);
  addToList(23,36,25);
  addToList(23,36,35);
  addToList(23,37,-37);
  //outgoing W+
  addToList(24,-37,25);
  addToList(24,-37,35);
  addToList(24,-37,36);
  //outgoing W-
  addToList(-24,37,25);
  addToList(-24,37,35);
  addToList(-24,37,36);
  VSSVertex::doinit();
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSWHHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  theSw = sqrt(sin2ThetaW());
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
}

void SSWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << theSw << theS2w << theC2w << thesbma << thecbma;
}

void SSWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> theSw >> theS2w >> theC2w >> thesbma >> thecbma;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSWHHVertex,VSSVertex>
describeHerwigSSWHHVertex("Herwig::SSWHHVertex", "HwSusy.so");

void SSWHHVertex::Init() {

  static ClassDocumentation<SSWHHVertex> documentation
    ("The coupling of a pair of higgs bosons and a SM gauge boson");

}

void SSWHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
			      tcPDPtr particle2, tcPDPtr particle3) {
  long gboson = particle1->id();
  assert(     gboson  == ParticleID::Z0    ||
	      gboson  == ParticleID::gamma || 
	  abs(gboson) == ParticleID::Wplus );
  long h1ID = particle2->id();
  long h2ID = particle3->id();
  Complex coup(0.);
  if( gboson == ParticleID::Z0 ) {
    if( abs(h1ID) == ParticleID::Hplus ) {
      coup = theC2w/theS2w;
      if(h1ID<0) coup *= -1.;
    }
    else if( h1ID == ParticleID::h0 ||  
	     h2ID == ParticleID::h0 ) {
      coup = Complex(0., 1.)*thecbma/theS2w;
    }
    else {
      coup =-Complex(0., 1.)*thesbma/theS2w;
    }
    if(h2ID==ParticleID::A0) coup *=-1.;
  }
  else if( gboson == ParticleID::gamma ) {
    coup = 1.;
    if(h1ID<0) coup *= -1.;
  }
  else {
    long higgs = abs(h1ID) == ParticleID::Hplus ? h2ID : h1ID;
    if( higgs == ParticleID::h0 ) {
      coup =  0.5*thecbma/theSw;
    }
    else if( higgs == ParticleID::H0) 
      coup = -0.5*thesbma/theSw;
    else 
      coup = -Complex(0., 0.5)/theSw;
    if(abs(h2ID) == ParticleID::Hplus ) coup *= -1.;
    if(gboson<0&&higgs!=ParticleID::A0) coup *= -1.;
  }
  if( q2 != theq2last || theElast==0.) {
    theq2last = q2;
    theElast = electroMagneticCoupling(q2);
  }
  norm(theElast*coup);
}
