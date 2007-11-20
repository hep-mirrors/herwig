// -*- C++ -*-
//
// SSFFHVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSFFHVertex class.
//

#include "SSFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSFFHVertex::SSFFHVertex() : thetanb(0.0), theMw(0.*MeV), 
			     theSw(0.0), theSa(0.0), theSb(0.0),
			     theCa(0.0), theCb(0.0), theCoupLast(0.0), 
			     theLLast(0.0), theRLast(0.0), theHLast(0), 
			     theFLast(0), theGlast(0.), theq2last() {
  vector<int> first, second, third;
  //neutral higgs
  int higgs = 25;
  for(unsigned int h = 0; h < 3; ++h) {
    if( h == 1 ) higgs = 35;
    if( h == 2 ) higgs = 36;
    //3rd generation quarks
    for(unsigned int i = 5; i < 7; ++i) {
      first.push_back(-i);
      second.push_back(i);
      third.push_back(higgs);
    }
    //tau lepton
    first.push_back(-15);
    second.push_back(15);
    third.push_back(higgs);
  }
  //outgoing H+
  first.push_back(-6);
  second.push_back(5);
  third.push_back(37);
  first.push_back(-16);
  second.push_back(15);
  third.push_back(37);
  //outgoing H-
  first.push_back(-5);
  second.push_back(6);
  third.push_back(-37);
  first.push_back(-15);
  second.push_back(16);
  third.push_back(-37);

  setList(first, second, third);
}

void SSFFHVertex::doinit() throw(InitException) {
  FFSVertex::doinit();
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSFFHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  
  theMw = getParticleData(ParticleID::Wplus)->mass();
  theSw = sqrt(theMSSM->sin2ThetaW());
  thetanb = theMSSM->tanBeta();
  theSb = thetanb/sqrt(1. + sqr(thetanb));
  theCb = sqrt( 1. - sqr(theSb) );
  theSa = sin(theMSSM->higgsMixingAngle());
  theCa = sqrt(1. - sqr(theSa));
  
  orderInGem(1);
  orderInGs(0);
}

void SSFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << theMSSM  << thetanb << ounit(theMw,GeV) << theSw << theSa
     << theSb << theCa << theCb;
}

void SSFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM  >> thetanb >> iunit(theMw,GeV) >> theSw >> theSa
     >> theSb >> theCa >> theCb;
}

ClassDescription<SSFFHVertex> SSFFHVertex::initSSFFHVertex;
// Definition of the static class description member.

void SSFFHVertex::Init() {

  static ClassDocumentation<SSFFHVertex> documentation
    ("The coupling of the higgs bosons to SM fermions in the MSSM");

}

void SSFFHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			      tcPDPtr particle3, int) {
  long id1(abs(particle1->id())), id2(abs(particle2->id())), 
    id3(abs(particle3->id())), higgsID(0), f1ID(0), f2ID(0);
  if( id1 == ParticleID::h0 || id1 == ParticleID::H0 || 
      id1 == ParticleID::A0 || id1 == ParticleID::Hplus ) {
    higgsID = id1;
    f1ID = id2;
    f2ID = id3;
  }
  else if( id2 == ParticleID::h0 || id2 == ParticleID::H0 || 
	   id2 == ParticleID::A0 || id2 == ParticleID::Hplus  ) {
    higgsID = id2;
    f1ID = id1;
    f2ID = id3;
  }
  else if( id3 == ParticleID::h0 || id3 == ParticleID::H0 || 
	   id3 == ParticleID::A0 || id3 == ParticleID::Hplus ) {
    higgsID = id3;
    f1ID = id1;
    f2ID = id2;
  }
  else {
    throw HelicityConsistencyError() 
      << "SSFFHVertex::setCoupling - There is no higgs particle in "
      << "this vertex. Particles: " << id1 << " " << id2 << " " << id3
      << Exception::warning;
    return;
  }
  
  if( ((f1ID > 6 && f1ID < 11) || f1ID > 16 ) ||
      ((f2ID > 6 && f1ID < 11) || f2ID > 16 ) ) {
    throw HelicityConsistencyError() 
      << "SSFFHVertex::setCoupling - The particles in this vertex " 
      << " are not SM fermions " << f1ID << " " << f2ID 
      << Exception::warning;
    setNorm(0.);
    setLeft(0.);
    setRight(0.);
    return;
  }
  if( q2 != theq2last ) {
    theGlast = sqrt(4.*Constants::pi*theMSSM->alphaEM(q2))/theSw;
    theq2last = q2;
  }

  if(higgsID == theHLast && f1ID == theFLast) {
    setNorm(theGlast*theCoupLast);
    setLeft(theLLast);
    setRight(theRLast);
    return;
  }
  theHLast = higgsID;
  if( higgsID == ParticleID::h0 || higgsID == ParticleID::H0 ||
      higgsID == ParticleID::A0 ) {
    assert( f1ID == f2ID);
    theFLast = f1ID;
    theCoupLast = getParticleData(f1ID)->mass()/2./theMw;
    
    if( higgsID == ParticleID::h0 ) {
      if( f1ID % 2 == 0 )
	theCoupLast *= -theCa/theSb;
      else
	theCoupLast *= theSa/theCb;
      theLLast = 1.;
      theRLast = 1.;
    }
    else if( higgsID == ParticleID::H0 ) {
      if( f1ID % 2 == 0 )
	theCoupLast *= -theSa/theSb;
      else
	theCoupLast *= -theCa/theCb;
      theLLast = 1.;
      theRLast = 1.;
    }
    else {
      if( f1ID % 2 == 0 )
	theCoupLast /= thetanb; 
      else
	theCoupLast *= thetanb;
      theCoupLast *= Complex(0.,1.);
      theLLast = 1.;
      theRLast = -1.;
    }
  }
  //H+
  else {
    if( id1 % 2 != 0 ) swap(f1ID, f2ID);
    theFLast = f1ID;
    theCoupLast = 1./sqrt(2);
    theLLast = getParticleData(f1ID)->mass()/thetanb/theMw;
    theRLast = getParticleData(f2ID)->mass()*thetanb/theMw;
  }
  setNorm(theGlast*theCoupLast);
  setLeft(theLLast);
  setRight(theRLast);

}

