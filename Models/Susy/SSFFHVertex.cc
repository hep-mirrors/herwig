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

SSFFHVertex::SSFFHVertex() : thetanb(0.0), theMw(ZERO), 
			     theSa(0.0), theSb(0.0),
			     theCa(0.0), theCb(0.0), theCoupLast(0.0), 
			     theLLast(0.0), theRLast(0.0), theHLast(0), 
			     theFLast(0), theGlast(0.), theq2last() {
  int higgs[] = { 25, 35, 36 };
  for ( long h = 0; h < 3; ++h ) {
    //neutral higgs
    //3rd generation quarks
    addToList(-5,5,higgs[h]);
    addToList(-6,6,higgs[h]);
    //tau lepton
    addToList(-15,15,higgs[h]);
  }
  //outgoing H+
  addToList(-6 , 5, 37);
  addToList(-16,15, 37);
  //outgoing H-
  addToList(-5 ,6 ,-37);
  addToList(-15,16,-37);
}

void SSFFHVertex::doinit() {
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSFFHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  
  theMw = getParticleData(ParticleID::Wplus)->mass();
  thetanb = theMSSM->tanBeta();
  theSb = thetanb/sqrt(1. + sqr(thetanb));
  theCb = sqrt( 1. - sqr(theSb) );
  theSa = sin(theMSSM->higgsMixingAngle());
  theCa = sqrt(1. - sqr(theSa));
  
  orderInGem(1);
  orderInGs(0);
  FFSVertex::doinit();
}

void SSFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << thetanb << ounit(theMw,GeV) << theSa
     << theSb << theCa << theCb;
}

void SSFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> thetanb >> iunit(theMw,GeV) >> theSa
     >> theSb >> theCa >> theCb;
}

ClassDescription<SSFFHVertex> SSFFHVertex::initSSFFHVertex;
// Definition of the static class description member.

void SSFFHVertex::Init() {

  static ClassDocumentation<SSFFHVertex> documentation
    ("The coupling of the higgs bosons to SM fermions in the MSSM");

}

void SSFFHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
			      tcPDPtr particle2,tcPDPtr particle3) {
  long f1ID(abs(particle1->id())), f2ID(abs(particle2->id())), 
    higgsID(particle3->id());
  // check higgs
  assert( higgsID == ParticleID::h0 ||     higgsID  == ParticleID::H0 || 
	  higgsID == ParticleID::A0 || abs(higgsID) == ParticleID::Hplus );
  // check fermions
  assert(!( ((f1ID > 6 && f1ID < 11) || f1ID > 16 ) ||
	    ((f2ID > 6 && f1ID < 11) || f2ID > 16 ) ));
  if( q2 != theq2last || theGlast==0.) {
    theGlast = weakCoupling(q2);
    theq2last = q2;
  }
  if(higgsID == theHLast && f1ID == theFLast) {
    norm(theGlast*theCoupLast);
    left(theLLast);
    right(theRLast);
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
    if( f1ID % 2 != 0 ) swap(f1ID, f2ID);
    theFLast = f1ID;
    theCoupLast = 1./sqrt(2);
    theLLast = getParticleData(f1ID)->mass()/thetanb/theMw;
    theRLast = getParticleData(f2ID)->mass()*thetanb/theMw;
    if( higgsID < 0 ) {
      Complex tmp = theLLast;
      theLLast = conj(theRLast);
      theRLast = conj(tmp);
    }
  }
  norm(theGlast*theCoupLast);
  left(theLLast);
  right(theRLast);

}

