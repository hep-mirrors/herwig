// -*- C++ -*-
//
// SSFFHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSFFHVertex class.
//

#include "SSFFHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSFFHVertex::SSFFHVertex() : thetanb(0.0), theMw(ZERO), 
			     theSa(0.0), theSb(0.0),
			     theCa(0.0), theCb(0.0),
			     theFLast(make_pair(0,0)), theGlast(0.),
			     theq2last(), theMassLast(make_pair(ZERO,ZERO)) {
  orderInGem(1);
  orderInGs(0);
}

void SSFFHVertex::doinit() {
  int higgs[] = { 25, 35, 36 };
  for ( long h = 0; h < 3; ++h ) {
    //neutral higgs
    // quarks
    for(long ix=1;ix<7;++ix)
      addToList(-ix,ix,higgs[h]);
    // charged leptons
    for(long ix=11;ix<16;ix+=2)
      addToList(-ix,ix,higgs[h]);
  }
  for(long ix=1;ix<6;ix+=2) {
    //outgoing H+
    addToList(-ix-1,  ix, 37);
    //outgoing H-
    addToList(-ix  ,ix+1,-37);
  }
  for(long ix=11;ix<16;ix+=2) {
    //outgoing H+
    addToList(-ix-1,  ix, 37);
    //outgoing H-
    addToList(-ix  ,ix+1,-37);
  }
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
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
  
  FFSVertex::doinit();
}

void SSFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << theMSSM << thetanb << ounit(theMw,GeV) << theSa
     << theSb << theCa << theCb;
}

void SSFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM >> thetanb >> iunit(theMw,GeV) >> theSa
     >> theSb >> theCa >> theCb;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSFFHVertex,FFSVertex>
describeHerwigSSFFHVertex("Herwig::SSFFHVertex", "HwSusy.so");

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
  }
  if( q2 != theq2last || theFLast.first  != f1ID) {
    theMassLast.first  = theMSSM->mass(q2,particle1);
    theFLast.first  = f1ID;
  }
  if( q2 != theq2last || theFLast.second != f2ID) {
    theMassLast.second = theMSSM->mass(q2,particle2);
    theFLast.second = f2ID;
  }
  theq2last = q2;
  Complex coup(0.);
  Complex lcoup(1.),rcoup(1.);
  if( higgsID == ParticleID::h0 || higgsID == ParticleID::H0 ||
      higgsID == ParticleID::A0 ) {
    coup = 0.5*theMassLast.first/theMw;
    if( higgsID == ParticleID::h0 ) {
      if( f1ID % 2 == 0 )
	coup *= -theCa/theSb;
      else
	coup *=  theSa/theCb;
    }
    else if( higgsID == ParticleID::H0 ) {
      if( f1ID % 2 == 0 )
	coup *= -theSa/theSb;
      else
	coup *= -theCa/theCb;
    }
    else {
      if( f1ID % 2 == 0 )
	coup /= thetanb; 
      else
	coup *= thetanb;
      coup *= Complex(0.,-1.);
      rcoup = -1.;
    }
  }
  //H+
  else {
    if( f1ID % 2 == 0 ) {
      lcoup = theMassLast.first /thetanb/theMw;
      rcoup = theMassLast.second*thetanb/theMw;
    }
    else {
      lcoup = theMassLast.second/thetanb/theMw;
      rcoup = theMassLast.first *thetanb/theMw;
    }
    coup = sqrt(0.5);
    if( higgsID > 0 ) swap(lcoup,rcoup);
  }
  norm(theGlast*coup);
  left (lcoup);
  right(rcoup);
}

