// -*- C++ -*-
//
// SSHSFSFVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHSFSFVertex class.
//

#include "SSHSFSFVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHSFSFVertex::SSHSFSFVertex() : theMix(3), theTriC(9, complex<Energy>(ZERO)), 
				 theSinA(0.0),
				 theCosA(0.0), theSinB(0.0), theCosB(0.0),
				 theTanB(0.0), theSinAB(0.0), theCosAB(0.0),
				 theMw(ZERO), theMz(ZERO), theMu(ZERO), 
				 theSw(0.0), theCw(0.0), theCoupLast(ZERO),
				 theq2Last(ZERO), theHLast(0), theSF1Last(0),
				 theSF2Last(0) {
  vector<long> first,second,third;
  int higgs = 25;
  //h0,H0
  for(unsigned int i = 0; i < 2; ++i) {
    if( i == 1 ) higgs = 35;
    //quarks
    for(unsigned int j = 1; j < 7; ++j) {
      long lj = 1000000 + j;
      long rj = 2000000 + j;
      //LLbar
      first.push_back(higgs);
      second.push_back(lj);
      third.push_back(-lj);
      //RRbar
      first.push_back(higgs);
      second.push_back(rj);
      third.push_back(-rj);
      //LRbar
      first.push_back(higgs);
      second.push_back(lj);
      third.push_back(-rj);
      //RLbar
      first.push_back(higgs);
      second.push_back(rj);
      third.push_back(-lj);
    }
    for(unsigned int j = 11; j < 17; ++j) {
      long lj = 1000000 + j;
      long rj = 2000000 + j;
      first.push_back(higgs);
      second.push_back(lj);
      third.push_back(-lj);
      if( j % 2 != 0) {
	first.push_back(higgs);
	second.push_back(rj);
	third.push_back(-rj);
	//LRbar
	first.push_back(higgs);
	second.push_back(lj);
	third.push_back(-rj);
	//RLbar
	first.push_back(higgs);
	second.push_back(rj);
	third.push_back(-lj);
      }
    }
  }
  //A0
  for(unsigned int j = 1; j < 7; ++j) {
    long lj = 1000000 + j;
    long rj = 2000000 + j;
    //LRbar
    first.push_back(36);
    second.push_back(lj);
    third.push_back(-rj);
    //RLbar
    first.push_back(36);
    second.push_back(rj);
    third.push_back(-lj);
  }
  for(unsigned int j = 11; j < 17; j += 2) {
    long lj = 1000000 + j;
    long rj = 2000000 + j;
    first.push_back(36);
    second.push_back(lj);
    third.push_back(-rj);
    first.push_back(36);
    second.push_back(rj);
    third.push_back(-lj);
  }
  //outgoing H+
  for(long ii = 2; ii < 7; ii += 2) {
    //LL
    first.push_back(37);
    second.push_back( 999999 + ii);
    third.push_back(-1000000 - ii);
    //RR
    first.push_back(37);
    second.push_back(1999999 + ii);
    third.push_back(-2000000 - ii);
    //RL
    first.push_back(37);
    second.push_back(1999999 + ii);
    third.push_back(-1000000 - ii);
    //LR
    first.push_back(37);
    second.push_back( 999999 + ii);
    third.push_back(-2000000 - ii);
  }
  for(long ii = 11; ii < 17; ii += 2) {
    first.push_back(37);
    second.push_back(1000000 + ii);
    third.push_back(-1000001 - ii);
    first.push_back(37);
    second.push_back(2000000 + ii);
    third.push_back(-1000001 - ii);
  }
  //outgoing H-
  for(long ii = 2; ii < 7; ii += 2) {
    //LL
    first.push_back(-37);
    second.push_back(1000000 + ii);
    third.push_back(- 999999 - ii);
    //RR
    first.push_back(-37);
    second.push_back(2000000 + ii);
    third.push_back(-1999999 - ii);
    //RL
    first.push_back(-37);
    second.push_back(1000000 + ii);
    third.push_back(-1999999 - ii);
    //LR
    first.push_back(-37);
    second.push_back(2000000 + ii);
    third.push_back(- 999999 - ii);
  }
  for(long ii = 11; ii < 17; ii += 2) {
    first.push_back(-37);
    second.push_back(1000001 + ii);
    third.push_back(-1000000 - ii);
    first.push_back(-37);
    second.push_back(1000001 + ii);
    third.push_back(-2000000 - ii);
  }
  setList(first, second, third);
}

void SSHSFSFVertex::doinit() {
  SSSVertex::doinit();
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() << "SSHSFSFVertex::doinit - A problem occurred"
			  << "while trying to cast the SM pointer to "
			  << "a Susy one." << Exception::abortnow;
  //mixing vector should have been sized correctly already
  assert( theMix.size() == 3 );
  theMix[0] = theMSSM->stopMix();
  theMix[1] = theMSSM->sbottomMix();
  theMix[2] = theMSSM->stauMix();

  if(!theMix[0] || !theMix[1] || !theMix[2])
    throw InitException() << "SSHSFSFVertex::doinit -  "
			  << "A null mixing matrix pointer. stop: " << theMix[0] 
			  << " sbottom: " << theMix[1] << " stau: " << theMix[2]
			  << Exception::abortnow;
  //trilinear vector should have been sized correctly already
  assert( theTriC.size() == 9 );
  //vector has been zeroed in constructor
  theTriC[4]=theMSSM->bottomTrilinear().real();
  theTriC[5]=theMSSM->topTrilinear().real();
  theTriC[8]=theMSSM->tauTrilinear().real();
  //Masses
  theMw = getParticleData(ParticleID::Wplus)->mass();
  theMz = getParticleData(ParticleID::Z0)->mass();
  //parameters
  theSinA = sin(theMSSM->higgsMixingAngle());
  theCosA = sqrt(1. - sqr(theSinA));
  theTanB = theMSSM->tanBeta();
  theMu = theMSSM->muParameter();
  theSinB = theTanB/sqrt(1. + sqr(theTanB));
  theCosB = sqrt( 1. - sqr(theSinB) );
  theSinAB = theSinA*theCosB + theCosA*theSinB;
  theCosAB = theCosA*theCosB - theSinA*theSinB;

  theSw = sqrt(sin2ThetaW());
  theCw = sqrt(1. - sin2ThetaW());

  orderInGem(1);
  orderInGs(0);
}

void SSHSFSFVertex::persistentOutput(PersistentOStream & os) const {
  os << theMix << theSinA << theCosA << theSinB
     << theCosB << theTanB << ounit(theMu, GeV) << theSinAB << theCosAB 
     << ounit(theMw,GeV) << ounit(theMz,GeV) << theSw << theCw 
     << ounit(theTriC,GeV);

}

void SSHSFSFVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMix >>  theSinA >> theCosA >> theSinB
     >> theCosB >> theTanB >> iunit(theMu, GeV) >> theSinAB >> theCosAB 
     >> iunit(theMw,GeV) >> iunit(theMz,GeV) >> theSw >> theCw
     >> iunit(theTriC,GeV);
}

ClassDescription<SSHSFSFVertex> SSHSFSFVertex::initSSHSFSFVertex;
// Definition of the static class description member.

void SSHSFSFVertex::Init() {

  static ClassDocumentation<SSHSFSFVertex> documentation
    ("The coupling of an MSSM Higgs to a pair of sfermions.");

}

void SSHSFSFVertex::setCoupling(Energy2 q2, tcPDPtr particle1, 
				tcPDPtr particle2, tcPDPtr particle3) {
  long id1(abs(particle1->id())), id2(abs(particle2->id())),
    id3(abs(particle3->id())), higgsID(0), sq1ID(0), sq2ID(0);

  if( id1 == ParticleID::h0 || id1 == ParticleID::H0 || 
      id1 == ParticleID::A0 || id1 == ParticleID::Hplus ) {
    higgsID = id1;
    sq1ID = id2;
    sq2ID = id3;
  }
  else if( id2 == ParticleID::h0 || id2 == ParticleID::H0 || 
	   id2 == ParticleID::A0 || id2 == ParticleID::Hplus  ) {
    higgsID = id2;
    sq1ID = id1;
    sq2ID = id3;
  }
  else if( id3 == ParticleID::h0 || id3 == ParticleID::H0 || 
	   id3 == ParticleID::A0 || id3 == ParticleID::Hplus ) {
    higgsID = id3;
    sq1ID = id1;
    sq2ID = id2;
  }
  else {
    throw HelicityConsistencyError() 
      << "SSHSFSFVertex::setCoupling - There is no higgs particle in "
      << "this vertex. Particles: " << id1 << " " << id2 << " " << id3
      << Exception::warning;
    return;
  }
  assert( higgsID != 0 && sq1ID != 0 && sq2ID != 0);
  
  if( q2 != theq2Last || thegLast==0.) {
    thegLast = weakCoupling(q2);
      theq2Last = q2;
    }
  
  if( higgsID == theHLast && sq1ID == theSF1Last && sq2ID == theSF2Last) {
    setNorm(thegLast*theCoupLast*UnitRemoval::InvE);
    return;
  }
  theHLast = higgsID;
  theSF1Last = sq1ID;
  theSF2Last = sq2ID;

  if( higgsID == ParticleID::Hplus ) 
    chargedHiggs(sq1ID, sq2ID);
  else {
    long smID(0);
    unsigned int alpha(sq1ID/1000000 - 1), beta(sq2ID/1000000 - 1);
    if( sq1ID/1000000 == 2 )
      smID = sq1ID - 2000000;
    else
      smID = sq1ID - 1000000;
    if( smID < 7  ) {
      if( smID % 2 == 0 ) 
	upSF(higgsID, smID, alpha, beta);
      else 
	downSF(higgsID, smID, alpha, beta);
    }
    else 
      leptonSF(higgsID, smID, alpha, beta);
  
  }
 
  setNorm(thegLast*theCoupLast*UnitRemoval::InvE);
}

 void SSHSFSFVertex::downSF(long higgs, long smID, 
			    unsigned int alpha, unsigned int beta) {
   assert( smID < 7 && smID % 2 != 0);
   Energy fmass = getParticleData(smID)->mass();
   double mfacta = 0.5*fmass/theMw;
   
   if( higgs == ParticleID::A0 ){
     theCoupLast = -Complex(1.,0.)*mfacta*(theTriC[smID - 1]*theTanB + theMu);
     return;
   }
   Energy mfactb = sqr(fmass)/theMw/theCosB;
   Energy facta = theMz/theCw;
   double factb = generator()->standardModel()->ed()*theSw*theSw;
   //mixing parameters
   Complex q1a(0.), q1b(0.), q2a(0.), q2b(0.);  
   if( smID == 1 || smID == 3) {
     q1a = (alpha == 0) ? 1.0 : 0.0;
     q1b = (beta == 0) ? 1.0 : 0.0;
     q2a = (alpha == 0) ? 0.0 : 1.0;
     q2b = (beta == 0) ? 0.0 : 1.0;
   }
   else {
     q1a = (*theMix[1])(alpha, 0);
     q1b = (*theMix[1])(beta, 0);
     q2a = (*theMix[1])(alpha, 1);
     q2b = (*theMix[1])(beta, 1);
   }
   Complex fbrac = (q1a*q1b*(0.5 + factb) - factb*q2a*q2b);
   Complex sbrac = (q1a*q1b + q2a*q2b);
   Complex tbrac = (q2a*q1b + q1a*q2b);
   if( higgs == ParticleID::h0) {
     theCoupLast = -facta*theSinAB*fbrac + mfactb*theSinA*sbrac
       + mfacta*(theTriC[smID - 1]*theSinA + theMu*theCosA)*tbrac/theCosB;
   }
   else if( higgs == ParticleID::H0 ) {
     theCoupLast = facta*theCosAB*fbrac - mfactb*theCosA*sbrac 
       + mfacta*(theMu*theSinA - theTriC[smID - 1]*theCosA)*tbrac/theCosB;
   }
   else
     throw HelicityConsistencyError() 
      << "SSHSFSFVertex::downSF - Unrecognised higgs particle " 
      << higgs << Exception::warning;
   
 }

void SSHSFSFVertex::upSF(long higgs, long smID, 
			 unsigned int alpha, unsigned int beta) {
  assert( smID < 7 && smID % 2 == 0);
  Energy fmass = getParticleData(smID)->mass();
  double mfacta = 0.5*fmass/theMw;
  if( higgs == ParticleID::A0 ){
    theCoupLast = -Complex(1.,0.)*mfacta*(theTriC[smID - 1]/theTanB + theMu);
    return;
  }
  Energy mfactb = sqr(fmass)/theMw/theSinB;
  Energy facta = theMz/theCw;
  double factb = generator()->standardModel()->eu()*theSw*theSw;
  //mixing parameters
  Complex q1a(0.), q1b(0.), q2a(0.), q2b(0.);  
  if( smID == 2 || smID == 4) {
    q1a = (alpha == 0) ? 1.0 : 0.0;
    q1b = (beta == 0) ? 1.0 : 0.0;
    q2a = (alpha == 0) ? 0.0 : 1.0;
    q2b = (beta == 0) ? 0.0 : 1.0;
  }
  else {
    q1a = (*theMix[2])(alpha, 0);
    q1b = (*theMix[2])(beta, 0);
    q2a = (*theMix[2])(alpha, 1);
    q2b = (*theMix[2])(beta, 1);
  }
  Complex fbrac = (q1a*q1b*(0.5 - factb) - factb*q2a*q2b);
  Complex sbrac = (q1a*q1b + q2a*q2b);
  Complex tbrac = (q2a*q1b + q1a*q2b);
  if( higgs == ParticleID::h0) {
    theCoupLast = facta*theSinAB*fbrac - mfactb*theCosA*sbrac
      - mfacta*(theTriC[smID - 1]*theCosA + theMu*theSinA)*tbrac/theSinB;
  }
  else if( higgs == ParticleID::H0 ) {
    theCoupLast = -facta*theCosAB*fbrac - mfactb*theSinA*sbrac 
      - mfacta*(theTriC[smID - 1]*theSinA - theMu*theCosA )*tbrac/theSinB;
  }
  else
    throw HelicityConsistencyError() 
      << "SSHSFSFVertex::upSF - Unrecognised higgs particle " 
      << higgs << Exception::warning;
}

void SSHSFSFVertex::leptonSF(long higgs, long smID, 
			     unsigned int alpha, unsigned int beta) {
  assert( smID >= 11 && smID <= 16 ); 
  Energy facta = theMz/theCw;
  if( smID % 2 == 0 ) {
    theCoupLast = complex<Energy>(facta/2.);
    if( higgs == ParticleID::h0) 
      theCoupLast *= theSinAB;
    else
      theCoupLast *= -theCosAB;
    return;
  }
  Energy fmass = getParticleData(smID)->mass();
  double mfacta = fmass/2./theMw;
  if( higgs == ParticleID::A0 ) {
    theCoupLast = -Complex(0.,1.)*mfacta*(theTriC[(smID + 1)/2]*theTanB + theMu);
    return;
  }
  Energy mfactb = fmass*fmass/theMw/theCosB;
  double factb = theSw*theSw;
   //mixing parameters
   Complex l1a(0.), l1b(0.), l2a(0.), l2b(0.);  
   if( smID == 11 || smID == 13) {
     l1a = (alpha == 0) ? 1.0 : 0.0;
     l1b = (beta == 0) ? 1.0 : 0.0;
     l2a = (alpha == 0) ? 0.0 : 1.0;
     l2b = (beta == 0) ? 0.0 : 1.0;
   }
   else {
     l1a = (*theMix[2])(alpha, 0);
     l1b = (*theMix[2])(beta, 0);
     l2a = (*theMix[2])(alpha, 1);
     l2b = (*theMix[2])(beta, 1);
   }
   Complex fbrac = (l1a*l1b*(0.5 - factb) + factb*l2a*l2b);
   Complex sbrac = (l1a*l1b + l2a*l2b);
   Complex tbrac = (l2a*l1b + l1a*l2b);
   if( higgs == ParticleID::h0) {
     theCoupLast = -facta*theSinAB*fbrac + mfactb*theSinA*sbrac
       + mfacta*(theTriC[(smID + 1)/2]*theSinA + theMu*theCosA)*tbrac/theCosB;
   }
   else if( higgs == ParticleID::H0 ) {
     theCoupLast = facta*theCosAB*fbrac - mfactb*theCosA*sbrac 
       + mfacta*(theMu*theSinA - theTriC[(smID + 1)/2]*theCosA)*tbrac/theCosB;
   }
   else
     throw HelicityConsistencyError() 
      << "SSHSFSFVertex::leptonSF - Unrecognised higgs particle " 
      << higgs << Exception::warning;
  
  
}

void SSHSFSFVertex::chargedHiggs(long id1, long id2) {
  //have id1 as up-type
  if( id1 % 2 != 0)
    swap(id1, id2);
  unsigned int beta(0);
  if( id2/1000000 == 2 )
    beta = 1;
  else 
    beta = 0;

  long smdID = (beta == 0) ? id2 - 1000000 : id2 - 2000000;
  Energy mfd = getParticleData(smdID)->mass();
  Energy2 facta = sqr(theMw)*2.*theSinB*theCosB;
  if( smdID == 11 || smdID == 13 || smdID == 15) {
    Complex l1b = (beta == 0) ? 1.0 : 0.0;
    Complex l2b = (beta == 0) ? 0.0 : 1.0;
    if( smdID == 15 ) {
      l1b = (*theMix[2])(beta, 0);
      l2b = (*theMix[2])(beta, 1);
    }
    theCoupLast = ( l1b*(mfd*mfd*theTanB - facta) 
		    + l2b*mfd*(theTriC[(smdID + 1)/2]*theTanB + theMu)
		   )/theMw/sqrt(2.);
  }
  else {
    unsigned int alpha(0);
    if( id1/1000000 == 2 )
    alpha = 1;
  else 
    alpha = 0;

  long smuID = (alpha == 0) ? id1 - 1000000 : id1 - 2000000;
  Energy mfu = getParticleData(smuID)->mass();
  Complex q1a(0.0), q1b(0.0), q2a(0.0), q2b(0.0);
  if( smuID == 2 || smuID == 4 ) {
    q1a = (alpha == 0) ? 1.0 : 0.0;
    q2a = (alpha == 0) ? 0.0 : 1.0;
  }
  else {
    q1a = (*theMix[0])(alpha, 0);
    q2a = (*theMix[0])(1, alpha);
  }
  if( smdID == 1 || smdID == 3 ) {
    q1b = (beta == 0) ? 1.0 : 0.0;
    q2b = (beta == 0) ? 0.0 : 1.0;
  }
  else {
    q1b = (*theMix[1])(0, beta);
    q2b = (*theMix[1])(1, beta);
  }
  
  theCoupLast = ( q1a*q1b*(mfd*mfd*theTanB + mfu*mfu/theTanB - facta)
		  + q2a*q2b*mfu*mfd*(theTanB + (1./theTanB))
		  + q1a*q1b*mfd*(theTriC[smdID - 1]*theTanB + theMu)
		  + q2a*q1b*mfu*(theMu + theTriC[(smuID + 1)/2]/theTanB)
		 )/theMw/sqrt(2.);
  }
}
