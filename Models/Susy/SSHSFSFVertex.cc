// -*- C++ -*-
//
// SSHSFSFVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
  orderInGem(1);
  orderInGs(0);
}

void SSHSFSFVertex::doinit() {
  int higgs = 25;
  //h0,H0
  for(unsigned int i = 0; i < 2; ++i) {
    if( i == 1 ) higgs = 35;
    //quarks
    for(unsigned int j = 1; j < 7; ++j) {
      long lj = 1000000 + j;
      long rj = 2000000 + j;
      //LLbar
      addToList(higgs,lj,-lj);
      //RRbar
      addToList(higgs,rj,-rj);
      //LRbar
      addToList(higgs,lj,-rj);
      //RLbar
      addToList(higgs,rj,-lj);
    }
    for(unsigned int j = 11; j < 17; ++j) {
      long lj = 1000000 + j;
      long rj = 2000000 + j;
      addToList(higgs,lj,-lj);
      if( j % 2 != 0) {
	addToList(higgs,rj,-rj);
	//LRbar
	addToList(higgs,lj,-rj);
	//RLbar
	addToList(higgs,rj,-lj);
      }
    }
  }
  //A0
  for(unsigned int j = 1; j < 7; ++j) {
    long lj = 1000000 + j;
    long rj = 2000000 + j;
    //LRbar
    addToList(36,lj,-rj);
    //RLbar
    addToList(36,rj,-lj);
  }
  for(unsigned int j = 11; j < 17; j += 2) {
    long lj = 1000000 + j;
    long rj = 2000000 + j;
    addToList(36,lj,-rj);
    addToList(36,rj,-lj);
  }
  //outgoing H+
  for(long ii = 2; ii < 7; ii += 2) {
    //LL
    addToList(37, 999999 + ii, -1000000 - ii);
    //RR
    addToList(37, 1999999 + ii, -2000000 - ii);
    //RL
    addToList(37, 1999999 + ii, -1000000 - ii);
    //LR
    addToList(37, 999999 + ii, -2000000 - ii);
  }
  for(long ii = 11; ii < 17; ii += 2) {
    addToList(37, 1000000 + ii, -1000001 - ii);
    addToList(37, 2000000 + ii, -1000001 - ii);
  }
  //outgoing H-
  for(long ii = 2; ii < 7; ii += 2) {
    //LL
    addToList(-37, 1000000 + ii, -999999 - ii);
    //RR
    addToList(-37, 2000000 + ii, -1999999 - ii);
    //RL
    addToList(-37, 1000000 + ii, -1999999 - ii);
    //LR
    addToList(-37, 2000000 + ii, -999999 - ii);
  }
  for(long ii = 11; ii < 17; ii += 2) {
    addToList(-37, 1000001 + ii, -1000000 - ii);
    addToList(-37, 1000001 + ii, -2000000 - ii);}
  SSSVertex::doinit();
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
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
  theCosA = cos(theMSSM->higgsMixingAngle());
  theTanB = theMSSM->tanBeta();
  theMu = theMSSM->muParameter();
  theSinB = theTanB/sqrt(1. + sqr(theTanB));
  theCosB = sqrt( 1. - sqr(theSinB) );
  theSinAB = theSinA*theCosB + theCosA*theSinB;
  theCosAB = theCosA*theCosB - theSinA*theSinB;

  theSw = sqrt(sin2ThetaW());
  theCw = sqrt(1. - sin2ThetaW());
}

void SSHSFSFVertex::persistentOutput(PersistentOStream & os) const {
  os << theMix << theSinA << theCosA << theSinB << theMSSM
     << theCosB << theTanB << ounit(theMu, GeV) << theSinAB << theCosAB 
     << ounit(theMw,GeV) << ounit(theMz,GeV) << theSw << theCw 
     << ounit(theTriC,GeV);

}

void SSHSFSFVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMix >>  theSinA >> theCosA >> theSinB >> theMSSM
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

void SSHSFSFVertex::setCoupling(Energy2 q2, tcPDPtr part1, 
				tcPDPtr part2, tcPDPtr part3) {
  // extract particle ids
  long higgs(part1->id()), sq1(part2->id()), sq2(part3->id());
  // higgs first
  if(abs(sq1)<100) swap(higgs,sq1);
  if(abs(sq2)<100) swap(higgs,sq2);
  // squark second, antisquark third
  if(sq1<0) swap(sq1,sq2);
  assert( higgs == 25 || higgs == 35 || 
	  higgs == 36 || abs(higgs) == 37 );
  sq2 *=-1; 
  assert(sq1>0&&sq2>0);
  
  // running coupling
  if( q2 != theq2Last || thegLast==0.) {
    thegLast = weakCoupling(q2);
    theq2Last = q2;
  }
  
  if( higgs == theHLast && sq1 == theSF1Last && sq2 == theSF2Last) {
    norm(thegLast*theCoupLast*UnitRemoval::InvE);
    return;
  }
  theHLast = higgs;
  theSF1Last = sq1;
  theSF2Last = sq2;
  if( abs(higgs) == ParticleID::Hplus ) 
    chargedHiggs(q2,sq1, sq2);
  else {
    long sm(0);
    unsigned int alpha(sq1/1000000 - 1), beta(sq2/1000000 - 1);
    if( sq1/1000000 == 2 )
      sm = sq1 - 2000000;
    else
      sm = sq1 - 1000000;
    if( sm < 7  ) {
      if( sm % 2 == 0 ) 
	upSF(q2,higgs, sm, alpha, beta);
      else 
	downSF(q2,higgs, sm, alpha, beta);
    }
    else 
      leptonSF(q2,higgs, sm, alpha, beta);
  
  }
 
  norm(thegLast*theCoupLast*UnitRemoval::InvE);
}

void SSHSFSFVertex::downSF(Energy2 q2, long higgs, long sm, 
			   unsigned int alpha, unsigned int beta) {
   assert( sm < 7 && sm % 2 != 0);
   Energy fmass = theMSSM->mass(q2,getParticleData(sm));
   double mfacta = 0.5*fmass/theMw;
   
   if( higgs == ParticleID::A0 ) {
     theCoupLast = -Complex(0.,1.)*mfacta*(theTriC[sm - 1]*theTanB + theMu);
     if(sm == 5) {
      // do not revert to *=, breaks with XCode 5.1
       theCoupLast = theCoupLast * (
	 (*theMix[1])(alpha, 1) * (*theMix[1])(beta , 0) -
	 (*theMix[1])(alpha, 0) * (*theMix[1])(beta , 1) );
     }
     else if(alpha<beta) theCoupLast *= -1.;
     return;
   }
   Energy mfactb = sqr(fmass)/theMw/theCosB;
   Energy facta = theMz/theCw;
   double factb = theMSSM->ed()*theSw*theSw;
   //mixing parameters
   Complex q1a(0.), q1b(0.), q2a(0.), q2b(0.);  
   if( sm == 1 || sm == 3) {
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
       + mfacta*(theTriC[sm - 1]*theSinA + theMu*theCosA)*tbrac/theCosB;
   }
   else if( higgs == ParticleID::H0 ) {
     theCoupLast = facta*theCosAB*fbrac - mfactb*theCosA*sbrac 
       + mfacta*(theMu*theSinA - theTriC[sm - 1]*theCosA)*tbrac/theCosB;
   }
   else
     throw HelicityConsistencyError() 
      << "SSHSFSFVertex::downSF - Unrecognised higgs particle " 
      << higgs << Exception::warning;
   
 }

void SSHSFSFVertex::upSF(Energy2 q2, long higgs, long sm, 
			 unsigned int alpha, unsigned int beta) {
  assert( sm < 7 && sm % 2 == 0);
  Energy fmass = theMSSM->mass(q2,getParticleData(sm));
  double mfacta = 0.5*fmass/theMw;
  if( higgs == ParticleID::A0 ){
    theCoupLast = -Complex(0.,1.)*mfacta*(theTriC[sm - 1]/theTanB + theMu);
    if(sm == 6) {
    // do not revert to *=, breaks with XCode 5.1
      theCoupLast = theCoupLast * (
	(*theMix[0])(alpha, 1) * (*theMix[0])(beta , 0) -
	(*theMix[0])(alpha, 0) * (*theMix[0])(beta , 1) );
    }
    else if(alpha<beta) theCoupLast *= -1.;
    return;
  }
  Energy mfactb = sqr(fmass)/theMw/theSinB;
  Energy facta = theMz/theCw;
  double factb = theMSSM->eu()*sqr(theSw);
  //mixing parameters
  Complex q1a(0.), q1b(0.), q2a(0.), q2b(0.);  
  if( sm == 2 || sm == 4) {
    q1a = (alpha == 0) ? 1.0 : 0.0;
    q1b = (beta  == 0) ? 1.0 : 0.0;
    q2a = (alpha == 0) ? 0.0 : 1.0;
    q2b = (beta  == 0) ? 0.0 : 1.0;
  }
  else {
    q1a = (*theMix[0])(alpha, 0);
    q1b = (*theMix[0])(beta , 0);
    q2a = (*theMix[0])(alpha, 1);
    q2b = (*theMix[0])(beta , 1);
  }
  Complex fbrac = (q1a*q1b*(0.5 - factb) + factb*q2a*q2b);
  Complex sbrac = (q1a*q1b + q2a*q2b);
  Complex tbrac = (q2a*q1b + q1a*q2b);
  if( higgs == ParticleID::h0) {
    theCoupLast = facta*theSinAB*fbrac - mfactb*theCosA*sbrac
      - mfacta*(theTriC[sm - 1]*theCosA + theMu*theSinA)*tbrac/theSinB;

  }
  else if( higgs == ParticleID::H0 ) {
    theCoupLast = -facta*theCosAB*fbrac - mfactb*theSinA*sbrac 
      - mfacta*(theTriC[sm - 1]*theSinA - theMu*theCosA )*tbrac/theSinB;

  }
  else
    assert(false);
}

void SSHSFSFVertex::leptonSF(Energy2 q2, long higgs, long sm, 
			     unsigned int alpha, unsigned int beta) {
  assert( sm >= 11 && sm <= 16 ); 
  Energy facta = theMz/theCw;
  if( sm % 2 == 0 ) {
    theCoupLast = complex<Energy>(facta/2.);
    if( higgs == ParticleID::h0) 
      theCoupLast *=  theSinAB;
    else
      theCoupLast *= -theCosAB;
    return;
  }
  Energy fmass = theMSSM->mass(q2,getParticleData(sm));
  double mfacta = fmass/2./theMw;
  if( higgs == ParticleID::A0 ) {
    theCoupLast = -Complex(0.,1.)*mfacta*(theTriC[(sm + 1)/2]*theTanB + theMu);
    if(sm == 15) {
      // do not revert to *=, breaks with XCode 5.1
      theCoupLast = theCoupLast * (
	(*theMix[2])(alpha, 1) * (*theMix[2])(beta , 0) -
	(*theMix[2])(alpha, 0) * (*theMix[2])(beta , 1) );
    }
    else if(alpha<beta) theCoupLast *= -1.;
    return;
  }
  Energy mfactb = fmass*fmass/theMw/theCosB;
  double factb = theSw*theSw;
   //mixing parameters
   Complex l1a(0.), l1b(0.), l2a(0.), l2b(0.);  
   if( sm == 11 || sm == 13) {
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
       + mfacta*(theTriC[(sm + 1)/2]*theSinA + theMu*theCosA)*tbrac/theCosB;
   }
   else if( higgs == ParticleID::H0 ) {
     theCoupLast = facta*theCosAB*fbrac - mfactb*theCosA*sbrac 
       + mfacta*(theMu*theSinA - theTriC[(sm + 1)/2]*theCosA)*tbrac/theCosB;
   }
   else
     throw HelicityConsistencyError() 
      << "SSHSFSFVertex::leptonSF - Unrecognised higgs particle " 
      << higgs << Exception::warning;
  
  
}

void SSHSFSFVertex::chargedHiggs(Energy2 q2, long id1, long id2) {
  //have id1 as up-type
  if( id1 % 2 != 0)
    swap(id1, id2);
  unsigned int beta(0);
  if( id2/1000000 == 2 )
    beta = 1;
  else 
    beta = 0;

  long smd = (beta == 0) ? id2 - 1000000 : id2 - 2000000;
  Energy mfd = theMSSM->mass(q2,getParticleData(smd));
  Energy2 facta = sqr(theMw)*2.*theSinB*theCosB;
  if( smd == 11 || smd == 13 || smd == 15) {
    Complex l1b = (beta == 0) ? 1.0 : 0.0;
    Complex l2b = (beta == 0) ? 0.0 : 1.0;
    if( smd == 15 ) {
      l1b = (*theMix[2])(beta, 0);
      l2b = (*theMix[2])(beta, 1);
    }
    theCoupLast = ( l1b*(mfd*mfd*theTanB - facta) 
		    + l2b*mfd*(theTriC[(smd + 1)/2]*theTanB + theMu)
		   )/theMw/sqrt(2.);
  }
  else {
    unsigned int alpha(0);
    if( id1/1000000 == 2 )
    alpha = 1;
  else 
    alpha = 0;

  long smu = (alpha == 0) ? id1 - 1000000 : id1 - 2000000;
  Energy mfu = theMSSM->mass(q2,getParticleData(smu));
  Complex q1a(0.0), q1b(0.0), q2a(0.0), q2b(0.0);
  if( smu == 2 || smu == 4 ) {
    q1a = (alpha == 0) ? 1.0 : 0.0;
    q2a = (alpha == 0) ? 0.0 : 1.0;
  }
  else {
    q1a = (*theMix[0])(alpha,0);
    q2a = (*theMix[0])(alpha,1);
  }
  if( smd == 1 || smd == 3 ) {
    q1b = (beta == 0) ? 1.0 : 0.0;
    q2b = (beta == 0) ? 0.0 : 1.0;
  }
  else {
    q1b = (*theMix[1])(beta,0);
    q2b = (*theMix[1])(beta,1);
  }
  
  theCoupLast = ( q1a*q1b*(mfd*mfd*theTanB + mfu*mfu/theTanB - facta)
		  + q2a*q2b*mfu*mfd*(theTanB + (1./theTanB))
		  + q1a*q2b*mfd*(theTriC[smd - 1]*theTanB + theMu)
		  + q2a*q1b*mfu*(theMu + theTriC[smu-1]/theTanB)
		 )/theMw/sqrt(2.);
  }
}
