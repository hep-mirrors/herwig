// -*- C++ -*-
//
// SSGOGOHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGOGOHVertex class.
//

#include "SSGOGOHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGOGOHVertex::SSGOGOHVertex() : theMw(), theSij(2, {0., 0.}),
				 theQij(2, {0., 0.}),
				 theQijLp(4, {0., 0.}),
				 theQijRp(4, {0., 0.}),
				 theSijdp(4, {0., 0., 0., 0.}),
				 theQijdp(4, {0., 0., 0., 0.}),
				 theSa(0.0), theSb(0.0),
				 theCa(0.0), theCb(0.0), theCoupLast(0.0),
				 theLLast(0.0), theRLast(0.0), theHLast(0),
				 theID1Last(0), theID2Last(0), theq2last() {
  orderInGem(1);
  orderInGs(0);
}

void SSGOGOHVertex::doinit() {
  long neu[4] = {1000022, 1000023, 1000025, 1000035};
  long chg[2] = {1000024, 1000037};
  long higgs[3] =  {25, 35, 36};
  for(unsigned int i = 0; i < 4; ++i) {
    //neutralinos
    for(unsigned int j = 0; j < 4; ++j) {
      for(unsigned int k = 0; k < 4; ++k) {
	if( i < 3 ) {
	  if(k<=j)
	    addToList(neu[j], neu[k], higgs[i]);
	  //charginos
	  if( j < 2 && k < 2 ) {
	    addToList(-chg[j], chg[k], higgs[i]);
	  }
	} 
	else {
	  if( k == 2 ) break;
	  addToList(-chg[k], neu[j], 37);
	  addToList( chg[k], neu[j],-37);
	}
      }
    }
  }
  FFSVertex::doinit();
  
  tMSSMPtr theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSGOGOHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  
  theMw = getParticleData(ParticleID::Wplus)->mass();
  double theSw = sqrt(sin2ThetaW());
  double tw = theSw/sqrt(1. - theSw*theSw);
  double tanb = theMSSM->tanBeta();
  theSb = tanb/sqrt(1. + sqr(tanb));
  theCb = sqrt( 1. - sqr(theSb) );
  theSa = sin(theMSSM->higgsMixingAngle());
  theCa = sqrt(1. - sqr(theSa));
  MixingMatrix nmix = *theMSSM->neutralinoMix();
  MixingMatrix umix = *theMSSM->charginoUMix();
  MixingMatrix vmix = *theMSSM->charginoVMix();

  for(unsigned int i = 0; i < 4; ++i) {
    for(unsigned int j = 0; j < 4; ++j) {
      if( i < 2 && j < 2 ) { 
	theQij[i][j] = vmix(i,0)*umix(j,1)/sqrt(2);
	theSij[i][j] = vmix(i,1)*umix(j,0)/sqrt(2);
      }
      if( j < 2 ) {
	theQijLp[i][j] = conj(nmix(i, 3)*vmix(j,0) 
			      + (nmix(i,1) + nmix(i,0)*tw)*vmix(j,1)/sqrt(2));
	theQijRp[i][j] = nmix(i, 2)*umix(j,0) 
	  - (nmix(i,1) + nmix(i,0)*tw)*umix(j,1)/sqrt(2);
      }
      theQijdp[i][j] = 0.5*( nmix(i,2)*( nmix(j,1) - tw*nmix(j,0) )
			     + nmix(j,2)*( nmix(i,1) - tw*nmix(i,0) ) );
      theSijdp[i][j] = 0.5*( nmix(i,3)*( nmix(j,1) - tw*nmix(j,0) )
			     + nmix(j,3)*( nmix(i,1) - tw*nmix(i,0) ) );
    }
  }
}

void SSGOGOHVertex::persistentOutput(PersistentOStream & os) const {
  os << theSij << theQij << theQijLp << theQijRp << theSijdp
     << theQijdp << ounit(theMw,GeV) << theSa << theSb << theCa 
     << theCb;
}

void SSGOGOHVertex::persistentInput(PersistentIStream & is, int) {
  is >> theSij >> theQij >> theQijLp >> theQijRp >> theSijdp
     >> theQijdp >> iunit(theMw,GeV) >> theSa >> theSb >> theCa 
     >> theCb;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSGOGOHVertex,FFSVertex>
describeHerwigSSGOGOHVertex("Herwig::SSGOGOHVertex", "HwSusy.so");

void SSGOGOHVertex::Init() {

  static ClassDocumentation<SSGOGOHVertex> documentation
    ("The coupling of the higgs bosons to SM fermions in the MSSM");

}

/// \todo fixme
void SSGOGOHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, 
				tcPDPtr particle2,tcPDPtr particle3) {
  long f1ID(particle1->id()), f2ID(particle2->id()), higgsID(particle3->id());
  assert(higgsID == ParticleID::h0 ||     higgsID  == ParticleID::H0 ||
	 higgsID == ParticleID::A0 || abs(higgsID) == ParticleID::Hplus);
  if( f1ID < 0 ) swap(f1ID, f2ID);
  
  if( q2 != theq2last || theCoupLast == 0. ) {
    theCoupLast = weakCoupling(q2);
    theq2last = q2;
  }
  if( higgsID == theHLast && f1ID == theID1Last && f2ID == theID2Last) {
    norm(theCoupLast);
    left(theLLast);
    right(theRLast);
    return;
  }
  theHLast = higgsID;
  theID1Last = f1ID;
  theID2Last = f2ID;
  
  if( higgsID == ParticleID::h0 ) {
    //charginos
    if( abs(f2ID) == ParticleID::SUSY_chi_1plus ||
	abs(f2ID) == ParticleID::SUSY_chi_2plus ) {
      unsigned int ei = (abs(f1ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      unsigned int ej = (abs(f2ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      theLLast =  conj(theQij[ej][ei])*theSa - conj(theSij[ej][ei])*theCa;
      theRLast = theQij[ei][ej]*theSa - theSij[ei][ej]*theCa;
    }
    //neutralinos
    else {
      unsigned int ei(f1ID - ParticleID::SUSY_chi_10), 
	ej(f2ID - ParticleID::SUSY_chi_10);
      if( ei > 1 )
	ei = ( ei == 13 ) ? 3 : 2;
      if( ej > 1 )
	ej = ( ej == 13 ) ? 3 : 2;
      theLLast = conj(theQijdp[ej][ei])*theSa + conj(theSijdp[ej][ei])*theCa;
      theRLast = theQijdp[ei][ej]*theSa + theSijdp[ei][ej]*theCa ;
    }
    
  }
  else if( higgsID == ParticleID::H0 ) {
    //charginos
    if( abs(f2ID) == ParticleID::SUSY_chi_1plus ||
	abs(f2ID) == ParticleID::SUSY_chi_2plus ) {
      unsigned int ei = (abs(f1ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      unsigned int ej = (abs(f2ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      theLLast =  -conj(theQij[ej][ei])*theCa - conj(theSij[ej][ei])*theSa;
      theRLast = -theQij[ei][ej]*theCa - theSij[ei][ej]*theSa;
    }
    //neutralinos
    else {
      unsigned int ei(f1ID - ParticleID::SUSY_chi_10), 
	ej(f2ID - ParticleID::SUSY_chi_10);
      if( ei > 1 )
	ei = ( ei == 13 ) ? 3 : 2;
      if( ej > 1 )
	ej = ( ej == 13 ) ? 3 : 2;
      
      theLLast = -conj(theQijdp[ej][ei])*theCa + conj(theSijdp[ej][ei])*theSa;
      theRLast = -theQijdp[ei][ej]*theCa + theSijdp[ei][ej]*theSa;
    }
  }
  else if( higgsID == ParticleID::A0 ) {
    if( abs(f2ID) == ParticleID::SUSY_chi_1plus ||
	abs(f2ID) == ParticleID::SUSY_chi_2plus ) {
      unsigned int ei = (abs(f1ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      unsigned int ej = (abs(f2ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;

      theLLast = Complex(0.,1.)*( conj(theQij[ej][ei])*theSb 
				  + conj(theSij[ej][ei])*theCb );
      theRLast = -Complex(0.,1.)*(theQij[ei][ej]*theSb + theSij[ei][ej]*theCb);
    }
    //neutralinos
    else {
      unsigned int ei(f1ID - ParticleID::SUSY_chi_10), 
	ej(f2ID - ParticleID::SUSY_chi_10);
      if( ei > 1 )
	ei = ( ei == 13 ) ? 3 : 2;
      if( ej > 1 )
	ej = ( ej == 13 ) ? 3 : 2;

      theLLast = Complex(0.,1.)*( conj(theQijdp[ej][ei])*theSb 
				  - conj(theSijdp[ej][ei])*theCb );
      theRLast = -Complex(0.,1.)*(theQijdp[ei][ej]*theSb - theSijdp[ei][ej]*theCb);
    }
  }
  //charged higgs
  else {
    unsigned int ei(0),ej(0);
    long chg(f2ID), neu(f1ID);
    if( abs(neu) == ParticleID::SUSY_chi_1plus || 
	abs(neu) == ParticleID::SUSY_chi_2plus ) swap(chg, neu);
    ej = ( abs(chg) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
    ei = neu - ParticleID::SUSY_chi_10;
    if( ei > 1 ) ei = ( ei == 13 ) ? 3 : 2;
    theLLast = -theQijLp[ei][ej]*theCb;
    theRLast = -theQijRp[ei][ej]*theSb;
    if( higgsID < 0 ) {
      Complex tmp = theLLast;
      theLLast = conj(theRLast);
      theRLast = conj(tmp);
    }
  }
  norm(theCoupLast);
  left(theLLast);
  right(theRLast);
}

