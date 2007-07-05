// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGOGOHVertex class.
//

#include "SSGOGOHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace Herwig::Helicity;

SSGOGOHVertex::SSGOGOHVertex() : theMw(0.*MeV), theRij(2, vector<Complex>(2,0.0)),
				 theQij(2, vector<Complex>(2,0.0)),
				 theQijLp(4, vector<Complex>(2,0.0)),
				 theQijRp(4, vector<Complex>(2,0.0)),
				 theRijdp(4, vector<Complex>(4,0.0)),
				 theQijdp(4, vector<Complex>(4,0.0)),
				 theSw(0.0), theSa(0.0), theSb(0.0),
				 theCa(0.0), theCb(0.0),theC2b(0.0),
				 theSba(0.), theCba(0.0), theCoupLast(0.0), 
				 theLLast(0.0), theRLast(0.0) {
  vector<int> first, second, third;
  int neu[4] = {1000022, 1000023, 1000025, 1000035};
  int chg[2] = {1000024, 1000037};
  int higgs[3] =  {25, 35, 36};
  for(unsigned int i = 0; i < 4; ++i) {
    //neutralinos
    for(unsigned int j = 0; j < 4; ++j) {
      for(unsigned int k = 0; k < 4; ++k) {
	if( i < 3 ) {
	  first.push_back(neu[j]);
	  second.push_back(neu[k]);
	  third.push_back(higgs[i]);
	  //charginos
	  if( j < 2 && k < 2 ) {
	    first.push_back(-chg[j]);
	    second.push_back(chg[k]);
	    third.push_back(higgs[i]);
	  }
	} 
	else {
	  if( k == 2 ) break;
	  first.push_back(-chg[j]);
	  second.push_back(neu[k]);
	  third.push_back(37);
	  first.push_back(chg[j]);
	  second.push_back(neu[k]);
	  third.push_back(-37);
	}
      }
    }
  }
  setList(first, second, third);
}

SSGOGOHVertex::~SSGOGOHVertex() {}

void SSGOGOHVertex::doinit() throw(InitException) {
  FFSVertex::doinit();
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM )
    throw InitException() 
      << "SSGOGOHVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  
  theMw = getParticleData(ParticleID::Wplus)->mass();
  theSw = sqrt(theMSSM->sin2ThetaW());
  double tw = theSw/(1. - theSw*theSw);
  double tanb = theMSSM->tanBeta();
  Energy mu = theMSSM->muParameter();
  theSb = tanb/sqrt(1. + sqr(tanb));
  theCb = sqrt( 1. - sqr(tanb) );
  theSa = sin(theMSSM->higgsMixingAngle());
  theCa = sqrt(1. - sqr(theSa));
  theC2b = theCb*theCb - theSb*theSb;
  
  MixingMatrix nmix = *theMSSM->neutralinoMix();
  MixingMatrix umix = *theMSSM->charginoUMix();
  MixingMatrix vmix = *theMSSM->charginoVMix();

  Energy mOne = theMSSM->softMOne();
  Energy mTwo = theMSSM->softMTwo();

  for(unsigned int i = 0; i < 4; ++i) {
    for(unsigned int j = 0; j < 4; ++j) {
      if( i < 2 && j < 2 ) { 
	theRij[i][j] = 
	  ( mTwo*umix(i,0)*vmix(j,0) + mu*umix(i,1)*vmix(j,1))/2./theMw;

	theQij[i][j] = umix(i,1)*vmix(j,0)/sqrt(2);
      }
      if( j < 2 ) {
	theQijLp[i][j] = theCb*( nmix(i, 3)*vmix(j,0) 
				 + (nmix(i,1) + nmix(i,0)*tw)*vmix(j,1)/sqrt(2));
	theQijRp[i][j] = theSb*( nmix(i, 3)*umix(j,0) 
				 - (nmix(i,1) + nmix(i,0)*tw)*umix(j,1)/sqrt(2));
      }

      theRijdp[i][j] = 
	( mTwo*nmix(i,1)*nmix(j,1) + mOne*nmix(i,0)*nmix(j,0)
	  - mu*(nmix(i,2)*nmix(j,3) + nmix(i,3)*nmix(j,2)) )/2./theMw;
      
      theQijdp[i][j] = 0.5*( nmix(i,2)*( nmix(j,1) - theSw*nmix(j,0) )
			     + nmix(j,2)*( nmix(i,1) - theSw*nmix(i,0) ) );
    }
  }
  
  orderInGem(1);
  orderInGs(0);
}

void SSGOGOHVertex::persistentOutput(PersistentOStream & os) const {
  os << theMSSM  << theRij << theQij << theQijLp << theQijRp << theRijdp
     << theQijdp << ounit(theMw,GeV) << theSw << theSa << theSb << theCa << theCb 
     << theC2b << theSba << theCba;
}

void SSGOGOHVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM  >> theRij >> theQij >> theQijLp >> theQijRp >> theRijdp
     >> theQijdp >> iunit(theMw,GeV) >> theSw >> theSa >> theSb >> theCa >> theCb 
     >> theC2b >> theSba >> theCba;
  theCoupLast = 0.0;
  theLLast = 0.0;
  theRLast = 0.0;
  theHLast = 0;
  theID1Last = 0;
  theID2Last = 0;
}

ClassDescription<SSGOGOHVertex> SSGOGOHVertex::initSSGOGOHVertex;
// Definition of the static class description member.

void SSGOGOHVertex::Init() {

  static ClassDocumentation<SSGOGOHVertex> documentation
    ("The coupling of the higgs bosons to SM fermions in the MSSM");

}

void SSGOGOHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			      tcPDPtr particle3, int) {
  long id1(particle1->id()), id2(particle2->id()), 
    id3(particle3->id()), higgsID(0), f1ID(0), f2ID(0);
  if( abs(id1) == ParticleID::h0 || abs(id1) == ParticleID::H0 || 
      abs(id1) == ParticleID::A0 || abs(id1) == ParticleID::Hplus ) {
    higgsID = abs(id1);
    f1ID = id2;
    f2ID = id3;
  }
  else if( abs(id2) == ParticleID::h0 || abs(id2) == ParticleID::H0 || 
	   abs(id2) == ParticleID::A0 || abs(id2) == ParticleID::Hplus  ) {
    higgsID = abs(id2);
    f1ID = id1;
    f2ID = id3;
  }
  else if( abs(id3) == ParticleID::h0 || abs(id3) == ParticleID::H0 || 
	   abs(id3) == ParticleID::A0 || abs(id3) == ParticleID::Hplus ) {
    higgsID = abs(id3);
    f1ID = id1;
    f2ID = id2;
  }
  else {
    throw HelicityConsistencyError() 
      << "SSGOGOHVertex::setCoupling - There is no higgs particle in "
      << "this vertex. Particles: " << id1 << " " << id2 << " " << id3
      << Exception::warning;
    return;
  }
  if( f1ID < 0 ) swap(f1ID, f2ID);

  if( higgsID == theHLast && f1ID == theID1Last && f2ID == theID2Last) {
    setNorm(theCoupLast);
    setLeft(theLLast);
    setRight(theRLast);
    return;
  }
  theHLast = higgsID;
  theID1Last = f1ID;
  theID2Last = f2ID;
  
  theCoupLast = -sqrt(4.*Constants::pi*theMSSM->alphaEM(q2))/theSw/theSb;

  if( higgsID == ParticleID::h0 ) {
    //charginos
    if( f2ID < 0 ) {
      unsigned int ei = (f1ID == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      unsigned int ej = (abs(f2ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      theLLast = -( conj(theQij[ei][ej])*theCba 
		    + conj(theRij[ei][ej])*theCa );
      theRLast = -( theQij[ej][ei]*theCba 
		    + theRij[ej][ei]*theCa );
      if( ei == ej ) {
	double delta = theCa*getParticleData(f1ID)->mass()/theMw/2.;
	theLLast += delta;
	theRLast += delta;
      }
    }
    //neutralinos
    else {
      unsigned int ei(f1ID - ParticleID::SUSY_chi_10), 
	ej(f2ID - ParticleID::SUSY_chi_10);
      if( ei > 1 )
	ei = ( ei == 13 ) ? 3 : 2;
      if( ej > 1 )
	ej = ( ej == 13 ) ? 3 : 2;
      theRLast = -( theQijdp[ei][ej]*theCba 
		    + theRijdp[ei][ej]*theCa ) ;
      theLLast = conj(theRLast);
      if( ei == ej ) {
	double delta = theCa*getParticleData(f1ID)->mass()/theMw/2.;
	theLLast += delta;
	theRLast += delta;
      }
    }
    
  }
  else if( higgsID == ParticleID::H0 ) {
    //charginos
    if( f2ID < 0 ) {
      unsigned int ei = (f1ID == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      unsigned int ej = (abs(f2ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      theLLast = conj(theQij[ei][ej])*theSba 
	- conj(theRij[ei][ej])*theSa;
      theRLast = theQij[ej][ei]*theSba 
	- theRij[ej][ei]*theSa ;
      if( ei == ej ) {
	double delta = theSa*getParticleData(f1ID)->mass()/theMw/2.;
	theLLast += delta;
	theRLast += delta;
      }
    }
    //neutralinos
    else {
      unsigned int ei(f1ID - ParticleID::SUSY_chi_10), 
	ej(f2ID - ParticleID::SUSY_chi_10);
      if( ei > 1 )
	ei = ( ei == 13 ) ? 3 : 2;
      if( ej > 1 )
	ej = ( ej == 13 ) ? 3 : 2;
      theRLast = theQijdp[ei][ej]*theSba 
	- theRijdp[ei][ej]*theSa;
      theLLast = conj(theRLast);
      if( ei == ej ) {
	double delta = theSa*getParticleData(f1ID)->mass()/theMw/2.;
	theLLast += delta;
	theRLast += delta;
      }
    }
  }
  else if( higgsID == ParticleID::A0 ) {
    theCoupLast *= Complex(0., 1.);
    if( f2ID < 0 ) {
      unsigned int ei = (f1ID == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      unsigned int ej = (abs(f2ID) == ParticleID::SUSY_chi_1plus) ? 0 : 1;
      theLLast = conj(theQij[ei][ej])*theC2b + conj(theRij[ei][ej])*theCb;
      theRLast = -theQij[ej][ei]*theC2b - theRij[ej][ei]*theCb;
      if( ei == ej ) {
	double delta = theCb*getParticleData(f1ID)->mass()/theMw/2.;
	theLLast += delta;
	theRLast -= delta;
      }
    }
    //neutralinos
    else {
      unsigned int ei(f1ID - ParticleID::SUSY_chi_10), 
	ej(f2ID - ParticleID::SUSY_chi_10);
      if( ei > 1 )
	ei = ( ei == 13 ) ? 3 : 2;
      if( ej > 1 )
	ej = ( ej == 13 ) ? 3 : 2;
      
      theRLast = -theQijdp[ei][ej]*theC2b - theRijdp[ei][ej]*theCb;
      theLLast = -conj(theRLast);
      if( ei == ej ) {
	double delta = theCb*getParticleData(f1ID)->mass()/theMw/2.;
	theLLast += delta;
	theRLast -= delta;
      }
    }
  }
  //charged higgs
  else {
    unsigned int ei(0),ej(0);
    long chg(f2ID), neu(f1ID);
    if( abs(neu) == ParticleID::SUSY_chi_1plus || 
	abs(neu) == ParticleID::SUSY_chi_2plus ) swap(chg, neu);
    ej = ( chg == ParticleID::SUSY_chi_1plus) ? 0 : 1;
    ei = neu - ParticleID::SUSY_chi_10;
    if( ei > 1 )
      ei = ( ei == 13 ) ? 3 : 2;
    
    theCoupLast *= theSb;
    theLLast = theQijLp[ei][ej];
    theRLast = theQijRp[ei][ej];
    if( chg < 0 ) {
      Complex tmp = theLLast;
      theLLast = conj(theRLast);
      theRLast = conj(tmp);
    }
  }
  setNorm(theCoupLast);
  setLeft(theLLast);
  setRight(theRLast);
}

