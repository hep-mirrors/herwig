// -*- C++ -*-
//
// SSHGGVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHGGVertex class.
//

#include "SSHGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHGGVertex::SSHGGVertex() : theSw(0.), theMw(), theZfact(),
			     theQt11(0.), theQt12(0.), theQt21(0.), theQt22(0.),
			     theQb11(0.), theQb12(0.), theQb21(0.), theQb22(0.),
			     theSqmass(4,ZERO),
			     theTanB(0.),theSinA(0.), 
			     theCosA(0.), theSinB(0.), theCosB(0.), 
			     theSinApB(0.), theCosApB(0.), theCouplast(0.), 
			     theq2last(), theHaveCoeff(false), theLastID(0) {
  //PDG codes for particles at vertices
  addToList(21,21,25);
  addToList(21,21,35);
  addToList(21,21,36);
}

void SSHGGVertex::doinit() {
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM ) 
    throw InitException()
      << "SSHGGVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  theMw = getParticleData(ParticleID::Wplus)->mass();
  thetop = getParticleData(ParticleID::t);
  thebot = getParticleData(ParticleID::b);
  theSw = sqrt(sin2ThetaW());
  theZfact = getParticleData(ParticleID::Z0)->mass()/sqrt(1. - sqr(theSw));
  
  theSinA = sin(theMSSM->higgsMixingAngle());
  theCosA = sqrt(1. - sqr(theSinA));
  theTanB = theMSSM->tanBeta();
  theSinB = theTanB/sqrt(1. + sqr(theTanB));
  theCosB = sqrt( 1. - sqr(theSinB) );
  theSinApB = theSinA*theCosB + theCosA*theSinB;
  theCosApB = theCosA*theCosB - theSinA*theSinB;
  
  MixingMatrix stop = *theMSSM->stopMix();
  MixingMatrix sbot = *theMSSM->sbottomMix();
  theQt11 = stop(0,0)*stop(0,0);
  theQt12 = stop(0,1)*stop(0,1);
  theQt21 = stop(1,0)*stop(1,0);
  theQt22 = stop(1,1)*stop(1,1);
  theQb11 = sbot(0,0)*sbot(0,0);
  theQb12 = sbot(0,1)*sbot(0,1);
  theQb21 = sbot(1,0)*sbot(1,0);
  theQb22 = sbot(1,1)*sbot(1,1);

  assert( theSqmass.size() == 4 );
  theSqmass[0] = getParticleData(ParticleID::SUSY_b_1)->mass();
  theSqmass[1] = getParticleData(ParticleID::SUSY_t_1)->mass();
  theSqmass[2] = getParticleData(ParticleID::SUSY_b_2)->mass();
  theSqmass[3] = getParticleData(ParticleID::SUSY_t_2)->mass();

  orderInGs(2);
  orderInGem(1);
  VVSLoopVertex::doinit();
}


void SSHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << theMSSM << theSw << ounit(theMw,GeV) << ounit(theZfact,GeV) 
     << theQt11 << theQt12 << theQt21 << theQt22
     << theQb11 << theQb12 << theQb21 << theQb22
     << thetop << thebot << theTanB
     << theSinA << theCosA << theSinB << theCosB << theSinApB << theCosApB
     << ounit(theSqmass, GeV);
}

void SSHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM >> theSw >> iunit(theMw,GeV) >> iunit(theZfact,GeV) 
     >> theQt11 >> theQt12 >> theQt21 >> theQt22 
     >> theQb11 >> theQb12 >> theQb21 >> theQb22 
     >> thetop >> thebot >> theTanB
     >> theSinA >> theCosA >> theSinB >> theCosB >> theSinApB >> theCosApB
     >> iunit(theSqmass, GeV);
}

ClassDescription<SSHGGVertex> SSHGGVertex::initSSHGGVertex;
// Definition of the static class description member.

void SSHGGVertex::Init() {
  
  static ClassDocumentation<SSHGGVertex> documentation
    ("This class implements the higgs-gluon-gluon effective "
     "vertex in the MSSM including stop, sbottom and top quarks " 
     "loops.");

}
  
void SSHGGVertex::setCoupling(Energy2 q2, tcPDPtr particle2,
			      tcPDPtr particle3, tcPDPtr particle1) {
  long higgs(abs(particle1->id()));
  if( higgs != ParticleID::h0 && 
      higgs != ParticleID::H0 &&
      higgs != ParticleID::A0 &&
      particle2->id() != ParticleID::g && particle3->id() != ParticleID::g )
    throw HelicityConsistencyError() 
      << "SSHGGVertex::setCoupling(): This vertex has the incorrect "
      << "particle content in it. " << higgs << " " 
      << particle2->id() << " " << particle3->id();
  if( q2 != theq2last || theCouplast == 0. || higgs != theLastID ) {
    theCouplast = sqr(strongCoupling(q2))*weakCoupling(q2);
    Energy mt = theMSSM->mass(q2, thetop);    
    Energy mb = theMSSM->mass(q2, thebot);
    masses.resize(0);
    type.resize(0);
    if( higgs == ParticleID::h0 || higgs == ParticleID::H0 ) {
      setNParticles(6);
      masses.insert(masses.begin(), theSqmass.begin(), theSqmass.end());
      masses.push_back(mt);
      masses.push_back(mb);
      type.resize(6, PDT::Spin0);
      type[4] = PDT::Spin1Half;
      type[5] = PDT::Spin1Half;
      couplings.resize(6, make_pair(0., 0.));  
      complex<Energy> brac1 = theZfact*(0.5 + theMSSM->ed()*sqr(theSw));
      complex<Energy> brac2 = theZfact*(0.5 - theMSSM->eu()*sqr(theSw));
      complex<Energy> brac3 = theZfact*theMSSM->ed()*sqr(theSw);
      complex<Energy> brac4 = theZfact*theMSSM->eu()*sqr(theSw);
      if( higgs == ParticleID::h0 ) {
	// lightest sbottom
	Complex coup = 0.5*UnitRemoval::InvE*
	  (theQb11*(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
	   theQb12*(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3));
	couplings[0] = make_pair(coup, coup);
	// lightest stop
	coup = 0.5*UnitRemoval::InvE*
	  (theQt11*( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
	   theQt12*( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4));
	couplings[1] = make_pair(coup, coup);
	// heavier sbottom
	coup = 0.5*UnitRemoval::InvE*
	  (theQb21*(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
	   theQb22*(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3));
	couplings[2] = make_pair(coup, coup);
	// heavier stop
	coup = 0.5*UnitRemoval::InvE*
	  (theQt21*( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
	   theQt22*( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4));
	couplings[3] = make_pair(coup, coup);
	// top
	coup = -0.5*(mt*theCosA/2./theMw/theSinB);
	couplings[4] = make_pair(coup, coup);
	// bottom
	coup = +0.5*(mb*theSinA/2./theMw/theCosB);
	couplings[5] = make_pair(coup, coup);
      }
      else {
	// lightest sbottom
	Complex coup = 0.5*UnitRemoval::InvE*
	  (theQb11*( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
	   theQb12*( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3));
	couplings[0] = make_pair(coup, coup);
	// lightest stop
	coup = 0.5*UnitRemoval::InvE*
	  (theQt11*( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
	   theQt12*( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4));
	couplings[1] = make_pair(coup, coup);
	// heavier sbottom
	coup = 0.5*UnitRemoval::InvE*
	  (theQb21*( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
	   theQb22*( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3)); 
	couplings[2] = make_pair(coup, coup);
	// heavier stop
	coup = 0.5*UnitRemoval::InvE*
	  (theQt21*( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
	   theQt22*( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4));
	couplings[3] = make_pair(coup, coup);
	// top
	coup = -0.5*mt*theSinA/2./theMw/theSinB;
	couplings[4] = make_pair(coup, coup);
	// bottom
	coup = -0.5*mb*theCosA/2./theMw/theCosB;
	couplings[5] = make_pair(coup, coup);
      }
    }
    else {
      setNParticles(2);
      masses.resize(2);
      couplings.resize(2);
      masses[0] = mt;
      masses[1] = mb;
      type.resize(2,PDT::Spin1Half);
      Complex coup = 0.5*Complex(0., 1.)*mt/2./theMw/theTanB;
      couplings[0] = make_pair(coup, -coup);
      coup = 0.5*Complex(0., 1.)*mb/2./theMw*theTanB;
      couplings[1] = make_pair(coup, -coup);
    }
    theq2last = q2;
    theLastID = higgs;
    theHaveCoeff = false;
  }
  norm(theCouplast);
  //calculate tensor coefficients
  if( !theHaveCoeff ) {
    VVSLoopVertex::setCoupling(q2, particle2, particle3, particle1);
    theHaveCoeff = true;
  }
}

 
