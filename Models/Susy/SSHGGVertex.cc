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
			     theQt11(0.), theQt22(0.), theQb11(0.),
			     theQb22(0.), theSqmass(4,ZERO),
			     theTanB(0.),theSinA(0.), 
			     theCosA(0.), theSinB(0.), theCosB(0.), 
			     theSinApB(0.), theCosApB(0.), theCouplast(0.), 
			     theq2last(), theHaveCoeff(false) {
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
  theQt22 = stop(1,1)*stop(1,1);
  theQb11 = sbot(0,0)*sbot(0,0);
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
  os << theMSSM << theSw << ounit(theMw,GeV) << ounit(theZfact,GeV) << theQt11 
     << theQt22 << theQb11 << theQb22 << thetop << thebot << theTanB
     << theSinA << theCosA << theSinB << theCosB << theSinApB << theCosApB
     << ounit(theSqmass, GeV);
}

void SSHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM >> theSw >> iunit(theMw,GeV) >> iunit(theZfact,GeV) >> theQt11 
     >> theQt22 >> theQb11 >> theQb22 >> thetop >> thebot >> theTanB
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
  if( q2 != theq2last || theCouplast == 0.)	{
    theCouplast = sqr(strongCoupling(q2))*weakCoupling(q2);
    Energy mt = theMSSM->mass(q2, thetop);
    if( higgs == ParticleID::h0 || higgs == ParticleID::H0 ) {
      setNParticles(5);
      masses.insert(masses.begin(), theSqmass.begin(), theSqmass.end());
      masses.push_back(mt);
      type.resize(5, PDT::Spin0);
      type[4] = PDT::Spin1Half;
      couplings.resize(5, make_pair(0., 0.));      
      Energy mb = theMSSM->mass(q2, thebot);
      complex<Energy> brac1 = theZfact*theQb11*(0.5 + theMSSM->ed()*theSw*theSw);
      complex<Energy> brac2 = theZfact*theQt11*(0.5 - theMSSM->eu()*theSw*theSw);
      complex<Energy> brac3 = theZfact*theQb22*theMSSM->ed()*theSw*theSw;
      complex<Energy> brac4 = theZfact*theQt22*theMSSM->eu()*theSw*theSw;
      if( higgs == ParticleID::h0 ) {
	Complex coup = 
	  (theQb11*mb*mb*theSinA/theMw/theCosB - theSinApB*brac1)*UnitRemoval::InvE;
	couplings[0] = make_pair(coup, coup);
	coup = (theSinApB*brac2 
		- theQt11*mt*mt*theCosA/theMw/theSinB)*UnitRemoval::InvE;
	couplings[1] = make_pair(coup, coup);
	coup = (theQb22*mb*mb*theSinA/theMw/theCosB 
		+ theSinApB*brac3)*UnitRemoval::InvE;
	couplings[2] = make_pair(coup, coup);
	coup = (theSinApB*brac4 
		- theQt22*mt*mt*theCosA/theMw/theSinB)*UnitRemoval::InvE;
	couplings[3] = make_pair(coup, coup);
	coup = -(mt*theCosA/2./theMw/theSinB);
	couplings[4] = make_pair(coup, coup);
      }
      else {
	Complex coup = (theCosApB*brac1 
			- theQb11*mb*mb*theCosA/theMw/theSinB)*UnitRemoval::InvE;
	couplings[0] = make_pair(coup, coup);
	coup = -(theCosApB*brac2 
		 + theQt11*mt*mt*theSinA/theMw/theSinB)*UnitRemoval::InvE;
	couplings[1] = make_pair(coup, coup);
	coup = -(theCosApB*brac3 
		 + theQb22*mb*mb*theCosA/theMw/theCosB)*UnitRemoval::InvE; 
	couplings[2] = make_pair(coup, coup);
	coup = -(theSinApB*brac4 
		 + theQt22*mt*mt*theSinA/theMw/theSinB)*UnitRemoval::InvE;
	couplings[3] = make_pair(coup, coup);
	coup = -mt*theSinA/2./theMw/theSinB;
	couplings[4] = make_pair(coup, coup);
      }
    }
    else {
      setNParticles(1);
      masses.resize(1, mt);
      type.resize(1, PDT::Spin1Half);
      Complex coup = Complex(0., 1.)*mt/2./theMw/theTanB;
      couplings.resize(1,  make_pair(coup, -coup));
    }
    theq2last = q2;
    theHaveCoeff = false;
  }
  norm(theCouplast);
  //calculate tensor coefficients
  if( !theHaveCoeff ) {
    VVSLoopVertex::setCoupling(q2, particle2, particle3, particle1);
    theHaveCoeff = true;
  }
}

 
