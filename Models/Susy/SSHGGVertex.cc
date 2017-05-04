// -*- C++ -*-
//
// SSHGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHGGVertex class.
//

#include "SSHGGVertex.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Looptools/clooptools.h"
#include <cassert>

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHGGVertex::SSHGGVertex() : theIncludeTriLinear(true),
			     thePseudoScalarTreatment(false),
			     theSw(0.), theMw(), theZfact(),
			     theQt1L(0.), theQt1R(0.), theQt1LR(0.), 
			     theQt2L(0.), theQt2R(0.), theQt2LR(0.),
			     theQb1L(0.), theQb1R(0.), theQb1LR(0.), 
			     theQb2L(0.), theQb2R(0.), theQb2LR(0.),
			     theSqmass(4,ZERO),
			     theTanB(0.),theSinA(0.), 
			     theCosA(0.), theSinB(0.), theCosB(0.), 
			     theSinApB(0.), theCosApB(0.), theCouplast(0.), 
			     theq2last(), theHaveCoeff(false), theLastID(0) {
  orderInGs(2);
  orderInGem(1);
}

void SSHGGVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(21,21,25);
  addToList(21,21,35);
  addToList(21,21,36);
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM ) 
    throw InitException()
      << "SSHGGVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  theMw = getParticleData(ParticleID::Wplus)->mass();
  thetop = getParticleData(ParticleID::t);
  thebot = getParticleData(ParticleID::b);
  theSw = sqrt(sin2ThetaW());
  theZfact = getParticleData(ParticleID::Wplus)->mass()/(1. - sqr(theSw));
  
  theSinA = sin(theMSSM->higgsMixingAngle());
  theCosA = sqrt(1. - sqr(theSinA));
  theTanB = theMSSM->tanBeta();
  theSinB = theTanB/sqrt(1. + sqr(theTanB));
  theCosB = sqrt( 1. - sqr(theSinB) );
  theSinApB = theSinA*theCosB + theCosA*theSinB;
  theCosApB = theCosA*theCosB - theSinA*theSinB;
  
  MixingMatrixPtr stop = theMSSM->stopMix();
  MixingMatrixPtr sbot = theMSSM->sbottomMix();
  theQt1L  = (*stop)(0,0)*(*stop)(0,0);
  theQt1R  = (*stop)(0,1)*(*stop)(0,1);
  theQt1LR = (*stop)(0,1)*(*stop)(0,0) + (*stop)(0,1)*(*stop)(0,0);
  theQt2L  = (*stop)(1,0)*(*stop)(1,0);
  theQt2R  = (*stop)(1,1)*(*stop)(1,1);
  theQt2LR = (*stop)(1,1)*(*stop)(1,0) + (*stop)(1,0)*(*stop)(1,1);
  theQb1L  = (*sbot)(0,0)*(*sbot)(0,0);
  theQb1R  = (*sbot)(0,1)*(*sbot)(0,1);
  theQb1LR = (*sbot)(0,1)*(*sbot)(0,0) + (*sbot)(0,1)*(*sbot)(0,0);
  theQb2L  = (*sbot)(1,0)*(*sbot)(1,0);
  theQb2R  = (*sbot)(1,1)*(*sbot)(1,1);
  theQb2LR = (*sbot)(1,1)*(*sbot)(1,0) + (*sbot)(1,0)*(*sbot)(1,1);
  
  assert( theSqmass.size() == 4 );
  theSqmass[0] = getParticleData(ParticleID::SUSY_b_1)->mass();
  theSqmass[1] = getParticleData(ParticleID::SUSY_t_1)->mass();
  theSqmass[2] = getParticleData(ParticleID::SUSY_b_2)->mass();
  theSqmass[3] = getParticleData(ParticleID::SUSY_t_2)->mass();

  VVSLoopVertex::doinit();
  if(loopToolsInitialized()) Looptools::ltexi();
}


void SSHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << theMSSM << theSw << ounit(theMw,GeV) << ounit(theZfact,GeV) 
     << theQt1L << theQt1R << theQt1LR << theQt2L << theQt2R << theQt2LR
     << theQb1L << theQb1R << theQb1LR << theQb2L << theQb2R << theQb2LR
     << thetop << thebot << theTanB << theIncludeTriLinear
     << theSinA << theCosA << theSinB << theCosB << theSinApB << theCosApB
     << ounit(theSqmass, GeV) << thePseudoScalarTreatment;
}

void SSHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM >> theSw >> iunit(theMw,GeV) >> iunit(theZfact,GeV) 
     >> theQt1L >> theQt1R >> theQt1LR >> theQt2L >> theQt2R >> theQt2LR 
     >> theQb1L >> theQb1R >> theQb1LR >> theQb2L >> theQb2R >> theQb2LR 
     >> thetop >> thebot >> theTanB >> theIncludeTriLinear
     >> theSinA >> theCosA >> theSinB >> theCosB >> theSinApB >> theCosApB
     >> iunit(theSqmass, GeV) >> thePseudoScalarTreatment;
}

ClassDescription<SSHGGVertex> SSHGGVertex::initSSHGGVertex;
// Definition of the static class description member.

void SSHGGVertex::Init() {
  
  static ClassDocumentation<SSHGGVertex> documentation
    ("This class implements the higgs-gluon-gluon effective "
     "vertex in the MSSM including stop, sbottom and top quarks " 
     "loops.");

  static Switch<SSHGGVertex,bool> interfaceIncludeTriLinear
    ("IncludeTriLinear",
     "Whether or not to include the A term squark trilinear couplings",
     &SSHGGVertex::theIncludeTriLinear, true, false, false);
  static SwitchOption interfaceIncludeTriLinearYes
    (interfaceIncludeTriLinear,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeTriLinearNo
    (interfaceIncludeTriLinear,
     "No",
     "Don't include them",
     false);

  static Switch<SSHGGVertex,bool> interfacePseudoScalarTreatment
    ("PseudoScalarTreatment",
     "Whether to treat the pseudoscalar as pseudoscalar or scalar, for testing only",
     &SSHGGVertex::thePseudoScalarTreatment, false, false, false);
  static SwitchOption interfacePseudoScalarTreatmentPseudoScalar
    (interfacePseudoScalarTreatment,
     "PseudoScalar",
     "Treat as a pseudoscalar, the right physics",
     false);
  static SwitchOption interfacePseudoScalarTreatmentScalar
    (interfacePseudoScalarTreatment,
     "Scalar",
     "Treat as a scalar for testing",
     true);
}
  
void SSHGGVertex::setCoupling(Energy2 q2, tcPDPtr particle2,
			      tcPDPtr particle3, tcPDPtr particle1) {
  long higgs(abs(particle1->id()));
  assert( higgs == ParticleID::h0 || higgs == ParticleID::H0 ||
	  higgs == ParticleID::A0 );
  assert(particle2->id() == ParticleID::g && particle3->id() == ParticleID::g );
  if( q2 != theq2last || theCouplast == 0. || higgs != theLastID ) {
    Looptools::clearcache();
    theCouplast = weakCoupling(q2)*sqr(strongCoupling(q2));
    Energy mt = theMSSM->mass(q2, thetop);    
    Energy mb = theMSSM->mass(q2, thebot);
    masses.clear();
    type.clear();
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
      Energy Trib=theMSSM->bottomTrilinear().real();
      Energy Trit=theMSSM->topTrilinear().real();
      Energy theMu = theMSSM->muParameter();
      if( higgs == ParticleID::h0 ) {
	// lightest sbottom
	complex<Energy> trilinear = theIncludeTriLinear ? 
	  theQb1LR*0.5*mb/theMw*(Trib*theSinA + theMu*theCosA)/theCosB : Energy();
	Complex coup = 0.5*UnitRemoval::InvE*
	  (theQb1L *(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
	   theQb1R *(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3)+trilinear);

	couplings[0] = make_pair(coup, coup);
	// lightest stop
	trilinear = theIncludeTriLinear ? 
	 -theQt1LR*0.5*mt/theMw*(Trit*theCosA + theMu*theSinA)/theSinB : Energy();
	coup = 0.5*UnitRemoval::InvE*
	  (theQt1L *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
	   theQt1R *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4)+trilinear);

	couplings[1] = make_pair(coup, coup);
	// heavier sbottom
	trilinear = theIncludeTriLinear ? 
	  theQb2LR*0.5*mb/theMw*(Trib*theSinA + theMu*theCosA)/theCosB : Energy();
	coup = 0.5*UnitRemoval::InvE*
	  (theQb2L *(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
	   theQb2R *(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3)+trilinear);

	couplings[2] = make_pair(coup, coup);
	// heavier stop
	trilinear = theIncludeTriLinear ? 
	  -theQt2LR*0.5*mt/theMw*(Trit*theCosA + theMu*theSinA)/theSinB : Energy();
	coup = 0.5*UnitRemoval::InvE*
	  (theQt2L *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
	   theQt2R *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4)+trilinear);
	   		
	couplings[3] = make_pair(coup, coup);
	// top
	coup = -0.25*(mt*theCosA/theMw/theSinB);

	couplings[4] = make_pair(coup, coup);
	// bottom
	coup = +0.25*(mb*theSinA/theMw/theCosB);

	couplings[5] = make_pair(coup, coup);
      }
      else {
	// lightest sbottom
	complex<Energy> trilinear = theIncludeTriLinear ? 
	   theQb1LR*0.5*mb/theMw*(theMu*theSinA - Trib*theCosA)/theCosB: Energy();
	Complex coup = 0.5*UnitRemoval::InvE*
	  (theQb1L *( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
	   theQb1R *( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3)+trilinear);

	couplings[0] = make_pair(coup, coup);
	// lightest stop
	trilinear = theIncludeTriLinear ? 
	   -theQt1LR*0.5*mt/theMw*(-theMu*theCosA + Trit*theSinA)/theSinB: Energy();
	coup = 0.5*UnitRemoval::InvE*
	  (theQt1L *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
	   theQt1R *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4)+trilinear);

	couplings[1] = make_pair(coup, coup);
	// heavier sbottom
	trilinear = theIncludeTriLinear ? 
	   theQb2LR*0.5*mb/theMw*(theMu*theSinA - Trib*theCosA)/theCosB: Energy();
	coup = 0.5*UnitRemoval::InvE*
	  (theQb2L *( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
	   theQb2R *( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3)+trilinear); 
	
	couplings[2] = make_pair(coup, coup);
	// heavier stop
	trilinear = theIncludeTriLinear ? 
	   -theQt2LR*0.5*mt/theMw*(-theMu*theCosA + Trit*theSinA)/theSinB: Energy();
	coup = 0.5*UnitRemoval::InvE*
	  (theQt2L *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
	   theQt2R *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4)+trilinear);
	
	couplings[3] = make_pair(coup, coup);
	// top
	coup = -0.25*mt*theSinA/theMw/theSinB;
	couplings[4] = make_pair(coup, coup);
	// bottom
	coup = -0.25*mb*theCosA/theMw/theCosB;
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
      Complex coup = 0.25*Complex(0., 1.)*mt/theMw/theTanB;
	  	
      couplings[0] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      coup = 0.25*Complex(0., 1.)*mb/theMw*theTanB;
	 
      couplings[1] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
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

 
