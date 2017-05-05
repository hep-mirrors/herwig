// -*- C++ -*-
//
// SSHPPVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSHPPVertex class.
//

#include "SSHPPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include <cassert>
#include "Herwig/Looptools/clooptools.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSHPPVertex::SSHPPVertex() : theIncludeTriLinear(true),
			     thePseudoScalarTreatment(false),
			     theSw(0.), theMw(), theZfact(),
			     theQt1L(0.), theQt1R(0.), theQt1LR(0.),
			     theQt2L(0.), theQt2R(0.), theQt2LR(0.),
			     theQb1L(0.), theQb1R(0.), theQb1LR(0.),
			     theQb2L(0.), theQb2R(0.), theQb2LR(0.),
			     theLt1L(0.), theLt1R(0.), theLt1LR(0.),
			     theLt2L(0.), theLt2R(0.), theLt2LR(0.),
			     theSfmass(6,ZERO),
			     theTanB(0.),theSinA(0.), 
			     theCosA(0.), theSinB(0.), theCosB(0.), 
			     theSinApB(0.), theCosApB(0.), 
			     theSinBmA(0.), theCosBmA(0.), 
			     theCouplast(0.), 
			     theq2last(), theHaveCoeff(false), theLastID(0) {
  orderInGs(0);
  orderInGem(3);
}

void SSHPPVertex::persistentOutput(PersistentOStream & os) const {
  os << theMSSM << theSw << ounit(theMw,GeV) << ounit(theZfact,GeV) 
     << theQt1L << theQt1R << theQt1LR << theQt2L << theQt2R << theQt2LR
     << theQb1L << theQb1R << theQb1LR << theQb2L << theQb2R << theQb2LR
     << theLt1L << theLt1R << theLt1LR << theLt2L << theLt2R << theLt2LR
     << thetop << thebot << thetau << theTanB << theIncludeTriLinear
     << theSinA << theCosA << theSinB << theCosB << theSinApB << theCosApB
     << theSinBmA << theCosBmA << thePseudoScalarTreatment
     << ounit(theSfmass, GeV) << theU << theV;
}

void SSHPPVertex::persistentInput(PersistentIStream & is, int) {
  is >> theMSSM >> theSw >> iunit(theMw,GeV) >> iunit(theZfact,GeV) 
     >> theQt1L >> theQt1R >> theQt1LR >> theQt2L >> theQt2R >> theQt2LR 
     >> theQb1L >> theQb1R >> theQb1LR >> theQb2L >> theQb2R >> theQb2LR
     >> theLt1L >> theLt1R >> theLt1LR >> theLt2L >> theLt2R >> theLt2LR
     >> thetop >> thebot >> thetau >> theTanB >> theIncludeTriLinear
     >> theSinA >> theCosA >> theSinB >> theCosB >> theSinApB >> theCosApB
     >> theSinBmA >> theCosBmA >> thePseudoScalarTreatment
     >> iunit(theSfmass, GeV) >> theU >> theV;
}

ClassDescription<SSHPPVertex> SSHPPVertex::initSSHPPVertex;
// Definition of the static class description member.

void SSHPPVertex::Init() {
  
  static ClassDocumentation<SSHPPVertex> documentation
    ("This class implements the higgs-gluon-gluon effective "
     "vertex in the MSSM including stop, sbottom and top quarks " 
     "loops.");

  static Switch<SSHPPVertex,bool> interfaceIncludeTriLinear
    ("IncludeTriLinear",
     "Whether or not to include the A term squark trilinear couplings",
     &SSHPPVertex::theIncludeTriLinear, true, false, false);
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

  static Switch<SSHPPVertex,bool> interfacePseudoScalarTreatment
    ("PseudoScalarTreatment",
     "Whether to treat the pseudoscalar as pseudoscalar or scalar, for testing only",
     &SSHPPVertex::thePseudoScalarTreatment, false, false, false);
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
  
void SSHPPVertex::setCoupling(Energy2 q2, tcPDPtr particle2,
			      tcPDPtr particle3, tcPDPtr particle1) {
  long higgs(abs(particle1->id()));
  // check allowed
  assert ( higgs == ParticleID::h0 || higgs == ParticleID::H0 ||
	   higgs == ParticleID::A0 );
  assert(particle2->id() == ParticleID::gamma && 
	 particle3->id() == ParticleID::gamma );
  // couplings
  if( q2 != theq2last || theCouplast == 0. || higgs != theLastID ) {
    Looptools::clearcache();
    theCouplast = sqr(electroMagneticCoupling(q2))*weakCoupling(q2);
    Energy mt   = theMSSM->mass(q2, thetop);    
    Energy mb   = theMSSM->mass(q2, thebot);
    Energy mtau = theMSSM->mass(q2, thetau);
    masses.clear();
    type.clear();
    if( higgs == ParticleID::h0 || higgs == ParticleID::H0 ) {
      setNParticles(13);
      masses.insert(masses.begin(), theSfmass.begin(), theSfmass.end());
      masses.push_back(mt);
      masses.push_back(mb);
      masses.push_back(mtau);
      masses.push_back(getParticleData(ParticleID::Hplus)->mass());
      masses.push_back(theMw);
      masses.push_back(getParticleData(ParticleID::SUSY_chi_1plus)->mass());
      masses.push_back(getParticleData(ParticleID::SUSY_chi_2plus)->mass());
      type.resize(13, PDT::Spin0);
      type[6] = PDT::Spin1Half;
      type[7] = PDT::Spin1Half;
      type[8] = PDT::Spin1Half;
      type[9] = PDT::Spin0;
      type[10] = PDT::Spin1;
      type[11] = PDT::Spin1Half;
      type[12] = PDT::Spin1Half;
      couplings.resize(13, make_pair(0., 0.));
      complex<Energy> brac1 = theZfact*(0.5 + theMSSM->ed()*sqr(theSw));
      complex<Energy> brac2 = theZfact*(0.5 - theMSSM->eu()*sqr(theSw));
      complex<Energy> brac3 = theZfact*theMSSM->ed()*sqr(theSw);
      complex<Energy> brac4 = theZfact*theMSSM->eu()*sqr(theSw);
      complex<Energy> brac5 = theZfact*(0.5 + theMSSM->ee()*sqr(theSw));
      complex<Energy> brac6 = theZfact*theMSSM->ee()*sqr(theSw);
      Energy Trib=theMSSM->bottomTrilinear().real();
      Energy Trit=theMSSM->topTrilinear().real();
      Energy Trita=theMSSM->tauTrilinear().real();
      Energy theMu = theMSSM->muParameter();
      if( higgs == ParticleID::h0 ) {
	// lightest sbottom
	complex<Energy> trilinear = theIncludeTriLinear ? 
	  theQb1LR*0.5*mb/theMw*(Trib*theSinA + theMu*theCosA)/theCosB : Energy();
	Complex coup = 3.*UnitRemoval::InvE*sqr(theMSSM->ed())*
	  (theQb1L *(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
	   theQb1R *(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3) +
	   trilinear);
	couplings[0] = make_pair(coup, coup);
	// lightest stop
	trilinear = theIncludeTriLinear ?
	  theQt1LR*0.5*mt/theMw*(Trit*theCosA + theMu*theSinA)/theSinB : Energy();
	coup = 3.*UnitRemoval::InvE*sqr(theMSSM->eu())*
	  (theQt1L *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
	   theQt1R *( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4) -
	   trilinear);
	couplings[1] = make_pair(coup, coup);
	// lightest stau
	trilinear = theIncludeTriLinear ?
	  theLt1LR*0.5*mtau/theMw*(Trita*theSinA + theMu*theCosA)/theCosB : Energy();
	coup = UnitRemoval::InvE*sqr(theMSSM->ee())*
	  (theLt1L *(   sqr(mtau)*theSinA/theMw/theCosB - theSinApB*brac5) +
	   theLt1R *(   sqr(mtau)*theSinA/theMw/theCosB + theSinApB*brac6) +
	   trilinear);
	couplings[2] = make_pair(coup, coup);
	// heavier sbottom
	trilinear = theIncludeTriLinear ? 
	   theQb2LR*0.5*mb/theMw*(Trib*theSinA + theMu*theCosA)/theCosB : Energy();
	coup = 3.*UnitRemoval::InvE*sqr(theMSSM->ed())*
	  (theQb2L *(   sqr(mb)*theSinA/theMw/theCosB - theSinApB*brac1) +
	   theQb2R *(   sqr(mb)*theSinA/theMw/theCosB + theSinApB*brac3) +
	   trilinear);
	couplings[3] = make_pair(coup, coup);
	// heavier stop
	trilinear = theIncludeTriLinear ? 
	  theQt2LR*0.5*mt/theMw*(Trit*theCosA + theMu*theSinA)/theSinB : Energy();
	coup = 3.*UnitRemoval::InvE*sqr(theMSSM->eu())*
	  (theQt2L*( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac2) +
	   theQt2R*( - sqr(mt)*theCosA/theMw/theSinB + theSinApB*brac4) -
	   trilinear);
	couplings[4] = make_pair(coup, coup);
	// heavier stau
	trilinear = theIncludeTriLinear ? 
	  theLt2LR*0.5*mtau/theMw*(Trita*theSinA + theMu*theCosA)/theCosB : Energy();
	coup = UnitRemoval::InvE*sqr(theMSSM->ee())*
	  (theLt2L *(   sqr(mtau)*theSinA/theMw/theCosB - theSinApB*brac5) +
	   theLt2R *(   sqr(mtau)*theSinA/theMw/theCosB + theSinApB*brac6)+
	   trilinear);
	couplings[5] = make_pair(coup, coup);
	// top
	coup = - 3.*mt*sqr(theMSSM->eu())*theCosA/2./theMw/theSinB;
	couplings[6] = make_pair(coup, coup);
	// bottom
	coup =   3.*mb*sqr(theMSSM->ed())*theSinA/2./theMw/theCosB;
	couplings[7] = make_pair(coup, coup);
	// tau
	coup =    mtau*sqr(theMSSM->ee())*theSinA/2./theMw/theCosB;
	couplings[8] = make_pair(coup, coup);
	// charged higgs
	coup = - UnitRemoval::InvE*theMw*(theSinBmA + 0.5/(1.-sqr(theSw))*
		       (sqr(theCosB)-sqr(theSinB))*theSinApB);
	couplings[9] = make_pair(coup, coup);
	// W boson
	coup = UnitRemoval::InvE*theMw*theSinBmA;
	couplings[10] = make_pair(coup, coup);
	// charginos
	for(unsigned int ix=0;ix<2;++ix) {
	  Complex Q = sqrt(0.5)*(*theV)(ix,0)*(*theU)(ix,1);
	  Complex S = sqrt(0.5)*(*theV)(ix,1)*(*theU)(ix,0);
	  coup = Q*theSinA-S*theCosA;
	  couplings[11+ix] = make_pair(conj(coup), coup);
	}
      }
      else {
	// lightest sbottom
	complex<Energy> trilinear = theIncludeTriLinear ? 
	  theQb1LR*0.5*mb/theMw*(theMu*theSinA - Trib*theCosA)/theCosB: Energy();
	Complex coup = 3.*UnitRemoval::InvE*sqr(theMSSM->ed())*
	  (theQb1L *( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
	   theQb1R *( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3)+trilinear);
	couplings[0] = make_pair(coup, coup);
	// lightest stop
	trilinear = theIncludeTriLinear ? 
	   -theQt1LR*0.5*mt/theMw*(-theMu*theCosA + Trit*theSinA)/theSinB: Energy();
	coup = 3.*UnitRemoval::InvE*sqr(theMSSM->eu())*
	  (theQt1L *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
	   theQt1R *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4)+trilinear);
	couplings[1] = make_pair(coup, coup);
	// lightest stau
	trilinear = theIncludeTriLinear ? 
	   theLt1LR*0.5*mtau/theMw*(theMu*theSinA - Trita*theCosA)/theCosB: Energy();
	coup = UnitRemoval::InvE*sqr(theMSSM->ee())*
	  (theLt1L *( - sqr(mtau)*theCosA/theMw/theCosB + theCosApB*brac5) +
	   theLt1R *( - sqr(mtau)*theCosA/theMw/theCosB - theCosApB*brac6)+trilinear);
	couplings[2] = make_pair(coup, coup);
	// heavier sbottom
	trilinear = theIncludeTriLinear ? 
	   theQb2LR*0.5*mb/theMw*(theMu*theSinA - Trib*theCosA)/theCosB: Energy();
	coup = 3.*UnitRemoval::InvE*sqr(theMSSM->ed())*
	  (theQb2L *( - sqr(mb)*theCosA/theMw/theCosB + theCosApB*brac1) +
	   theQb2R *( - sqr(mb)*theCosA/theMw/theCosB - theCosApB*brac3)+trilinear); 
	couplings[3] = make_pair(coup, coup);
	// heavier stop
	trilinear = theIncludeTriLinear ? 
	  -theQt2LR*0.5*mt/theMw*(-theMu*theCosA + Trit*theSinA)/theSinB: Energy();
	coup = 3.*UnitRemoval::InvE*sqr(theMSSM->eu())*
	  (theQt2L *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac2) +
	   theQt2R *( - sqr(mt)*theSinA/theMw/theSinB - theCosApB*brac4)+trilinear);
	couplings[4] = make_pair(coup, coup);
	// heavier stau
	trilinear = theIncludeTriLinear ? 
	   theLt2LR*0.5*mtau/theMw*(theMu*theSinA - Trita*theCosA)/theCosB: Energy();
	coup = UnitRemoval::InvE*sqr(theMSSM->ee())*
	  (theLt2L *( - sqr(mtau)*theCosA/theMw/theCosB + theCosApB*brac5) +
	   theLt2R *( - sqr(mtau)*theCosA/theMw/theCosB - theCosApB*brac6)+trilinear);
	couplings[5] = make_pair(coup, coup);
	// top
	coup = -3.*mt*sqr(theMSSM->eu())*theSinA/2./theMw/theSinB;
	couplings[6] = make_pair(coup, coup);
	// bottom
	coup = -3.*mb*sqr(theMSSM->ed())*theCosA/2./theMw/theCosB;
	couplings[7] = make_pair(coup, coup);
	// tau
	coup = -mtau*sqr(theMSSM->ee())*theCosA/2./theMw/theCosB;
	couplings[8] = make_pair(coup, coup);
	// charged higgs
	coup = - UnitRemoval::InvE*theMw*(theCosBmA - 0.5/(1.-sqr(theSw))*
					  (sqr(theCosB)-sqr(theSinB))*theCosApB);
	couplings[9] = make_pair(coup, coup);
	// W boson
	coup = UnitRemoval::InvE*theMw*theCosBmA;
	couplings[10] = make_pair(coup, coup);
	// charginos
	for(unsigned int ix=0;ix<2;++ix) {
	  Complex Q = sqrt(0.5)*(*theV)(ix,0)*(*theU)(ix,1);
	  Complex S = sqrt(0.5)*(*theV)(ix,1)*(*theU)(ix,0);
	  coup = -Q*theCosA-S*theSinA;
	  couplings[11+ix] = make_pair(conj(coup), coup);
	}
      }
    }
    else {
      setNParticles(5);
      masses.resize(5);
      couplings.resize(5);
      masses[0] = mt;
      masses[1] = mb;
      masses[2] = mtau;
      masses[3] = getParticleData(ParticleID::SUSY_chi_1plus)->mass();
      masses[4] = getParticleData(ParticleID::SUSY_chi_2plus)->mass();
      type.resize(5,PDT::Spin1Half);
      // top
      Complex coup = 3.*Complex(0., 1.)*sqr(theMSSM->eu())*mt/2./theMw/theTanB;
      couplings[0] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      // bottom
      coup = 3.*Complex(0., 1.)*sqr(theMSSM->ed())*mb/2./theMw*theTanB;
      couplings[1] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      // tau
      coup = Complex(0., 1.)*sqr(theMSSM->ee())*mtau/2./theMw*theTanB;
      couplings[2] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      // charginos
      for(unsigned int ix=0;ix<2;++ix) {
	Complex Q = sqrt(0.5)*(*theV)(ix,0)*(*theU)(ix,1);
	Complex S = sqrt(0.5)*(*theV)(ix,1)*(*theU)(ix,0);
	coup = - Complex(0., 1.)*(Q*theSinB+S*theCosB);
	couplings[3+ix] = make_pair(coup, thePseudoScalarTreatment ? coup : -coup);
      }
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

// functions for loops for testing
// namespace {

// Complex F0(double tau) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return tau*(1.-tau*ft);
// }

// Complex FHalf(double tau,double eta) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return -2.*tau*(eta+(1.-tau*eta)*ft);
// }

// Complex F1(double tau) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return 2.+3.*tau+3.*tau*(2.-tau)*ft;
// }
// }

void SSHPPVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(22,22,25);
  addToList(22,22,35);
  addToList(22,22,36);
  theMSSM = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !theMSSM ) 
    throw InitException()
      << "SSHPPVertex::doinit() - The pointer to the MSSM object is null!"
      << Exception::abortnow;
  theMw = getParticleData(ParticleID::Wplus)->mass();
  thetop = getParticleData(ParticleID::t);
  thebot = getParticleData(ParticleID::b);
  thetau = getParticleData(ParticleID::tauminus);
  theSw = sqrt(sin2ThetaW());
  theZfact = getParticleData(ParticleID::Wplus)->mass()/(1. - sqr(theSw));
  
  theSinA = sin(theMSSM->higgsMixingAngle());
  theCosA = sqrt(1. - sqr(theSinA));
  theTanB = theMSSM->tanBeta();
  theSinB = theTanB/sqrt(1. + sqr(theTanB));
  theCosB = sqrt( 1. - sqr(theSinB) );
  theSinApB = theSinA*theCosB + theCosA*theSinB;
  theCosApB = theCosA*theCosB - theSinA*theSinB;
  theSinBmA =-theSinA*theCosB + theCosA*theSinB;
  theCosBmA = theCosA*theCosB + theSinA*theSinB;
  
  MixingMatrix stop = *theMSSM->stopMix();
  MixingMatrix sbot = *theMSSM->sbottomMix();
  MixingMatrix stau = *theMSSM->stauMix();
  theQt1L  = stop(0,0)*stop(0,0);
  theQt1R  = stop(0,1)*stop(0,1);
  theQt1LR = stop(0,1)*stop(0,0) + stop(0,1)*stop(0,0);
  theQt2L  = stop(1,0)*stop(1,0);
  theQt2R  = stop(1,1)*stop(1,1);
  theQt2LR = stop(1,1)*stop(1,0) + stop(1,0)*stop(1,1);
  theQb1L  = sbot(0,0)*sbot(0,0);
  theQb1R  = sbot(0,1)*sbot(0,1);
  theQb1LR = sbot(0,1)*sbot(0,0) + sbot(0,1)*sbot(0,0);
  theQb2L  = sbot(1,0)*sbot(1,0);
  theQb2R  = sbot(1,1)*sbot(1,1);
  theQb2LR = sbot(1,1)*sbot(1,0) + sbot(1,0)*sbot(1,1);
  theLt1L  = stau(0,0)*stau(0,0);
  theLt1R  = stau(0,1)*stau(0,1);
  theLt1LR = stau(0,1)*stau(0,0) + stau(0,1)*stau(0,0);
  theLt2L  = stau(1,0)*stau(1,0);
  theLt2R  = stau(1,1)*stau(1,1);
  theLt2LR = stau(1,1)*stau(1,0) + stau(1,0)*stau(1,1);
  theU = theMSSM->charginoUMix();
  theV = theMSSM->charginoVMix();

  assert( theSfmass.size() == 6 );
  theSfmass[0] = getParticleData(ParticleID::SUSY_b_1)->mass();
  theSfmass[1] = getParticleData(ParticleID::SUSY_t_1)->mass();
  theSfmass[2] = getParticleData(ParticleID::SUSY_tau_1minus)->mass();
  theSfmass[3] = getParticleData(ParticleID::SUSY_b_2)->mass();
  theSfmass[4] = getParticleData(ParticleID::SUSY_t_2)->mass();
  theSfmass[5] = getParticleData(ParticleID::SUSY_tau_2minus)->mass();

  VVSLoopVertex::doinit();
  // test calc of the width
//   for(unsigned int ix=0;ix<2;++ix) {
//     Energy mh   = getParticleData(25+long(ix)*10)->mass();
//     Energy mt   = theMSSM->mass(sqr(mh  ), thetop);    
//     Energy mb   = theMSSM->mass(sqr(mh  ), thebot);    
//     Energy mtau = theMSSM->mass(sqr(mh  ), thetau);
//     Energy mhp  = getParticleData(ParticleID::Hplus)->mass();
//     Energy mc[2] = {getParticleData(ParticleID::SUSY_chi_1plus)->mass(),
// 		    getParticleData(ParticleID::SUSY_chi_2plus)->mass()};
//     // sbottom
//     Complex rsb1,rsb2;
//     if(ix==0) {
//       rsb1 = 
// 	+theQb1L*(-sqr(mb/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  -(-0.5+sqr(theSw)/3.)*theSinApB)
// 	+theQb1R*(-sqr(mb/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  +sqr(theSw)/3.*theSinApB);
//       rsb2 = 
// 	+theQb2L*(-sqr(mb/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  -(-0.5+sqr(theSw)/3.)*theSinApB)
// 	+theQb2R*(-sqr(mb/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  +sqr(theSw)/3.*theSinApB);
//     }
//     else {
//       rsb1 = 
// 	+theQb1L*(+sqr(mb/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  +(-0.5+sqr(theSw)/3.)*theCosApB)
// 	+theQb1R*(+sqr(mb/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  -sqr(theSw)/3.*theCosApB);
//       rsb2 = 
// 	+theQb2L*(+sqr(mb/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  +(-0.5+sqr(theSw)/3.)*theCosApB)
// 	+theQb2R*(+sqr(mb/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  -sqr(theSw)/3.*theCosApB);
//     }
//     Complex Isb1 = 3.*sqr(1./3.)*rsb1*sqr(theMw/theSfmass[0])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[0]/mh));
//     Complex Isb2 = 3.*sqr(1./3.)*rsb2*sqr(theMw/theSfmass[3])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[3]/mh));
//     // stop
//     Complex rst1,rst2;
//     if(ix==0) {
//       rst1 = 
// 	+theQt1L*(+sqr(mt/theMw)*(1.-sqr(theSw))*theCosA/theSinB
// 		  -(+0.5-2.*sqr(theSw)/3.)*theSinApB)
// 	+theQt1R*(+sqr(mt/theMw)*(1.-sqr(theSw))*theCosA/theSinB
// 		  -2.*sqr(theSw)/3.*theSinApB);
//       rst2 = 
// 	+theQt2L*(+sqr(mt/theMw)*(1.-sqr(theSw))*theCosA/theSinB
// 		  -(+0.5-2.*sqr(theSw)/3.)*theSinApB)
// 	+theQt2R*(+sqr(mt/theMw)*(1.-sqr(theSw))*theCosA/theSinB
// 		  -2.*sqr(theSw)/3.*theSinApB);
//     }
//     else {
//       rst1 = 
// 	+theQt1L*(+sqr(mt/theMw)*(1.-sqr(theSw))*theSinA/theSinB
// 		  +(+0.5-2.*sqr(theSw)/3.)*theCosApB)
// 	+theQt1R*(+sqr(mt/theMw)*(1.-sqr(theSw))*theSinA/theSinB
// 		  +2.*sqr(theSw)/3.*theCosApB);
//       rst2 = 
// 	+theQt2L*(+sqr(mt/theMw)*(1.-sqr(theSw))*theSinA/theSinB
// 		  +(+0.5-2.*sqr(theSw)/3.)*theCosApB)
// 	+theQt2R*(+sqr(mt/theMw)*(1.-sqr(theSw))*theSinA/theSinB
// 		  +2.*sqr(theSw)/3.*theCosApB);
//     }
//     Complex Ist1 = 3.*sqr(2./3.)*rst1*sqr(theMw/theSfmass[1])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[1]/mh));
//     Complex Ist2 = 3.*sqr(2./3.)*rst2*sqr(theMw/theSfmass[4])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[4]/mh));

//     // stau
//     Complex rstau1,rstau2;
//     if(ix==0) {
//       rstau1 = 
// 	+theLt1L*(-sqr(mtau/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  -(-0.5+sqr(theSw))*theSinApB)
// 	+theLt1R*(-sqr(mtau/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  +sqr(theSw)*theSinApB);
//       rstau2 = 
// 	+theLt2L*(-sqr(mtau/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  -(-0.5+sqr(theSw))*theSinApB)
// 	+theLt2R*(-sqr(mtau/theMw)*(1.-sqr(theSw))*theSinA/theCosB
// 		  +sqr(theSw)*theSinApB);
//     }
//     else {
//       rstau1 = 
// 	+theLt1L*(+sqr(mtau/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  +(-0.5+sqr(theSw))*theCosApB)
// 	+theLt1R*(+sqr(mtau/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  -sqr(theSw)*theCosApB);
//       rstau2 = 
// 	+theLt2L*(+sqr(mtau/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  +(-0.5+sqr(theSw))*theCosApB)
// 	+theLt2R*(+sqr(mtau/theMw)*(1.-sqr(theSw))*theCosA/theCosB
// 		  -sqr(theSw)*theCosApB);
//     }
//     Complex Istau1 = rstau1*sqr(theMw/theSfmass[2])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[2]/mh));
//     Complex Istau2 = rstau2*sqr(theMw/theSfmass[5])/(1.-sqr(theSw))
//       *F0(sqr(2.*theSfmass[5]/mh));
//     // charged higgs
//     Complex rh;
//     if(ix==0) {
//       rh = theSinBmA+0.5*(sqr(theCosB)-sqr(theSinB))*theSinApB/(1.-sqr(theSw));
//     }
//     else {
//       rh = theCosBmA-0.5*(sqr(theCosB)-sqr(theSinB))*theCosApB/(1.-sqr(theSw));
//     }
//     Complex Ih = rh*sqr(theMw/mhp)*F0(sqr(2.*mhp/mh));
//     // W
//     Complex rw;
//     if(ix==0) {
//       rw = theSinBmA;
//     }
//     else {
//       rw = theCosBmA;
//     }
//     Complex IW = rw*F1(sqr(2.*theMw/mh));
//     // top
//     Complex rt;
//     if(ix==0) {
//       rt = theCosA/theSinB;
//     }
//     else {
//       rt = theSinA/theSinB;
//     }
//     Complex Itop = 3.*sqr(2./3.)*rt*FHalf(sqr(2.*mt/mh),1.);
//     // bottom
//     Complex rb;
//     if(ix==0) {
//       rb =-theSinA/theCosB;
//     }
//     else {
//       rb = theCosA/theCosB;
//     }
//     Complex Ibot = 3.*sqr(1./3.)*rb*FHalf(sqr(2.*mb/mh),1.);
//     // tau
//     Complex rtau;
//     if(ix==0) {
//       rtau =-theSinA/theCosB;
//     }
//     else {
//       rtau = theCosA/theCosB;
//     }
//     Complex Itau = rtau*FHalf(sqr(2.*mtau/mh),1.);
//     // charginos
//     Complex rc[2],IC[2];
//     for(unsigned int ic=0;ic<2;++ic) {
//       Complex Q = sqrt(0.5)*(*theV)(ic,0)*(*theU)(ic,1);
//       Complex S = sqrt(0.5)*(*theV)(ic,1)*(*theU)(ic,0);
//       if(ix==0) {
// 	rc[ic] = 2.*(S*theCosA-Q*theSinA);
//       }
//       else {
// 	rc[ic] = 2.*(S*theSinA+Q*theCosA);
//       }
//       IC[ic] = rc[ic]*FHalf(sqr(2.*mc[ic]/mh),1.)*theMw/mc[ic];
//     }
//     Energy pre = sqr(mh/theMw)*mh/1024./pow(Constants::pi,3)
//       *sqr(weakCoupling(sqr(mh))*sqr(electroMagneticCoupling(sqr(mh)))/4./Constants::pi);
//     cerr << "testing lighter sbottom" << ix << " " 
// 	 << pre*std::norm(Isb1)/GeV << "\n";
//     cerr << "testing heavier sbottom" << ix << " " 
// 	 << pre*std::norm(Istau2)/GeV << "\n";
//     cerr << "testing lighter stop" << ix << " " 
// 	 << pre*std::norm(Ist1)/GeV << "\n";
//     cerr << "testing heavier stop" << ix << " " 
// 	 << pre*std::norm(Ist2)/GeV << "\n";
//     cerr << "testing lighter stau" << ix << " " 
// 	 << pre*std::norm(Istau1)/GeV << "\n";
//     cerr << "testing heavier stau" << ix << " " 
// 	 << pre*std::norm(Isb2)/GeV << "\n";
//     cerr << "testing top " << ix << " " 
// 	 << pre*std::norm(Itop)/GeV << "\n";
//     cerr << "testing bottom " << ix << " " 
// 	 << pre*std::norm(Ibot)/GeV << "\n";
//     cerr << "testing tau " << ix << " " 
// 	 << pre*std::norm(Itau)/GeV << "\n";
//     cerr << "testing higgs " << ix << " " 
// 	 << pre*std::norm(Ih)/GeV << "\n";
//     cerr << "testing W " << ix << " " 
// 	 << pre*std::norm(IW)/GeV << "\n";
//     cerr << "testing chi1 " << ix << " " 
// 	 << pre*std::norm(IC[0])/GeV << "\n";
//     cerr << "testing chi2 " << ix << " " 
// 	 << pre*std::norm(IC[1])/GeV << "\n";
//     cerr << "testing chi " << ix << " " 
// 	 << pre*std::norm(IC[0]+IC[1])/GeV << "\n";
//     cerr << "testing higgs width " << ix << " " 
// 	 << pre*std::norm(Isb1+Isb2+Ist1+Ist2+Istau1+Istau2+
// 			  Itop+Ibot+Itau+Ih+IW+IC[0]+IC[1])/GeV << "\n";
//   }
  if(loopToolsInitialized()) Looptools::ltexi();
}
 
