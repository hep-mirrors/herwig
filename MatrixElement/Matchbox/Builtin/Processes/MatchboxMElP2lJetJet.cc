// -*- C++ -*-
//
// MatchboxMElP2lJetJet.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMElP2lJetJet class.

      
#include "MatchboxMElP2lJetJet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxMElP2lJetJet::MatchboxMElP2lJetJet() 
  : MatchboxMEBase(), 
    theZMass2(0.0*GeV2), theZWidth2(0.0*GeV2),
    theUserScale(0.0*GeV) {}

MatchboxMElP2lJetJet::~MatchboxMElP2lJetJet() {}

bool MatchboxMElP2lJetJet::generateKinematics(const double * r) {

  if ( phasespace() )
    return MatchboxMEBase::generateKinematics(r);

  if ( theZMass2 == 0.0*GeV2 ) {
    theZMass2 = sqr(getParticleData(ParticleID::Z0)->mass());
    theZWidth2 = sqr(getParticleData(ParticleID::Z0)->width());
  }
  
  double x;
  Lorentz5Momentum nloMomenta[5];
  Energy LeptonMass;
  if(fabs(mePartonData()[1]->id()) != 11 ||
     fabs(mePartonData()[1]->id()) != 13 ||
     fabs(mePartonData()[1]->id()) != 15){
    x = lastX2();
    nloMomenta[0] = meMomenta()[0];
    nloMomenta[1] = meMomenta()[1];
    LeptonMass = mePartonData()[0]->mass();  
      } else {
    x = lastX1();
    nloMomenta[1] = meMomenta()[0];
    nloMomenta[0] = meMomenta()[1];
    LeptonMass = mePartonData()[1]->mass();
  }

  // Defining variables
  double QZmin = 2.;
  double QZmax = sqrt(lastSHat())/ThePEG::GeV;

  Energy2 Q2=exp((log(sqr(QZmax))-log(sqr(QZmin)))*r[0]+log(sqr(QZmin)))*ThePEG::GeV2;

  double x0 = .01;

  double xB = x;
  double xi1Min = -log(1.+x0-xB);        
  double xi1Max = -log(x0);            
   
  double xi2Min = 1./(1.+2.*x0)*log(x0/(1.+x0));                            
  double xi2Max = 1./(1.+2.*x0)*log((x0+1.)/x0);  
                                                             
  double xi1 = xi1Min + r[1]*(xi1Max-xi1Min);  
  double xi2 = xi2Min + r[2]*(xi2Max-xi2Min); 
      
  double xp = (1. + x0 - exp(-xi1));                                        
  double zp = ((1.+x0)*exp(xi2*(1.+2.*x0))-x0)/(1.+exp(xi2*(1.+2.*x0))); 

  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2; 
  
  double PhiParton = 2.*Constants::pi*r[3];
  double PhiLepton = 2.*Constants::pi*r[4];
 
  // Partonic bit
  Energy pTParton = sqrt(Q2*(1.-xp)*(1.-zp)*zp/xp);
  Energy Parton2Energy = sqrt(Q2*sqr(x2)/4.+sqr(pTParton));
  Energy Parton4Energy = sqrt(Q2*sqr(x3)/4.+sqr(pTParton));
  
  // Parton momenta in the Breit frame (BF)
  nloMomenta[3] =
    Lorentz5Momentum(pTParton*std::cos(PhiParton),
		     pTParton*std::sin(PhiParton),
		     -0.5*sqrt(Q2)*x2,Parton2Energy);
  nloMomenta[4] =
    Lorentz5Momentum(-pTParton*std::cos(PhiParton),
		     -pTParton*std::sin(PhiParton),
		     -0.5*sqrt(Q2)*x3,Parton4Energy);

  Lorentz5Momentum px = (nloMomenta[3]+nloMomenta[4]);

  if(px.m2()> lastSHat()+sqr(LeptonMass)) {
    jacobian(0.0);
    return false; 
  }
       
  // Leptonic bit
  Energy2 OutLeptonMagnitude2 =
    sqr(lastSHat()+sqr(LeptonMass)-px.m2())/
    4./lastSHat()-sqr(LeptonMass);

  if ( OutLeptonMagnitude2 < ZERO ) {
    jacobian(0.0);
    return false;
  }

  Energy OutLeptonMagnitude = 
    sqrt(OutLeptonMagnitude2);
  Energy OutLeptonEnergy = 
    sqrt(sqr(OutLeptonMagnitude)+sqr(LeptonMass));
 
  if(OutLeptonEnergy >= sqrt(lastSHat())) {
    jacobian(0.0);
    return false;
  }

  Energy InLeptonMagnitude =
    nloMomenta[0].vect().mag();

  double cThetaLepton = 
    (nloMomenta[0].e()*OutLeptonEnergy-0.5*Q2-sqr(LeptonMass))/
    InLeptonMagnitude/OutLeptonMagnitude;

  if(cThetaLepton < -1. || cThetaLepton > 1.) {
    jacobian(0.0);
    return false;
  }

  double sThetaLepton = sqrt(1.-sqr(cThetaLepton));

  // Lepton momentum in the CM frame
  nloMomenta[2] =   
  Lorentz5Momentum(OutLeptonMagnitude*sThetaLepton*std::cos(PhiLepton),
		   OutLeptonMagnitude*sThetaLepton*std::sin(PhiLepton),
		   OutLeptonMagnitude*cThetaLepton,OutLeptonEnergy);

  double betaZ = (px.z()*px.e()+OutLeptonMagnitude*
		  (sqrt(lastSHat())-OutLeptonEnergy))/
                  (sqr(sqrt(lastSHat())-OutLeptonEnergy)
                   +sqr(px.z()));
  
  Boost betaBFtoCM(0,0,-betaZ); 
  if(betaZ>=1.) {
    jacobian(0.0);
    return false;
  }

  LorentzRotation rotX;
  double Theta = acos(cThetaLepton);
  rotX = LorentzRotation();
  rotX.setRotate(-Theta,Axis(1.,0.,0.));

  LorentzRotation rotZ;
  rotZ = LorentzRotation();
  rotZ.setRotate(PhiLepton-0.5*Constants::pi,
		Axis(0.,0.,1.));

  // Boosting partonic momenta
  // from BF to CM
  nloMomenta[3].boost(betaBFtoCM);
  nloMomenta[3] *= rotX;
  nloMomenta[3] *= rotZ;
 
  nloMomenta[4].boost(betaBFtoCM);
  nloMomenta[4] *= rotX;
  nloMomenta[4] *= rotZ;
  
  // Definig Weight
  double psWeight = 1./(64.*pow(2.*Constants::pi,3));
  double jac = Q2/(InLeptonMagnitude*sqrt(lastSHat()))/sqr(xp);
  double  InvQ2Weight = Q2*(log(sqr(QZmax))-log(sqr(QZmin)))/lastSHat();
  
  // Mapping to flatten the Jacobian in the singular domain
  double mapX = (xi1Max-xi1Min)*(xi2Max-xi2Min)*(zp+x0)*(1.+x0-zp)*(1.+x0-xp);

  for(int i = 2; i < 5; i++)
    meMomenta()[i] = nloMomenta[i];
  
  jacobian(mapX*psWeight*InvQ2Weight*jac);
  
  setScale();
  logGenerateKinematics(r);
  return true;

}

Energy2 MatchboxMElP2lJetJet::factorizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return -(meMomenta()[0]-meMomenta()[2]).m2();

}

Energy2 MatchboxMElP2lJetJet::renormalizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return -(meMomenta()[0]-meMomenta()[2]).m2();

}

double MatchboxMElP2lJetJet::colourCorrelatedME2(pair<int,int> ij) const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::colourCorrelatedME2(ij);

  generator()->logWarning(Exception() 
			  << "The matrix element '" << name() << "' "
			  << "is not capable of calculating colour- or spin correlated "
			  << "matrix element squares."
			  << Exception::warning);
  
  lastME2(0.0);
  return lastME2();

}

double MatchboxMElP2lJetJet::spinColourCorrelatedME2(pair<int,int> ij,
						     const SpinCorrelationTensor& c) const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::spinColourCorrelatedME2(ij,c);

  generator()->logWarning(Exception() 
			  << "The matrix element '" << name() << "' "
			  << "is not capable of calculating colour- or spin correlated "
			  << "matrix element squares."
			  << Exception::warning);

  lastME2(0.0);
  return lastME2();

}


void MatchboxMElP2lJetJet::doinit() {
  MatchboxMEBase::doinit();
  MatchboxMEllbarqqbarg::doinit(*this);
  for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	q != theQuarkFlavours.end(); ++q )
    if ( (**q).mass() != ZERO )
      Throw<InitException>() << "The matrix element '"
			     << name() << "' is only capable of "
			     << "producing massless quarks.";
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxMElP2lJetJet::persistentOutput(PersistentOStream & os) const {
  MatchboxMEllbarqqbarg::persistentOutput(os);
  os << theLeptonFlavours << theQuarkFlavours << ounit(theUserScale,GeV);
}

void MatchboxMElP2lJetJet::persistentInput(PersistentIStream & is, int) {
  MatchboxMEllbarqqbarg::persistentInput(is);
  is >> theLeptonFlavours >> theQuarkFlavours >> iunit(theUserScale,GeV);
}

AbstractClassDescription<MatchboxMElP2lJetJet> MatchboxMElP2lJetJet::initMatchboxMElP2lJetJet;
// Definition of the static class description member.

void MatchboxMElP2lJetJet::Init() {

  static ClassDocumentation<MatchboxMElP2lJetJet> documentation
    ("MatchboxMElP2lJetJet");

  static RefVector<MatchboxMElP2lJetJet,ParticleData> interfaceLeptonFlavours
    ("LeptonFlavours",
     "The lepton flavours for this matrix element.",
     &MatchboxMElP2lJetJet::theLeptonFlavours, -1, false, false, true, true, false);

  static RefVector<MatchboxMElP2lJetJet,ParticleData> interfaceQuarkFlavours
    ("QuarkFlavours",
     "The quark flavours for this matrix element.",
     &MatchboxMElP2lJetJet::theQuarkFlavours, -1, false, false, true, true, false);


  static Parameter<MatchboxMElP2lJetJet,Energy> interfaceUserScale
    ("UserScale",
     "A user defined renormalization scale.",
     &MatchboxMElP2lJetJet::theUserScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

}

