// -*- C++ -*-
//
// MatchboxMEPP2llbarJet.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMEPP2llbarJet class.
//

#include "MatchboxMEPP2llbarJet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxMEPP2llbarJet::MatchboxMEPP2llbarJet() 
  : MatchboxMEBase(), 
    theZMass2(0.0*GeV2), theZWidth2(0.0*GeV2),
    theUserScale(0.0*GeV) {}

MatchboxMEPP2llbarJet::~MatchboxMEPP2llbarJet() {}

bool MatchboxMEPP2llbarJet::generateKinematics(const double * r) {

  if ( phasespace() )
    return MatchboxMEBase::generateKinematics(r);

  if ( theZMass2 == 0.0*GeV2 ) {
    theZMass2 = sqr(getParticleData(ParticleID::Z0)->mass());
    theZWidth2 = sqr(getParticleData(ParticleID::Z0)->width());
  }

  Energy Q = sqrt(lastSHat());
  Energy m1 = mePartonData()[2]->mass();
  Energy m2 = mePartonData()[3]->mass();

  // get the z mass and jacobian to map out the BW

  Energy QZmin = m1 + m2;
  Energy QZmax = Q;

  double muZ = theZMass2/sqr(Q);
  double gammaZ = theZWidth2/sqr(Q);

  double xqMin = atan((sqr(QZmin/Q)-muZ)/sqrt(muZ*gammaZ))/sqrt(muZ*gammaZ);
  double xqMax = atan((sqr(QZmax/Q)-muZ)/sqrt(muZ*gammaZ))/sqrt(muZ*gammaZ);
  double xq = xqMin+r[0]*(xqMax-xqMin);

  double qZ = muZ + sqrt(muZ*gammaZ)*tan(xq*sqrt(muZ*gammaZ));
  Energy QZ = Q*sqrt(qZ);

  double bwWeight = sqr(qZ-muZ)+muZ*gammaZ;
  bwWeight *= xqMax-xqMin;

  double zMassWeight = bwWeight/(2.*Constants::pi);

  // the parton parton -> Z parton phasespace
  // map out the 1/pt^2

  Energy ptMin = 1.0*GeV;
  Energy ptMax = (Q/2.)*(1.-sqr(QZ/Q));

  double xtMin = 2.*log(ptMin/Q);
  double xtMax = log(sqr(ptMax/Q)+sqr(ptMin/Q));
  double xt = xtMin + r[1]*(xtMax-xtMin);

  Energy pt = Q*sqrt(exp(xt)-sqr(ptMin/Q));

  double mapPt = (sqr(pt)+sqr(ptMin))/sqr(Q);
  mapPt *= xtMax-xtMin;

  double phiJet = 2.*Constants::pi*r[2];
  double cPhiJet = cos(phiJet);
  double sPhiJet =
    phiJet < Constants::pi ? 
    sqrt(1.-sqr(cPhiJet)) : 
    - sqrt(1.-sqr(cPhiJet));

  Energy eJet = (sqr(Q)-sqr(QZ))/(2.*Q);
  double hemisphere = rnd() < 0.5 ? 1.0 : -1.0;

  Energy pJet3z = hemisphere*sqrt(sqr(eJet)-sqr(pt));

  if ( abs((eJet-pt)/GeV) < Constants::epsilon )
    pJet3z = ZERO;

  ThreeVector<Energy> pJet3 (pt*cPhiJet,pt*sPhiJet,pJet3z);

  if ( sqr(1.-sqr(QZ/Q))-4.*sqr(pt/Q) <= 0.0 ) {
    jacobian(0.0);
    return false;
  }

  double twoBodyWeight = mapPt*(1./sqrt(sqr(1.-sqr(QZ/Q))-4.*sqr(pt/Q)))/(8.*Constants::pi);

  Lorentz5Momentum pJet(0.0*GeV,pJet3);
  Lorentz5Momentum pZ(QZ,-pJet3);

  // decay the leptons in the Z rest frame

  double phiLepton = 2.*Constants::pi*r[3];
  double cPhiLepton = cos(phiLepton);
  double sPhiLepton =
    phiLepton < Constants::pi ? 
    sqrt(1.-sqr(cPhiLepton)) : 
    - sqrt(1.-sqr(cPhiLepton));

  double ctLepton = 2.*r[4]-1.;
  double stLepton = sqrt(1.-sqr(ctLepton));

  Energy qLepton = sqrt((sqr(QZ)-sqr(m1-m2))*(sqr(QZ)-sqr(m1+m2)))/(2.*QZ);

  double decayWeight = (qLepton/QZ)/(4.*Constants::pi);

  Axis nLepton(cPhiLepton*stLepton,sPhiLepton*stLepton,ctLepton);

  Lorentz5Momentum pLepton(m1,qLepton*nLepton);
  Lorentz5Momentum pLeptonBar(m2,-qLepton*nLepton);

  Boost leptonBoost = pZ.boostVector();
  pLepton.boost(leptonBoost);
  pLeptonBar.boost(leptonBoost);

  meMomenta()[2] = pLepton;
  meMomenta()[3] = pLeptonBar;
  meMomenta()[4] = pJet;

  jacobian(zMassWeight*twoBodyWeight*decayWeight);

  setScale();
  logGenerateKinematics(r);
  return true;

}

Energy2 MatchboxMEPP2llbarJet::factorizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return (meMomenta()[2]+meMomenta()[3]).m2();

}

Energy2 MatchboxMEPP2llbarJet::renormalizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return (meMomenta()[2]+meMomenta()[3]).m2();

}

double MatchboxMEPP2llbarJet::colourCorrelatedME2(pair<int,int> ij) const {

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

double MatchboxMEPP2llbarJet::spinColourCorrelatedME2(pair<int,int> ij,
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


void MatchboxMEPP2llbarJet::doinit() {
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


void MatchboxMEPP2llbarJet::persistentOutput(PersistentOStream & os) const {
  MatchboxMEllbarqqbarg::persistentOutput(os);
  os << theLeptonFlavours << theQuarkFlavours << ounit(theUserScale,GeV);
}

void MatchboxMEPP2llbarJet::persistentInput(PersistentIStream & is, int) {
  MatchboxMEllbarqqbarg::persistentInput(is);
  is >> theLeptonFlavours >> theQuarkFlavours >> iunit(theUserScale,GeV);
}

AbstractClassDescription<MatchboxMEPP2llbarJet> MatchboxMEPP2llbarJet::initMatchboxMEPP2llbarJet;
// Definition of the static class description member.

void MatchboxMEPP2llbarJet::Init() {

  static ClassDocumentation<MatchboxMEPP2llbarJet> documentation
    ("MatchboxMEPP2llbarJet");

  static RefVector<MatchboxMEPP2llbarJet,ParticleData> interfaceLeptonFlavours
    ("LeptonFlavours",
     "The lepton flavours for this matrix element.",
     &MatchboxMEPP2llbarJet::theLeptonFlavours, -1, false, false, true, true, false);

  static RefVector<MatchboxMEPP2llbarJet,ParticleData> interfaceQuarkFlavours
    ("QuarkFlavours",
     "The quark flavours for this matrix element.",
     &MatchboxMEPP2llbarJet::theQuarkFlavours, -1, false, false, true, true, false);


  static Parameter<MatchboxMEPP2llbarJet,Energy> interfaceUserScale
    ("UserScale",
     "A user defined renormalization scale.",
     &MatchboxMEPP2llbarJet::theUserScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

}

