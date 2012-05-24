// -*- C++ -*-
//
// MatchboxMEllbar2qqbar.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMEllbar2qqbar class.
//

#include "MatchboxMEllbar2qqbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <gsl/gsl_sf_dilog.h>

using namespace Herwig;

MatchboxMEllbar2qqbar::MatchboxMEllbar2qqbar() 
  : MatchboxMEBase(), theUserScale(0.0*GeV) {}

MatchboxMEllbar2qqbar::~MatchboxMEllbar2qqbar() {}

IBPtr MatchboxMEllbar2qqbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMEllbar2qqbar::fullclone() const {
  return new_ptr(*this);
}

double MatchboxMEllbar2qqbar::me2() const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::me2();

  double res;
  if ( !calculateME2(res) )
    return res;

  Lorentz5Momentum pq = 
    mePartonData()[2]->id() > 0 ? 
    meMomenta()[2] : 
    meMomenta()[3];

  Lorentz5Momentum pqbar = 
    mePartonData()[2]->id() > 0 ? 
    meMomenta()[3] : 
    meMomenta()[2];

  Lorentz5Momentum pl = 
    mePartonData()[0]->id() > 0 ? 
    - meMomenta()[1] : 
    - meMomenta()[0];

  Lorentz5Momentum plbar = 
    mePartonData()[0]->id() > 0 ? 
    - meMomenta()[0] : 
    - meMomenta()[1];

  prepare(pl,plbar,pq,pqbar,lastSHat(),mePartonData()[0],mePartonData()[2]);

  lastME2(evaluateME2()*me2Norm());
  cacheME2(lastME2());

  logME2();

  return lastME2();

}

double MatchboxMEllbar2qqbar::oneLoopInterference() const {
  
  // case of massless quarks stolen from MatchboxMEllbarqqbarVirtual::me2()
  if ( mass(2)==0 && mass(3)==0 ) {
    return 
      (lastAlphaS()/(2.*Constants::pi)) *
      ((SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc())) *
      ( sqr(Constants::pi) - 8. ) * 
      me2();
  }
  
  // massive quarks (same mass) & massless leptons
  // taken from Jersak,Laermann,Zerwas, Phys.Rev.D 25, pp.1218-1228
  double v = sqrt( 1.-4*sqr(mass(2)) );
  double realf1 = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc()) *
    lastAlphaS()/(2*Constants::pi) *
    ( -2 - (1.+2.*v*v)/(2*v)*log((1.-v)/(1.+v)) +
      (1.+v*v)/v * ( gsl_sf_dilog((1.-v)/(1.+v)) + sqr(Constants::pi)/3. -
		     1./4.*sqr(log((1.-v)/(1.+v))) +
		     log((1.-v)/(1.+v))*log(2.*v/(1.+v)) )
    // last term from different rausziehen of eps terms (right scaling of mass???)
//       -(1.+(1.+v*v)/(2.*v)*log((1.-v)/(1.+v))) * log(1./sqr(mass(2))) );
//       +(1.+(1.+v*v)/(2.*v)*log((1.-v)/(1.+v))) * log(s/(4.*Constants::pi*renScale*renScale)) );
      -(1.+(1.+v*v)/(2.*v)*log((1.-v)/(1.+v))) * log(4./(1.-v*v)) );
  double realf2 = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc()) *
    lastAlphaS()/(2*Constants::pi) *
    (1.-v*v)/(2.*v) * log((1.-v)/(1.+v));
  double imagf2 = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc()) *
    lastAlphaS()/(2*Constants::pi) *
    (1.-v*v)/(2*v) * Constants::pi;
    
  double xQ2 = (momentum(0)+momentum(1)).m2();
  double M2 = sqr(BMass/amplitudeScale());
  double G2 = sqr(BWidth/amplitudeScale());
  
  double sVV = 32./3.*SM().Nc()/v *
    ( M2*(G2+M2)*sqr(el)*sqr(eq) - 2.*M2*el*eq*(vl*vq+el*eq)*xQ2 +
      (sqr(al)*sqr(vq)+sqr(vl*vq+el*eq))*sqr(xQ2) ) /
    ( sqr(xQ2-M2) + M2*G2 );
  double sAA = 32./3.*SM().Nc()/v *
    sqr(aq)*(sqr(al)*sqr(vl))*sqr(xQ2) / ( sqr(xQ2-M2) + M2*G2 );
  double sVA = 64./3.*SM().Nc()/v *
    ( 2.*vl*vq*xQ2 + el*eq*(xQ2-M2) ) / ( sqr(xQ2-M2) + M2*G2 );
  double sVAprime = 64./3.*SM().Nc()/v *
    al*aq*el*eq*sqrt(G2)*sqrt(M2) / ( sqr(xQ2-M2) + M2*G2 );
  
  double cost = ( -1. + 2.*invariant(0,2) ) /
    sqrt( (1.-4.*sqr(mass(0))) * (1.-4.*sqr(mass(2))) );
  
  return 2.*realf1*me2() +
    3./8.*(1.+sqr(cost)) * 2.*v*realf2*(sVV-v*v*sAA) +
    3./4.*(1.-sqr(cost)) * v*realf2*sVV +
    3./4.*cost*v*v*sVA * (-2)*v*v*sVAprime*imagf2;
}

Energy2 MatchboxMEllbar2qqbar::factorizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return lastSHat();

}

Energy2 MatchboxMEllbar2qqbar::renormalizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return lastSHat();

}

double MatchboxMEllbar2qqbar::colourCorrelatedME2(pair<int,int> ij) const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::colourCorrelatedME2(ij);

  if ( ij.first == ij.second ||
       ij.first < 2 || ij.first > 3 ||
       ij.second < 2 || ij.second > 3 ) {
    generator()->logWarning(Exception() 
			    << "A non-exisiting colour correlation was requested "
			    << "from the matrix element '" << name() << "'."
			    << Exception::warning);
    lastME2(0.0);
    return lastME2();  
  }

  return -me2();

}

double MatchboxMEllbar2qqbar::spinColourCorrelatedME2(pair<int,int>,
						      const SpinCorrelationTensor&) const {

  generator()->logWarning(Exception() 
			  << "A non-exisiting spin correlation was requested "
			  << "from the matrix element '" << name() << "'."
			  << Exception::warning);

  lastME2(0.0);
  return lastME2();  

}

void MatchboxMEllbar2qqbar::getDiagrams() const {

  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);

  for ( PDVector::const_iterator l = theLeptonFlavours.begin();
	l != theLeptonFlavours.end(); ++l )
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ) {
      add(new_ptr((Tree2toNDiagram(2), *l, (**l).CC(), 1, gamma, 3, *q, 3, (**q).CC(), -1)));
      add(new_ptr((Tree2toNDiagram(2), *l, (**l).CC(), 1, Z0, 3, *q, 3, (**q).CC(), -2)));
    }

}

Selector<MEBase::DiagramIndex> 
MatchboxMEllbar2qqbar::diagrams(const DiagramVector &) const {
  Selector<MEBase::DiagramIndex> sel;
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  double wGamma = 
    sqr( sqr(generator()->maximumCMEnergy())/lastSHat() );
  double wZ = 
    sqr( sqr(generator()->maximumCMEnergy()) ) /
    ( sqr(lastSHat()-sqr(Z0->mass())) + sqr(Z0->mass())*sqr(Z0->width()) );
  sel.insert(wGamma,0);
  sel.insert(wZ,1);
  return sel;
}

Selector<const ColourLines *>
MatchboxMEllbar2qqbar::colourGeometries(tcDiagPtr) const {
  static const ColourLines c("4 -5");
  static const ColourLines cbar("-4 5");
  Selector<const ColourLines *> sel;
  if ( mePartonData()[2]->id() > 0 )
    sel.insert(1.0, &c);
  else
    sel.insert(1.0, &cbar);
  return sel;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxMEllbar2qqbar::doinit() {
  MatchboxMEBase::doinit();
  MatchboxMEllbarqqbar::doinit(*this);
}

void MatchboxMEllbar2qqbar::persistentOutput(PersistentOStream & os) const {
  MatchboxMEllbarqqbar::persistentOutput(os);
  os << theLeptonFlavours << theQuarkFlavours << ounit(theUserScale,GeV);
}

void MatchboxMEllbar2qqbar::persistentInput(PersistentIStream & is, int) {
  MatchboxMEllbarqqbar::persistentInput(is);
  is >> theLeptonFlavours >> theQuarkFlavours >> iunit(theUserScale,GeV);
}

ClassDescription<MatchboxMEllbar2qqbar> MatchboxMEllbar2qqbar::initMatchboxMEllbar2qqbar;
// Definition of the static class description member.

void MatchboxMEllbar2qqbar::Init() {

  static ClassDocumentation<MatchboxMEllbar2qqbar> documentation
    ("MatchboxMEllbar2qqbar");

  static RefVector<MatchboxMEllbar2qqbar,ParticleData> interfaceLeptonFlavours
    ("LeptonFlavours",
     "The lepton flavours for this matrix element.",
     &MatchboxMEllbar2qqbar::theLeptonFlavours, -1, false, false, true, true, false);

  static RefVector<MatchboxMEllbar2qqbar,ParticleData> interfaceQuarkFlavours
    ("QuarkFlavours",
     "The quark flavours for this matrix element.",
     &MatchboxMEllbar2qqbar::theQuarkFlavours, -1, false, false, true, true, false);


  static Parameter<MatchboxMEllbar2qqbar,Energy> interfaceUserScale
    ("UserScale",
     "A user defined renormalization scale.",
     &MatchboxMEllbar2qqbar::theUserScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


}

