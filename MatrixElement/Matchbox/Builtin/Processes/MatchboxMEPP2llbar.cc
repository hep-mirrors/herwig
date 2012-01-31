// -*- C++ -*-
//
// MatchboxMEPP2llbar.cc.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMEPP2llbar class.
//

#include "MatchboxMEPP2llbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxMEPP2llbar::MatchboxMEPP2llbar() 
  : MatchboxMEBase(), theUserScale(0.0*GeV) {}

MatchboxMEPP2llbar::~MatchboxMEPP2llbar() {}

IBPtr MatchboxMEPP2llbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMEPP2llbar::fullclone() const {
  return new_ptr(*this);
}

double MatchboxMEPP2llbar::me2() const {

  double res;
  if ( !calculateME2(res) )
    return res;

  Lorentz5Momentum pq = - (
    mePartonData()[0]->id() > 0 ? 
    meMomenta()[1] : 
    meMomenta()[0] );

  Lorentz5Momentum pqbar = - (
    mePartonData()[1]->id() > 0 ? 
    meMomenta()[1] : 
    meMomenta()[0] );

  Lorentz5Momentum pl = 
    mePartonData()[2]->id() > 0 ? 
    meMomenta()[2] : 
    meMomenta()[3];

  Lorentz5Momentum plbar = 
    mePartonData()[2]->id() > 0 ? 
    meMomenta()[3] : 
    meMomenta()[2];

  prepare(pl,plbar,pq,pqbar,lastSHat(),mePartonData()[2],mePartonData()[0]);

  lastME2(evaluateME2()*me2Norm());
  cacheME2(lastME2());

  logME2();

  return lastME2();

}

Energy2 MatchboxMEPP2llbar::factorizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return lastSHat();

}

Energy2 MatchboxMEPP2llbar::renormalizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return lastSHat();

}

double MatchboxMEPP2llbar::colourCorrelatedME2(pair<int,int> ij) const {

  if ( ij.first == ij.second ||
       ij.first > 1 ||
       ij.second > 1 ) {
    generator()->logWarning(Exception() 
			    << "A non-exisiting colour correlation was requested "
			    << "from the matrix element '" << name() << "'."
			    << Exception::warning);
    lastME2(0.0);
    return lastME2();
  }

  return -me2();

}

double MatchboxMEPP2llbar::spinColourCorrelatedME2(pair<int,int>,
						   const SpinCorrelationTensor&) const {

  generator()->logWarning(Exception() 
			  << "A non-exisiting spin correlation was requested "
			  << "from the matrix element '" << name() << "'."
			  << Exception::warning);

  lastME2(0.0);
  return lastME2();

}

void MatchboxMEPP2llbar::getDiagrams() const {

  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);

  for ( PDVector::const_iterator l = theLeptonFlavours.begin();
	l != theLeptonFlavours.end(); ++l )
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ) {
      add(new_ptr((Tree2toNDiagram(2), *q, (**q).CC(), 1, gamma, 3, *l, 3, (**l).CC(), -1)));
      add(new_ptr((Tree2toNDiagram(2), *q, (**q).CC(), 1,    Z0, 3, *l, 3, (**l).CC(), -2)));
    }

}

Selector<MEBase::DiagramIndex> 
MatchboxMEPP2llbar::diagrams(const DiagramVector & diags) const {
  Selector<MEBase::DiagramIndex> sel;
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  double wGamma = 
    sqr( sqr(generator()->maximumCMEnergy())/lastSHat() );
  double wZ = 
    sqr( sqr(generator()->maximumCMEnergy()) ) /
    ( sqr(lastSHat()-sqr(Z0->mass())) + sqr(Z0->mass())*sqr(Z0->width()) );
  assert(diags.size() == 2);
  sel.insert(wGamma,0);
  sel.insert(wZ,1);
  return sel;
}


Selector<const ColourLines *>
MatchboxMEPP2llbar::colourGeometries(tcDiagPtr) const {
  static const ColourLines c("1 -2");
  static const ColourLines cbar("-1 2");
  Selector<const ColourLines *> sel;
  if ( mePartonData()[0]->id() > 0 )
    sel.insert(1.0, &c);
  else
    sel.insert(1.0, &cbar);
  return sel;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxMEPP2llbar::doinit() {
  MatchboxMEBase::doinit();
  MatchboxMEllbarqqbar::doinit(*this);
}

void MatchboxMEPP2llbar::persistentOutput(PersistentOStream & os) const {
  MatchboxMEllbarqqbar::persistentOutput(os);
  os << theLeptonFlavours << theQuarkFlavours << ounit(theUserScale,GeV);
}

void MatchboxMEPP2llbar::persistentInput(PersistentIStream & is, int) {
  MatchboxMEllbarqqbar::persistentInput(is);
  is >> theLeptonFlavours >> theQuarkFlavours >> iunit(theUserScale,GeV);
}

ClassDescription<MatchboxMEPP2llbar> MatchboxMEPP2llbar::initMatchboxMEPP2llbar;
// Definition of the static class description member.

void MatchboxMEPP2llbar::Init() {

  static ClassDocumentation<MatchboxMEPP2llbar> documentation
    ("MatchboxMEPP2llbar");

  static RefVector<MatchboxMEPP2llbar,ParticleData> interfaceLeptonFlavours
    ("LeptonFlavours",
     "The lepton flavours for this matrix element.",
     &MatchboxMEPP2llbar::theLeptonFlavours, -1, false, false, true, true, false);

  static RefVector<MatchboxMEPP2llbar,ParticleData> interfaceQuarkFlavours
    ("QuarkFlavours",
     "The quark flavours for this matrix element.",
     &MatchboxMEPP2llbar::theQuarkFlavours, -1, false, false, true, true, false);


  static Parameter<MatchboxMEPP2llbar,Energy> interfaceUserScale
    ("UserScale",
     "A user defined renormalization scale.",
     &MatchboxMEPP2llbar::theUserScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


}

