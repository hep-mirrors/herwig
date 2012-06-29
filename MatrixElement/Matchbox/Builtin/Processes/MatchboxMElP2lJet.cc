// -*- C++ -*-
//
// MatchboxMElP2lJet.cc.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMElP2lJet class.
//

#include "MatchboxMElP2lJet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxMElP2lJet::MatchboxMElP2lJet() 
  : MatchboxMEBase(), theUserScale(0.0*GeV) {}

MatchboxMElP2lJet::~MatchboxMElP2lJet() {}

IBPtr MatchboxMElP2lJet::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMElP2lJet::fullclone() const {
  return new_ptr(*this);
}

double MatchboxMElP2lJet::me2() const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::me2();

  double res;
  if ( !calculateME2(res) )
    return res;

  Lorentz5Momentum pq = 
    mePartonData()[1]->id() > 0 ? 
    meMomenta()[3] : 
    meMomenta()[1];

  if ( mePartonData()[1]->id() < 0 )
    pq = -pq;

  Lorentz5Momentum pqbar = 
    mePartonData()[1]->id() < 0 ? 
    meMomenta()[3] : 
    meMomenta()[1];

  if ( mePartonData()[1]->id() > 0 )
    pqbar = -pqbar;

  Lorentz5Momentum pl = 
    mePartonData()[0]->id() > 0 ? 
    meMomenta()[2] : 
    meMomenta()[0];

  if ( mePartonData()[0]->id() < 0 )
    pl = -pl;

  Lorentz5Momentum plbar = 
    mePartonData()[0]->id() < 0 ? 
    meMomenta()[2] : 
    meMomenta()[0];

  if ( mePartonData()[0]->id() > 0 )
    plbar = -plbar;

  prepare(pl,plbar,pq,pqbar,lastSHat(),mePartonData()[0],mePartonData()[1]);

  lastME2(evaluateME2()*me2Norm());
  cacheME2(lastME2());

  logME2();

  return lastME2();

}

Energy2 MatchboxMElP2lJet::factorizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return -(meMomenta()[0]-meMomenta()[2]).m2();

}

Energy2 MatchboxMElP2lJet::renormalizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return -(meMomenta()[0]-meMomenta()[2]).m2();

}

double MatchboxMElP2lJet::colourCorrelatedME2(pair<int,int> ij) const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::colourCorrelatedME2(ij);

  if ( !((ij.first == 1 && ij.second == 3) || 
	 (ij.second == 1 && ij.first == 3) ) ) { 
    generator()->logWarning(Exception()  
			    << "A non-exisiting colour correlation was requested " 
			    << "from the matrix element '" << name() << "'." 
			    << Exception::warning); 
    lastME2(0.0); 
    return lastME2(); 
  } 
	 	 
  return -me2(); 

}

double MatchboxMElP2lJet::spinColourCorrelatedME2(pair<int,int>,
						  const SpinCorrelationTensor&) const {

  generator()->logWarning(Exception() 
			  << "A non-exisiting spin correlation was requested "
			  << "from the matrix element '" << name() << "'."
			  << Exception::warning);

  lastME2(0.0);
  return lastME2();

}

void MatchboxMElP2lJet::getDiagrams() const {

  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);

  for ( PDVector::const_iterator l = theLeptonFlavours.begin();
	l != theLeptonFlavours.end(); ++l )
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ) {
      add(new_ptr((Tree2toNDiagram(3), *l, gamma, *q, 1, *l, 2, *q, -1)));
      add(new_ptr((Tree2toNDiagram(3), *l, Z0, *q, 1, *l, 2, *q, -2)));
    }

}

Selector<const ColourLines *>
MatchboxMElP2lJet::colourGeometries(tcDiagPtr) const {
  static const ColourLines c("3 5");
  static const ColourLines cbar("-5 -3");
  Selector<const ColourLines *> sel;
  if ( mePartonData()[1]->id() > 0 )
    sel.insert(1.0, &c);
  else
    sel.insert(1.0, &cbar);
  return sel;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxMElP2lJet::doinit() {
  MatchboxMEBase::doinit();
  MatchboxMEllbarqqbar::doinit(*this);
}

void MatchboxMElP2lJet::persistentOutput(PersistentOStream & os) const {
  MatchboxMEllbarqqbar::persistentOutput(os);
  os << theLeptonFlavours << theQuarkFlavours << ounit(theUserScale,GeV);
}

void MatchboxMElP2lJet::persistentInput(PersistentIStream & is, int) {
  MatchboxMEllbarqqbar::persistentInput(is);
  is >> theLeptonFlavours >> theQuarkFlavours >> iunit(theUserScale,GeV);
}

ClassDescription<MatchboxMElP2lJet> MatchboxMElP2lJet::initMatchboxMElP2lJet;
// Definition of the static class description member.

void MatchboxMElP2lJet::Init() {

  static ClassDocumentation<MatchboxMElP2lJet> documentation
    ("MatchboxMElP2lJet");

  static RefVector<MatchboxMElP2lJet,ParticleData> interfaceLeptonFlavours
    ("LeptonFlavours",
     "The lepton flavours for this matrix element.",
     &MatchboxMElP2lJet::theLeptonFlavours, -1, false, false, true, true, false);

  static RefVector<MatchboxMElP2lJet,ParticleData> interfaceQuarkFlavours
    ("QuarkFlavours",
     "The quark flavours for this matrix element.",
     &MatchboxMElP2lJet::theQuarkFlavours, -1, false, false, true, true, false);


  static Parameter<MatchboxMElP2lJet,Energy> interfaceUserScale
    ("UserScale",
     "A user defined renormalization scale.",
     &MatchboxMElP2lJet::theUserScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


}

