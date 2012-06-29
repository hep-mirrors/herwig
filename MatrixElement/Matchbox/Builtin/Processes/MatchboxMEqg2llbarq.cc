// -*- C++ -*-
//
// MatchboxMEqg2llbarq.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMEqg2llbarq class.
//

#include "MatchboxMEqg2llbarq.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxMEqg2llbarq::MatchboxMEqg2llbarq() 
  : MatchboxMEPP2llbarJet(), theWhichGluon(0) {}

MatchboxMEqg2llbarq::~MatchboxMEqg2llbarq() {}

IBPtr MatchboxMEqg2llbarq::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMEqg2llbarq::fullclone() const {
  return new_ptr(*this);
}

double MatchboxMEqg2llbarq::me2() const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::me2();

  double res;
  if ( !calculateME2(res) )
    return res;

  cPDPtr qid =
    theWhichGluon == 0 ?
    mePartonData()[1] :
    mePartonData()[0];

  Lorentz5Momentum pq;

  if ( theWhichGluon == 0 ) {
    if ( mePartonData()[1]->id() < 0 )
      pq = -meMomenta()[1];
    else
      pq = meMomenta()[4];
  } else {
    if ( mePartonData()[0]->id() < 0 )
      pq = -meMomenta()[0];
    else
      pq = meMomenta()[4];
  }

  Lorentz5Momentum pqbar;

  if ( theWhichGluon == 0 ) {
    if ( mePartonData()[1]->id() > 0 )
      pqbar = -meMomenta()[1];
    else
      pqbar = meMomenta()[4];
  } else {
    if ( mePartonData()[0]->id() > 0 )
      pqbar = -meMomenta()[0];
    else
      pqbar = meMomenta()[4];
  }

  Lorentz5Momentum pl = 
    mePartonData()[2]->id() > 0 ? 
    meMomenta()[2] : 
    meMomenta()[3];

  Lorentz5Momentum plbar = 
    mePartonData()[2]->id() < 0 ? 
    meMomenta()[2] : 
    meMomenta()[3];

  Lorentz5Momentum pg = - (
    theWhichGluon == 0 ?
    meMomenta()[0] :
    meMomenta()[1] );

  prepare(pl,plbar,pq,pqbar,pg,lastSHat(),mePartonData()[2],qid);

  lastME2(-MatchboxMEllbarqqbarg::evaluateME2()*me2Norm());
  cacheME2(lastME2());

  logME2();

  return lastME2();

}

void MatchboxMEqg2llbarq::getDiagrams() const {

  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);

  for ( PDVector::const_iterator l = theLeptonFlavours.begin();
	l != theLeptonFlavours.end(); ++l )
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ) {
      if ( theWhichGluon == 0 ) {
	add(new_ptr((Tree2toNDiagram(3), g, (**q).CC(), *q, 2, gamma, 4, *l, 4, (**l).CC(), 1, *q, -1)));
	add(new_ptr((Tree2toNDiagram(2), g, *q, 1, *q, 3, gamma, 4, *l, 4, (**l).CC(), 3, *q, -2)));
	add(new_ptr((Tree2toNDiagram(3), g, (**q).CC(), *q, 2, Z0, 4, *l, 4, (**l).CC(), 1, *q, -3)));
	add(new_ptr((Tree2toNDiagram(2), g, *q, 1, *q, 3, Z0, 4, *l, 4, (**l).CC(), 3, *q, -4)));
      } else {
	add(new_ptr((Tree2toNDiagram(3), *q, *q, g, 1, gamma, 4, *l, 4, (**l).CC(), 2, *q, -1)));
	add(new_ptr((Tree2toNDiagram(2), *q, g, 1, *q, 3, gamma, 4, *l, 4, (**l).CC(), 3, *q, -2)));
	add(new_ptr((Tree2toNDiagram(3), *q, *q, g, 1, Z0, 4, *l, 4, (**l).CC(), 2, *q, -3)));
	add(new_ptr((Tree2toNDiagram(2), *q, g, 1, *q, 3, Z0, 4, *l, 4, (**l).CC(), 3, *q, -4)));
      }
    }

}

Selector<MEBase::DiagramIndex> 
MatchboxMEqg2llbarq::diagrams(const DiagramVector &) const {
  Selector<MEBase::DiagramIndex> sel;
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  double wGamma = 
    sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[2]*meMomenta()[3])) );
  double wZ = 
    sqr( sqr(generator()->maximumCMEnergy()) ) /
    ( sqr(2.*(meMomenta()[2]*meMomenta()[3])-sqr(Z0->mass())) + sqr(Z0->mass())*sqr(Z0->width()) );
  Energy2 wt = 
    theWhichGluon == 0 ?
    meMomenta()[0]*meMomenta()[4] : meMomenta()[1]*meMomenta()[4];
  if ( wt < lastSHat() ) {
    sel.insert(wGamma,0);
    sel.insert(wZ,2);
  } else {
    sel.insert(wGamma,1);
    sel.insert(wZ,3);
  }
  return sel;
}

Selector<const ColourLines *>
MatchboxMEqg2llbarq::colourGeometries(tcDiagPtr d) const {
  Selector<const ColourLines *> sel;
  if ( theWhichGluon == 0 ) {
    static const ColourLines ct("1 7, -1 -2 3"); 
    static const ColourLines cbart("-1 -7, 1 2 -3"); 
    static const ColourLines cs("-1 2, 1 3 7"); 
    static const ColourLines cbars("1 -2, -1 -3 -7"); 
    if ( mePartonData()[1]->id() > 0 ) {
      if ( abs(d->id()) % 2 != 0 )
	sel.insert(1.0, &ct);
      else
	sel.insert(1.0, &cs);
    } else {
      if ( abs(d->id()) % 2 != 0 )
	sel.insert(1.0, &cbart);
      else
	sel.insert(1.0, &cbars);
    }
  } else {
    static const ColourLines ct("1 2 -3, 3 7"); 
    static const ColourLines cbart("-1 -2 3, -3 -7"); 
    static const ColourLines cs("1 -2, 2 3 7"); 
    static const ColourLines cbars("-1 2, -2 -3 -7"); 
    if ( mePartonData()[0]->id() > 0 ) {
      if ( abs(d->id()) % 2 != 0 )
	sel.insert(1.0, &ct);
      else
	sel.insert(1.0, &cs);
    } else {
      if ( abs(d->id()) % 2 != 0 )
	sel.insert(1.0, &cbart);
      else
	sel.insert(1.0, &cbars);
    }
  }
  return sel;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxMEqg2llbarq::persistentOutput(PersistentOStream & os) const {
  os << theWhichGluon;
}

void MatchboxMEqg2llbarq::persistentInput(PersistentIStream & is, int) {
  is >> theWhichGluon;
}

ClassDescription<MatchboxMEqg2llbarq> MatchboxMEqg2llbarq::initMatchboxMEqg2llbarq;
// Definition of the static class description member.

void MatchboxMEqg2llbarq::Init() {

  static ClassDocumentation<MatchboxMEqg2llbarq> documentation
    ("MatchboxMEqg2llbarq");

  static Switch<MatchboxMEqg2llbarq,int> interfaceWhichGluon
    ("WhichGluon",
     "Set the position of the incoming gluon.",
     &MatchboxMEqg2llbarq::theWhichGluon, 0, false, false);
  static SwitchOption interfaceWhichGluonFirst
    (interfaceWhichGluon,
     "First",
     "From first incoming hadron.",
     0);
  static SwitchOption interfaceWhichGluonSecond
    (interfaceWhichGluon,
     "Second",
     "From second incoming hadron.",
     1);

}

