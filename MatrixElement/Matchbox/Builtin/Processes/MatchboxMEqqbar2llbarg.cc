// -*- C++ -*-
//
// MatchboxMEqqbar2llbarg.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMEqqbar2llbarg class.
//

#include "MatchboxMEqqbar2llbarg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxMEqqbar2llbarg::MatchboxMEqqbar2llbarg() 
  : MatchboxMEPP2llbarJet() {}

MatchboxMEqqbar2llbarg::~MatchboxMEqqbar2llbarg() {}

IBPtr MatchboxMEqqbar2llbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMEqqbar2llbarg::fullclone() const {
  return new_ptr(*this);
}

double MatchboxMEqqbar2llbarg::me2() const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::me2();

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

  Lorentz5Momentum pg = meMomenta()[4];

  prepare(pl,plbar,pq,pqbar,pg,lastSHat(),mePartonData()[2],mePartonData()[0]);

  lastME2(MatchboxMEllbarqqbarg::evaluateME2()*me2Norm());
  cacheME2(lastME2());

  logME2();

  return lastME2();

}

void MatchboxMEqqbar2llbarg::getDiagrams() const {

  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);

  for ( PDVector::const_iterator l = theLeptonFlavours.begin();
	l != theLeptonFlavours.end(); ++l )
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ) {
      add(new_ptr((Tree2toNDiagram(3), *q, *q, (**q).CC(), 2, gamma, 4, *l, 4, (**l).CC(), 1, g, -1)));
      add(new_ptr((Tree2toNDiagram(3), *q, *q, (**q).CC(), 1, gamma, 4, *l, 4, (**l).CC(), 2, g, -2)));
      add(new_ptr((Tree2toNDiagram(3), *q, *q, (**q).CC(), 2, Z0   , 4, *l, 4, (**l).CC(), 1, g, -3)));
      add(new_ptr((Tree2toNDiagram(3), *q, *q, (**q).CC(), 1, Z0   , 4, *l, 4, (**l).CC(), 2, g, -4)));
    }

}

Selector<MEBase::DiagramIndex> 
MatchboxMEqqbar2llbarg::diagrams(const DiagramVector &) const {
  Selector<MEBase::DiagramIndex> sel;
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  double wGamma = 
    sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[2]*meMomenta()[3])) );
  double wZ = 
    sqr( sqr(generator()->maximumCMEnergy()) ) /
    ( sqr(2.*(meMomenta()[2]*meMomenta()[3])-sqr(Z0->mass())) + sqr(Z0->mass())*sqr(Z0->width()) );
  if ( meMomenta()[0]*meMomenta()[4] < meMomenta()[1]*meMomenta()[4] ) {
    sel.insert(wGamma,0);
    sel.insert(wZ,2);
  } else {
    sel.insert(wGamma,1);
    sel.insert(wZ,3);
  }
  return sel;
}

Selector<const ColourLines *>
MatchboxMEqqbar2llbarg::colourGeometries(tcDiagPtr d) const {
  static const ColourLines cq("1 7, -7 2 -3");
  static const ColourLines cqbar("1 2 7, -3 -7");
  static const ColourLines cbarq("-1 -7, 7 -2 3");
  static const ColourLines cbarqbar("-1 -2 -7, 3 7");
  Selector<const ColourLines *> sel;
  if ( mePartonData()[0]->id() > 0 ) {
    if ( abs(d->id()) % 2 != 0 )
      sel.insert(1.0, &cq);
    else
      sel.insert(1.0, &cqbar);
  } else {
    if ( abs(d->id()) % 2 != 0 )
      sel.insert(1.0, &cbarq);
    else
      sel.insert(1.0, &cbarqbar);
  }
  return sel;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxMEqqbar2llbarg::persistentOutput(PersistentOStream &) const {
}

void MatchboxMEqqbar2llbarg::persistentInput(PersistentIStream &, int) {
}

ClassDescription<MatchboxMEqqbar2llbarg> MatchboxMEqqbar2llbarg::initMatchboxMEqqbar2llbarg;
// Definition of the static class description member.

void MatchboxMEqqbar2llbarg::Init() {

  static ClassDocumentation<MatchboxMEqqbar2llbarg> documentation
    ("MatchboxMEqqbar2llbarg");

}

