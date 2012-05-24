// -*- C++ -*-
//
// MatchboxMElg2lqqbar.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMElg2lqqbar class.
//

#include "MatchboxMElg2lqqbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxMElg2lqqbar::MatchboxMElg2lqqbar() 
  : MatchboxMElP2lJetJet() {}

MatchboxMElg2lqqbar::~MatchboxMElg2lqqbar() {}

IBPtr MatchboxMElg2lqqbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMElg2lqqbar::fullclone() const {
  return new_ptr(*this);
}

double MatchboxMElg2lqqbar::me2() const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::me2();

  double res;
  if ( !calculateME2(res) )
    return res;

  Lorentz5Momentum pq = 
    mePartonData()[3]->id() > 0 ? 
    meMomenta()[3] : 
    meMomenta()[4];

  Lorentz5Momentum pqbar = 
    mePartonData()[3]->id() < 0 ? 
    meMomenta()[3] : 
    meMomenta()[4];

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

  Lorentz5Momentum pg = -meMomenta()[1];

  prepare(pl,plbar,pq,pqbar,pg,lastSHat(),mePartonData()[0],mePartonData()[3]);

  lastME2(-MatchboxMEllbarqqbarg::evaluateME2()*me2Norm());
  cacheME2(lastME2());

  logME2();

  return lastME2();

}

void MatchboxMElg2lqqbar::getDiagrams() const {

  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);

  for ( PDVector::const_iterator l = theLeptonFlavours.begin();
	l != theLeptonFlavours.end(); ++l )
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ) {
      add(new_ptr((Tree2toNDiagram(4), *l, gamma, (**q).CC(), g, 1, *l, 2, *q, 3, (**q).CC(), -1)));
      add(new_ptr((Tree2toNDiagram(4), *l, gamma, *q, g, 1, *l, 3, *q, 2, (**q).CC(), -2)));
      add(new_ptr((Tree2toNDiagram(4), *l, Z0, (**q).CC(), g, 1, *l, 2, *q, 3, (**q).CC(), -1)));
      add(new_ptr((Tree2toNDiagram(4), *l, Z0, *q, g, 1, *l, 3, *q, 2, (**q).CC(), -2)));
    }

}

Selector<const ColourLines *>
MatchboxMElg2lqqbar::colourGeometries(tcDiagPtr d) const {
  static const ColourLines ct("4 -3 6, -4 -7");
  static const ColourLines cbart("-4 3 -6, 4 7");
  static const ColourLines cu("-4 3 -7, 4 6");
  static const ColourLines cbaru("4 -3 7, -4 -6");
  Selector<const ColourLines *> sel;
  if ( mePartonData()[3]->id() > 0 ) {
    if ( abs(d->id()) % 2 != 0 )
      sel.insert(1.0, &ct);
    else
      sel.insert(1.0, &cu);
  } else {
    if ( abs(d->id()) % 2 != 0 )
      sel.insert(1.0, &cbart);
    else
      sel.insert(1.0, &cbaru);
  }
  return sel;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxMElg2lqqbar::persistentOutput(PersistentOStream &) const {
}

void MatchboxMElg2lqqbar::persistentInput(PersistentIStream &, int) {
}

ClassDescription<MatchboxMElg2lqqbar> MatchboxMElg2lqqbar::initMatchboxMElg2lqqbar;
// Definition of the static class description member.

void MatchboxMElg2lqqbar::Init() {

  static ClassDocumentation<MatchboxMElg2lqqbar> documentation
    ("MatchboxMElg2lqqbar");

}

