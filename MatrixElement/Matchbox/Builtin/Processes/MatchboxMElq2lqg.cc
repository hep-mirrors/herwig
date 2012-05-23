// -*- C++ -*-
//
// MatchboxMElq2lqg.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMElq2lqg class.
//

#include "MatchboxMElq2lqg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxMElq2lqg::MatchboxMElq2lqg() 
  : MatchboxMElP2lJetJet() {}

MatchboxMElq2lqg::~MatchboxMElq2lqg() {}

IBPtr MatchboxMElq2lqg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMElq2lqg::fullclone() const {
  return new_ptr(*this);
}

double MatchboxMElq2lqg::me2() const {

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

  Lorentz5Momentum pg = meMomenta()[4];

  prepare(pl,plbar,pq,pqbar,pg,lastSHat(),mePartonData()[0],mePartonData()[1]);

  lastME2(MatchboxMEllbarqqbarg::evaluateME2()*me2Norm());
  cacheME2(lastME2());

  logME2();

  return lastME2();

}

void MatchboxMElq2lqg::getDiagrams() const {

  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);

  for ( PDVector::const_iterator l = theLeptonFlavours.begin();
	l != theLeptonFlavours.end(); ++l )
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ) {
      add(new_ptr((Tree2toNDiagram(4), *l, gamma, (**q).CC(), *q, 1, *l, 2, *q, 3, g, -1)));
      add(new_ptr((Tree2toNDiagram(3), *l, gamma, *q, 1, *l, 2, *q, 5, *q, 5, g, -2)));
      add(new_ptr((Tree2toNDiagram(4), *l, Z0, (**q).CC(), *q, 1, *l, 2, *q, 3, g, -3)));
      add(new_ptr((Tree2toNDiagram(3), *l, Z0, *q, 1, *l, 2, *q, 5, *q, 5, g, -4)));

    }

}

Selector<const ColourLines *>
MatchboxMElq2lqg::colourGeometries(tcDiagPtr d) const {
  static const ColourLines ci("4 7, -7 -3 6");
  static const ColourLines cbari("-4 -7, 7 3 -6");
  static const ColourLines cf("3 5 7, -7 6");
  static const ColourLines cbarf("-3 -5 -7, 7 -6");
  Selector<const ColourLines *> sel;
  if ( mePartonData()[1]->id() > 0 ) {
    if ( abs(d->id()) % 2 != 0 )
      sel.insert(1.0, &ci);
    else
      sel.insert(1.0, &cf);
  } else {
    if ( abs(d->id()) % 2 != 0 )
      sel.insert(1.0, &cbari);
    else
      sel.insert(1.0, &cbarf);
  }
  return sel;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxMElq2lqg::persistentOutput(PersistentOStream &) const {
}

void MatchboxMElq2lqg::persistentInput(PersistentIStream &, int) {
}

ClassDescription<MatchboxMElq2lqg> MatchboxMElq2lqg::initMatchboxMElq2lqg;
// Definition of the static class description member.

void MatchboxMElq2lqg::Init() {

  static ClassDocumentation<MatchboxMElq2lqg> documentation
    ("MatchboxMElq2lqg");

}

