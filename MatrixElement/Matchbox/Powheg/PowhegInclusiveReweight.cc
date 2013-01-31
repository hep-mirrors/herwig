// -*- C++ -*-
//
// PowhegInclusiveReweight.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegInclusiveReweight class.
//

#include "PowhegInclusiveReweight.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PowhegInclusiveReweight::PowhegInclusiveReweight()
  : theBornScreening(true) {}

PowhegInclusiveReweight::~PowhegInclusiveReweight() {}

IBPtr PowhegInclusiveReweight::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegInclusiveReweight::fullclone() const {
  return new_ptr(*this);
}

double PowhegInclusiveReweight::evaluate() const {

  if ( projectionDipole()->verbose() )
    generator()->log() << "'" << name() << "' evaluating inclusive reweight\n";

  double sratio;
  double ratio = ME2byDipoles::evaluate(sratio,true);

  if ( bornScreening() ) {
    if ( !projectionDipole()->underlyingBornME()->
	 lastXCombPtr()->willPassCuts() )
      return 0.;
    double born = scaledBorn();
    double screen = scaledBornScreen();
    ratio *= born / ( born + screen );
  }

  if ( projectionDipole()->verbose() )
    generator()->log() << "'" << name() << "' done evaluating inclusive reweight\n";

  return ratio - sratio;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void PowhegInclusiveReweight::persistentOutput(PersistentOStream & os) const {
  os << theBornScreening;
}

void PowhegInclusiveReweight::persistentInput(PersistentIStream & is, int) {
  is >> theBornScreening;
}

void PowhegInclusiveReweight::Init() {

  static ClassDocumentation<PowhegInclusiveReweight> documentation
    ("PowhegInclusiveReweight");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PowhegInclusiveReweight,ME2byDipoles>
describeHerwigPowhegInclusiveReweight("Herwig::PowhegInclusiveReweight", "HwMatchbox.so");
