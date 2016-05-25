// -*- C++ -*-
//
// DipoleEvolutionOrdering.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleEvolutionOrdering class.
//

#include "DipoleEvolutionOrdering.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DipoleEvolutionOrdering::DipoleEvolutionOrdering() 
  : HandlerBase() {}

DipoleEvolutionOrdering::~DipoleEvolutionOrdering() {}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleEvolutionOrdering::persistentOutput(PersistentOStream &) const {
}

void DipoleEvolutionOrdering::persistentInput(PersistentIStream &, int) {
}

AbstractClassDescription<DipoleEvolutionOrdering> DipoleEvolutionOrdering::initDipoleEvolutionOrdering;
// Definition of the static class description member.

void DipoleEvolutionOrdering::Init() {

  static ClassDocumentation<DipoleEvolutionOrdering> documentation
    ("DipoleEvolutionOrdering defines a particular evolution "
     "algortihm for the dipole shower.");

}

