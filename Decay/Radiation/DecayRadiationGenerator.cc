// -*- C++ -*-
//
// DecayRadiationGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayRadiationGenerator class.
//

#include "DecayRadiationGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

AbstractNoPIOClassDescription<DecayRadiationGenerator> DecayRadiationGenerator::initDecayRadiationGenerator;
// Definition of the static class description member.

void DecayRadiationGenerator::Init() {

  static ClassDocumentation<DecayRadiationGenerator> documentation
    ("The DecayRadiationGenerator class is the base class for the implementation of"
     "QED radiation in particle decays.");

}

