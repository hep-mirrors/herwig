// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardestEmissionGenerator class.
//

#include "HardestEmissionGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

AbstractNoPIOClassDescription<HardestEmissionGenerator> HardestEmissionGenerator::initHardestEmissionGenerator;
// Definition of the static class description member.

void HardestEmissionGenerator::Init() {

  static ClassDocumentation<HardestEmissionGenerator> documentation
    ("The HardestEmissionGenerator class is the base class for the generation"
     "of the hardest emission in the Nason shower approach");

}

