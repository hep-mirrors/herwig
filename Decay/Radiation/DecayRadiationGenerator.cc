// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayRadiationGenerator class.
//

#include "DecayRadiationGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DecayRadiationGenerator.tcc"
#endif

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

