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

DecayRadiationGenerator::~DecayRadiationGenerator() {}

void DecayRadiationGenerator::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void DecayRadiationGenerator::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

AbstractClassDescription<DecayRadiationGenerator> DecayRadiationGenerator::initDecayRadiationGenerator;
// Definition of the static class description member.

void DecayRadiationGenerator::Init() {

  static ClassDocumentation<DecayRadiationGenerator> documentation
    ("There is no documentation for the DecayRadiationGenerator class");

}

