// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DefaultEmissionGenerator class.
//

#include "DefaultEmissionGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void DefaultEmissionGenerator::persistentOutput(PersistentOStream & os) const {
}

void DefaultEmissionGenerator::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<DefaultEmissionGenerator> DefaultEmissionGenerator::initDefaultEmissionGenerator;
// Definition of the static class description member.

void DefaultEmissionGenerator::Init() {

  static ClassDocumentation<DefaultEmissionGenerator> documentation
    ("The DefaultEmissionGenerator class uses the shower approach to generate"
     " the hardest emission.");

}

void DefaultEmissionGenerator::generateHardest(ShowerTreePtr) {
}

bool DefaultEmissionGenerator::canHandle(ShowerTreePtr) {
  return true;
}
