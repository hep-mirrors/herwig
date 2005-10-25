// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ModelGenerator class.
//

#include "ModelGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ModelGenerator.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ModelGenerator::~ModelGenerator() {}

void ModelGenerator::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void ModelGenerator::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<ModelGenerator> ModelGenerator::initModelGenerator;
// Definition of the static class description member.

void ModelGenerator::Init() {

  static ClassDocumentation<ModelGenerator> documentation
    ("There is no documentation for the ModelGenerator class");

}

