// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQED class.
//

#include "ShowerAlphaQED.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerAlphaQED.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ShowerAlphaQED::~ShowerAlphaQED() {}

void ShowerAlphaQED::persistentOutput(PersistentOStream & os) const {

}

void ShowerAlphaQED::persistentInput(PersistentIStream & is, int) {

}

ClassDescription<ShowerAlphaQED> ShowerAlphaQED::initShowerAlphaQED;
// Definition of the static class description member.

void ShowerAlphaQED::Init() {

  static ClassDocumentation<ShowerAlphaQED> documentation
    ("This (concrete) class describes the QED alpha running.");

}

double ShowerAlphaQED::value(const Energy2 scale) {
  return scaleFactor()*_alpha;
}

