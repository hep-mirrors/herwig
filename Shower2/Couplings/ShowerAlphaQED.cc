// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQED class.
//

#include "ShowerAlphaQED.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerAlphaQED.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ShowerAlphaQED::~ShowerAlphaQED() {}

void ShowerAlphaQED::persistentOutput(PersistentOStream & os) const {
  os << _alpha;
}

void ShowerAlphaQED::persistentInput(PersistentIStream & is, int) {
  is >> _alpha;
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

double ShowerAlphaQED::overestimateValue() {
  return scaleFactor()*_alpha;
}

double ShowerAlphaQED::ratio(const Energy2 scale) {
  return 1.;
}

void ShowerAlphaQED::doinit() throw(InitException) {
  ShowerAlpha::doinit();
  _alpha=generator()->standardModel()->alphaEM();
}
