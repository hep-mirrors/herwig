// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQED class.
//

#include "ShowerAlphaQED.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


ShowerAlphaQED::~ShowerAlphaQED() {}


ClassDescription<ShowerAlphaQED> ShowerAlphaQED::initShowerAlphaQED;
// Definition of the static class description member.


void ShowerAlphaQED::Init() {

  static ClassDocumentation<ShowerAlphaQED> documentation
    ("This (concrete) class describes the QED alpha running.");

}


double ShowerAlphaQED::value(const Energy2 scale) {
  return scaleFactor()*_alpha;
}
