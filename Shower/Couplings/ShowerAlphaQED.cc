// -*- C++ -*-
//
// ShowerAlphaQED.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQED class.
//

#include "ShowerAlphaQED.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

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

double ShowerAlphaQED::value(const Energy2) const {
  return scaleFactor()*_alpha;
}

double ShowerAlphaQED::overestimateValue() const {
  return scaleFactor()*_alpha;
}

double ShowerAlphaQED::ratio(const Energy2,double ) const {
  return 1.;
}

void ShowerAlphaQED::doinit() throw(InitException) {
  ShowerAlpha::doinit();
  _alpha=generator()->standardModel()->alphaEM();
}
