// -*- C++ -*-
//
// ShowerAlphaQED.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQED class.
//

#include "ShowerAlphaQED.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void ShowerAlphaQED::persistentOutput(PersistentOStream & os) const {
  os << _alpha << couplingSource_;
}

void ShowerAlphaQED::persistentInput(PersistentIStream & is, int) {
  is >> _alpha >> couplingSource_;
}

ClassDescription<ShowerAlphaQED> ShowerAlphaQED::initShowerAlphaQED;
// Definition of the static class description member.

void ShowerAlphaQED::Init() {

  static ClassDocumentation<ShowerAlphaQED> documentation
    ("This (concrete) class describes the QED alpha running.");

  static Parameter<ShowerAlphaQED,double> interfaceAlpha
    ("Alpha",
     "The value of alpha_EM",
     &ShowerAlphaQED::_alpha, 1./137., 0., 1.,
     false, false, Interface::limited);


  static Switch<ShowerAlphaQED,unsigned int> interfaceCouplingSource
    ("CouplingSource",
     "Where to get the coupling from",
     &ShowerAlphaQED::couplingSource_, 0, false, false);
  static SwitchOption interfaceCouplingSourceLocal
    (interfaceCouplingSource,
     "Local",
     "Use the local value",
     0);
  static SwitchOption interfaceCouplingSourceThompson
    (interfaceCouplingSource,
     "Thompson",
     "Use the Thompson value from the StamdardModel object",
     1);
  static SwitchOption interfaceCouplingSourceMZ
    (interfaceCouplingSource,
     "MZ",
     "Use the value at MZ from the StandardModel object",
     2);

}

double ShowerAlphaQED::value(const Energy2) const {
  return _alpha;
}

double ShowerAlphaQED::overestimateValue() const {
  return _alpha;
}

double ShowerAlphaQED::ratio(const Energy2,double ) const {
  return 1.;
}

void ShowerAlphaQED::doinit() {
  ShowerAlpha::doinit();
  if(couplingSource_==1)
    _alpha=generator()->standardModel()->alphaEM();
  else if(couplingSource_==2)
    _alpha=generator()->standardModel()->alphaEMMZ();
}
