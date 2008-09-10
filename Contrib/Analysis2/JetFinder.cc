// -*- C++ -*-
//
// JetFinder.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the JetFinder class.
//

#include "JetFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "JetFinder.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

JetFinder::~JetFinder() {}

void JetFinder::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << ounit(_dCut,GeV) << _yCut << _R << _lastEvent;
}

void JetFinder::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> iunit(_dCut,GeV) >> _yCut >> _R >> _lastEvent;
}

AbstractClassDescription<JetFinder> JetFinder::initJetFinder;
// Definition of the static class description member.

void JetFinder::Init() {

  static ClassDocumentation<JetFinder> documentation
    ("JetFinder is the base class for all jet finders "
     "used by the Herwig++ analysis code.");


  static Parameter<JetFinder,Energy> interfacedCut
    ("dCut",
     "Set the dimensionful resolution parameter",
     &JetFinder::_dCut, GeV, 10.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


  static Parameter<JetFinder,double> interfaceyCut
    ("yCut",
     "Set the dimensionless resolution parameter",
     &JetFinder::_yCut, 0.1, 0.0, 1.0,
     false, false, Interface::limited);


  static Parameter<JetFinder,double> interfaceR
    ("R",
     "Set the optional cone radius.",
     &JetFinder::_R, 0.7, 0.0, 0,
     false, false, Interface::lowerlim);

}

