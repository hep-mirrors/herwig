// -*- C++ -*-
//
// ShowerVeto.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerVeto class.
//

#include "ShowerVeto.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerVeto.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ShowerVeto::~ShowerVeto() {}

void ShowerVeto::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _vetoType;
}

void ShowerVeto::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _vetoType;
}

AbstractClassDescription<ShowerVeto> ShowerVeto::initShowerVeto;
// Definition of the static class description member.

void ShowerVeto::Init() {

  static ClassDocumentation<ShowerVeto> documentation
    ("ShowerVeto is the base class for all vetoes on showering other than "
     "from matrix element corrections.");

}

