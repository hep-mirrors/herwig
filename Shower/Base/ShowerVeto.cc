// -*- C++ -*-
//
// ShowerVeto.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
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
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"

using namespace Herwig;

void ShowerVeto::persistentOutput(PersistentOStream & os) const {
  os << oenum(_vetoType);
}

void ShowerVeto::persistentInput(PersistentIStream & is, int) {
  is >> ienum(_vetoType);
}

AbstractClassDescription<ShowerVeto> ShowerVeto::initShowerVeto;
// Definition of the static class description member.

void ShowerVeto::Init() {

  static ClassDocumentation<ShowerVeto> documentation
    ("ShowerVeto is the base class for all vetoes on showering other than "
     "from matrix element corrections.");

}

