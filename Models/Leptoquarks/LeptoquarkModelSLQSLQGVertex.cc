// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModelSLQSLQGVertex class.
//

#include "LeptoquarkModelSLQSLQGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


using namespace Herwig;
using namespace ThePEG;

LeptoquarkModelSLQSLQGVertex::LeptoquarkModelSLQSLQGVertex() {
  addToList(21,9911561,-9911561);
}

void LeptoquarkModelSLQSLQGVertex::doinit() {
  VSSVertex::doinit();
  //  _theModel = generator()->standardModel();
  // tcHwLeptoquarkPtr hwLeptoquark=dynamic_ptr_cast<tcHwLeptoquarkPtr>(_theModel);
}

void LeptoquarkModelSLQSLQGVertex::persistentOutput(PersistentOStream & os) const {
  // os << _theModel;
}

void LeptoquarkModelSLQSLQGVertex::persistentInput(PersistentIStream & is, int) {
  // is >> _theModel;
}

ClassDescription<LeptoquarkModelSLQSLQGVertex> LeptoquarkModelSLQSLQGVertex::initLeptoquarkModelSLQSLQGVertex;
// Definition of the static class description member.

void LeptoquarkModelSLQSLQGVertex::Init() {
  static ClassDocumentation<LeptoquarkModelSLQSLQGVertex> documentation
    ("The LeptoquarkModelSLQSLQGVertex class is the implementation of"
     " the LeptoquarkModel scalar LQ-scalar LQ-gluon vertex");
  
}

void LeptoquarkModelSLQSLQGVertex::setCoupling(Energy2 q2,tcPDPtr ,tcPDPtr ,tcPDPtr ) {
  //  norm(sqrt(4.* Constants::pi*strongCoupling(q2)));
  norm(strongCoupling(q2));
}

