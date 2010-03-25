// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGGVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModelSLQSLQGGVertex class.
//

#include "LeptoquarkModelSLQSLQGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


using namespace Herwig;
using namespace ThePEG;

LeptoquarkModelSLQSLQGGVertex::LeptoquarkModelSLQSLQGGVertex() :  _couplast(0.),_q2last(ZERO) {
  addToList(21,21,9911561,-9911561);
  addToList(21,21,9911551,-9911551);
  addToList(21,21,9921561,-9921561);
  addToList(21,21,9921551,-9921551);
  addToList(21,21,9921661,-9921661);
}

void LeptoquarkModelSLQSLQGGVertex::doinit() {
  orderInGs(2);
  orderInGem(0);
  VVSSVertex::doinit();
  //  _theModel = generator()->standardModel();
  // tcHwLeptoquarkPtr hwLeptoquark=dynamic_ptr_cast<tcHwLeptoquarkPtr>(_theModel);
}

void LeptoquarkModelSLQSLQGGVertex::persistentOutput(PersistentOStream & os) const {
   os << _theModel;
}

void LeptoquarkModelSLQSLQGGVertex::persistentInput(PersistentIStream & is, int) {
   is >> _theModel;
}

ClassDescription<LeptoquarkModelSLQSLQGGVertex> LeptoquarkModelSLQSLQGGVertex::initLeptoquarkModelSLQSLQGGVertex;
// Definition of the static class description member.

void LeptoquarkModelSLQSLQGGVertex::Init() {
  static ClassDocumentation<LeptoquarkModelSLQSLQGGVertex> documentation
    ("The LeptoquarkModelSLQSLQGGVertex class is the implementation of"
     " the LeptoquarkModel scalar LQ-scalar LQ-gluon-gluon vertex");
  
}

void LeptoquarkModelSLQSLQGGVertex::setCoupling(Energy2 q2,tcPDPtr ,tcPDPtr ,tcPDPtr, tcPDPtr ) {
  if(q2 != _q2last || _couplast == 0.) {
    _couplast = sqr(strongCoupling(q2));  
    _q2last = q2;
  }
  norm(_couplast);
}

