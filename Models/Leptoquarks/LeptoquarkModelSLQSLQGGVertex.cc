// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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

LeptoquarkModelSLQSLQGGVertex::LeptoquarkModelSLQSLQGGVertex() : _q2last(ZERO),
								 _couplast(0.) {
  orderInGs(2);
  orderInGem(0);
}

void LeptoquarkModelSLQSLQGGVertex::doinit() {
  addToList(21,21,9941551,-9941551);
  addToList(21,21,9911561,-9911561);
  addToList(21,21,9921551,-9921551);
  addToList(21,21,9931561,-9931561);
  addToList(21,21,9931551,-9931551);
  addToList(21,21,9931661,-9931661);
  addToList(21,21,9941561,-9941561);

  addToList(21,21,9951551,-9951551);
  addToList(21,21,9951651,-9951651);
  addToList(21,21,9961551,-9961551);

  addToList(21,21,9971561,-9971561);
  addToList(21,21,9981561,-9981561);
  addToList(21,21,9981551,-9981551);
  addToList(21,21,9981651,-9981651);

  addToList(21,21,9991551,-9991551);
  addToList(21,21,9991561,-9991561);
  addToList(21,21,9901561,-9901561);
  addToList(21,21,9901661,-9901661);
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

