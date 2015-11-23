// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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

LeptoquarkModelSLQSLQGVertex::LeptoquarkModelSLQSLQGVertex() : _q2last(ZERO),
							       _couplast(0.) {
  orderInGs(1);
  orderInGem(0);
}

void LeptoquarkModelSLQSLQGVertex::doinit() {
  addToList(21,9941551,-9941551);
  addToList(21,9911561,-9911561);
  addToList(21,9921551,-9921551);
  addToList(21,9931561,-9931561);
  addToList(21,9931551,-9931551);
  addToList(21,9931661,-9931661);
  addToList(21,9941561,-9941561);
  addToList(21,9951551,-9951551);
  addToList(21,9951651,-9951651);
  addToList(21,9961551,-9961551);
  addToList(21,9971561,-9971561);
  addToList(21,9981561,-9981561);
  addToList(21,9981551,-9981551);
  addToList(21,9981651,-9981651);
  addToList(21,9991551,-9991551);
  addToList(21,9991561,-9991561);
  addToList(21,9901561,-9901561);
  addToList(21,9901661,-9901661);
  VSSVertex::doinit();
}

void LeptoquarkModelSLQSLQGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theModel;
}

void LeptoquarkModelSLQSLQGVertex::persistentInput(PersistentIStream & is, int) {
   is >> _theModel;
}

ClassDescription<LeptoquarkModelSLQSLQGVertex> LeptoquarkModelSLQSLQGVertex::initLeptoquarkModelSLQSLQGVertex;
// Definition of the static class description member.

void LeptoquarkModelSLQSLQGVertex::Init() {
  static ClassDocumentation<LeptoquarkModelSLQSLQGVertex> documentation
    ("The LeptoquarkModelSLQSLQGVertex class is the implementation of"
     " the LeptoquarkModel scalar LQ-scalar LQ-gluon vertex");
  
}

void LeptoquarkModelSLQSLQGVertex::setCoupling(Energy2 q2,tcPDPtr ,tcPDPtr ,tcPDPtr ) { 
   if(q2 != _q2last || _couplast == 0.) {
     _couplast = strongCoupling(q2);
     _q2last = q2;
   }
  norm(_couplast);
}

