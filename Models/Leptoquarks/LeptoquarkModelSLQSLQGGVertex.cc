// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModelSLQSLQGGVertex class.
//

#include "LeptoquarkModelSLQSLQGGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"


using namespace Herwig;
using namespace ThePEG;

LeptoquarkModelSLQSLQGGVertex::LeptoquarkModelSLQSLQGGVertex() : q2last_(ZERO),
								 couplast_(0.) {
  orderInGs(2);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TTFUNDS);
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
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<LeptoquarkModelSLQSLQGGVertex,VVSSVertex>
describeHerwigLeptoquarkModelSLQSLQGGVertex("Herwig::LeptoquarkModelSLQSLQGGVertex", "HwLeptoquarkModel.so");

void LeptoquarkModelSLQSLQGGVertex::Init() {
  static ClassDocumentation<LeptoquarkModelSLQSLQGGVertex> documentation
    ("The LeptoquarkModelSLQSLQGGVertex class is the implementation of"
     " the LeptoquarkModel scalar LQ-scalar LQ-gluon-gluon vertex");
  
}

void LeptoquarkModelSLQSLQGGVertex::setCoupling(Energy2 q2,tcPDPtr ,tcPDPtr ,tcPDPtr, tcPDPtr ) {
  if(q2 != q2last_ || couplast_ == 0.) {
    couplast_ = sqr(strongCoupling(q2));  
    q2last_ = q2;
  }
  norm(couplast_);
}
