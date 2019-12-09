// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModelSLQSLQGVertex class.
//

#include "LeptoquarkModelSLQSLQGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

LeptoquarkModelSLQSLQGVertex::LeptoquarkModelSLQSLQGVertex() : q2last_(ZERO),
							       couplast_(0.) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
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

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<LeptoquarkModelSLQSLQGVertex,VSSVertex>
describeHerwigLeptoquarkModelSLQSLQGVertex("Herwig::LeptoquarkModelSLQSLQGVertex", "HwLeptoquarkModel.so");

void LeptoquarkModelSLQSLQGVertex::Init() {
  static ClassDocumentation<LeptoquarkModelSLQSLQGVertex> documentation
    ("The LeptoquarkModelSLQSLQGVertex class is the implementation of"
     " the LeptoquarkModel scalar LQ-scalar LQ-gluon vertex");
  
}

void LeptoquarkModelSLQSLQGVertex::setCoupling(Energy2 q2,tcPDPtr ,tcPDPtr p2,tcPDPtr ) { 
  if(q2 != q2last_ || couplast_ == 0.) {
    couplast_ = strongCoupling(q2);
    q2last_ = q2;
  }
  if(p2->id()<0)
    norm( couplast_);
  else
    norm(-couplast_);
}

