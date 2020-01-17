// -*- C++ -*-
//
// SMGGGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMGGGGVertex class.
//

#include "SMGGGGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SMGGGGVertex::SMGGGGVertex() : _couplast(0.),_q2last() {
  orderInGs(2);
  orderInGem(0);
  colourStructure(ColourStructure::SU3FF);
}

void SMGGGGVertex::doinit() {
  // particles
  addToList(21,21,21,21);
  VVVVVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SMGGGGVertex,VVVVVertex>
describeHerwigSMGGGGVertex("Herwig::SMGGGGVertex", "Herwig.so");

void SMGGGGVertex::Init() {
  static ClassDocumentation<SMGGGGVertex> documentation
    ("The SMGGGGVertex class is the implementation of the"
     " Standard Model quartic gluon coupling");
  
}


// couplings for the GGGG vertex
void SMGGGGVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
			       tcPDPtr,tcPDPtr) {
  // set the order and type
  setType(1);
  setOrder(0,1,2,3);
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = sqr(strongCoupling(q2));
    _q2last=q2;
  }
  norm(_couplast);
}

