// -*- C++ -*-
//
// SMGGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMGGGVertex class.
//

#include "SMGGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SMGGGVertex::SMGGGVertex() : _couplast(0.), _q2last(0.*GeV2) {
  orderInGs(1);
  orderInGem(0);
}

void SMGGGVertex::doinit() {
  // the particles
  addToList(21,21,21);
  VVVVertex::doinit();
}

NoPIOClassDescription<SMGGGVertex> SMGGGVertex::initSMGGGVertex;
// Definition of the static class description member.

void SMGGGVertex::Init() {
 static ClassDocumentation<SMGGGVertex> documentation
    ("The SMGGGVertex class is the implementation"
     " of the Standard Model triple gluon vertex.");
  
}

// couplings for the GGG vertex
void SMGGGVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr, tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = strongCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
}
