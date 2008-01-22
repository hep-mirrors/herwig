// -*- C++ -*-
//
// SMGGGVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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

SMGGGVertex::SMGGGVertex() : _couplast(0.), _q2last() {
  // the particles
  vector<int> first,second,third;
  first.push_back(21);
  second.push_back(21);
  third.push_back(21);
  setList(first,second,third);
}

void SMGGGVertex::doinit() throw(InitException) {
  orderInGs(1);
  orderInGem(0);
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
  if(q2!=_q2last) {
    _couplast = strongCoupling(q2);
    _q2last=q2;
  }
  setNorm(_couplast);
}
