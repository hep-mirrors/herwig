// -*- C++ -*-
//
// SMGGGGVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMGGGGVertex class.
//

#include "SMGGGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SMGGGGVertex::SMGGGGVertex() : _couplast(0.),_q2last() {
  // particles
  vector<long> first,second,third,fourth;
  first.push_back(21);
  second.push_back(21);
  third.push_back(21);
  fourth.push_back(21);
  setList(first,second,third,fourth);
}

void SMGGGGVertex::doinit() {
  orderInGs(2);
  orderInGem(0);
  VVVVVertex::doinit();
}

NoPIOClassDescription<SMGGGGVertex>
SMGGGGVertex::initSMGGGGVertex;
// Definition of the static class description member.

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
  setNorm(_couplast);
}

