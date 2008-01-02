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

namespace Herwig {
using namespace ThePEG;

void SMGGGGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM;
}

void SMGGGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM;
}

ClassDescription<SMGGGGVertex>
SMGGGGVertex::initSMGGGGVertex;
// Definition of the static class description member.

void SMGGGGVertex::Init() {
  static ClassDocumentation<SMGGGGVertex> documentation
    ("The SMGGGGVertex class is the implementation of the"
     " Standard Model quartic gluon coupling");
  
}


// couplings for the GGGG vertex
void SMGGGGVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
    				      tcPDPtr,tcPDPtr)
{
  // set the order and type
  setType(1);setOrder(0,1,2,3);
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaS(q2);
      _couplast = 4.0*Constants::pi*alpha;
      _q2last=q2;
    }
  setNorm(_couplast);
}
}

