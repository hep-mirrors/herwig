// -*- C++ -*-
//
// SMFFPVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModelFFPVertex class.
//

#include "SMFFPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void SMFFPVertex::persistentOutput(PersistentOStream & os) const {
  os << _charge <<  _theSM;
}

void SMFFPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _charge >> _theSM;
}

ClassDescription<SMFFPVertex> 
SMFFPVertex::initSMFFPVertex;
// Definition of the static class description member.

void SMFFPVertex::Init() {
 static ClassDocumentation<SMFFPVertex> documentation
    ("The SMFFPVertex class is the implementation of"
     "the coupling of the photon to the Standard Model fermions");
}

// coupling for FFP vertex
void SMFFPVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,tcPDPtr)
{
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaEM(q2);
      _couplast = -sqrt(4.0*3.14159265*alpha);
      _q2last=q2;
    }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16))
    {
      setLeft(_charge[iferm]);
      setRight(_charge[iferm]);
    }
  else
    {
      throw HelicityConsistencyError() << "SMGFFPVertex::setCoupling "
				       << "Unknown particle in photon vertex" 
				       << Exception::warning;
      setLeft(0.);setRight(0.);
    }
}
  
}


