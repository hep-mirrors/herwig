// -*- C++ -*-
//
// SMFFGVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModelFFGVertex class.
//

#include "SMFFGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using namespace ThePEG;

NoPIOClassDescription<SMFFGVertex> 
SMFFGVertex::initSMFFGVertex;
// Definition of the static class description member.

void SMFFGVertex::Init() {

  static ClassDocumentation<SMFFGVertex> documentation
    ("The SMFFGVertex class is the implementation of"
     "the coupling of the gluon to the Standard Model fermions");
  
}

// coupling for FFG vertex
void SMFFGVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = -strongCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm=abs(aa->id());
  if(iferm>=1 && iferm<=6) {
    left(1.);
    right(1.);
  }
  else
    throw HelicityConsistencyError() << "SMFFGVertex::setCoupling" 
				     << "Unknown particle in gluon vertex" 
				     << Exception::runerror;
}

SMFFGVertex::SMFFGVertex() : _couplast(0.), _q2last(ZERO) {
  // PDG codes for the particles
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,21);
  }

}
  
void SMFFGVertex::doinit() {
  orderInGs(1);
  orderInGem(0);
  FFVVertex::doinit();
}
