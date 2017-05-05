// -*- C++ -*-
//
// RSModelVVGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelVVGRVertex class.
//

#include "RSModelVVGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void RSModelVVGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV);
}

void RSModelVVGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV);
}

ClassDescription<RSModelVVGRVertex> RSModelVVGRVertex::initRSModelVVGRVertex;
// Definition of the static class description member.

void RSModelVVGRVertex::Init() {
 static ClassDocumentation<RSModelVVGRVertex> documentation
    ("The RSModelVVGRVertex class is the implementation"
     " of the RSModel vector-vector-graviton vertex");
  
}
  
void RSModelVVGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(UnitRemoval::E * kappa_));
}

RSModelVVGRVertex::RSModelVVGRVertex() : kappa_(ZERO) {
  orderInGem(1);
  orderInGs (0);
}

void RSModelVVGRVertex::doinit() {
  addToList(23,23,39);
  addToList(22,22,39);
  addToList(24,-24,39);
  addToList(21,21,39);
  VVTVertex::doinit();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS)
    throw Exception() << "Must be RSModel in RSModelVVGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwRS->lambda_pi();
}
