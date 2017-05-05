// -*- C++ -*-
//
// RSModelSSGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelSSGRVertex class.
//

#include "RSModelSSGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

RSModelSSGRVertex::RSModelSSGRVertex() : kappa_(ZERO) {
  orderInGem(1);
  orderInGs (0);
}

void RSModelSSGRVertex::doinit() {
  addToList(25,25,39);
  SSTVertex::doinit();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) 
    throw Exception() << "Must have RSModel in RSModelSSGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwRS->lambda_pi();
}

void RSModelSSGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV);
}

void RSModelSSGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV);
}

ClassDescription<RSModelSSGRVertex> RSModelSSGRVertex::initRSModelSSGRVertex;
// Definition of the static class description member.

void RSModelSSGRVertex::Init() {
  static ClassDocumentation<RSModelSSGRVertex> documentation
    ("The RSModelSSGRVertex class is the implementation of"
     " the RSModel scalar-scalar-graviton vertex");
  
}

void RSModelSSGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
    norm(Complex(kappa_ * UnitRemoval::E));
}
