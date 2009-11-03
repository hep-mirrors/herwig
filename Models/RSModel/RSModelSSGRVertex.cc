// -*- C++ -*-
//
// RSModelSSGRVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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

RSModelSSGRVertex::RSModelSSGRVertex() {
  addToList(25,25,39);
  _theKappa=InvEnergy();
}

void RSModelSSGRVertex::doinit() {
  SSTVertex::doinit();
  _theModel = generator()->standardModel();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(_theModel);
  if(hwRS){_theKappa=2./hwRS->lambda_pi();}
  else{throw InitException();}
}

void RSModelSSGRVertex::persistentOutput(PersistentOStream & os) const {
  os << _theModel << ounit(_theKappa,InvGeV);
}

void RSModelSSGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> iunit(_theKappa,InvGeV);
}

ClassDescription<RSModelSSGRVertex> RSModelSSGRVertex::initRSModelSSGRVertex;
// Definition of the static class description member.

void RSModelSSGRVertex::Init() {
  static ClassDocumentation<RSModelSSGRVertex> documentation
    ("The RSModelSSGRVertex class is the implementation of"
     " the RSModel scalar-scalar-graviton vertex");
  
}

void RSModelSSGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
    norm(Complex(_theKappa * UnitRemoval::E));
}
