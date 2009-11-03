// -*- C++ -*-
//
// RSModelFFGRVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFGRVertex class.
//

#include "RSModelFFGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void RSModelFFGRVertex::persistentOutput(PersistentOStream & os) const {
  os << _theModel << ounit(_theKappa,InvGeV);
}

void RSModelFFGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> iunit(_theKappa,InvGeV);
}

ClassDescription<RSModelFFGRVertex> RSModelFFGRVertex::initRSModelFFGRVertex;
// Definition of the static class description member.

void RSModelFFGRVertex::Init() {
  static ClassDocumentation<RSModelFFGRVertex> documentation
    ("The RSModelFFGRVertex class is the RSModel calculation"
     " of the fermion-antifermion-graviton vertex");
  
}
  
void RSModelFFGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(_theKappa * UnitRemoval::E));
}

RSModelFFGRVertex::RSModelFFGRVertex() {
  // PDG codes for the particles
  // the quarks
  for (int ix=1;ix<7;++ix) addToList(-ix,ix,39);
  // the leptons
  for (int ix=11;ix<17;++ix) addToList(-ix,ix,39);
  _theKappa=InvEnergy();
}

void RSModelFFGRVertex::doinit() {
  FFTVertex::doinit();
  _theModel = generator()->standardModel();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(_theModel);
  if(hwRS){_theKappa=2./hwRS->lambda_pi();}
  else{throw InitException();}
}
