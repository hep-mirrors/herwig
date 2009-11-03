// -*- C++ -*-
//
// RSModelVVGRVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
  os << _theModel << ounit(_theKappa,InvGeV);
}

void RSModelVVGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theModel >> iunit(_theKappa,InvGeV);
}

ClassDescription<RSModelVVGRVertex> RSModelVVGRVertex::initRSModelVVGRVertex;
// Definition of the static class description member.

void RSModelVVGRVertex::Init() {
 static ClassDocumentation<RSModelVVGRVertex> documentation
    ("The RSModelVVGRVertex class is the implementation"
     " of the RSModel vector-vector-graviton vertex");
  
}
  
void RSModelVVGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  setNorm(Complex(UnitRemoval::E * _theKappa));
}

RSModelVVGRVertex::RSModelVVGRVertex() {
  vector<long> first,second,third;
  first.push_back(23);second.push_back( 23);third.push_back(39);
  first.push_back(22);second.push_back( 22);third.push_back(39);
  first.push_back(24);second.push_back(-24);third.push_back(39);
  first.push_back(21);second.push_back( 21);third.push_back(39);
  setList(first,second,third);
  _theKappa=InvEnergy();
}

void RSModelVVGRVertex::doinit() {
  VVTVertex::doinit();
  _theModel = generator()->standardModel();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(_theModel);
  if(hwRS){_theKappa=2./hwRS->lambda_pi();}
  else{throw InitException();}
}
