// -*- C++ -*-
//
// SSGGSQSQVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGGSQSQVertex class.
//

#include "SSGGSQSQVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGGSQSQVertex::SSGGSQSQVertex() : _q2last(),_couplast(0.) {
  //L-L squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(21,21,ix,-ix);
  }
  //R-R squarks
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(21,21,ix,-ix);
  }
}

NoPIOClassDescription<SSGGSQSQVertex> SSGGSQSQVertex::initSSGGSQSQVertex;
// Definition of the static class description member.

void SSGGSQSQVertex::Init() {

  static ClassDocumentation<SSGGSQSQVertex> documentation
    ("This implements the gluon-gluon-squark-squark vertex.");

}

void SSGGSQSQVertex::setCoupling(Energy2 q2, tcPDPtr, tcPDPtr, tcPDPtr,
				 tcPDPtr) { 
  if(q2 != _q2last || _couplast == 0.) {
    _couplast = sqr(strongCoupling(q2));
    _q2last = q2;
  }
  norm(_couplast);
}

void SSGGSQSQVertex::doinit() {
  orderInGs(2);
  orderInGem(0);
  VVSSVertex::doinit();
}


