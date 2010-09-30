// -*- C++ -*-
//
// ShowerVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerVertex class.
//

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/SpinInfo.h"
#include "ShowerVertex.h"

using namespace Herwig;
using namespace Herwig::Helicity;
using namespace ThePEG;

NoPIOClassDescription<ShowerVertex> ShowerVertex::initShowerVertex;
// Definition of the static class description member.

void ShowerVertex::Init() {

  static ClassDocumentation<ShowerVertex> documentation
    ("The ShowerVertex class is the implementation of a "
     "vertex for a shower for the Herwig++ spin correlation algorithm");

}

// method to get the rho matrix for a given outgoing particle
RhoDMatrix ShowerVertex::getRhoMatrix(int i, bool) const {
  // get the rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout;
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix) {
    if(int(ix)!=i)
      rhoout.push_back(outgoing()[ix]->DMatrix());
  }
  // calculate the spin density matrix
  RhoDMatrix input=incoming()[0]->rhoMatrix();
  RhoDMatrix temp=_matrixelement.calculateRhoMatrix(i,input,rhoout);
  return temp;
}

// method to get the D matrix for an incoming particle
RhoDMatrix ShowerVertex::getDMatrix(int) const {
  // get the decay matrices for the outgoing particles
  vector<RhoDMatrix> Dout;
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix) {
    Dout.push_back(outgoing()[ix]->DMatrix());
  }
  // calculate the spin density matrix and return the answer
  RhoDMatrix temp = _matrixelement.calculateDMatrix(Dout);
  return temp;
}
