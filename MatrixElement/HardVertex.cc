// -*- C++ -*-
//
// HardVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardVertex class.
//
// Author: Peter Richardson
//

#include "ThePEG/EventRecord/SpinInfo.h"
#include "HardVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

NoPIOClassDescription<HardVertex> HardVertex::initHardVertex;
  // Definition of the static class description member.
    
void HardVertex::Init() {
  
  static ClassDocumentation<HardVertex> documentation
    ("The HardVertex class implements the vertex for a hard "
     "interaction for the Herwig spin correlation algorithm");
  
}
 
// method to get the rho matrix for a given outgoing particle
RhoDMatrix HardVertex::getRhoMatrix(int i,bool) const {
  // get the rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout(outgoing().size()-1);
  for(int ix=0,N=outgoing().size();ix<N;++ix) {
    if(ix<i)      rhoout[ix  ] = outgoing()[ix]->DMatrix();
    else if(ix>i) rhoout[ix-1] = outgoing()[ix]->DMatrix();
  }
  // calculate the spin density matrix
  return _matrixelement.
    calculateRhoMatrix(i,
		       incoming()[0]->rhoMatrix(),
		       incoming()[1]->rhoMatrix(),rhoout);
}

// method to get the D matrix for an incoming particle
RhoDMatrix HardVertex::getDMatrix(int i) const {
  // get rho rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout(outgoing().size());
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix)
    rhoout[ix] = outgoing()[ix]->DMatrix();
  unsigned int iother = i==0 ? 1 : 0;
  // calculate the decay matrix
  return _matrixelement.
    calculateDMatrix(i,incoming()[iother]->rhoMatrix(),rhoout);
}
