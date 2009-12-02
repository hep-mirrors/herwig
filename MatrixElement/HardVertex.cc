// -*- C++ -*-
//
// HardVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardVertex class.
//
// Author: Peter Richardson
//

#include "HardVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Helicity/SpinInfo.h"

using namespace Herwig;
using ThePEG::Helicity::SpinInfo;
using ThePEG::Helicity::tcSpinfoPtr;
  
NoPIOClassDescription<HardVertex> HardVertex::initHardVertex;
  // Definition of the static class description member.
    
void HardVertex::Init() {
  
  static ClassDocumentation<HardVertex> documentation
    ("The HardVertex class implements the vertex for a hard "
     "interaction for the Herwig++ spin correlation algorithm");
  
}
 
// method to get the rho matrix for a given outgoing particle
RhoDMatrix HardVertex::getRhoMatrix(int i,bool) const {
  // get the rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout(outgoing().size()-1);
  for(int ix=0,N=outgoing().size();ix<N;++ix) {
    if(ix<i)      rhoout[ix  ] = 
      dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix();
    else if(ix>i) rhoout[ix-1] = 
      dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix();
  }
  // calculate the spin density matrix
  return _matrixelement.
    calculateRhoMatrix(i,dynamic_ptr_cast<tcSpinfoPtr>(incoming()[0])->DMatrix(),
		       dynamic_ptr_cast<tcSpinfoPtr>(incoming()[1])->DMatrix(),rhoout);
}

// method to get the D matrix for an incoming particle
RhoDMatrix HardVertex::getDMatrix(int i) const {
  // get rho rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout(outgoing().size());
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix)
    rhoout[ix] = dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix();
  // calculate the decay matrix
  return _matrixelement.
    calculateDMatrix(i,dynamic_ptr_cast<tcSpinfoPtr>(incoming()[1])->DMatrix(),rhoout);
}


