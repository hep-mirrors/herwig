// -*- C++ -*-
//
// DecayVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayVertex class.
//
//  Author: Peter Richardson
//

#include "DecayVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Helicity/SpinInfo.h"

using namespace Herwig;
using ThePEG::Helicity::SpinInfo;
using ThePEG::Helicity::tcSpinfoPtr;

using namespace ThePEG;

NoPIOClassDescription<DecayVertex> DecayVertex::initDecayVertex;
  // Definition of the static class description member.
    
void DecayVertex::Init() {
  
  static ClassDocumentation<DecayVertex> documentation
    ("The DecayVertex class is the implementation of a "
     "vertex for a decay for the Herwig++ spin correlation algorithm");
  
}

// method to get the rho matrix for a given outgoing particle
RhoDMatrix DecayVertex::getRhoMatrix(int i) {
  // get the rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout(outgoing().size()-1);
  for(int ix=0,N=outgoing().size();ix<N;++ix) {
    if(ix<i)      rhoout[ix] = 
      dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix();
    else if(ix>i) rhoout[ix-1] = 
      dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix();
  }
  // calculate the spin density matrix
  RhoDMatrix input=dynamic_ptr_cast<tcSpinfoPtr>(incoming()[0])->rhoMatrix();
  return _matrixelement.calculateRhoMatrix(i,input,rhoout);
}

// method to get the D matrix for an incoming particle
RhoDMatrix DecayVertex::getDMatrix(int) {
  // get the decay matrices for the outgoing particles
  vector<RhoDMatrix> Dout(outgoing().size());
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix)
    Dout[ix] = dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix();
  // calculate the spin density matrix and return the answer
  return _matrixelement.calculateDMatrix(Dout);
}
