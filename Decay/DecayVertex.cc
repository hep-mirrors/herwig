// -*- C++ -*-
//
// DecayVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayVertex class.
//
//  Author: Peter Richardson
//

#include <ThePEG/EventRecord/SpinInfo.h>
#include "DecayVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;
using namespace ThePEG;

// Static variable needed for the type description system in ThePEG.
DescribeNoPIOClass<DecayVertex,HelicityVertex>
describeHerwigDecayVertex("Herwig::DecayVertex", "Herwig.so");
    
void DecayVertex::Init() {
  
  static ClassDocumentation<DecayVertex> documentation
    ("The DecayVertex class is the implementation of a "
     "vertex for a decay for the Herwig spin correlation algorithm");
  
}

// method to get the rho matrix for a given outgoing particle
RhoDMatrix DecayVertex::getRhoMatrix(int i,bool recursive) const {
  // get the rho matrix of the decaying particle
  RhoDMatrix input;
  tcSpinPtr inspin = incoming()[0];
  assert(inspin);
  if(recursive&&inspin->productionVertex()&&
     inspin->iSpin()!=PDT::Spin0) {
    input = inspin->productionVertex()->
      getRhoMatrix(inspin->productionLocation(),true);
    inspin->rhoMatrix() = input;
    inspin->needsUpdate();
  }
  else {
    input = inspin->rhoMatrix();
  }
  // get the D matrices for the outgoing particles
  vector<RhoDMatrix> rhoout(outgoing().size()-1);
  for(int ix=0,N=outgoing().size();ix<N;++ix) {
    if(ix<i)      
      rhoout[ix]   = outgoing()[ix]->DMatrix();
    else if(ix>i) 
      rhoout[ix-1] = outgoing()[ix]->DMatrix();
  }
  // calculate the spin density matrix
  return matrixElement_->calculateRhoMatrix(i,input,rhoout);
}

// method to get the D matrix for an incoming particle
RhoDMatrix DecayVertex::getDMatrix(int) const {
  tcSpinPtr inspin = incoming()[0];
  if(inspin->developed()==SpinInfo::Developed) 
    return inspin->DMatrix();
  // get the decay matrices for the outgoing particles
  vector<RhoDMatrix> Dout(outgoing().size());
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix) {
    tcSpinPtr hwspin = outgoing()[ix];
    if(hwspin->developed()!=SpinInfo::Developed) 
      hwspin->develop();
    Dout[ix] = hwspin->DMatrix();
  }
  // calculate the spin density matrix and return the answer
  return matrixElement_->calculateDMatrix(Dout);
}
