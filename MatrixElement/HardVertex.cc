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
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HardVertex.tcc"
#endif

namespace Herwig {

using ThePEG::Helicity::SpinInfo;
using ThePEG::Helicity::tcSpinfoPtr;

using namespace ThePEG;
  
NoPIOClassDescription<HardVertex> HardVertex::initHardVertex;
  // Definition of the static class description member.
    
void HardVertex::Init() {
  
  static ClassDocumentation<HardVertex> documentation
    ("The HardVertex class implements the vertex for a hard "
     "interaction for the Herwig++ spin correlation algorithm");
  
}
 
// method to get the rho matrix for a given outgoing particle
RhoDMatrix HardVertex::getRhoMatrix(int i)
{
  // get the rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout;
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix)
    {
      if(int(ix)!=i)
	{rhoout.push_back(dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix());}
    }
  // calculate the spin density matrix
  RhoDMatrix temp=_matrixelement.calculateRhoMatrix(i,dynamic_ptr_cast<tcSpinfoPtr>(incoming()[0])->DMatrix(),dynamic_ptr_cast<tcSpinfoPtr>(incoming()[1])->DMatrix(),rhoout);
  return temp;
}

// method to get the D matrix for an incoming particle
RhoDMatrix HardVertex::getDMatrix(int i)
{
  // get rho rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout;
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix)
    {rhoout.push_back(dynamic_ptr_cast<tcSpinfoPtr>(outgoing()[ix])->DMatrix());}
  // calculate the decay matrix
  int j=0;if(i==0){j=1;}
  RhoDMatrix temp=_matrixelement.calculateDMatrix(i,dynamic_ptr_cast<tcSpinfoPtr>(incoming()[1])->DMatrix(),rhoout);
  return temp;
}
}

