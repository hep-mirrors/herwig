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

#include "ThePEG/EventRecord/SpinInfo.h"
#include "HardVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"

using namespace Herwig;

NoPIOClassDescription<HardVertex> HardVertex::initHardVertex;
  // Definition of the static class description member.
    
void HardVertex::Init() {
  
  static ClassDocumentation<HardVertex> documentation
    ("The HardVertex class implements the vertex for a hard "
     "interaction for the Herwig++ spin correlation algorithm");
  
}
 
// method to get the rho matrix for a given outgoing particle
RhoDMatrix HardVertex::getRhoMatrix(int i,bool) const {
  bool debug_basis_states(false);
  if(debug_basis_states) {
    // Debug basis states:
    cout << "HardVertex:" << endl;
    for(unsigned int ix=0;ix<incoming().size();++ix)
      cout << "in " << ix << " = " 
	   << incoming()[ix]->currentMomentum()/GeV << "   ";
    cout << endl;
    for(unsigned int ix=0;ix<outgoing().size();++ix)
      cout << "out " << ix <<" = " 
	   << outgoing()[ix]->currentMomentum()/GeV << "   ";
    cout << endl;
    tcFermionSpinPtr in0=dynamic_ptr_cast<tcFermionSpinPtr>(incoming()[0]);
    tcFermionSpinPtr in1=dynamic_ptr_cast<tcFermionSpinPtr>(incoming()[1]);
    tcVectorSpinPtr  out0=dynamic_ptr_cast<tcVectorSpinPtr>(outgoing()[0]);
    tcVectorSpinPtr  out1=dynamic_ptr_cast<tcVectorSpinPtr>(outgoing()[1]);
    cout << "in 0  : p.vtx " << in0->productionVertex() << "  ";
    cout << " hel -: ";
    for(unsigned int ix=0;ix<4;ix++) 
      cout << in0->getProductionBasisState(0)[ix]/UnitRemoval::SqrtE << " ";
    cout << " hel +: ";
    for(unsigned int ix=0;ix<4;ix++) 
      cout << in0->getProductionBasisState(1)[ix]/UnitRemoval::SqrtE << " ";
    cout << endl;
    cout << "in 1  : p.vtx " << in1->productionVertex() << "  ";
    cout << " hel -: ";
    for(unsigned int ix=0;ix<4;ix++) 
      cout << in1->getProductionBasisState(0)[ix]/UnitRemoval::SqrtE << " ";
    cout << " hel +: ";
    for(unsigned int ix=0;ix<4;ix++) 
      cout << in1->getProductionBasisState(1)[ix]/UnitRemoval::SqrtE << " ";
    cout << endl;
    cout << "out 0 : p.vtx " << out0->productionVertex() << "  ";
    cout << " hel -: ";
    cout << out0->getProductionBasisState(0).t() << " "
	 << out0->getProductionBasisState(0).x() << " "	   
	 << out0->getProductionBasisState(0).y() << " "	   
	 << out0->getProductionBasisState(0).z() << " ";
    cout << " hel +: ";
    cout << out0->getProductionBasisState(2).t() << ", "
	 << out0->getProductionBasisState(2).x() << ", "	   
	 << out0->getProductionBasisState(2).y() << ", "	   
	 << out0->getProductionBasisState(2).z() << endl;
    cout << "out 1 : p.vtx " << out1->productionVertex() << "  ";
    cout << " hel -: ";
    cout << out1->getProductionBasisState(0).t() << " "
	 << out1->getProductionBasisState(0).x() << " "	   
	 << out1->getProductionBasisState(0).y() << " "	   
	 << out1->getProductionBasisState(0).z() << " ";
    cout << " hel +: ";
    cout << out1->getProductionBasisState(2).t() << ", "
	 << out1->getProductionBasisState(2).x() << ", "	   
	 << out1->getProductionBasisState(2).y() << ", "	   
	 << out1->getProductionBasisState(2).z() << endl;
    if(outgoing().size()>2) {
      tcVectorSpinPtr out2=dynamic_ptr_cast<tcVectorSpinPtr>(outgoing()[2]);
      cout << "out 2 : p.vtx " << out2->productionVertex() << "  ";
      cout << " hel -: ";
      cout << out2->getProductionBasisState(0).t() << " "
	   << out2->getProductionBasisState(0).x() << " "	   
	   << out2->getProductionBasisState(0).y() << " "	   
	   << out2->getProductionBasisState(0).z() << " ";
      cout << " hel +: ";
      cout << out2->getProductionBasisState(2).t() << ", "
	   << out2->getProductionBasisState(2).x() << ", "	   
	   << out2->getProductionBasisState(2).y() << ", "	   
	   << out2->getProductionBasisState(2).z() << endl;
    }
  }
  // get the rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout(outgoing().size()-1);
  for(int ix=0,N=outgoing().size();ix<N;++ix) {
    if(ix<i)      rhoout[ix  ] = 
		    outgoing()[ix]->DMatrix();
    else if(ix>i) rhoout[ix-1] = 
		    outgoing()[ix]->DMatrix();
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
  // calculate the decay matrix
  return _matrixelement.
    calculateDMatrix(i,incoming()[1]->rhoMatrix(),rhoout);
}
