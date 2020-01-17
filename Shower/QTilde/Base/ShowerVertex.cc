// -*- C++ -*-
//
// ShowerVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerVertex class.
//

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/EventRecord/SpinInfo.h"
#include "ShowerVertex.h"

using namespace Herwig;
using namespace Herwig::Helicity;
using namespace ThePEG;

DescribeNoPIOClass<ShowerVertex,HelicityVertex>
describeShowerVertex ("Herwig::ShowerVertex","");

void ShowerVertex::Init() {

  static ClassDocumentation<ShowerVertex> documentation
    ("The ShowerVertex class is the implementation of a "
     "vertex for a shower for the Herwig spin correlation algorithm");

}

// method to get the rho matrix for a given outgoing particle
RhoDMatrix ShowerVertex::getRhoMatrix(int i, bool) const {
  assert(matrixElement_->nOut()==2);
  // calculate the incoming spin density matrix
  RhoDMatrix input=incoming()[0]->rhoMatrix();
  if(convertIn_) input = mapIncoming(input);
  // get the rho matrices for the outgoing particles
  vector<RhoDMatrix> rhoout;
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix) {
    if(int(ix)!=i)
      rhoout.push_back(outgoing()[ix]->DMatrix());
  }
  // calculate the spin density matrix
  return matrixElement_->calculateRhoMatrix(i,input,rhoout);
}

// method to get the D matrix for an incoming particle
RhoDMatrix ShowerVertex::getDMatrix(int) const {
  assert(matrixElement_->nOut()==2);
  // get the decay matrices for the outgoing particles
  vector<RhoDMatrix> Dout;
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix) {
    Dout.push_back(outgoing()[ix]->DMatrix());
  }
  // calculate the spin density matrix 
  RhoDMatrix rho = matrixElement_->calculateDMatrix(Dout);
  // map if needed
  if(convertIn_) {
    RhoDMatrix rhop(rho.iSpin(),false);
    for(int ixb=0;ixb<rho.iSpin();++ixb) {
      for(int iyb=0;iyb<rho.iSpin();++iyb) {
	if(inMatrix_(iyb,ixb)==0.)continue;
	for(int iya=0;iya<rho.iSpin();++iya) {
	  if(rho(iya,iyb)==0.)continue;
	  for(int ixa=0;ixa<rho.iSpin();++ixa) {
	    rhop(ixa,ixb) += rho(iya,iyb)*inMatrix_(ixa,iya)*conj(inMatrix_(ixb,iyb));
	  }
	}
      }
    }
    rhop.normalize();
    rho = rhop;
  }
  // return the answer
  return rho;
}

RhoDMatrix ShowerVertex::mapIncoming(RhoDMatrix rho) const {
  RhoDMatrix output(rho.iSpin());
  for(int ixa=0;ixa<rho.iSpin();++ixa) {
    for(int ixb=0;ixb<rho.iSpin();++ixb) {
      for(int iya=0;iya<rho.iSpin();++iya) {
	for(int iyb=0;iyb<rho.iSpin();++iyb) {
	  output(ixa,ixb) += rho(iya,iyb)*inMatrix_(iya,ixa)*conj(inMatrix_(iyb,ixb));
	}
      }
    }
  }
  output.normalize();
  return output;
}

