// -*- C++ -*-
//
// GeneralDecayMatrixElement.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralDecayMatrixElement class.
//
// Author: Peter Richardson
//

#include "GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG;
    
// calculate the decay matrix for this decay
RhoDMatrix GeneralDecayMatrixElement::calculateDMatrix(const vector<RhoDMatrix> & rhoout) const {
  // vectors for the helicities
  vector<int> ihel1(outspin().size()+1),ihel2(outspin().size()+1);
  // rhomatrix to be returned
  RhoDMatrix output(inspin(), false);
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  unsigned int ix,iy,iz;
  int ixa,iya;
  for(ix=0;ix<matrixElement_.size();++ix) {
    // map the vector index to the helicities
    for(ixa=outspin().size();ixa>=0;--ixa) 
      ihel1[ixa]=(ix%constants_[ixa])/constants_[ixa+1];
    // inner loop
    for(iy=0;iy<matrixElement_.size();++iy) {
      // map the vector index to the helicities	   
      for(iya=outspin().size();iya>=0;--iya)
	ihel2[iya]=(iy%constants_[iya])/constants_[iya+1];
      // matrix element piece
      temp=matrixElement_[ix]*conj(matrixElement_[iy]);
      // spin density matrices for the outgoing particles
      for(iz=0;iz<outspin().size();++iz)
	temp*=rhoout[iz](ihel1[iz+1],ihel2[iz+1]);
      output(ihel1[0],ihel2[0])+=temp;
    }
  }
  // ensure unit trace for the matrix
  output.normalize();
  // return the answer
  return output;
}

// calculate the rho matrix for a given outgoing particle
RhoDMatrix GeneralDecayMatrixElement::
calculateRhoMatrix(int id,const RhoDMatrix & rhoin,
		   const vector<RhoDMatrix> & rhoout) const {
  // vectors for the helicities
  vector<int> ihel1(outspin().size()+1),ihel2(outspin().size()+1);
  // rhomatrix to be returned
  RhoDMatrix output(outspin()[id], false);
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  unsigned int ix,iy,iz;
  int ixa,iya;
  for(ix=0;ix<matrixElement_.size();++ix) {
    // map the vector index to the helicities
    for(ixa=outspin().size();ixa>=0;--ixa)
      ihel1[ixa]=(ix%constants_[ixa])/constants_[ixa+1];
    // inner loop
    for(iy=0;iy<matrixElement_.size();++iy) {
      // map the vector index to the helicities	   
      for(iya=outspin().size();iya>=0;--iya)
	ihel2[iya]=(iy%constants_[iya])/constants_[iya+1];
      // matrix element piece
      temp=matrixElement_[ix]*conj(matrixElement_[iy]);
      // spin denisty matrix for the incoming particle
      temp *= rhoin(ihel1[0],ihel2[0]);
      // spin density matrix for the outgoing particles
      for(iz=0;iz<outspin().size()-1;++iz) {
	if(int(iz)<id) temp*=rhoout[iz](ihel1[iz+1],ihel2[iz+1]);
	else           temp*=rhoout[iz](ihel1[iz+2],ihel2[iz+2]);
      }
      // add to the rho matrix
      output(ihel1[id+1],ihel2[id+1])+=temp;
    }
  }
  // return the answer
  output.normalize();
  return output;
}

// contract the matrix element with the rho matrix of the incoming particle
Complex GeneralDecayMatrixElement::contract(const RhoDMatrix & in) const {
  unsigned int ispin(abs(int(inspin())));
  Complex me=0.;
  for(unsigned int ix=0;ix<constants_[1];++ix) {
    for(unsigned int inhel1=0;inhel1<ispin;++inhel1) {
      for(unsigned int inhel2=0;inhel2<ispin;++inhel2) {
	// compute the term
	me+=matrixElement_[inhel1*constants_[1]+ix]*
	  conj(matrixElement_[inhel2*constants_[1]+ix])*in(inhel1,inhel2);
      }
    }
  }
  return me;
}

// contract the matrix element with the rho matrix of the incoming particle
Complex GeneralDecayMatrixElement::contract(const GeneralDecayMatrixElement & con, 
					    const RhoDMatrix & in) {
  unsigned int ispin(abs(int(inspin())));
  Complex me=0.;
  unsigned int ix,inhel1,inhel2;
  for(ix=0;ix<constants_[1];++ix) {
    for(inhel1=0;inhel1<ispin;++inhel1) {
      for(inhel2=0;inhel2<ispin;++inhel2) {
	// compute the term
	me+=matrixElement_[inhel1*constants_[1]+ix]*
	  conj(con.matrixElement_[inhel2*constants_[1]+ix])*in(inhel1,inhel2);
      }
    }
  }
  return me;
}
