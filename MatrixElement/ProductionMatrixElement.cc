// -*- C++ -*-
//
// ProductionMatrixElement.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ProductionMatrixElement class.
//
// Author: Peter Richardson
//
#include "ProductionMatrixElement.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;
  
// calculate a decay matrix for one of the incoming particles
RhoDMatrix ProductionMatrixElement::
calculateDMatrix(int id, const RhoDMatrix & rhoin,
		 const vector<RhoDMatrix> & rhoout) const {
  // vectors for the helicities
  vector<unsigned int> ihel1(_outspin.size()+2),ihel2(_outspin.size()+2);
  // rhomatrix to be returned
  RhoDMatrix output(_inspin[id], false);
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  unsigned int ix,iy;
  int ixa,iya;
  for(ix=0;ix<_matrixelement.size();++ix)
    {
      // map the vector index to the helicities
      for(ixa=_outspin.size()+1;ixa>=0;--ixa)
	{ihel1[ixa]=(ix%_constants[ixa])/_constants[ixa+1];}
      // inner loop
      for(iy=0;iy<_matrixelement.size();++iy)
	{
	  // map the vector index to the helicities	   
	  for(iya=_outspin.size()+1;iya>=0;--iya)
	    {ihel2[iya]=(iy%_constants[iya])/_constants[iya+1];}
	  // matrix element piece
	  temp=_matrixelement[ix]*conj(_matrixelement[iy]);
	  // spin density matrices for the outgoing particles
	  for(unsigned int iz=0;iz<_outspin.size();++iz)
	    {temp*=rhoout[iz](ihel1[iz+2],ihel2[iz+2]);}
	  // construct the spin density matrix
	  if(id==0)
	    {
	      temp*=rhoin(ihel1[1],ihel2[1]);
	      output(ihel1[0],ihel2[0])+=temp;
	    }
	  else
	    {
	      temp*=rhoin(ihel1[0],ihel2[0]);
	      output(ihel1[1],ihel2[1])+=temp;
	    }
	}
    }
  // return the answer
  return output;
}

// calculate the rho matrix for a given outgoing particle
RhoDMatrix ProductionMatrixElement::
calculateRhoMatrix(int id,const RhoDMatrix & rhoin0,
		   const RhoDMatrix & rhoin1,
		   const vector<RhoDMatrix> & rhoout) const {
  unsigned int ix,iy;
  int ixa,iya;
  // vectors for the helicities
  vector<unsigned int> ihel1(_outspin.size()+2),ihel2(_outspin.size()+2);
  // rhomatrix to be returned
  RhoDMatrix output(_outspin[id], false);
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  for(ix=0;ix<_matrixelement.size();++ix) {
    // map the vector index to the helicities
    for(ixa=_outspin.size()+1;ixa>=0;--ixa)
      ihel1[ixa]=(ix%_constants[ixa])/_constants[ixa+1];
    // inner loop
    for(iy=0;iy<_matrixelement.size();++iy) {
      // map the vector index to the helicities	   
      for(iya=_outspin.size()+1;iya>=0;--iya)
	ihel2[iya] = (iy%_constants[iya])/_constants[iya+1];
      // matrix element piece
      temp = _matrixelement[ix]*conj(_matrixelement[iy]);
      // spin denisty matrix for the incoming particles
      temp *= rhoin0(ihel1[0],ihel2[0]);
      temp *= rhoin1(ihel1[1],ihel2[1]);
      // spin density matrix for the outgoing particles
      for(unsigned int iz=0;iz<_outspin.size()-1;++iz) {
	if(int(iz)<id) temp *= rhoout[iz](ihel1[iz+2],ihel2[iz+2]);
	else           temp *= rhoout[iz](ihel1[iz+3],ihel2[iz+3]);
      }
      output(ihel1[id+2],ihel2[id+2])+=temp;
    }
  }      
  // normalise the matrix so it has unit trace
  output.normalize();
  // return the answer
  return output;
}

double ProductionMatrixElement::average() const {
  double output(0.);
  for(unsigned int ix=0;ix<_matrixelement.size();++ix) {
    output += norm(_matrixelement[ix]);
  }
  return output;
}
 
double ProductionMatrixElement::average(const RhoDMatrix & in1, 
					const RhoDMatrix & in2) const {
  Complex output(0.);
  for( int ihel1=0;ihel1<int(_inspin[0]);++ihel1) {
    for( int ihel2=0;ihel2<int(_inspin[1]);++ihel2) {
      int loc1 = ihel1*_constants[1] + ihel2*_constants[2];
      for( int jhel1=0;jhel1<int(_inspin[0]);++jhel1) {
	for( int jhel2=0;jhel2<int(_inspin[1]);++jhel2) {
	  int loc2 = jhel1*_constants[1] + jhel2*_constants[2];
	  Complex fact = in1(ihel1,jhel1)*in2(ihel2,jhel2);
	  for(int ohel=0;ohel<_constants[2];++ohel) {
	    output += fact
	      *_matrixelement[loc1+ohel]*conj(_matrixelement[loc2+ohel]);
	  }
	}
      }
    }
  }
  return real(output);
}

Complex ProductionMatrixElement::average(const ProductionMatrixElement & me2,
					 const RhoDMatrix & in1, 
					 const RhoDMatrix & in2) const {
  Complex output(0.);
  for( int ihel1=0;ihel1<int(_inspin[0]);++ihel1) {
    for( int ihel2=0;ihel2<int(_inspin[1]);++ihel2) {
      int loc1 = ihel1*_constants[1] + ihel2*_constants[2];
      for( int jhel1=0;jhel1<int(_inspin[0]);++jhel1) {
	for( int jhel2=0;jhel2<int(_inspin[1]);++jhel2) {
	  int loc2 = jhel1*_constants[1] + jhel2*_constants[2];
	  Complex fact = in1(ihel1,jhel1)*in2(ihel2,jhel2);
	  for(int ohel=0;ohel<_constants[2];++ohel) {
	    output += fact
	      *_matrixelement[loc1+ohel]*conj(me2._matrixelement[loc2+ohel]);
	  }
	}
      }
    }
  }
  return output;
}
