// -*- C++ -*-
//
// ProductionMatrixElement.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ProductionMatrixElement class.
//
// Author: Peter Richardson
//
#include "ProductionMatrixElement.h"

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
  for(ix=0;ix<_matrixelement.size();++ix) {
    // map the vector index to the helicities
    for(ixa=_outspin.size()+1;ixa>=0;--ixa)
      {ihel1[ixa]=(ix%_constants[ixa])/_constants[ixa+1];}
    // inner loop
    for(iy=0;iy<_matrixelement.size();++iy) {
      // map the vector index to the helicities
      for(iya=_outspin.size()+1;iya>=0;--iya)
	{ihel2[iya]=(iy%_constants[iya])/_constants[iya+1];}
      // matrix element piece
      temp=_matrixelement[ix]*conj(_matrixelement[iy]);
      // spin density matrices for the outgoing particles
      for(unsigned int iz=0;iz<_outspin.size();++iz)
	{temp*=rhoout[iz](ihel1[iz+2],ihel2[iz+2]);}
      // construct the spin density matrix
      if(id==0) {
	temp*=rhoin(ihel1[1],ihel2[1]);
	output(ihel1[0],ihel2[0])+=temp;
      }
      else {
	temp*=rhoin(ihel1[0],ihel2[0]);
	output(ihel1[1],ihel2[1])+=temp;
      }
    }
  }
  // normalise the matrix so it has unit trace
  output.normalize();
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

ProductionMatrixElement::ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out) {
  _nout=2;
  _inspin.resize(2);
  _inspin[0]=in1;
  _inspin[1]=in2;
  _outspin.push_back(out);
  setMESize();
}

ProductionMatrixElement::ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
						 PDT::Spin out2) {
  _nout=2;
  _inspin.resize(2);
  _inspin[0]=in1; 
  _inspin[1]=in2;
  _outspin.push_back(out1);
  _outspin.push_back(out2);
  setMESize();
}
  
ProductionMatrixElement::ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
						 PDT::Spin out2,PDT::Spin out3) {
  _inspin.resize(2);
  _nout=3;
  _inspin[0]=in1;
  _inspin[1]=in2;
  _outspin.push_back(out1);
  _outspin.push_back(out2);
  _outspin.push_back(out3);
  setMESize();
}

ProductionMatrixElement::ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
						 PDT::Spin out2,PDT::Spin out3, PDT::Spin out4) {
  _nout=4;
  _inspin.resize(2);
  _inspin[0]=in1;
  _inspin[1]=in2;
  _outspin.push_back(out1);
  _outspin.push_back(out2);
  _outspin.push_back(out3);
  _outspin.push_back(out4);
  setMESize();
}

ProductionMatrixElement::ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
						 PDT::Spin out2,PDT::Spin out3, PDT::Spin out4,
						 PDT::Spin out5) {
  _nout=5;
  _inspin.resize(2);
  _inspin[0]=in1;
  _inspin[1]=in2;
  _outspin.push_back(out1);
  _outspin.push_back(out2);
  _outspin.push_back(out3);
  _outspin.push_back(out4);
  _outspin.push_back(out5);
  setMESize();
}

ProductionMatrixElement::ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
						 PDT::Spin out2,PDT::Spin out3, PDT::Spin out4,
						 PDT::Spin out5, PDT::Spin out6) {
  _nout=6;
  _inspin.resize(2);
  _inspin[0]=in1;
  _inspin[1]=in2;
  _outspin.push_back(out1);
  _outspin.push_back(out2);
  _outspin.push_back(out3);
  _outspin.push_back(out4);
  _outspin.push_back(out5);
  _outspin.push_back(out6);
  setMESize();
}
  
ProductionMatrixElement::ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,vector<PDT::Spin> out) {
  _inspin.resize(2);
  _nout=out.size(); 
  _inspin[0]=in1;
  _inspin[1]=in2;
  _outspin=out;
  setMESize();
}
  
Complex   ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out) const {
  assert(_outspin.size()==1);
  unsigned int iloc = in1*_constants[1] + in2*_constants[2] + out*_constants[3];
  assert(iloc<_matrixelement.size());
  return _matrixelement[iloc];
}
  
Complex & ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out) {
  assert(_outspin.size()==1);
  unsigned int iloc = in1*_constants[1] + in2*_constants[2] + out*_constants[3];
  assert(iloc<_matrixelement.size());
  return _matrixelement[iloc];
}

Complex   ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2) const {
  assert(_outspin.size()==2);
  unsigned int iloc = in1*_constants[1] + in2*_constants[2] +
    out1*_constants[3] + out2*_constants[4];
  assert(iloc<_matrixelement.size());
  return _matrixelement[iloc];
}

Complex & ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2) {
  assert(_outspin.size()==2);
  unsigned int iloc = in1*_constants[1] + in2*_constants[2] +
    out1*_constants[3] + out2*_constants[4];
  assert(iloc<_matrixelement.size());
  return _matrixelement[iloc];
}

Complex   ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2,
						unsigned int out3) const {
  assert(_outspin.size()==3);
  vector<unsigned int> ivec(5);
  ivec[0]=in1;
  ivec[1]=in2;
  ivec[2]=out1;
  ivec[3]=out2;
  ivec[4]=out3;
  return (*this)(ivec);
}

Complex & ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2,
						unsigned int out3) {
  assert(_outspin.size()==3);
  vector<unsigned int> ivec(5);
  ivec[0]=in1;
  ivec[1]=in2;
  ivec[2]=out1;
  ivec[3]=out2;
  ivec[4]=out3;
  return (*this)(ivec);
}

Complex   ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2,
						unsigned int out3,unsigned int out4) const {
  assert(_outspin.size()==4);
  vector<unsigned int> ivec(6);
  ivec[0]=in1;
  ivec[1]=in2;
  ivec[2]=out1;
  ivec[3]=out2;
  ivec[4]=out3;
  ivec[5]=out4;
  return (*this)(ivec);
}
  
Complex & ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2,
						unsigned int out3, unsigned int out4) {
  assert(_outspin.size()==4);
  vector<unsigned int> ivec(6);
  ivec[0]=in1;
  ivec[1]=in2;
  ivec[2]=out1;
  ivec[3]=out2;
  ivec[4]=out3;
  ivec[5]=out4;
  return (*this)(ivec);
}

Complex   ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2,
						unsigned int out3,unsigned int out4,
						unsigned int out5) const {
  assert(_outspin.size()==5);
  vector<unsigned int> ivec(7);
  ivec[0]=in1;
  ivec[1]=in2;
  ivec[2]=out1;
  ivec[3]=out2;
  ivec[4]=out3;
  ivec[5]=out4;
  ivec[6]=out5;
  return (*this)(ivec);
}
  
Complex & ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2,
						unsigned int out3, unsigned int out4,
						unsigned int out5) {
  assert(_outspin.size()==5);
  vector<unsigned int> ivec(7);
  ivec[0]=in1;
  ivec[1]=in2;
  ivec[2]=out1;
  ivec[3]=out2;
  ivec[4]=out3;
  ivec[5]=out4;
  ivec[6]=out5;
  return (*this)(ivec);
}

Complex   ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2,
						unsigned int out3,unsigned int out4,
						unsigned int out5,unsigned int out6) const {
  assert(_outspin.size()==6);
  vector<unsigned int> ivec(8);
  ivec[0]=in1;
  ivec[1]=in2;
  ivec[2]=out1;
  ivec[3]=out2;
  ivec[4]=out3;
  ivec[5]=out4;
  ivec[6]=out5;
  ivec[7]=out6;
  return (*this)(ivec);
}

Complex & ProductionMatrixElement::operator () (unsigned int in1,unsigned int in2,
						unsigned int out1,unsigned int out2,
						unsigned int out3, unsigned int out4,
						unsigned int out5, unsigned int out6) {
  assert(_outspin.size()==6);
  vector<unsigned int> ivec(8);
  ivec[0]=in1;
  ivec[1]=in2;
  ivec[2]=out1;
  ivec[3]=out2;
  ivec[4]=out3;
  ivec[5]=out4;
  ivec[6]=out5;
  ivec[7]=out6;
  return (*this)(ivec);
}

Complex   ProductionMatrixElement::operator () (vector<unsigned int> hel) const {
  assert(_outspin.size() == hel.size()-2);
  unsigned int iloc(0),ix;
  // incoming and outgoing particles
  for(ix=0;ix<hel.size();++ix)
    iloc += hel[ix]*_constants[ix+1];
  assert(iloc<_matrixelement.size());
  return _matrixelement[iloc];
}
  
Complex & ProductionMatrixElement::operator () (vector<unsigned int> hel) {
  assert(_outspin.size() == hel.size()-2);
  unsigned int iloc=0,ix;
  // incoming particles
  for(ix=0;ix<hel.size();++ix)
    iloc += hel[ix]*_constants[ix+1];
  assert(iloc<_matrixelement.size());
  return _matrixelement[iloc];
}
//@}

void ProductionMatrixElement::reset(const ProductionMatrixElement & x) const {
  _nout = x._nout;
  _inspin = x._inspin;
  _outspin = x._outspin;
  _matrixelement = x._matrixelement;
  _constants     = x._constants;
}
  
void ProductionMatrixElement::setMESize() {
  unsigned int ix;
  int isize=_inspin[0]*_inspin[1];
  for(ix=0;ix<_outspin.size();++ix)
    isize*=_outspin[ix];
  // zero the matrix element
  _matrixelement.resize(isize,0.);
  // set up the constants for the mapping of helicity to vectro index
  _constants.resize(_outspin.size()+3);
  unsigned int temp=1;
  for(ix=_outspin.size()+1;ix>1;--ix) {
    temp*=_outspin[ix-2];
    _constants[ix]=temp;
  }
  temp*=_inspin[1];_constants[1]=temp;
  temp*=_inspin[0];_constants[0]=temp;
  _constants[_outspin.size()+2]=1;
}
  
