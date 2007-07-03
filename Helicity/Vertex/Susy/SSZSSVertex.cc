// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSZSSVertex class.
//

#include "SSZSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

inline SSZSSVertex::SSZSSVertex() : _sw(0.),_cw(0.),_zfact(0.),_couplast(0.),
				    _q2last(),_idlast(0) {
  vector<int> first,second,third;
  //LL-sleptons
  for(unsigned int ix=1000011;ix<1000017;++ix) {
    first.push_back(23);
    second.push_back(ix);
    third.push_back(-ix);
  }
  //RR-sleptons
  for(unsigned int ix=2000011;ix<2000016;ix+=2) {
    first.push_back(23);
    second.push_back(ix);
    third.push_back(-ix);
  }
  //L-Rbar stau
  first.push_back(23);
  second.push_back(1000015);
  third.push_back(-2000015);
  //Lbar-R stau
  first.push_back(23);
  second.push_back(-1000015);
  third.push_back(2000015);
   
  //LL squarks
  for(unsigned int ix=1000001;ix<1000007;++ix) {
    first.push_back(23);
    second.push_back(ix);
    third.push_back(-ix);
  }
  //RR squarks
  for(unsigned int ix=2000001;ix<2000007;++ix) {
    first.push_back(23);
    second.push_back(ix);
    third.push_back(-ix);
  }
 //L-Rbar stop
  first.push_back(23);
  second.push_back(1000006);
  third.push_back(-2000006);
  //Lbar-R stop
  first.push_back(23);
  second.push_back(-1000006);
  third.push_back(2000006);

  //L-Rbar sbottom
  first.push_back(23);
  second.push_back(1000005);
  third.push_back(-2000005);
  //Lbar-R sbottom
  first.push_back(23);
  second.push_back(-1000005);
  third.push_back(2000005);
  
  setList(first,second,third);
}

void SSZSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _sw << _cw << _stop << _sbottom << _stau;
}

void SSZSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS >> _sw >> _cw >> _stop >> _sbottom >> _stau;
  _couplast=0;
  _q2last=0*GeV2;
  _idlast=0;
}

ClassDescription<SSZSSVertex> SSZSSVertex::initSSZSSVertex;
// Definition of the static class description member.

void SSZSSVertex::Init() {

  static ClassDocumentation<SSZSSVertex> documentation
    ("The SSZSSVertex implements the coupling of a Z0 to "
     "2 sfermions");

}

void SSZSSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long isf1(0),isf2(0);
  if(part1->id()==23) {
    isf1 = abs(part2->id());
    isf2 = abs(part3->id());
  }
  else if(part2->id()==23) {
    isf1 = abs(part1->id());
    isf2 = abs(part3->id());
  }
  else {
    isf1 = abs(part1->id());
    isf2 = abs(part2->id());
  }
  unsigned int eig1(isf1/1000000 -1),eig2(isf2/1000000 -1);
  if(eig1 == eig2 || isf1 == 1000006 || isf1 == 2000006 || 
     isf1 == 1000005 || isf1 == 2000005|| isf1 == 1000015 || 
     isf1==2000015) {
    if(q2 != _q2last) {
      double alpha = _theSS->alphaEM(q2);
      _couplast = -sqrt(4.*Constants::pi*alpha)/_cw/_sw;
      _q2last=q2;
    }
    if(isf1!=_idlast) {
      if(isf1==1000006||isf1==2000006) {
	_zfact = -0.5*(*_stop)(0,eig1)*(*_stop)(0,eig2);
	if(eig1==eig2) {
	  _zfact += (_theSS->eu())*_sw*_sw;
	}
      }
      else if(isf1==1000005||isf1==2000005) {
	_zfact = 0.5*(*_sbottom)(0,eig1)*(*_sbottom)(0,eig2);
	if(eig1==eig2) {
	  _zfact += (_theSS->ed())*_sw*_sw;
	}
      }
      else if(isf1==1000015||isf1==2000015) {
	_zfact = 0.5*(*_stau)(0,eig1)*(*_stau)(0,eig2);
	if(eig1==eig2) {
	  _zfact -= _sw*_sw; 
	}
      }
      else {
	if(isf1==1000012||isf1==1000014||isf1==1000016) {
	  _zfact = -0.5;
	}
	else if(isf1==1000011||isf1==1000013||
		isf1==2000011||isf1==2000013) {
	  if(eig1==0 && eig2==0) {
	    _zfact = 0.5*(1.- 2.*_sw*_sw);
	  }
	  else {
	    _zfact = -_sw*_sw;
	  }
	}
	else if(isf1==1000001||isf1==1000003||
		isf1==2000001||isf1==2000003) {
	  if(eig1==0 && eig2==0) {
	    _zfact = 0.5*(1. + 2.*(_theSS->ed())*_sw*_sw);
	  }
	  else {
	    _zfact = _theSS->ed()*_sw*_sw;
	  }
	}
	else {
	  if(eig1==0 && eig2==0) {
	    _zfact = 0.5*(-1. + 2.*(_theSS->eu())*_sw*_sw);
	  }
	  else {
	    _zfact = _theSS->eu()*_sw*_sw;
	  }
	}
      }
     _idlast = isf1;
    }
    setNorm(_couplast*_zfact);
  }
  else {
    throw HelicityConsistencyError() << "Incorrect particle found in "
				     << "SSZSSVertex. "
				     << Exception::warning;
  }
}
