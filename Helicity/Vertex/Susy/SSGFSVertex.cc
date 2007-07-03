// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGFSVertex class.
//

#include "SSGFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

SSGFSVertex::SSGFSVertex() :_q2last(0.*sqr(MeV)),_couplast(0.), 
			    _id1last(0), _id2last(0) {
  vector<int> first,second,third;
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(1000021);
    second.push_back(ix);
    third.push_back(-(ix+1000000));
  }
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(1000021);
    second.push_back(ix);
    third.push_back(-(ix+2000000));
  }
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(1000021);
    second.push_back(-ix);
    third.push_back(ix+1000000);
  }
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(1000021);
    second.push_back(-ix);
    third.push_back(ix+2000000);
  }
  setList(first,second,third);
}

void SSGFSVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _stop << _sbottom;
}

void SSGFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS >> _stop >> _sbottom;
  _couplast = 0.;
  _q2last = 0.*sqr(MeV);
  _id1last = 0;
  _id2last = 0;
}

ClassDescription<SSGFSVertex> SSGFSVertex::initSSGFSVertex;
// Definition of the static class description member.

void SSGFSVertex::Init() {

  static ClassDocumentation<SSGFSVertex> documentation
    ("The SSGFSVertex implements coupling of the gluinos to the "
     "squarks and quarks");

}

void SSGFSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2, tcPDPtr part3, int iinc){
  long iferm(0),isc(0);
  if(abs(part1->id()) == 1000021) {
    if(part2->iSpin() == PDT::Spin1Half) {
      iferm = part2->id();
      isc = part3->id();
    }
    else {
      iferm = part3->id();
      isc = part2->id();
    }
  }
  if(abs(part2->id()) == 1000021) {
    if(part1->iSpin() == PDT::Spin1Half) {
      iferm = part1->id();
      isc = part3->id();
    }
    else {
      iferm = part3->id();
      isc = part1->id();
    }
  }
  if(abs(part3->id()) == 1000021) {
    if(part1->iSpin() == PDT::Spin1Half) {
      iferm = part1->id();
      isc = part2->id();
    }
    else {
      iferm = part2->id();
      isc = part1->id();
    }
  }
  if(abs(iferm) >=1 && abs(iferm) <=6) {
    if(q2 != _q2last) {
      double alphaStr = _theSS->alphaS(q2);
      _couplast = -2.*sqrt(2.*Constants::pi*alphaStr);
      _q2last = q2;
    }
    if(abs(iferm) != _id1last || abs(isc) != _id2last) { 
      _id1last = abs(iferm);
      _id2last = abs(isc);
      unsigned int eig = (abs(isc)/1000000) - 1;
      if(iferm == 6) {
	_leftlast = -(*_stop)(1,eig);
	_rightlast = (*_stop)(0,eig);
      }
      else if(iferm == 5){
	_leftlast = -(*_sbottom)(1,eig);
	_rightlast = (*_sbottom)(0,eig);
      }
      else {
	if(eig == 0) { 
	  _leftlast = 0.;
	  _rightlast = 1.;
	}
	else {
	  _leftlast = -1.;
	  _rightlast = 0.;
	}
      }
    }
    setNorm(_couplast);
    //arrange l/r couplings
    tcPDPtr incoming, fermion;
    //find incoming
    if(iinc == 1) {
      incoming = part1;
      if(part2->iSpin() == PDT::Spin1Half) fermion = part2;
      else fermion = part3;
    }
    else if(iinc == 2) {
      incoming = part2;
      if(part1->iSpin() == PDT::Spin1Half) fermion = part1;
      else fermion = part3;
    } 
    else {
      incoming = part3;
      if(part1->iSpin() == PDT::Spin1Half) fermion = part1;
      else fermion = part2;
    }
    //determine whether to flip couplings 
    if(incoming->iSpin() == PDT::Spin0) {
      if(incoming->id() > 0) {
	setLeft(_leftlast);
	setRight(_rightlast);
      }
      else {
	setLeft(conj(_rightlast));
	setRight(conj(_leftlast));
      }
    }
    else if(incoming->iSpin() == PDT::Spin1Half && incoming->CC()) {
      if(incoming->id() > 0) {
	setLeft(conj(_rightlast));
	setRight(conj(_leftlast));
      }
      else {
	setLeft(_leftlast);
	setRight(_rightlast);
      }
    }
    else {
      if(fermion->id() < 0) {
	setLeft(conj(_rightlast));
	setRight(conj(_leftlast));
      }
      else {
	setLeft(_leftlast);
	setRight(_rightlast);
      }
    }

  }
  else{
    throw HelicityConsistencyError() << "Incorrect particles detected in "
				     << "SSGFSVertex!"
				     << Exception::warning;
    setNorm(0.);
    setLeft(0.);
    setRight(0.);
  }
}
