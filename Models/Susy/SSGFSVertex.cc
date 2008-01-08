// -*- C++ -*-
//
// SSGFSVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGFSVertex class.
//

#include "SSGFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

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
  os << _stop << _sbottom;
}

void SSGFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _stop >> _sbottom;
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
      iferm = abs(part2->id());
      isc = abs(part3->id());
    }
    else {
      iferm = abs(part3->id());
      isc = abs(part2->id());
    }
  }
  else if(abs(part2->id()) == 1000021) {
    if(part1->iSpin() == PDT::Spin1Half) {
      iferm = abs(part1->id());
      isc = abs(part3->id());
    }
    else {
      iferm = abs(part3->id());
      isc = abs(part1->id());
    }
  }
  else if(abs(part3->id()) == 1000021) {
    if(part1->iSpin() == PDT::Spin1Half) {
      iferm = abs(part1->id());
      isc = abs(part2->id());
    }
    else {
      iferm = abs(part2->id());
      isc = abs(part1->id());
    }
  }
  else {
    throw HelicityConsistencyError()
      << "SSGFSVertex::setCoupling() - There is no gluino in this vertex!"
      << part1->id() << " " << part2->id() << " " << part3->id()
      << Exception::warning;
    setNorm(0.);
    setLeft(0.);
    setRight(0.);
    return;
  }    
  if(iferm >=1 && iferm <=6) {
    if(q2 != _q2last) {
      _couplast = -strongCoupling(q2)*sqrt(2.);
      _q2last = q2;
    }
    if(iferm != _id1last || isc != _id2last) { 
      _id1last = iferm;
      _id2last = isc;
      unsigned int eig = (isc/1000000) - 1;
      if(iferm == 6) {
	_leftlast = -(*_stop)(eig,1);
	_rightlast = (*_stop)(eig,0);
      }
      else if(iferm == 5){
	_leftlast = -(*_sbottom)(eig,1);
	_rightlast = (*_sbottom)(eig,0);
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
    throw HelicityConsistencyError() 
      << "SSGFSVertex::setCoupling() - There are unknown particles in "
      << "this vertex. " << part1->id() << " " << part2->id() << " " 
      << part3->id() << Exception::warning;
    setNorm(0.);
    setLeft(0.);
    setRight(0.);
  }
}
