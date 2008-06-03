// -*- C++ -*-
//
// RSModelFFVGRVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFVGRVertex class.
//

#include "RSModelFFVGRVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


using namespace Herwig;
using namespace ThePEG;

RSModelFFVGRVertex::RSModelFFVGRVertex() 
  : _charge(17,0.), _couplast(2,0.), _q2last(2,0.*GeV2) {
  vector<long> first,second,third,fourth;
  for(int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(22);
    fourth.push_back(39);
  }
  for(int ix=11;ix<17;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(22);
    fourth.push_back(39);
  }
  for(int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(21);
    fourth.push_back(39);
  }
  setList(first,second,third,fourth);
  _theKappa=InvEnergy();
}

void RSModelFFVGRVertex::doinit() throw(InitException) {
  orderInGem(1);
  FFVTVertex::doinit();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  for(int ix=1;ix<4;++ix) {
    _charge[2*ix-1]  = (generator()->standardModel()->ed());
    _charge[2*ix ]   = (generator()->standardModel()->eu());
    _charge[2*ix+9 ] = (generator()->standardModel()->ee());
    _charge[2*ix+10] = (generator()->standardModel()->enu());
  }
  if(hwRS) _theKappa=2./hwRS->lambda_pi();
}

void RSModelFFVGRVertex::persistentOutput(PersistentOStream & os) const {
  os << _charge << ounit(_theKappa,InvGeV);
}

void RSModelFFVGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> _charge >> iunit(_theKappa,InvGeV);
}

ClassDescription<RSModelFFVGRVertex> RSModelFFVGRVertex::initRSModelFFVGRVertex;
// Definition of the static class description member.

void RSModelFFVGRVertex::Init() {
  static ClassDocumentation<RSModelFFVGRVertex> documentation
    ("The RSModelFFVGRVertexxs class is the implementation"
     " of the two fermion vector coupling for the RS model.");
  
}
// FFVGR coupling
void RSModelFFVGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,
    				 tcPDPtr c, tcPDPtr) {
  // work out the particles
  int iferm=abs(a->id());
  int ibos =c->id();
  Complex norm;
  // overall factor
  // photon
  if(ibos==22) {
    // alpha
    if(_q2last[0]!=q2) {
      _couplast[0] = electroMagneticCoupling(q2);
      _q2last[0]=q2;
    }
    norm = UnitRemoval::E * _theKappa*_couplast[0];
    // _charge of particle
    if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16))
      {norm = norm*_charge[iferm];}
    else throw HelicityConsistencyError() << "RSModelFFVGRVertex::setCoupling " 
					  << "Unknown particle in FFVGR vertex" 
					  << Exception::runerror;
  }
  // gluon
  else if (ibos==21||ibos==9) {
    if(_q2last[1]!=q2) {
      _couplast[1] = strongCoupling(q2);
      _q2last[1]=q2;
    }
    norm = UnitRemoval::E * _theKappa*_couplast[1];
  }
  else throw HelicityConsistencyError() << "RSModelFFVGRVertex::setCoupling " 
					<< "Unknown boson in FFVGR vertex" 
					<< Exception::runerror;
  // set the coupling
  setNorm(norm);
}

