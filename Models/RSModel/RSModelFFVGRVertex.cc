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
  : _charge(17,0.), _couplast(2,0.), _q2last(2,ZERO) {
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,22,39);
  }
  for(int ix=11;ix<17;++ix) {
    addToList(-ix,ix,22,39);
  }
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,21,39);
  }
  _theKappa=InvEnergy();
}

void RSModelFFVGRVertex::doinit() {
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
void RSModelFFVGRVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,
				     tcPDPtr cc, tcPDPtr) {
  // work out the particles
  int iferm=abs(aa->id());
  int ibos =cc->id();
  Complex coup;
  // overall factor
  assert(ibos==22||ibos==21);
  // photon
  if(ibos==22) {
    // alpha
    if(_q2last[0]!=q2||_couplast[0]==0.) {
      _couplast[0] = electroMagneticCoupling(q2);
      _q2last[0]=q2;
    }
    coup = UnitRemoval::E * _theKappa*_couplast[0];
    // _charge of particle
    assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
    coup *= _charge[iferm];
  }
  // gluon
  else if (ibos==21) {
    if(_q2last[1]!=q2||_couplast[1]==0.) {
      _couplast[1] = strongCoupling(q2);
      _q2last[1]=q2;
    }
    coup = UnitRemoval::E * _theKappa*_couplast[1];
  }
  // set the coupling
  norm(coup);
}

