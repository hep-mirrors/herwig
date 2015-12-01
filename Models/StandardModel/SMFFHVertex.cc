// -*- C++ -*-
//
// SMFFHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFHVertex class.
//

#include "SMFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

IBPtr SMFFHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SMFFHVertex::fullclone() const {
  return new_ptr(*this);
}


SMFFHVertex::SMFFHVertex()  {
  // set up for the couplings
  _couplast=InvEnergy();
  _idlast=0;
  _q2last=ZERO;
  _masslast=ZERO;
  _mw=ZERO;
  _fermion = 0;
  orderInGem(1);
  orderInGs(0);
}

void SMFFHVertex::doinit() {
  if ( !_fermion ) {
    // PDG codes for the particles
    // the quarks
    for(int ix=1;ix<7;++ix) {
      addToList(-ix, ix, 25);
    }
    // the leptons
    for(int ix=11;ix<17;ix+=2) {
      addToList(-ix, ix, 25);
    }
  } else {
    addToList(-_fermion,_fermion,25);
  }
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if (!_theSM) 
    throw InitException();
  _mw= getParticleData(ThePEG::ParticleID::Wplus)->mass();
  FFSVertex::doinit();
}

void SMFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << _fermion;
}

void SMFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> _fermion;
}

ClassDescription<SMFFHVertex> 
SMFFHVertex::initSMFFHVertex;
// Definition of the static class description member.

void SMFFHVertex::Init() {

  static ClassDocumentation<SMFFHVertex> documentation
    ("The SMFFHVertex class is the implementation"
     " of the helicity amplitude calculation of the Standard Model Higgs"
     " fermion-antiferiom vertex.");

  static Parameter<SMFFHVertex,int> interfaceFermion
    ("Fermion",
     "The fermion to couple to. If not set all fermions are considered.",
     &SMFFHVertex::_fermion, 0, 0, 16,
     false, false, Interface::limited);
 
}

void SMFFHVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr, tcPDPtr) {
  int iferm=abs(aa->id());
  // left and right couplings set to one
  left (1.);
  right(1.);
  // first the overall normalisation
  if(q2!=_q2last||_couplast==complex<InvEnergy>()) {
    _couplast = -0.5*weakCoupling(q2)/_mw;
    _q2last=q2;
    _idlast=iferm;
    assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
    _masslast=_theSM->mass(q2,aa);
  }
  else if(iferm!=_idlast) {
    _idlast=iferm;
    assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
    _masslast=_theSM->mass(q2,aa);
  }
  norm(_couplast*_masslast);
}
