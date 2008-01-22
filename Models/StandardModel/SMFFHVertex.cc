// -*- C++ -*-
//
// SMFFHVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFHVertex class.
//

#include "SMFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

SMFFHVertex::SMFFHVertex()  {
  // PDG codes for the particles
  vector<int> first,second,third;
  // the quarks
  for(unsigned int ix=1;ix<7;++ix) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(25);
  }
  // the leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(25);
  }
  setList(first,second,third);
  // set up for the couplings
  _couplast=InvEnergy();
  _idlast=0;
  _q2last=0.*GeV2;
  _masslast=0.*GeV;
  _mw=0.*GeV;
}

void SMFFHVertex::doinit() throw(InitException) {
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if (!_theSM) 
    throw InitException();
  _mw= getParticleData(ThePEG::ParticleID::Wplus)->mass();
  orderInGem(1);
  orderInGs(0);
  FFSVertex::doinit();
}

void SMFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV);
}

void SMFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV);
}

ClassDescription<SMFFHVertex> 
SMFFHVertex::initSMFFHVertex;
// Definition of the static class description member.

void SMFFHVertex::Init() {

  static ClassDocumentation<SMFFHVertex> documentation
    ("The SMFFHVertex class is the implementation"
     " of the helicity amplitude calculation of the Standard Model Higgs"
     " fermion-antiferiom vertex.");
  
}

void SMFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr, tcPDPtr, int) {
  int iferm=abs(a->id());
  // left and right couplings set to one
  setLeft(1.); setRight(1.);
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -0.5*weakCoupling(q2)/_mw;
    _q2last=q2;
    _idlast=iferm;
    if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
      _masslast=_theSM->mass(q2,a);
    }
    else {
      throw HelicityConsistencyError() << "SMFFHVertex::setCoupling " 
				       << "Unknown particle in Higgs vertex" 
				       << Exception::warning;
      _masslast = 0*MeV;
    }
  }
  else if(iferm!=_idlast) {
    _idlast=iferm;
    if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
      _masslast=_theSM->mass(q2,a);
    }
    else {
      throw HelicityConsistencyError() << "SMFFHVertex::setCoupling " 
				       << "Unknown particle in Higgs vertex" 
				       << Exception::warning;
      _masslast = 0*MeV;
    }
  }
  setNorm(_couplast*_masslast);
}
