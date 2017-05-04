// -*- C++ -*-
//
// TTbAModelZPQQVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModelZPQQVertex class.
//

#include "TTbAModelZPQQVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;



IBPtr TTbAModelZPQQVertex::clone() const {
  return new_ptr(*this);
}

IBPtr TTbAModelZPQQVertex::fullclone() const {
  return new_ptr(*this);
}

TTbAModelZPQQVertex::TTbAModelZPQQVertex()  {
  addToList(-2,6,32);
  addToList(-6,2,32);
  addToList(-2,2,32);
  addToList(-4,4,32);
  orderInGem(1);
  orderInGs(1);
}

void TTbAModelZPQQVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwTTbAPtr hwTTbA=dynamic_ptr_cast<tcHwTTbAPtr>(_theModel);
  if(hwTTbA) {
    _cZPTU_R =hwTTbA->_cZPTU_right();
    _cZPTU_L =hwTTbA->_cZPTU_left();
    _cZPUU_R =hwTTbA->_cZPUU_right();
    _cZPUU_L =hwTTbA->_cZPUU_left();  
    _cZPCC_R =hwTTbA->_cZPCC_right();
    _cZPCC_L =hwTTbA->_cZPCC_left();
    _models =hwTTbA->_model();

  }
  FFVVertex::doinit();
}

void TTbAModelZPQQVertex::persistentOutput(PersistentOStream & os) const {
  os << _cZPTU_R << _cZPTU_L << _cZPUU_R << _cZPUU_L << _cZPCC_R << _cZPCC_L << _models;
}

void TTbAModelZPQQVertex::persistentInput(PersistentIStream & is, int) {
  is >> _cZPTU_R >> _cZPTU_L >> _cZPUU_R >> _cZPUU_L >> _cZPCC_R >> _cZPCC_L >> _models;
}

ClassDescription<TTbAModelZPQQVertex> 
TTbAModelZPQQVertex::initTTbAModelZPQQVertex;
// Definition of the static class description member.


void TTbAModelZPQQVertex::Init() {
  
  static ClassDocumentation<TTbAModelZPQQVertex> documentation
    ("The TTbAModelZPQQVertex class is the implementation"
     " of the helicity amplitude calculation of the TTbA"
     " Z prime Quark-antiQuark vertex.");
}

void TTbAModelZPQQVertex::setCoupling(Energy2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  double _cR = 0, _cL = 0;
  if( abs(aa->id()) == 6 || abs(bb->id()) == 6 || abs(cc->id()) == 6) { 
    _cR = _cZPTU_R; 
    _cL = _cZPTU_L; 
  } else {
    if( abs(aa->id()) != 4 && abs(bb->id()) != 4 && abs(cc->id()) != 4) { 
      _cR = _cZPUU_R; 
      _cL = _cZPUU_L;
    }
    if( abs(aa->id()) == 4 || abs(bb->id()) == 4 || abs(cc->id()) == 4) { 
      _cR = _cZPCC_R; 
      _cL = _cZPCC_L;
    }
  }
  
  if(_models!=1) { _cL = 1E-10; _cR = 1E-10; }
  right(_cR);
  left(_cL);

  norm(1.0);
}
