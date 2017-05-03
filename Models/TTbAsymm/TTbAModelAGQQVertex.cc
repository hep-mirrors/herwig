// -*- C++ -*-
//
// TTbAModelAGQQVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModelAGQQVertex class.
//

#include "TTbAModelAGQQVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

IBPtr TTbAModelAGQQVertex::clone() const {
  return new_ptr(*this);
}

IBPtr TTbAModelAGQQVertex::fullclone() const {
  return new_ptr(*this);
}

TTbAModelAGQQVertex::TTbAModelAGQQVertex()  {
  orderInGem(1);
  orderInGs(1);
  addToList(-1,1,63);
  addToList(-2,2,63);
  addToList(-3,3,63);
  addToList(-4,4,63);
  addToList(-5,5,63);
  addToList(-6,6,63);

  

}

void TTbAModelAGQQVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwTTbAPtr hwTTbA=dynamic_ptr_cast<tcHwTTbAPtr>(_theModel);
  if(hwTTbA) {
    _cAGQQ_R =hwTTbA->_cAGQQ_right();
    _cAGQQ_L =hwTTbA->_cAGQQ_left();
    _cAGTT_R =hwTTbA->_cAGTT_right();
    _cAGTT_L =hwTTbA->_cAGTT_left();
    _models = hwTTbA->_model();
  }
  FFVVertex::doinit();
}

void TTbAModelAGQQVertex::persistentOutput(PersistentOStream & os) const {
  os << _cAGQQ_R << _cAGQQ_L << _cAGTT_R << _cAGTT_L << _models;
}

void TTbAModelAGQQVertex::persistentInput(PersistentIStream & is, int) {
  is >> _cAGQQ_R >> _cAGQQ_L >>_cAGTT_R >> _cAGTT_L >> _models;
}

ClassDescription<TTbAModelAGQQVertex> 
TTbAModelAGQQVertex::initTTbAModelAGQQVertex;
// Definition of the static class description member.


void TTbAModelAGQQVertex::Init() {
  
  static ClassDocumentation<TTbAModelAGQQVertex> documentation
    ("The TTbAModelAGQQVertex class is the implementation"
     " of the helicity amplitude calculation of the TTbA"
     " quark-lepton vertex.");
}


void TTbAModelAGQQVertex::setCoupling(Energy2 q2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  
  double _cL = 0, _cR = 0;
  double gstrong = 1.0;
  gstrong = strongCoupling(q2);

  if(abs(aa->id()) == 63 || abs(bb->id()) == 63 || abs(cc->id()) == 63) {
    if(abs(aa->id()) !=6 && abs(bb->id()) !=6 && abs(cc->id()) != 6) {
      _cR = _cAGQQ_R; 
      _cL = _cAGQQ_L; 
    } else { 
      _cR = _cAGTT_R; 
      _cL = _cAGTT_L; 
    }
    
  }

  if(_models!=2) { _cL = 1E-10; _cR = 1E-10; }
  left(_cL);
  right(_cR);
  
  norm(gstrong);
}
