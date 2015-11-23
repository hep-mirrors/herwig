// -*- C++ -*-
//
// TTbAModelWPTDVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModelWPTDVertex class.
//

#include "TTbAModelWPTDVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

IBPtr TTbAModelWPTDVertex::clone() const {
  return new_ptr(*this);
}

IBPtr TTbAModelWPTDVertex::fullclone() const {
  return new_ptr(*this);
}

TTbAModelWPTDVertex::TTbAModelWPTDVertex()  {
  addToList(-1,6,-34);
  addToList(-6,1,34);
  orderInGem(1);
  orderInGs(1);
}

void TTbAModelWPTDVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwTTbAPtr hwTTbA=dynamic_ptr_cast<tcHwTTbAPtr>(_theModel);
  if(hwTTbA) {
    _cWPTD_R =hwTTbA->_cWPTD_right();
    _cWPTD_L =hwTTbA->_cWPTD_left();
    _models =hwTTbA->_model();
  }
  FFVVertex::doinit();
}

void TTbAModelWPTDVertex::persistentOutput(PersistentOStream & os) const {
  os << _cWPTD_R << _cWPTD_L << _models;
}

void TTbAModelWPTDVertex::persistentInput(PersistentIStream & is, int) {
  is >> _cWPTD_R >> _cWPTD_L >> _models; 
}

ClassDescription<TTbAModelWPTDVertex> 
TTbAModelWPTDVertex::initTTbAModelWPTDVertex;
// Definition of the static class description member.


void TTbAModelWPTDVertex::Init() {
  
  static ClassDocumentation<TTbAModelWPTDVertex> documentation
    ("The TTbAModelWPTDVertex class is the implementation"
     " of the helicity amplitude calculation of the TTbA"
     " quark-lepton vertex.");
}


void TTbAModelWPTDVertex::setCoupling(Energy2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  
  double _cL = 0, _cR = 0;
  
  if(abs(aa->id()) == 34 || abs(bb->id()) == 34 || abs(cc->id()) == 34) {
    _cR = _cWPTD_R; 
    _cL = _cWPTD_L; 
  }
  if(_models!=0) { _cL = 1E-10; _cR = 1E-10; }
  left(_cL);
  right(_cR);
  
  norm(1.0);
}
