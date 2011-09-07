// -*- C++ -*-
//
// ZprimeModelZPQQVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZprimeModelZPQQVertex class.
//

#include "ZprimeModelZPQQVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;



IBPtr ZprimeModelZPQQVertex::clone() const {
  return new_ptr(*this);
}

IBPtr ZprimeModelZPQQVertex::fullclone() const {
  return new_ptr(*this);
}

ZprimeModelZPQQVertex::ZprimeModelZPQQVertex()  {
  addToList(-6,6,32);
  addToList(-5,5,32);
  addToList(-4,4,32);
  addToList(-3,3,32);
  addToList(-2,2,32);
  addToList(-1,1,32);
  orderInGem(0);
  orderInGs(0);
}

void ZprimeModelZPQQVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwZprimePtr hwZprime=dynamic_ptr_cast<tcHwZprimePtr>(_theModel);
  if(hwZprime) {
    _cZPTT_R =hwZprime->_cZPTT_right();
    _cZPTT_L =hwZprime->_cZPTT_left();
    _cZPUU_R =hwZprime->_cZPUU_right();
    _cZPUU_L =hwZprime->_cZPUU_left();  
    _cZPCC_R =hwZprime->_cZPCC_right();
    _cZPCC_L =hwZprime->_cZPCC_left();
    _cZPDD_R =hwZprime->_cZPDD_right();
    _cZPDD_L =hwZprime->_cZPDD_left();
    _cZPBB_R =hwZprime->_cZPBB_right();
    _cZPBB_L =hwZprime->_cZPBB_left();  
    _cZPSS_R =hwZprime->_cZPSS_right();
    _cZPSS_L =hwZprime->_cZPSS_left();
    _cZP_o =hwZprime->_cZPoverallCoup();

  }
  FFVVertex::doinit();
}

void ZprimeModelZPQQVertex::persistentOutput(PersistentOStream & os) const {
  os << _cZPTT_R << _cZPTT_L << _cZPUU_R << _cZPUU_L << _cZPCC_R << _cZPCC_L << _cZPDD_R << _cZPDD_L << _cZPSS_R << _cZPSS_L  << _cZPBB_R << _cZPBB_L << _cZP_o;
}

void ZprimeModelZPQQVertex::persistentInput(PersistentIStream & is, int) {

is >> _cZPTT_R >> _cZPTT_L >> _cZPUU_R >> _cZPUU_L >> _cZPCC_R >> _cZPCC_L >> _cZPDD_R >> _cZPDD_L >> _cZPSS_R >> _cZPSS_L  >> _cZPBB_R >> _cZPBB_L >> _cZP_o;

}

ClassDescription<ZprimeModelZPQQVertex> 
ZprimeModelZPQQVertex::initZprimeModelZPQQVertex;
// Definition of the static class description member.


void ZprimeModelZPQQVertex::Init() {
  
  static ClassDocumentation<ZprimeModelZPQQVertex> documentation
    ("The ZprimeModelZPQQVertex class is the implementation"
     " of the helicity amplitude calculation of the Zprime"
     " Z prime Quark-antiQuark vertex.");
}

void ZprimeModelZPQQVertex::setCoupling(Energy2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  double _cR = 1.0, _cL = 1.0;

  right(_cR);
  left(_cL);

  norm(_cZP_o);
}
