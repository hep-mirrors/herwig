// -*- C++ -*-
//
// LeptoquarkModelSLQFFVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModelSLQFFVertex class.
//

#include "LeptoquarkModelSLQFFVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

IBPtr LeptoquarkModelSLQFFVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LeptoquarkModelSLQFFVertex::fullclone() const {
  return new_ptr(*this);
}

LeptoquarkModelSLQFFVertex::LeptoquarkModelSLQFFVertex()  {
  //S0
  addToList(-15,-6,9911561);
  addToList(15,6,-9911561);
  addToList(-16,-5,9911561);
  addToList(16,5,-9911561);

  //~S0
  addToList(-15,-5,9911551);
  addToList(15,5,-9911551);

  //S1 triplet
  //S1p
  addToList(-15,-5,9921551);
  addToList(15,5,-9921551);
  //S1z
  addToList(-15,-6,9921561);
  addToList(15,6,-9921561);
  addToList(-16,-5,9921561);
  addToList(16,5,-9921561);
  //S1m
  addToList(-16,-6,99216611);
  addToList(16,6,-9921661);
  
  
  /*addToList(-11,-2,9911561);
  addToList(11,2,-9911561);
  addToList(-12,-1,9911561);
  addToList(12,1,-9911561);
  */
  _q2last=ZERO;
}

void LeptoquarkModelSLQFFVertex::doinit() {
  orderInGem(0);
  orderInGs(0);
  _theModel = generator()->standardModel();
  tcHwLeptoquarkPtr hwLeptoquark=dynamic_ptr_cast<tcHwLeptoquarkPtr>(_theModel);
  if(hwLeptoquark){

    _CFF=hwLeptoquark->_cfermion();
    _cL =hwLeptoquark->_cleft();
    _cR =hwLeptoquark->_cright();
    
  }
  FFSVertex::doinit();
}

void LeptoquarkModelSLQFFVertex::persistentOutput(PersistentOStream & os) const {
  os << _CFF << _cL << _cR;
}

void LeptoquarkModelSLQFFVertex::persistentInput(PersistentIStream & is, int) {
  is >> _CFF >> _cL >> _cR;
}

ClassDescription<LeptoquarkModelSLQFFVertex> 
LeptoquarkModelSLQFFVertex::initLeptoquarkModelSLQFFVertex;
// Definition of the static class description member.


void LeptoquarkModelSLQFFVertex::Init() {
  
  static ClassDocumentation<LeptoquarkModelSLQFFVertex> documentation
    ("The LeptoquarkModelSLQFFVertex class is the implementation"
     " of the helicity amplitude calculation of the Leptoquark"
     " quark-lepton vertex.");
}


void LeptoquarkModelSLQFFVertex::setCoupling(Energy2 q2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  
  long isc(cc->id()), ism(aa->id()), 
    ichg(bb->id());
  
 
  for(int nl = 0; nl < 3; nl++) {
    if( isc == -(12+2*nl) || ism == -(12+2*nl) || ichg == -(12+2*nl) ) { 
      left(0.0); right(_cR);
    }
    
    
    if( isc == (12+2*nl) || ism == (12+2*nl) || ichg == (12+2*nl) ) { 
      left(_cL); right(0.0);
    }
    if( fabs(isc) != (12+2*nl) &&  fabs(ism) != (12+2*nl) && fabs(ichg) != (12+2*nl) ) {
      if( isc == (11+2*nl) || ichg == (11+2*nl) || ism == (11+2*nl)) {
	left(_cL);
	right(_cR);
      }
      if( isc == -(11+2*nl) || ichg == -(11+2*nl) || ism == -(11+2*nl)) {
	left(_cR);
	right(_cL);
      }
    }
  }
  norm(_CFF);
}
