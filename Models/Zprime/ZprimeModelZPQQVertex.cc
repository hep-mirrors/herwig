// -*- C++ -*-
//
// ZprimeModelZPQQVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
  addToList(-2,6,32);
  addToList(-6,2,32);
  addToList(-6,6,32);
  addToList(-5,5,32);
  addToList(-4,4,32);
  addToList(-3,3,32);
  addToList(-2,2,32);
  addToList(-1,1,32);
  addToList(-11,11,32);
  addToList(-13,13,32);
  addToList(-15,15,32);
  addToList(-12,12,32);
  addToList(-14,14,32);
  addToList(-16,16,32);
  orderInGem(1);
  orderInGs(0);
}

void ZprimeModelZPQQVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwZprimePtr hwZprime=dynamic_ptr_cast<tcHwZprimePtr>(_theModel);
  if(hwZprime) {

    _cZPTU_R =hwZprime->_cZPTU_right();
    _cZPTU_L =hwZprime->_cZPTU_left();

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

    _cZPee_R =hwZprime->_cZPee_right();
    _cZPee_L =hwZprime->_cZPee_left();
    _cZPmm_R =hwZprime->_cZPmm_right();
    _cZPmm_L =hwZprime->_cZPmm_left();  
    _cZPtt_R =hwZprime->_cZPtt_right();
    _cZPtt_L =hwZprime->_cZPtt_left();

    _cZPnuenue_R =hwZprime->_cZPnuenue_right();
    _cZPnuenue_L =hwZprime->_cZPnuenue_left();
    _cZPnuenue_R =hwZprime->_cZPnumnum_right();
    _cZPnumnum_L =hwZprime->_cZPnumnum_left();  
    _cZPnutnut_R =hwZprime->_cZPnutnut_right();
    _cZPnutnut_L =hwZprime->_cZPnutnut_left();

    _cZP_o =hwZprime->_cZPoverallCoup();

  }
  FFVVertex::doinit();
}

void ZprimeModelZPQQVertex::persistentOutput(PersistentOStream & os) const {
  os << _cZPTU_R << _cZPTU_L  << _cZPTT_R << _cZPTT_L << _cZPUU_R << _cZPUU_L << _cZPCC_R << _cZPCC_L << _cZPDD_R << _cZPDD_L << _cZPSS_R << _cZPSS_L  << _cZPBB_R << _cZPBB_L << _cZPee_R << _cZPee_L << _cZPmm_R << _cZPmm_L  << _cZPtt_R << _cZPtt_L << _cZPnuenue_R << _cZPnuenue_L << _cZPnumnum_R << _cZPnumnum_L  << _cZPnutnut_R << _cZPnutnut_L <<  _cZP_o;
}

void ZprimeModelZPQQVertex::persistentInput(PersistentIStream & is, int) {

  is  >> _cZPTU_R >> _cZPTU_L  >> _cZPTT_R >> _cZPTT_L >> _cZPUU_R >> _cZPUU_L >> _cZPCC_R >> _cZPCC_L >> _cZPDD_R >> _cZPDD_L >> _cZPSS_R >> _cZPSS_L  >> _cZPBB_R >> _cZPBB_L >> _cZPee_R >> _cZPee_L >> _cZPmm_R >> _cZPmm_L >> _cZPtt_R >> _cZPtt_L >> _cZPnuenue_R >> _cZPnuenue_L >> _cZPnumnum_R >> _cZPnumnum_L >> _cZPnutnut_R >> _cZPnutnut_L >> _cZP_o;

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

  long ccc(cc->id()), aaa(aa->id()), bbb(bb->id());

  if( abs(aaa) == 6 || abs(bbb) == 6  || abs(ccc) == 6 ) {
    if( abs(aaa) !=2 && abs(bbb) !=2  && abs(ccc) != 2 ) {
    _cL = _cZPTT_L; _cR = _cZPTT_R;
    } else if( abs(aaa) ==2 || abs(bbb) ==2  || abs(ccc) == 2 ) {
      _cL = _cZPTU_L; _cR = _cZPTU_R;
    }
  }
  
  if( abs(aaa) == 5 || abs(bbb) == 5  || abs(ccc) == 5 ) {
    _cL = _cZPBB_L; _cR = _cZPBB_R;
  }
    
  if( abs(aaa) == 4 || abs(bbb) == 4  || abs(ccc) == 4 ) {
    _cL = _cZPCC_L; _cR = _cZPCC_R;
  }
  
  
  if( abs(aaa) == 3 || abs(bbb) == 3  || abs(ccc) == 3 ) {
    _cL = _cZPSS_L; _cR = _cZPSS_R;
  }
  
  
  if( (abs(aaa) == 2 || abs(bbb) == 2  || abs(ccc) == 2) &&  (abs(aaa) !=6 && abs(bbb) !=6  && abs(ccc) != 6)) {
    _cL = _cZPUU_L; _cR = _cZPUU_R;
  }

 
  if( abs(aaa) == 1 || abs(bbb) == 1  || abs(ccc) == 1 ) {
    _cL = _cZPDD_L; _cR = _cZPDD_R;
  }
  
   if( abs(aaa) == 11 || abs(bbb) == 11  || abs(ccc) == 11 ) {
    _cL = _cZPee_L; _cR = _cZPee_R;
  }
  
  
  if( abs(aaa) == 13 || abs(bbb) == 13  || abs(ccc) == 13 ) {
    _cL = _cZPmm_L; _cR = _cZPmm_R;
  }

 
  if( abs(aaa) == 15 || abs(bbb) == 15  || abs(ccc) == 15 ) {
    _cL = _cZPtt_L; _cR = _cZPtt_R;
  }


   if( abs(aaa) == 12 || abs(bbb) == 12  || abs(ccc) == 12 ) {
    _cL = _cZPnuenue_L; _cR = _cZPnuenue_R;
  }
  
  
  if( abs(aaa) == 14 || abs(bbb) == 14 || abs(ccc) == 14 ) {
    _cL = _cZPnumnum_L; _cR = _cZPnumnum_R;
  }

 
  if( abs(aaa) == 16 || abs(bbb) == 16  || abs(ccc) == 16 ) {
    _cL = _cZPnutnut_L; _cR = _cZPnutnut_R;
  }

  right(_cR);
  left(_cL);

  norm(_cZP_o);
}
