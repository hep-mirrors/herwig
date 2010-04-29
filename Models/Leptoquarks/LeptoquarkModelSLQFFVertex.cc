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
  addToList(-16,-6,9921661);
  addToList(16,6,-9921661);

  //S1/2 doublet
  addToList(15,-6,9941561);
  addToList(-15,6,-9941561);
  
  addToList(-15,5,-9941551);
  addToList(-16,6,-9941551);
  addToList(15,-5,9941551);
  addToList(16,-6,9941551);



  //S1/2 tilde doublet
  addToList(16,-5,9951651);
  addToList(15,-5,9951551);

  addToList(-16,5,-9951651);
  addToList(-15,5,-9951551);

}

void LeptoquarkModelSLQFFVertex::doinit() {
  orderInGem(0);
  orderInGs(0);
  _theModel = generator()->standardModel();
  tcHwLeptoquarkPtr hwLeptoquark=dynamic_ptr_cast<tcHwLeptoquarkPtr>(_theModel);
  if(hwLeptoquark){
    _CFF=hwLeptoquark->_cfermion();
    _cL0 =hwLeptoquark->_cleft();
    _cR0 =hwLeptoquark->_cright();
    _cR0t = hwLeptoquark->_crighttilde();
    _cL1 =hwLeptoquark->_cleft1(); 
    _cL12 =hwLeptoquark->_cleft12(); 
    _cR12 =hwLeptoquark->_cright12(); 
    _cL12t =hwLeptoquark->_cleft12tilde(); 
    
  }
  FFSVertex::doinit();
}

void LeptoquarkModelSLQFFVertex::persistentOutput(PersistentOStream & os) const {
  os << _CFF << _cL0 << _cR0 << _cR0t << _cL1 << _cL12 << _cR12 << _cL12t;
}

void LeptoquarkModelSLQFFVertex::persistentInput(PersistentIStream & is, int) {
  is >> _CFF >> _cL0 >> _cR0 >> _cR0t >> _cL1 >> _cL12 >> _cR12 >> _cL12t;
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


void LeptoquarkModelSLQFFVertex::setCoupling(Energy2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  
  long isc(cc->id()), ism(aa->id()), 
    ichg(bb->id());
  
  //set the couplings to left and right 
  //S0
  if( fabs(isc) == 9911561 || fabs(ism) == 9911561 || fabs(ichg) == 9911561 ) {
    _cL = _cL0; _cR = _cR0;
  }
  //~S0
  if( fabs(isc) == 9911551 || fabs(ism) == 9911551 || fabs(ichg) == 9911551 ) {
    _cL = 0; _cR = _cR0t;
  }
  
  //S1 triplet
  //Q = + 4/3
  if( fabs(isc) == 9921551 || fabs(ism) == 9921551 || fabs(ichg) == 9921551 ) {
    _cL = - sqrt(2.)* _cL1; _cR = 0;
  }
  //Q = + 1/3
  if( fabs(isc) == 9921561 || fabs(ism) == 9921561 || fabs(ichg) == 9921561 ) {
    _cL = - _cL1; _cR = 0;
  }
  //Q = - 2/3
  if( fabs(isc) == 9921661 || fabs(ism) == 9921661 || fabs(ichg) == 9921661 ) {
    _cL = sqrt(2.) * _cL1; _cR = 0;
  }
  
  //S1/2 doublet

  //Q = + 5/3 
  if( fabs(isc) == 9941561 || fabs(ism) == 9941561 || fabs(ichg) == 9941561 ) {
    _cL = _cR12; _cR = _cL12;
  }
  
  
  //Q = + 2/3 
  if( fabs(isc) == 9941551 || fabs(ism) == 9941551 || fabs(ichg) == 9941551 ) {
    if(isc == 5 || ism == 5 || ichg == 5) { 
      _cL = - _cR12; _cR = 0.;
    }
    if(isc == 6 || ism == 6 || ichg == 6) { 
      _cL = 0.; _cR = _cL12;
    }
  }

  //S1/2 tilde doublet

  //Q = + 2/3 
  if( fabs(isc) == 9951551 || fabs(ism) == 9951551 || fabs(ichg) == 9951551 ) {
    _cL = 0.; _cR = _cL12t;
  }
  
  
  //Q = - 1/3 
  if( fabs(isc) == 9951651 || fabs(ism) == 9951651 || fabs(ichg) == 9951651 ) {
    _cL = 0; _cR = _cL12t;
  }
  
  left(_cL); right(_cR);
	 
  if( fabs(isc) != 9951551 && fabs(ism) != 9951551 && fabs(ichg) != 9951551 ) {
    //loop over generations (currently only third)
    for(int nl = 0; nl < 3; nl++) {
      //no right-handed neutrino (or left-handed anti-neutrino)
      //neutrino
      if( isc == (12+2*nl) || ism == (12+2*nl) || ichg == (12+2*nl) ) { 
	left(_cL); right(0.0);
      }
      //anti-neutrino      
      if( isc == -(12+2*nl) || ism == -(12+2*nl) || ichg == -(12+2*nl) ) { 
	left(0.0); right(_cL);
      }
    }
  }
  //swap couplings
  if( fabs(isc) > 9900000 || fabs(ism) > 9900000 || fabs(ichg) > 9900000 ) {
    if( isc < 0 || ism < 0 || ichg < 0 ) {
      left(_cR); right(_cL);
    }
  }
  
  
  norm(_CFF);
}
