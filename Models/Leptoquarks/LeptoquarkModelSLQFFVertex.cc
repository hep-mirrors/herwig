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
  addToList(-15,-5,9921551);
  addToList(15,5,-9921551);

  //S1 triplet
  //S1p
  addToList(-15,-5,9931551);
  addToList(15,5,-9931551);
  //S1z
  addToList(-15,-6,9931561);
  addToList(15,6,-9931561);
  addToList(-16,-5,9931561);
  addToList(16,5,-9931561);
  //S1m
  addToList(-16,-6,9931661);
  addToList(16,6,-9931661);

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

  //dS0
  addToList(-15,-5,9961551);
  addToList(15,5,-9961551);

  addToList(-16,-6,9961551);
  addToList(16,6,-9961551);

  //~dS0
  addToList(-15,-6,9971561);
  addToList(15,6,-9971561);


  //dS1 triplet

  //dS1p
  addToList(-15,-6,9981561);
  addToList(15,6,-9981561);

  //dS1z
  addToList(-16,-6,9981551);
  addToList(16,6,-9981551);

  addToList(-15,-5,9981551);
  addToList(15,5,-9981551);

  //dS1m
  addToList(-16,-5,9981651);
  addToList(16,5,-9981651);

  //dS1/2 doublet
  addToList(-15,-5,9991551);
  addToList(15,5,-9991551);

  addToList(-15,-6,9991561);
  addToList(15,6,-9991561);

  addToList(-16,-5,9991561);
  addToList(16,5,-9991561);

  //dS1/2 tilde doublet
  addToList(-15,-6,9901561);
  addToList(15,6,-9901561);

  addToList(-16,-6,9901661);
  addToList(16,6,-9901661);

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

    
    _derivscale = hwLeptoquark->_fscale();
    _dcL0 =hwLeptoquark->_dcleft();
    _dcR0 =hwLeptoquark->_dcright();
    _dcR0t = hwLeptoquark->_dcrighttilde();
    _dcL1 =hwLeptoquark->_dcleft1(); 
    _dcL12 =hwLeptoquark->_dcleft12(); 
    _dcR12 =hwLeptoquark->_dcright12(); 
    _dcL12t =hwLeptoquark->_dcleft12tilde(); 
    
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
  double mtop = 174.2;
  double mbot = 4.2;
  double mtau = 1.77699;
  int lqtype1[9] = { 9911561, 9921551, 9931551, 9931561, 9931661, 9991551, 9991561, 9901561, 9901661 };
  int lqtype2[9] = { 9941561, 9941551, 9951551, 9951651, 9961651, 9971561, 9981561, 9981551, 9981651 };


  //set the couplings to left and right 
  //S0
  if( fabs(isc) == 9911561 || fabs(ism) == 9911561 || fabs(ichg) == 9911561 ) {
    if(fabs(isc) == 5 || fabs(ism) == 5 || fabs(ichg) == 5) { 
      _cL = -_cL0; _cR = 0.;
    }
    if(fabs(isc) == 6 || fabs(ism) == 6 || fabs(ichg) == 6) { 
      _cL = _cL0; _cR = _cR0;
    }
  }
  //~S0
  if( fabs(isc) == 9921551 || fabs(ism) == 9921551 || fabs(ichg) == 9921551 ) {
    _cL = 0; _cR = _cR0t;
  }
  
  //S1 triplet
  //Q = + 4/3
  if( fabs(isc) == 9931551 || fabs(ism) == 9931551 || fabs(ichg) == 9931551 ) {
    _cL = - sqrt(2.)* _cL1; _cR = 0;
  }
  //Q = + 1/3
  if( fabs(isc) == 9931561 || fabs(ism) == 9931561 || fabs(ichg) == 9931561 ) {
    _cL = - _cL1; _cR = 0;
  }
  //Q = - 2/3
  if( fabs(isc) == 9931661 || fabs(ism) == 9931661 || fabs(ichg) == 9931661 ) {
    _cL = sqrt(2.) * _cL1; _cR = 0;
  }
  
  //S1/2 doublet

  //Q = + 5/3 
  if( fabs(isc) == 9941561 || fabs(ism) == 9941561 || fabs(ichg) == 9941561 ) {
    _cL = _cL12; _cR = _cR12;
  }
  
  
  //Q = + 2/3 
  if( fabs(isc) == 9941551 || fabs(ism) == 9941551 || fabs(ichg) == 9941551 ) {
    if(fabs(isc) == 5 || fabs(ism) == 5 || fabs(ichg) == 5) { 
      _cL = 0.; _cR = - _cR12;
    }
    if(fabs(isc) == 6 || fabs(ism) == 6 || fabs(ichg) == 6) { 
      _cL = 0.; _cR = _cL12;
    }
  }

  //S1/2 tilde doublet

  //Q = + 2/3 
  if( fabs(isc) == 9951551 || fabs(ism) == 9951551 || fabs(ichg) == 9951551 ) {
    _cL = _cL12t; _cR = 0.;
  }
  
  
  //Q = - 1/3 
  if( fabs(isc) == 9951651 || fabs(ism) == 9951651 || fabs(ichg) == 9951651 ) {
    _cL = _cL12t; _cR = 0.;
  }

  //dS0
  if( fabs(isc) == 9961551 || fabs(ism) == 9961551 || fabs(ichg) == 9961551) {
    if(fabs(isc) == 5 || fabs(ism) == 5 || fabs(ichg) == 5) { 
      _cL = _dcL0 * mbot +_dcR0 * mtau; _cR = _dcR0 * mbot + _dcL0 * mtau;
    }
    if(fabs(isc) == 6 || fabs(ism) == 6 || fabs(ichg) == 6) { 
      _cL = _dcR0 * mtop; _cR = 0;
    }
    _cL /= sqrt(2.) * _derivscale; 
    _cR /= sqrt(2.) * _derivscale; 
  }

  //d~S0
  if( fabs(isc) == 9971561 || fabs(ism) == 9971561 || fabs(ichg) ==  9971561) {
    _cL = _dcR0t * mtau / (sqrt(2.) * _derivscale); 
    _cR = _dcR0t * mtop / (sqrt(2.) * _derivscale); 
  }

  //dS1 triplet
  if( fabs(isc) == 9981561 || fabs(ism) == 9981561 || fabs(ichg) ==  9981561) {
    _cL = sqrt(2.) * _dcL1 * mtop / (sqrt(2.) * _derivscale);
    _cR = sqrt(2.) * _dcL1 * mtau / (sqrt(2.) * _derivscale);
  }
  if( fabs(isc) == 9981551 || fabs(ism) == 9981551 || fabs(ichg) ==  9981551) {
    if(fabs(isc) == 5 || fabs(ism) == 5 || fabs(ichg) == 5) { 
      _cL = -_dcL1 * mbot; _cR = -_dcL1 * mtau;
    }
    if(fabs(isc) == 6 || fabs(ism) == 6 || fabs(ichg) == 6) { 
      _cL = _dcL1 * mtop; _cR = 0.;
    }
    _cL /= sqrt(2.) * _derivscale; 
    _cR /= sqrt(2.) * _derivscale; 
  }

  if( fabs(isc) == 9981651 || fabs(ism) == 9981651 || fabs(ichg) ==  9981651) {
    _cL = sqrt(2.) * _dcL1 * mbot / (sqrt(2.) * _derivscale);
    _cR = 0.;
  }
  
  
  //dS1/2 doublet
  if( fabs(isc) == 9991551 || fabs(ism) == 9991551 || fabs(ichg) == 9991551 ) {
    _cL = _dcL12 * mbot + _dcR12 * mtau;
    _cR = _dcR12 * mbot + _dcL12 * mtau;
    _cL /= sqrt(2.) * _derivscale; 
    _cR /= sqrt(2.) * _derivscale; 
  }

  if( fabs(isc) == 9991561 || fabs(ism) == 9991561 || fabs(ichg) == 9991561 ) {
    if(fabs(isc) == 6 || fabs(ism) == 6 || fabs(ichg) == 6) { 
      _cL = _dcR12 * mtau; 
      _cR = _dcR12 * mtop;
    }
    if(fabs(isc) == 5 || fabs(ism) == 5 || fabs(ichg) == 5) { 
      _cL = _dcL12 * mbot;
    }      
    _cL /= sqrt(2.) * _derivscale; 
    _cR /= sqrt(2.) * _derivscale; 
  }

  //dS1/2 tilde doublet
  if( fabs(isc) == 9901561 || fabs(ism) == 9901561  || fabs(ichg) == 9901561 ) {
    _cL = _dcL12t * mtau  / (sqrt(2.) * _derivscale); _cR = _dcL12t * mtop / (sqrt(2.) * _derivscale);
  }
  
  if( fabs(isc) == 9901661 || fabs(ism) == 9901661  || fabs(ichg) == 9901661 ) {
    _cL = _dcL12t * mtop / (sqrt(2.) * _derivscale); _cR = 0.;
  }
  
  

  //this coupling selection is for the case of anti-LQ's (as dicated by the Lagrangian couplings to the quarks and leptons)
  left(_cL); right(_cR);

  int idt1, idt2;
  for(int lt = 0; lt < 9; lt++) {

    idt1 = lqtype1[lt];
    idt2 = lqtype2[lt];
    if( isc == -idt1 || ism == -idt1 || ichg == -idt1  ) {
      	left(_cR); right(_cL);
    }
 
    if( isc == idt2 || ism == idt2 || ichg == idt2  ) {
      	left(_cR); right(_cL);
    }
  }

  norm(_CFF);
}
