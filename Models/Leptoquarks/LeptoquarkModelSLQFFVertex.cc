// -*- C++ -*-
//
// LeptoquarkModelSLQFFVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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

LeptoquarkModelSLQFFVertex::LeptoquarkModelSLQFFVertex() {
  orderInGem(1);
  orderInGs(0);
}

void LeptoquarkModelSLQFFVertex::doinit() {
  //S0
  addToList( 15, 6,-9911561);
  addToList(-15,-6, 9911561);
  
  addToList(-16,-5, 9911561);
  addToList( 16, 5,-9911561);

  //~S0
  addToList(-15,-5, 9921551);
  addToList( 15, 5,-9921551);

  //S1 triplet
  //S1p
  addToList(-15,-5, 9931551);
  addToList( 15, 5,-9931551);
  //S1z
  addToList(-15,-6, 9931561);
  addToList( 15, 6,-9931561);
  addToList(-16,-5, 9931561);
  addToList( 16, 5,-9931561);
  //S1m
  addToList(-16,-6, 9931661);
  addToList( 16, 6,-9931661);

  //S1/2 doublet
  addToList( 15,-6, 9941561);
  addToList(-15, 6,-9941561);
  
  addToList(-15, 5,-9941551);
  addToList(-16, 6,-9941551);
  addToList( 15,-5, 9941551);
  addToList( 16,-6, 9941551);


  //S1/2 tilde doublet
  addToList( 5,-16,-9951651);
  addToList(-5, 16, 9951651);

  addToList(-5, 15, 9951551);
  addToList( 5,-15,-9951551);


  //dS0
  addToList( 15,-5, 9961551);
  addToList(-15, 5,-9961551);

  addToList( 16,-6, 9961551);
  addToList(-16, 6,-9961551);

  //~dS0
  addToList( 15,-6, 9971561);
  addToList(-15, 6,-9971561);


  //dS1 triplet

  //dS1p
  addToList( 15,-6, 9981561);
  addToList(-15, 6,-9981561);

  //dS1z
  addToList( 16,-6, 9981551);
  addToList(-16, 6,-9981551);

  addToList( 15,-5, 9981551);
  addToList(-15, 5,-9981551);

  //dS1m
  addToList( 16,-5, 9981651);
  addToList(-16, 5,-9981651);

  //dS1/2 doublet
  addToList(-15,-5, 9991551);
  addToList( 15, 5,-9991551);

  addToList(-15,-6, 9991561);
  addToList( 15, 6,-9991561);

  addToList(-16,-5, 9991561);
  addToList( 16, 5,-9991561);

  //dS1/2 tilde doublet
  addToList(-15,-6, 9901561);
  addToList( 15, 6,-9901561);

  addToList(-16,-6, 9901661);
  addToList( 16, 6,-9901661);


  _theModel = generator()->standardModel();
  tcHwLeptoquarkPtr hwLeptoquark=dynamic_ptr_cast<tcHwLeptoquarkPtr>(_theModel);
  if(hwLeptoquark){
    _CFF=hwLeptoquark->cfermion();
    _cL0 =hwLeptoquark->cleft();
    _cR0 =hwLeptoquark->cright();
    _cR0t = hwLeptoquark->crighttilde();
    _cL1 =hwLeptoquark->cleft1(); 
    _cL12 =hwLeptoquark->cleft12(); 
    _cR12 =hwLeptoquark->cright12(); 
    _cL12t =hwLeptoquark->cleft12tilde(); 
    
    
    _derivscale = hwLeptoquark->fscale();
    _dcL0 =hwLeptoquark->dcleft();
    _dcR0 =hwLeptoquark->dcright();
    _dcR0t = hwLeptoquark->dcrighttilde();
    _dcL1 =hwLeptoquark->dcleft1(); 
    _dcL12 =hwLeptoquark->dcleft12(); 
    _dcR12 =hwLeptoquark->dcright12(); 
    _dcL12t =hwLeptoquark->dcleft12tilde(); 
    
  }
  FFSVertex::doinit();
}

void LeptoquarkModelSLQFFVertex::persistentOutput(PersistentOStream & os) const {
  os << _CFF << _cL0 << _cR0 << _cR0t 
     << _cL1 << _cL12 << _cR12 << _cL12t 
     << _dcL0 << _dcR0 << _dcR0t 
     << _dcL1 << _dcL12 << _dcR12 << _dcL12t 
     << ounit(_derivscale,GeV);
}

void LeptoquarkModelSLQFFVertex::persistentInput(PersistentIStream & is, int) {
  is >> _CFF >> _cL0 >> _cR0 >> _cR0t 
     >> _cL1 >> _cL12 >> _cR12 >> _cL12t 
     >>_dcL0 >> _dcR0 >> _dcR0t 
     >> _dcL1 >> _dcL12 >> _dcR12 >> _dcL12t 
     >> iunit(_derivscale,GeV);
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
  long lqid = isc;
  long smid_1 = ism;
  long smid_2 = ichg;
  if(abs(lqid) < 9900000) { 
    lqid = ism; 
    smid_1 = ichg;
    smid_2 = isc; 
  }
  if(abs(lqid) < 9900000) {
    smid_1 = ism;
    smid_2 = isc;
  }
  if( abs(smid_1) > abs(smid_2) ) { swap(smid_1, smid_2); }

  const Energy denom = sqrt(2.) * _derivscale;
  const double mtop = getParticleData(ParticleID::t)->mass()        / denom;
  const double mbot = getParticleData(ParticleID::b)->mass()        / denom;
  const double mtau = getParticleData(ParticleID::tauminus)->mass() / denom;

  //set the couplings to left and right 
  //S0
  if( abs(isc) == 9911561 || abs(ism) == 9911561 || abs(ichg) == 9911561 ) {
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cL = -_cL0; _cR = Complex(0.);
    }
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cL = _cL0;
      _cR = _cR0;
    }
  }
  //~S0
  if( abs(isc) == 9921551 || abs(ism) == 9921551 || abs(ichg) == 9921551 ) {
    _cL = Complex(0.); _cR = _cR0t;
  }
  
  //S1 triplet
  //Q = + 4/3
  if( abs(isc) == 9931551 || abs(ism) == 9931551 || abs(ichg) == 9931551 ) {
    _cL = sqrt(2.)* _cL1; _cR = Complex(0.);
  }
  //Q = + 1/3
  if( abs(isc) == 9931561 || abs(ism) == 9931561 || abs(ichg) == 9931561 ) {
    _cL = - _cL1; _cR = Complex(0.);
  }
  //Q = - 2/3
  if( abs(isc) == 9931661 || abs(ism) == 9931661 || abs(ichg) == 9931661 ) {
    _cL = sqrt(2.) * _cL1; _cR = Complex(0.);
  }
  
  //S1/2 doublet

  //Q = + 5/3 
  if( abs(isc) == 9941561 || abs(ism) == 9941561 || abs(ichg) == 9941561 ) {
    _cR = _cL12; _cL = _cR12;
  }
  
  
  //Q = + 2/3 
  if( abs(isc) == 9941551 || abs(ism) == 9941551 || abs(ichg) == 9941551 ) {
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cR = Complex(0.); _cL = - _cR12;
    }
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cL = Complex(0.); _cR = _cL12;
    }
  }

  //S1/2 tilde doublet

  //Q = + 2/3 
  if( abs(isc) == 9951551 || abs(ism) == 9951551 || abs(ichg) == 9951551 ) {
    _cR = _cL12t; _cL = Complex(0.);
  }
  
  
  //Q = - 1/3 
  if( abs(isc) == 9951651 || abs(ism) == 9951651 || abs(ichg) == 9951651 ) {
    _cR = _cL12t; _cL = Complex(0.);
  }



  //dS0
  if( abs(isc) == 9961551 || abs(ism) == 9961551 || abs(ichg) == 9961551) {
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cR = _dcL0 * mbot +_dcR0 * mtau; 
      _cL = _dcR0 * mbot + _dcL0 * mtau;
    }
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cR = _dcL0 * mtop; 
      _cL = Complex(0.);
    }
  }

  //d~S0
  if( abs(isc) == 9971561 || abs(ism) == 9971561 || abs(ichg) ==  9971561) {
    _cR = _dcR0t * mtau;
    _cL = _dcR0t * mtop;
  }

  //dS1 triplet
  if( abs(isc) == 9981561 || abs(ism) == 9981561 || abs(ichg) ==  9981561) {
    _cR = sqrt(2.) * _dcL1 * mtop;
    _cL = sqrt(2.) * _dcL1 * mtau;
  }
  if( abs(isc) == 9981551 || abs(ism) == 9981551 || abs(ichg) ==  9981551) {
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cR = -_dcL1 * mbot; 
      _cL = -_dcL1 * mtau;
    }
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cR = _dcL1 * mtop;
      _cL = Complex(0.);
    }
  }

  if( abs(isc) == 9981651 || abs(ism) == 9981651 || abs(ichg) ==  9981651) {
    _cL = sqrt(2.) * _dcL1 * mbot;
    _cR = Complex(0.);
  }
  
  
  //dS1/2 doublet
  if( abs(isc) == 9991551 || abs(ism) == 9991551 || abs(ichg) == 9991551 ) {
    _cL = _dcL12 * mbot + _dcR12 * mtau;
    _cR = _dcR12 * mbot + _dcL12 * mtau;
  }

  if( abs(isc) == 9991561 || abs(ism) == 9991561 || abs(ichg) == 9991561 ) {
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cL = _dcR12 * mtau; 
      _cR = _dcR12 * mtop;
    }
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cL = _dcL12 * mbot;
    }    
  }

  //dS1/2 tilde doublet
  if( abs(isc) == 9901561 || abs(ism) == 9901561  || abs(ichg) == 9901561 ) {
    _cL = _dcL12t * mtop; 
    _cR = _dcL12t * mtau;
  }
  
  if( abs(isc) == 9901661 || abs(ism) == 9901661  || abs(ichg) == 9901661 ) {
    _cL = _dcL12t * mtop;
    _cR = Complex(0.);
  }


  if(smid_1 > 0) { 
    left(conj(_cR)); 
    right(conj(_cL));
  } 
  else { 
    left(_cL); 
    right(_cR); 
  } 

  norm(_CFF);
}
