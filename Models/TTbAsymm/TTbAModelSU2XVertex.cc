// -*- C++ -*-
//
// TTbAModelSU2XVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModelSU2XVertex class.
//

#include "TTbAModelSU2XVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;



IBPtr TTbAModelSU2XVertex::clone() const {
  return new_ptr(*this);
}

IBPtr TTbAModelSU2XVertex::fullclone() const {
  return new_ptr(*this);
}

TTbAModelSU2XVertex::TTbAModelSU2XVertex()  {
  orderInGem(1);
  orderInGs(1);
  
  addToList(-6,6,70);
  addToList(-2,2,70);
  addToList(-2,6,70);
  addToList(2,-6,70);

  addToList(-6,6,71);
  addToList(-2,2,71);
  addToList(-2,6,71);
  addToList(2,-6,71);

   
  addToList(-6,6,-71);
  addToList(-2,2,-71);
  addToList(-2,6,-71);
  addToList(2,-6,-71);


}

void TTbAModelSU2XVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwTTbAPtr hwTTbA=dynamic_ptr_cast<tcHwTTbAPtr>(_theModel);
  if(hwTTbA) {
    _alphaX =hwTTbA->_alphaX_value();
    _costhetaX =hwTTbA->_costhetaX_value();
    _models =hwTTbA->_model();

  }
  FFVVertex::doinit();
}

void TTbAModelSU2XVertex::persistentOutput(PersistentOStream & os) const {
  os << _alphaX << _costhetaX << _models;
}

void TTbAModelSU2XVertex::persistentInput(PersistentIStream & is, int) {
  is >> _alphaX >> _costhetaX >> _models;
}

ClassDescription<TTbAModelSU2XVertex> 
TTbAModelSU2XVertex::initTTbAModelSU2XVertex;
// Definition of the static class description member.


void TTbAModelSU2XVertex::Init() {
  
  static ClassDocumentation<TTbAModelSU2XVertex> documentation
    ("The TTbAModelSU2XVertex class is the implementation"
     " of the helicity amplitude calculation of the TTbA"
     " SU(2)_X vertex.");
}

void TTbAModelSU2XVertex::setCoupling(Energy2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {

  double _cR = 0, _fac = 0;
  _gX = sqrt( 4 * Constants::pi * _alphaX ); 
  double ct = _costhetaX;
  double st = sqrt(1 - pow(ct,2)); 

  //Vz
  if( abs(aa->id()) == 70 || abs(bb->id()) == 70 || abs(cc->id()) == 70) { 
    _fac = _gX / 2.0; 
    if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR = pow(ct,2) - pow(st,2);
      }
    }
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
       if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
	_cR = pow(st,2) - pow(ct,2);
      }
    }
    if( abs(aa->id()) == 2 || abs(bb->id()) == 2 || abs(cc->id()) == 2) {
      if( abs(aa->id()) == 6 || abs(bb->id()) == 6 || abs(cc->id()) == 6) {
	_cR = 2 * ct * st;
      }
    }
  }
  
  //Ym
  if( aa->id() == -71 || bb->id() == -71 || cc->id() == -71 ) { 
    
    _fac = _gX / sqrt(2.0); 
    
    if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR = - ct * st;
      }
    }
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
      if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
	_cR = ct * st;
      }
    }    
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR = - pow(st,2); 
      }
    }
    if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
      if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
	_cR = pow(ct,2);
      }
    }   
  }


  //Yp 
  if( aa->id() == 71 || bb->id() == 71 || cc->id() == 71 ) { 

    _fac = _gX / sqrt(2.0); 

    if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR =  - ct * st;
      }
    }
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
       if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
	_cR = ct * st;
      }
    }
    
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR = pow(ct,2);
      }
    }
    if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
      if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
	_cR = - pow(st,2);
      }
    } 
  }

  //normalise according to Lagrangian factor
  _cR *= _fac;


  //If this model is not selected set coupling to zero.
  if(_models!=3) { _cR = 1E-10; }
 

  right(_cR);
  left(0.);
 
  norm(1.0);

}
