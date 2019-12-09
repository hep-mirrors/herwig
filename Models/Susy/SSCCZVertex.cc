// -*- C++ -*-
//
// SSCCZVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCCZVertex class.
//

#include "SSCCZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSCCZVertex::SSCCZVertex() : _sw2(0.), _cw(0.), _couplast(0.),
			     _q2last(), _id1last(0), _id2last(0),
			     _leftlast(0.), _rightlast(0.), _gblast(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

void SSCCZVertex::doinit() {
  addToList(-1000024, 1000024, 23);
  addToList(-1000024, 1000037, 23);

  addToList(-1000037, 1000024, 23);
  addToList(-1000037, 1000037, 23);

  //photon
  addToList(-1000024, 1000024, 22);
  addToList(-1000037, 1000037, 22);
  FFVVertex::doinit();
  tSusyBasePtr theSS = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!theSS) 
    throw InitException() << "SSCCZVertex::doinit - The model pointer "
				     << "is null! "
				     << Exception::abortnow;
  _sw2 = sin2ThetaW();
  _cw = sqrt(1. - _sw2);
  _theU = theSS->charginoUMix();
  _theV = theSS->charginoVMix();
  if(!_theU || !_theV)
    throw InitException() << "SSCCZVertex::doinit - "
			  << "A mixing matrix pointer is null.  U: " 
			  << _theU << "  V: " << _theV
			  << Exception::abortnow;
}

void SSCCZVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw2 << _cw << _theU << _theV;
}

void SSCCZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw2 >> _cw >> _theU >> _theV;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSCCZVertex,Helicity::FFVVertex>
describeSSCCZVertex("Herwig::SSCCZVertex", "HwSusy.so");

void SSCCZVertex::Init() {

  static ClassDocumentation<SSCCZVertex> documentation
    ("This class implements the coupling of a Z/gamma to a pair of"
     " charginos. ");

}

void SSCCZVertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
			      tcPDPtr part3) {
  long ichar1(part1->id()), ichar2(part2->id()), boson(part3->id());
  assert( boson == ParticleID::gamma || boson == ParticleID::Z0);
  assert( abs(ichar1) == 1000024 || abs(ichar1) == 1000037); 
  assert( abs(ichar2) == 1000024 || abs(ichar2) == 1000037); 
  if(_q2last != q2||_couplast==0.) {
    _q2last = q2;
    _couplast = electroMagneticCoupling(q2);
  }
  norm(_couplast);
  if(boson != _gblast || ichar1 != _id1last || ichar2 != _id2last) {
    _gblast = boson;
    _id1last = ichar1;
    _id2last = ichar2;
    if( boson == ParticleID::Z0 ) {
      unsigned int ic1(0), ic2(0);
      if(abs(ichar1) == 1000037) ic1 = 1;
      if(abs(ichar2) == 1000037) ic2 = 1;
      _leftlast = -(*_theV)(ic1, 0)*conj((*_theV)(ic2, 0)) - 
	0.5*(*_theV)(ic1, 1)*conj((*_theV)(ic2, 1));
      _rightlast = -conj((*_theU)(ic1, 0))*(*_theU)(ic2, 0) - 
	0.5*conj((*_theU)(ic1, 1))*(*_theU)(ic2, 1);
      if(abs(ichar1) == abs(ichar2)) {
	_leftlast += _sw2;
	_rightlast += _sw2;
      }
      _leftlast /= sqrt(_sw2)*_cw;
      _rightlast /= sqrt(_sw2)*_cw;
    }
    else {
      if(abs(ichar1) == abs(ichar2)) {
	_leftlast  = -1.;
	_rightlast = -1.;
      }
      else {
	_leftlast  = 0.;
	_rightlast = 0.;
      }
    }
    if(ichar1>0) {
      Complex temp = _leftlast;
      _leftlast  = -_rightlast;
      _rightlast = -temp;
    }
  }
  left(_leftlast);
  right(_rightlast);
}
