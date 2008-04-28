// -*- C++ -*-
//
// SSCCZVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCCZVertex class.
//

#include "SSCCZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSCCZVertex::SSCCZVertex() : _sw2(0.), _cw(0.), _couplast(0.),
			     _q2last(), _id1last(0), _id2last(0),
			     _leftlast(0.), _rightlast(0.), _gblast(0){
  vector<int> first, second, third;
  for(unsigned int ix = 0; ix < 2; ++ix) {
    int ic1(1000024);
    if(ix == 1) ic1 = 1000037;
    for(unsigned int iy = 0; iy < 2; ++iy) {
      int ic2(1000024);
      if(iy == 1) ic2 = 1000037;
      first.push_back(-ic1);
      second.push_back(ic2);
      third.push_back(23);
    }
  }
  //photon
  first.push_back(-1000024);
  second.push_back(1000024);
  third.push_back(22);
  first.push_back(-1000037);
  second.push_back(1000037);
  third.push_back(22);
  setList(first, second, third);
}

void SSCCZVertex::doinit() throw(InitException) {
  FFVVertex::doinit();
  tSusyBasePtr theSS = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!theSS) 
    throw InitException() << "SSCCZVertex::doinit - The model pointer "
				     << "is null! "
				     << Exception::abortnow;
  _sw2 = theSS->sin2ThetaW();
  _cw = sqrt(1. - _sw2);
  _theU = theSS->charginoUMix();
  _theV = theSS->charginoVMix();
  if(!_theU || !_theV)
    throw InitException() << "SSCCZVertex::doinit - "
			  << "A mixing matrix pointer is null.  U: " 
			  << _theU << "  V: " << _theV
			  << Exception::abortnow;
  orderInGs(0);
  orderInGem(1);
}

void SSCCZVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw2 << _cw << _theU << _theV;
}

void SSCCZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw2 >> _cw >> _theU >> _theV;
}

ClassDescription<SSCCZVertex> SSCCZVertex::initSSCCZVertex;
// Definition of the static class description member.

void SSCCZVertex::Init() {

  static ClassDocumentation<SSCCZVertex> documentation
    ("This class implements the coupling of a Z/gamma to a pair of"
     " charginos. ");

}

void SSCCZVertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
			      tcPDPtr part3) {
  long ichar1(abs(part1->id())), ichar2(abs(part2->id())),
    boson(part3->id());
  if( (boson != ParticleID::gamma && boson != ParticleID::Z0) ||
      (ichar1 != 1000024 && ichar1 != 1000037) ||
      (ichar2 != 1000024 && ichar2 != 1000037) ) {
    throw HelicityConsistencyError() 
      << "SSCCZVertex::setCoupling() - An incorrect particle has been found. "
      << part1->id() << " " << part2->id() << " " << part3->id()
      << Exception::warning;
    setNorm(0.); setLeft(0.), setRight(0.);
    return;
  }
  if(_q2last != q2) {
    _q2last = q2;
    _couplast = electroMagneticCoupling(q2);
  }
  setNorm(_couplast);
  if(boson != _gblast || ichar1 != _id1last || ichar2 != _id2last) {
    _gblast = boson;
    _id1last = ichar1;
    _id2last = ichar2;
    if( boson == ParticleID::Z0 ){
      unsigned int ic1(0), ic2(0);
      if(ichar1 == 1000037) ic1 = 1;
      if(ichar2 == 1000037) ic2 = 1;
      _leftlast = -(*_theV)(ic1, 0)*conj((*_theV)(ic2, 0)) - 
	0.5*(*_theV)(ic1, 1)*conj((*_theV)(ic2, 1));
      _rightlast = -conj((*_theU)(ic1, 0))*(*_theU)(ic2, 0) - 
	0.5*conj((*_theU)(ic1, 1))*(*_theU)(ic2, 1);
      if(ichar1 == ichar2) {
	_leftlast += _sw2;
	_rightlast += _sw2;
      }
      _leftlast /= sqrt(_sw2)*_cw;
      _rightlast /= sqrt(_sw2)*_cw;
    }
    else {
      _leftlast = -1.;
      _rightlast = -1.;
    }
  }
  setLeft(_leftlast);
  setRight(_rightlast);
}
