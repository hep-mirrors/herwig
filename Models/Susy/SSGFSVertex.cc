// -*- C++ -*-
//
// SSGFSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGFSVertex class.
//

#include "SSGFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGFSVertex::SSGFSVertex() :_q2last(0.*GeV2),_couplast(0.), 
			    _id1last(0), _id2last(0) {
  orderInGs(1);
  orderInGem(0);
}

void SSGFSVertex::persistentOutput(PersistentOStream & os) const {
  os << _stop << _sbottom << gluinoPhase_;
}

void SSGFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _stop >> _sbottom >> gluinoPhase_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSGFSVertex,FFSVertex>
describeHerwigSSGFSVertex("Herwig::SSGFSVertex", "HwSusy.so");

void SSGFSVertex::Init() {

  static ClassDocumentation<SSGFSVertex> documentation
    ("The SSGFSVertex implements coupling of the gluinos to the "
     "squarks and quarks");

}

void SSGFSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2, tcPDPtr part3) {
  tcPDPtr ferm;
  long isc(0);
  if(abs(part1->id()) == 1000021) {
    if(part2->iSpin() == PDT::Spin1Half) {
      ferm = part2;
      isc = abs(part3->id());
    }
    else {
      ferm = part3;
      isc = abs(part2->id());
    }
  }
  else if(abs(part2->id()) == 1000021) {
    if(part1->iSpin() == PDT::Spin1Half) {
      ferm = part1;
      isc = abs(part3->id());
    }
    else {
      ferm = part3;
      isc = abs(part1->id());
    }
  }
  else if(abs(part3->id()) == 1000021) {
    if(part1->iSpin() == PDT::Spin1Half) {
      ferm = part1;
      isc = abs(part2->id());
    }
    else {
      ferm = part2;
      isc = abs(part1->id());
    }
  }
  else throw HelicityConsistencyError()
    << "SSGFSVertex::setCoupling() - There is no gluino in this vertex!"
    << part1->id() << " " << part2->id() << " " << part3->id()
    << Exception::runerror;
  long iferm = abs(ferm->id());
  assert(iferm >=1 && iferm <=6);

  if(q2 != _q2last  || _couplast==0.) {
    _couplast = -strongCoupling(q2)*sqrt(2.);
    _q2last = q2;
  }
  if(iferm != _id1last || isc != _id2last) { 
    _id1last = iferm;
    _id2last = isc;
    unsigned int eig = (isc/1000000) - 1;
    if(iferm == 6) {
      _leftlast  = -(*_stop)(eig,1)*conj(gluinoPhase_);
      _rightlast =  (*_stop)(eig,0)*     gluinoPhase_ ;
    }
    else if(iferm == 5){
      _leftlast  = -(*_sbottom)(eig,1)*conj(gluinoPhase_);
      _rightlast =  (*_sbottom)(eig,0)*     gluinoPhase_ ;
    }
    else {
      if(eig == 0) { 
	_leftlast  =  0.;
	_rightlast =  gluinoPhase_;
      }
      else {
	_leftlast  = -conj(gluinoPhase_);
	_rightlast =  0.;
      }
    }
  }
  norm(_couplast);
  //arrange l/r couplings
  if(ferm->id() < 0) {
    left(conj(_rightlast));
    right(conj(_leftlast));
  }
  else {
    left(_leftlast);
    right(_rightlast);
  }
}

void SSGFSVertex::doinit() {
  for(long ix=1;ix<7;++ix) {
    addToList(1000021, ix, -(ix+1000000));
    addToList(1000021, ix, -(ix+2000000));
    addToList(1000021, -ix, (ix+1000000));
    addToList(1000021, -ix, (ix+2000000));
  }
  FFSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());

  _stop = model->stopMix();
  _sbottom = model->sbottomMix();
  gluinoPhase_ = model->gluinoPhase();
  if(!_stop || !_sbottom)
    throw InitException() << "SSGFSVertex::doinit() - "
			  << "There is a null mixing matrix pointer. "
			  << "stop: " << _stop << " sbottom: " << _sbottom 
			  << Exception::abortnow;
}
