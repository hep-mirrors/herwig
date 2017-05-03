// -*- C++ -*-
//
// SSNNZVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNNZVertex class.
//

#include "SSNNZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "MixingMatrix.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSNNZVertex::SSNNZVertex() : _sw(0.), _cw(0.), _id1last(0), 
			     _id2last(0), _q2last(), _couplast(0.),
			     _leftlast(0.), _rightlast(0.) {
  orderInGem(1);
  orderInGs(0);
}

void SSNNZVertex::doinit() {
  long neu[] = { 1000022, 1000023, 1000025, 1000035, 1000045 };
  for(unsigned int i = 0; i < 5; ++i)
    for(unsigned int j = 0; j < 5; ++j)
      addToList(neu[i], neu[j], 23);
  FFVVertex::doinit();
  tSusyBasePtr theSS = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSNNZVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  
  _theN  = theSS->neutralinoMix();
  if(!_theN)
    throw InitException() << "SSNNZVertex::doinit - The neutralino "
			  << "mixing matrix pointer is null." 
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1 - _sw*_sw);
}

void SSNNZVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _cw << _theN;
}

void SSNNZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _theN;
  _id1last = 0;
  _id2last = 0;
  _q2last = ZERO;
  _couplast = 0.;
  _leftlast = 0.;
  _rightlast = 0.;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SSNNZVertex,Helicity::FFVVertex>
describeSSNNZVertex("Herwig::SSNNZVertex", "HwSusy.so");

void SSNNZVertex::Init() {

  static ClassDocumentation<SSNNZVertex> documentation
    ("The coupling of a Z-boson to a pair of neutralinos");

}

void SSNNZVertex::setCoupling(Energy2 q2,tcPDPtr part1,
#ifndef NDEBUG
			      tcPDPtr part2,tcPDPtr part3) {
#else
			      tcPDPtr part2,tcPDPtr) {
#endif
  assert(part3->id() == ParticleID::Z0);
  long ic1 = part2->id();
  long ic2 = part1->id();
  assert(ic1 == ParticleID::SUSY_chi_10 || ic1 == ParticleID::SUSY_chi_20 ||
	 ic1 == ParticleID::SUSY_chi_30 || ic1 == ParticleID::SUSY_chi_40 ||
	 ic1 == 1000045);
  assert(ic2 == ParticleID::SUSY_chi_10 || ic2 == ParticleID::SUSY_chi_20 ||
	 ic2 == ParticleID::SUSY_chi_30 || ic2 == ParticleID::SUSY_chi_40 ||
	 ic2 == 1000045);
  if(q2 != _q2last || _couplast==0.) {
    _q2last = q2;
    _couplast = weakCoupling(q2)/_cw;
  }
  if(ic1 != _id1last || ic2 != _id2last) {
    _id1last = ic1;
    _id2last = ic2;
    unsigned int neu1(ic1 - 1000022), neu2(ic2 - 1000022);
    if(neu1 > 1) {
      if(ic1 == 1000025)
	neu1 = 2;
      else if(ic1 == 1000035)
	neu1 = 3;
      else 
	neu1 = 4;
    }
    if(neu2 > 1) {
      if(ic2 == 1000025)
	neu2 = 2;
      else if(ic2 == 1000035)
	neu2 = 3;
      else
	neu2 = 4;
    }
    _leftlast = 0.5*( (*_theN)(neu1, 3)*conj((*_theN)(neu2, 3)) -
		      (*_theN)(neu1, 2)*conj((*_theN)(neu2, 2)) );
    _rightlast = -conj(_leftlast);
  }
  norm(_couplast);
  left(_leftlast);
  right(_rightlast);
}
