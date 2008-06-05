// -*- C++ -*-
//
// SSNNZVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNNZVertex class.
//

#include "SSNNZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/Susy/MixingMatrix.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSNNZVertex::SSNNZVertex() : _sw(0.), _cw(0.), _id1last(0), 
			     _id2last(0), _q2last(), _couplast(0.),
			     _leftlast(0.), _rightlast(0.) {
  vector<long> first, second, third(25, 23);
  for(unsigned int i = 0; i < 5; ++i) {
    long neu1;
    if(i == 0) neu1 = 1000022;
    else if(i == 1) neu1 = 1000023;
    else if(i == 2) neu1 = 1000025;
    else if(i == 3)  neu1 = 1000035;
    else neu1 = 1000045;
    for(unsigned int j = 0; j < 5; ++j) {
      long neu2;
      if(j == 0) neu2 = 1000022;
      else if(j == 1) neu2 = 1000023;
      else if(j == 2) neu2 = 1000025;
      else if(j == 3)  neu2 = 1000035;
      else neu2 = 1000045;
      first.push_back(neu1);
      second.push_back(neu2);		      
    }
  }
  setList(first, second, third);
}

void SSNNZVertex::doinit() throw(InitException) {
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
  _sw = sqrt(theSS->sin2ThetaW());
  _cw = sqrt(1 - _sw*_sw);
  orderInGem(1);
  orderInGs(0);
}

void SSNNZVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _cw << _theN;
}

void SSNNZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _theN;
  _id1last = 0;
  _id2last = 0;
  _q2last = 0.*GeV2;
  _couplast = 0.;
  _leftlast = 0.;
  _rightlast = 0.;
}

ClassDescription<SSNNZVertex> SSNNZVertex::initSSNNZVertex;
// Definition of the static class description member.

void SSNNZVertex::Init() {

  static ClassDocumentation<SSNNZVertex> documentation
    ("The coupling of a Z-boson to a pair of neutralinos");

}

void SSNNZVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long ic1(0), ic2(0);
  if(part1->id() == ParticleID::Z0) {
    ic1 = part2->id();
    ic2 = part3->id();
  }
  else if(part2->id() == ParticleID::Z0) {
    ic1 = part1->id();
    ic2 = part3->id();
  }
  else if(part3->id() == ParticleID::Z0) {
    ic1 = part1->id();
    ic2 = part2->id();
  }
  else
    throw HelicityConsistencyError() << "There is no Z0 in the ZNNVertex!"
				     << Exception::warning;
  if(ic1 == ParticleID::SUSY_chi_10 || ic1 == ParticleID::SUSY_chi_20 ||
     ic1 == ParticleID::SUSY_chi_30 || ic1 == ParticleID::SUSY_chi_40 || 
     ic2 == ParticleID::SUSY_chi_10 || ic2 == ParticleID::SUSY_chi_20 ||
     ic2 == ParticleID::SUSY_chi_30 || ic2 == ParticleID::SUSY_chi_40 ||
     ic1 == 1000045                 || ic2 == 1000045                ) {
    if(q2 != _q2last) {
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
    setNorm(_couplast);
    setLeft(_leftlast);
    setRight(_rightlast);
  }
  else 
    throw HelicityConsistencyError() << "A particle other than a Z0 " 
				     << "or a ~chi_i0 exists in the "
				     << "SSNNZVertex " << ic1 << " "
				     << ic2 << Exception::warning;
	
    
  
  
}

