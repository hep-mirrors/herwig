// -*- C++ -*-
//
// SSCNWVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCNWVertex class.
//

#include "SSCNWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSCNWVertex::SSCNWVertex() : _sw(0.),  _couplast(0.), _q2last(ZERO), 
			     _id1last(0), _id2last(0), _leftlast(0.),
			     _rightlast(0.) {
  orderInGs(0);
  orderInGem(1);
}

void SSCNWVertex::doinit() {
  long neu[] = { 1000022, 1000023, 1000025, 1000035, 1000045 };
  long cha[] = { 1000024, 1000037 };
  // sign == -1 outgoing W-, sign == +1 outgoing W+
  for(int sign = -1; sign < 2; sign += 2)
    for(unsigned int ine = 0; ine < 5; ++ine)
      for(unsigned int ic = 0; ic < 2; ++ic)
	addToList(-sign*cha[ic], neu[ine], sign*24);
  FFVVertex::doinit();
  tSusyBasePtr theSS = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSCNWVertex::doinit() - The model pointer is null!"
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  
  _theN = theSS->neutralinoMix();
  _theU = theSS->charginoUMix();
  _theV = theSS->charginoVMix();

  if(!_theN || !_theU || ! _theV)
    throw InitException() << "SSCNWVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " N: " << _theN << " U: " << _theU << " V: "
			  << _theV << Exception::abortnow;
}

void SSCNWVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _theN << _theU << _theV;
}

void SSCNWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _theN >> _theU >> _theV;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SSCNWVertex,Helicity::FFVVertex>
describeSSCNWVertex("Herwig::SSCNWVertex", "HwSusy.so");

void SSCNWVertex::Init() {

  static ClassDocumentation<SSCNWVertex> documentation
    ("This class implements the coupling of a W boson to a "
     "neutralino and a chargino");

}

void SSCNWVertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
#ifndef NDEBUG
			      tcPDPtr part3) {
#else
			      tcPDPtr) {
#endif
  assert(abs(part3->id()) == ParticleID::Wplus);
  long neu, cha;
  if(part1->charged()) {
    cha = part1->id();
    neu = part2->id();
  }
  else {
    cha = part2->id();
    neu = part1->id();
  }
  assert((abs(cha) == 1000024 || abs(cha) == 1000037) && 
    (neu == 1000022 || neu == 1000023 || 
     neu == 1000025 || neu == 1000035 || 
     neu == 1000045) );
  if(q2 != _q2last||_couplast==0.) {
    _q2last = q2;
    _couplast = weakCoupling(q2);
  }
  norm(_couplast);
  if(cha != _id1last || neu != _id2last) {
    _id1last = cha;
    _id2last = neu;
    unsigned int eigc = abs(cha) == 1000037 ? 1 : 0;
    unsigned int eign(0);
    if     (neu == 1000023) eign = 1;
    else if(neu == 1000025) eign = 2;
    else if(neu == 1000035) eign = 3;
    else if(neu == 1000045) eign = 4;
    _leftlast = (*_theN)(eign, 1)*conj((*_theV)(eigc, 0)) - 
      ( (*_theN)(eign, 3)*conj((*_theV)(eigc, 1))/sqrt(2));
    _rightlast = conj((*_theN)(eign, 1))*(*_theU)(eigc, 0) +
      ( conj((*_theN)(eign, 2))*(*_theU)(eigc, 1)/sqrt(2));
  }
  Complex ltemp = _leftlast;
  Complex rtemp = _rightlast;
  // conjugate if +ve chargino
  if(cha>0) {
    ltemp = conj(ltemp);
    rtemp = conj(rtemp);
  }
  if((part1->id()==cha&&cha>0)||(part2->id()==cha&&cha<0)) {
    Complex temp = ltemp;
    ltemp  = -rtemp;
    rtemp = -temp;
  }
  left (ltemp);
  right(rtemp);
}
