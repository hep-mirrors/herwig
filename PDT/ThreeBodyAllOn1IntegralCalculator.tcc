// -*- C++ -*-
//
// ThreeBodyAllOn1IntegralCalculator.tcc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the ThreeBodyAllOn1IntegralCalculator class.
//

namespace Herwig{
using namespace ThePEG;

template <class T> 
Energy ThreeBodyAllOn1IntegralCalculator<T>::operator() (double x) const {
  Energy2 scale;
  if(_intmass>ZERO) scale = _intmass*(_intmass+_intwidth*tan(x));
  else           scale = UnitRemoval::E2 * pow(x,1./(_intpower + 1.));
  InvEnergy output=_theDgamma.threeBodydGammads(_mode,_m2[0],scale,_m[1],_m[2],_m[3]);
  // the jacobian
  InvEnergy2 term;
  Energy2 rm2,rw2;
  if(_intmass>ZERO) {
    rm2 = sqr(_intmass);
    rw2 = sqr(_intwidth);
    term = _intmass*_intwidth / (sqr(scale-rm2) + rw2 * rm2);
  }
  else {
    term = UnitRemoval::InvE2 * (_intpower+1.)*pow(scale * UnitRemoval::InvE2, _intpower);
  }
  return output/term;
}

template <class T>
Energy ThreeBodyAllOn1IntegralCalculator<T>::partialWidth(Energy2 scale) const {
  _m2[0]=scale;
  _m[0]=sqrt(scale);
  // limits for the outer integral
  Energy2 upp=ZERO,low=ZERO;
  switch(_variabletype) {
  case 1:
    upp = (_m[0]-_m[3])*(_m[0]-_m[3]);
    low = (_m[1]+_m[2])*(_m[1]+_m[2]);
    break;
  case 2:
    upp = (_m[0]-_m[2])*(_m[0]-_m[2]);
    low = (_m[1]+_m[3])*(_m[1]+_m[3]);
    break;
  case 3:
    upp = (_m[0]-_m[1])*(_m[0]-_m[1]);
    low = (_m[2]+_m[3])*(_m[2]+_m[3]);
    break;
  }
  double rupp, rlow;
  // transform them
  if(_intmass>ZERO) {
    rupp = atan2((upp-_intmass*_intmass), _intmass*_intwidth);
    rlow = atan2((low-_intmass*_intmass), _intmass*_intwidth);
  }
  else {
    rupp = pow(upp * UnitRemoval::InvE2, _intpower+1.);
    rlow = pow(low * UnitRemoval::InvE2, _intpower+1.);
  }
  return _integrator.value(*this,rlow,rupp);
}

}

