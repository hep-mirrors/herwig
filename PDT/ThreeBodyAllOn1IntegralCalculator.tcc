// -*- C++ -*-
//
// This is the implementation of the non-inlined templated member
// functions of the ThreeBodyAllOn1IntegralCalculator class.
//

namespace Herwig{
using namespace ThePEG;

template <class T> 
double ThreeBodyAllOn1IntegralCalculator<T>::operator() (double x) const {
  Energy2 scale;
  if(_intmass>0) scale = _intmass*(_intmass+_intwidth*tan(x));
  else           scale = pow(x,1./(_intwidth+1.));
  Energy output=_theDgamma.threeBodydGammads(_mode,_m2[0],scale,_m[1],_m[2],_m[3]);
  // the jacobian
  double term;
  Energy2 rm2,rw2;
  if(_intmass>0) {
    rm2=_intmass*_intmass;
    rw2 = _intwidth*_intwidth;
    term = (scale-rm2)*(scale-rm2)+rw2*rm2;
    term = _intmass*_intwidth/term;
  }
  else {
    term = (_intwidth+1.)*pow(scale,_intwidth);
  }
  return output/term;
}

template <class T>
Energy ThreeBodyAllOn1IntegralCalculator<T>::partialWidth(Energy2 scale) const {
  _m2[0]=scale;
  _m[0]=sqrt(scale);
  // limits for the outer integral
  Energy2 upp=0.,low=0.;
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
  // transform them
  if(_intmass>0) {
    upp = atan((upp-_intmass*_intmass)/_intmass/_intwidth);
    low = atan((low-_intmass*_intmass)/_intmass/_intwidth);
  }
  else {
    upp = pow(upp,_intwidth+1.);
    low = pow(low,_intwidth+1.);
  }
  return _integrator.value(*this,low,upp);
}

}

