// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeBodyAllOn1IntegralCalculator class.
//

#include "ThreeBodyAllOn1IntegralCalculator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreeBodyAllOn1IntegralCalculator.tcc"
#endif

namespace Herwig{
using namespace ThePEG;

ThreeBodyAllOn1IntegralCalculator::~ThreeBodyAllOn1IntegralCalculator() 
{
  delete _theIntegrand;
  delete _Integrator;
}

Energy ThreeBodyAllOn1IntegralCalculator::integrand(double x) {
  Energy2 scale;
  if(_intmass>0)
    {scale = _intmass*(_intmass+_intwidth*tan(x));}
  else
    {scale = pow(x,1./(_intwidth+1.));}
  Genfun::Argument a(7); a[0]=_m2[0];a[1]=scale;a[2]=_m[1];a[3]=_m[2];a[4]=_m[3];
  Energy output=(*_theDgamma)(a);
  // the jacobian
  double term;
  Energy2 rm2,rw2;
  if(_intmass>0)
    {
      rm2=_intmass*_intmass;
      rw2 = _intwidth*_intwidth;
      term = (scale-rm2)*(scale-rm2)+rw2*rm2;
      term = _intmass*_intwidth/term;
    }
  else
    {term = (_intwidth+1.)*pow(scale,_intwidth);}
  return output/term;
}

Energy ThreeBodyAllOn1IntegralCalculator::partialWidth(Energy2 scale) const
 {
   _m2[0]=scale;_m[0]=sqrt(scale);
   // limits for the outer integral
   Energy2 upp=0.,low=0.;
   switch(_variabletype)
     {
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
   if(_intmass>0)
     {
       upp = atan((upp-_intmass*_intmass)/
		  _intmass/_intwidth);
       low =  atan((low-_intmass*_intmass)/
		   _intmass/_intwidth);
     }
   else
     {
       upp = pow(upp,_intwidth+1.);
       low = pow(low,_intwidth+1.);
     }
   _Integrator->resetLimits(low,upp);
  return (*_Integrator)[*_theIntegrand];
}
}

// the integrand
namespace Herwig {
using namespace ThePEG;
using namespace Genfun;

// inner integral
FUNCTION_OBJECT_IMP(ThreeBodyAllOn1IntegralOuter)
    
ThreeBodyAllOn1IntegralOuter::ThreeBodyAllOn1IntegralOuter(ThreeBodyAllOn1IntegralCalculatorPtr in)
{_theIntegrator=in;}

ThreeBodyAllOn1IntegralOuter::~ThreeBodyAllOn1IntegralOuter() {}
  
ThreeBodyAllOn1IntegralOuter::ThreeBodyAllOn1IntegralOuter(const ThreeBodyAllOn1IntegralOuter & right) 
{  }
  
double ThreeBodyAllOn1IntegralOuter::operator() (double x) const 
{
  // calculate the integrand
  double output=_theIntegrator->integrand(x);
  return output;
}

}
