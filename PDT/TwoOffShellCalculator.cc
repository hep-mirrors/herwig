// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoOffShellCalculator class.
//

#include "TwoOffShellCalculator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TwoOffShellCalculator.tcc"
#endif


namespace Herwig {
using namespace ThePEG;

TwoOffShellCalculator::~TwoOffShellCalculator() 
{
  delete _Integrator;
  delete _integrand;
}

NoPIOClassDescription<TwoOffShellCalculator> TwoOffShellCalculator::initTwoOffShellCalculator;
// Definition of the static class description member.

void TwoOffShellCalculator::Init() {

  static ClassDocumentation<TwoOffShellCalculator> documentation
    ("The TwoOffShellCalculator class performs the integral of the partial width"
     " when two decay products are off-shell");

}

// calculate the width for a given mass
Energy TwoOffShellCalculator::partialWidth(Energy2 q2) const
{
  _scale=q2;
  // the limits
  Energy upp=min(sqrt(q2)-_mother,_massptr->upperLimit());
  Energy low=max(_minmass,_massptr->lowerLimit());
  if(low>upp){return 0.;}
  // transform the limits of BW smoothing
  Energy2 mass2  =_massptr->nominalMass()*_massptr->nominalMass();
  Energy2 mwidth =_massptr->nominalMass()*_massptr->nominalWidth();
  double  rhomin=atan((low*low-mass2)/mwidth);
  double  rhomax=atan((upp*upp-mass2)/mwidth);
  _Integrator->resetLimits(rhomin,rhomax);
  return (*_Integrator)[*_integrand];
}
}

// integrand
namespace Herwig {
using namespace Genfun;
FUNCTION_OBJECT_IMP(TwoOffShellIntegrand)
    
TwoOffShellIntegrand::TwoOffShellIntegrand(TwoOffShellCalculatorPtr in, Energy2 m2,
					   Energy2 mw)
 {_integrand=in;_mass2=m2;_mwidth=mw;}

double TwoOffShellIntegrand::operator() (double x) const 
{return _integrand->dGamma(sqrt(_mass2+_mwidth*tan(x)));}
}
