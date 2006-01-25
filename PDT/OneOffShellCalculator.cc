// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneOffShellCalculator class.
//

#include "OneOffShellCalculator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "OneOffShellCalculator.tcc"
#endif

namespace Herwig {
using namespace ThePEG;

OneOffShellCalculator::~OneOffShellCalculator() 
{delete _Integrator;delete _integrand;}

NoPIOClassDescription<OneOffShellCalculator> OneOffShellCalculator::initOneOffShellCalculator;
// Definition of the static class description member.

void OneOffShellCalculator::Init() {

  static ClassDocumentation<OneOffShellCalculator> documentation
    ("The OneOffShellCalculator class performs the integration of the partial"
     " width when one particle outgoing particle is off-shell");

}

// calculate the width for a given mass
  Energy OneOffShellCalculator::partialWidth(Energy2 q2) const
{
  _scale=q2;
  // the limits
  Energy upp=min(sqrt(q2)-otherMass(_themass),_massptr->upperLimit());
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
FUNCTION_OBJECT_IMP(OneOffShellIntegrand)
    
  OneOffShellIntegrand::OneOffShellIntegrand(tOneOffShellCalculatorPtr in, Energy2 m2,
					     Energy2 mw)
 {_integrand=in;_mass2=m2;_mwidth=mw;}

OneOffShellIntegrand::~OneOffShellIntegrand() {}
  
OneOffShellIntegrand::OneOffShellIntegrand(const OneOffShellIntegrand & right) 
{  }
  
double OneOffShellIntegrand::operator() (double x) const 
{return _integrand->dGamma(sqrt(_mass2+_mwidth*tan(x)));}

}
