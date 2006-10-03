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


// calculate the width for a given mass
Energy OneOffShellCalculator::partialWidth(Energy2 q2) const {
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
  return _integrator.value(_integrand,rhomin,rhomax);
}

}


 

