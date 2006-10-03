// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoOffShellCalculator class.
//

#include "TwoOffShellCalculator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

namespace Herwig {
using namespace ThePEG;

// calculate the width for a given mass
Energy TwoOffShellCalculator::partialWidth(Energy2 q2) const {
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
  return _integrator.value(_integrand,rhomin,rhomax);
}
}

