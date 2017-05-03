// -*- C++ -*-
//
// OneOffShellCalculator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneOffShellCalculator class.
//

#include "OneOffShellCalculator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

// calculate the width for a given mass
Energy OneOffShellCalculator::partialWidth(Energy2 q2) const {
  OneOffShellIntegrand integrand(this,sqr(_massptr->nominalMass()),
				 _massptr->nominalWidth()*_massptr->nominalMass());
  _scale=q2;
  // the limits
  Energy upp=min(sqrt(q2)-otherMass(_themass),_massptr->upperLimit());
  Energy low=max(_minmass,_massptr->lowerLimit());
  if(low>upp){return ZERO;}
  // transform the limits of BW smoothing
  Energy2 mass2  =_massptr->nominalMass()*_massptr->nominalMass();
  Energy2 mwidth =_massptr->nominalMass()*_massptr->nominalWidth();
  double  rhomin=atan2((low*low-mass2),mwidth);
  double  rhomax=atan2((upp*upp-mass2),mwidth);
  return _integrator.value(integrand,rhomin,rhomax);
}
