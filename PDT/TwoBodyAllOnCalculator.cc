// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoBodyAllOnCalculator class.
//

#include "TwoBodyAllOnCalculator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "GenericWidthGenerator.h"

namespace Herwig {
using namespace ThePEG;

Energy TwoBodyAllOnCalculator::partialWidth(Energy2 scale) const {
  return _widthgen->partial2BodyWidth(_mode,sqrt(scale),_mass1,_mass2);
}

}
