// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoBodyAllOnCalculator class.
//

#include "TwoBodyAllOnCalculator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "GenericWidthGenerator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TwoBodyAllOnCalculator.tcc"
#endif

namespace Herwig {
using namespace ThePEG;


TwoBodyAllOnCalculator::~TwoBodyAllOnCalculator() {}

NoPIOClassDescription<TwoBodyAllOnCalculator> TwoBodyAllOnCalculator::initTwoBodyAllOnCalculator;
// Definition of the static class description member.

void TwoBodyAllOnCalculator::Init() {

  static ClassDocumentation<TwoBodyAllOnCalculator> documentation
    ("The TwoBodyAllOnCalculator class calculates the partial width for a"
     " two body decay where both decay products are on-shell");

}

Energy TwoBodyAllOnCalculator::partialWidth(Energy2 scale) const
 {return _widthgen->partial2BodyWidth(_mode,sqrt(scale),_mass1,_mass2);}

}
