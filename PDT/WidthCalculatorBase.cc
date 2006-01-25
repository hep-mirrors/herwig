// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WidthCalculatorBase class.
//

#include "WidthCalculatorBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "WidthCalculatorBase.tcc"
#endif

namespace Herwig {
using namespace ThePEG;

WidthCalculatorBase::~WidthCalculatorBase() {}

AbstractNoPIOClassDescription<WidthCalculatorBase> WidthCalculatorBase::initWidthCalculatorBase;
// Definition of the static class description member.

void WidthCalculatorBase::Init() {

  static ClassDocumentation<WidthCalculatorBase> documentation
    ("The WidthCalculatorBase class is the base class for the calculation of"
     " running widths in Herwig++");

}

}
