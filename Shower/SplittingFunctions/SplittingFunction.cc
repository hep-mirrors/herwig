// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplittingFunction class.
//

#include "SplittingFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SplittingFunction.tcc"
#endif


using namespace Herwig;

SplittingFunction::~SplittingFunction() {}

AbstractNoPIOClassDescription<SplittingFunction> SplittingFunction::initSplittingFunction;
// Definition of the static class description member.

void SplittingFunction::Init() {

  static ClassDocumentation<SplittingFunction> documentation
    ("There is no documentation for the SplittingFunction class");

}

