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

AbstractNoPIOClassDescription<SplittingFunction> SplittingFunction::initSplittingFunction;
// Definition of the static class description member.

void SplittingFunction::Init() {

  static ClassDocumentation<SplittingFunction> documentation
    ("The SplittingFunction class is the based class for 1->2 splitting functions"
     " in Herwig++");

}

double SplittingFunction::integOverPPDFFactor(const double) const {
  throw Exception() << "Base class SplittingFunction::integOverPPDFFactor() "
		    << "used you should not be using this method for this"
		    << "splitting function\n"
		    << Exception::runerror;

}

double SplittingFunction::invIntegOverPPDFFactor(const double) const {
  throw Exception() << "Base class SplittingFunction::invIntegOverPPDFFactor()"
		    << "used you should not be using this method for this"
		    << "splitting function\n"
		    << Exception::runerror;
}
