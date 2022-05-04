// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmnesFunction class.
//

#include "OmnesFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<OmnesFunction,Interfaced>
describeHerwigOmnesFunction("Herwig::OmnesFunction", "Herwig.so");

void OmnesFunction::Init() {

  static ClassDocumentation<OmnesFunction> documentation
    ("The OmnesFunction class provides a basis class for the implementation of the Omnes function.");

}

