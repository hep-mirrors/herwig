// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarAmplitude class.
//

#include "ScalarAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"



using namespace Herwig;

ScalarAmplitude::ScalarAmplitude() {}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<ScalarAmplitude,Interfaced>
describeHerwigScalarAmplitude("Herwig::ScalarAmplitude", "Herwig.so");

void ScalarAmplitude::Init() {

  static ClassDocumentation<ScalarAmplitude> documentation
    ("Base class for the implementation of scalar amplitudes");

}

