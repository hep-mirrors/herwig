// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SudakovCutOff class.
//

#include "SudakovCutOff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<SudakovCutOff,Interfaced>
describeHerwigSudakovCutOff("Herwig::SudakovCutOff", "HwShower.so");

void SudakovCutOff::Init() {

  static ClassDocumentation<SudakovCutOff> documentation
    ("Base class for cut-offs in the form factor");

}

