// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGammaAmplitude class.
//

#include "GammaGammaAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<GammaGammaAmplitude,Interfaced>
describeThePEGGammaGammaAmplitude("ThePEG::GammaGammaAmplitude", "HwMEGammaGamma.so");

void GammaGammaAmplitude::Init() {

  static ClassDocumentation<GammaGammaAmplitude> documentation
    ("The GammaGammaAmplitude class provides a base class for"
     " the implementation of gamma gamma -> X processes");

}

Selector<const ColourLines *>
GammaGammaAmplitude::colourGeometries(unsigned int, tcDiagPtr ) const {
  static ColourLines c("");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c);
  return sel;
}
