// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PerturbativeDecayer class.
//

#include "PerturbativeDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void PerturbativeDecayer::persistentOutput(PersistentOStream & os) const {
}

void PerturbativeDecayer::persistentInput(PersistentIStream & is, int) {
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<PerturbativeDecayer,DecayIntegrator>
describeHerwigPerturbativeDecayer("Herwig::PerturbativeDecayer", "Herwig.so HwPerturbativeDecay.so");

void PerturbativeDecayer::Init() {

  static ClassDocumentation<PerturbativeDecayer> documentation
    ("The PerturbativeDecayer class is the mase class for perturbative decays in Herwig");

}
