// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronSpectrum class.
//

#include "HadronSpectrum.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HadronSpectrum::HadronSpectrum() {}

HadronSpectrum::~HadronSpectrum() {}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HadronSpectrum::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void HadronSpectrum::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<HadronSpectrum,Interfaced>
  describeHerwigHadronSpectrum("Herwig::HadronSpectrum", "HwHadronization.so");

void HadronSpectrum::Init() {

  static ClassDocumentation<HadronSpectrum> documentation
    ("There is no documentation for the HadronSpectrum class");

}

