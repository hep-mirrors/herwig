// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModelHadronSpectrum class.
//

#include "StandardModelHadronSpectrum.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

StandardModelHadronSpectrum::StandardModelHadronSpectrum() {}

StandardModelHadronSpectrum::~StandardModelHadronSpectrum() {}

IBPtr StandardModelHadronSpectrum::clone() const {
  return new_ptr(*this);
}

IBPtr StandardModelHadronSpectrum::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void StandardModelHadronSpectrum::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void StandardModelHadronSpectrum::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<StandardModelHadronSpectrum,HadronSpectrum>
  describeHerwigStandardModelHadronSpectrum("Herwig::StandardModelHadronSpectrum", "HwHadronization.so");

void StandardModelHadronSpectrum::Init() {

  static ClassDocumentation<StandardModelHadronSpectrum> documentation
    ("There is no documentation for the StandardModelHadronSpectrum class");

}

