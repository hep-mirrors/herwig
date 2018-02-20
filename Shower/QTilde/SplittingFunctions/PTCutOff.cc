// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PTCutOff class.
//

#include "PTCutOff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr PTCutOff::clone() const {
  return new_ptr(*this);
}

IBPtr PTCutOff::fullclone() const {
  return new_ptr(*this);
}

void PTCutOff::persistentOutput(PersistentOStream & os) const {
}

void PTCutOff::persistentInput(PersistentIStream & is, int) {
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PTCutOff,SudakovCutOff>
describeHerwigPTCutOff("Herwig::PTCutOff", "HwShower.so");

void PTCutOff::Init() {

  static ClassDocumentation<PTCutOff> documentation
    ("There is no documentation for the PTCutOff class");

}

void PTCutOff::doinit() {
  SudakovCutOff::doinit();
}
