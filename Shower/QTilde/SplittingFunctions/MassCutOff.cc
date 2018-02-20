// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MassCutOff class.
//

#include "MassCutOff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr MassCutOff::clone() const {
  return new_ptr(*this);
}

IBPtr MassCutOff::fullclone() const {
  return new_ptr(*this);
}

void MassCutOff::persistentOutput(PersistentOStream & os) const {
}

void MassCutOff::persistentInput(PersistentIStream & is, int) {
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MassCutOff,SudakovCutOff>
describeHerwigMassCutOff("Herwig::MassCutOff", "HwShower.so");

void MassCutOff::Init() {

  static ClassDocumentation<MassCutOff> documentation
    ("There is no documentation for the MassCutOff class");

}

