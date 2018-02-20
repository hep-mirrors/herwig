// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VariableMassCutOff class.
//

#include "VariableMassCutOff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr VariableMassCutOff::clone() const {
  return new_ptr(*this);
}

IBPtr VariableMassCutOff::fullclone() const {
  return new_ptr(*this);
}

void VariableMassCutOff::persistentOutput(PersistentOStream & os) const {
}

void VariableMassCutOff::persistentInput(PersistentIStream & is, int) {
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VariableMassCutOff,SudakovCutOff>
describeHerwigVariableMassCutOff("Herwig::VariableMassCutOff", "HwShower.so");

void VariableMassCutOff::Init() {

  static ClassDocumentation<VariableMassCutOff> documentation
    ("There is no documentation for the VariableMassCutOff class");

}

