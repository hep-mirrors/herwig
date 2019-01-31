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
#include "ThePEG/Interface/Parameter.h"

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
  os << ounit(pTmin_,GeV) << ounit(pT2min_,GeV2);
}

void PTCutOff::persistentInput(PersistentIStream & is, int) {
  is >> iunit(pTmin_,GeV) >> iunit(pT2min_,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PTCutOff,SudakovCutOff>
describeHerwigPTCutOff("Herwig::PTCutOff", "HwShower.so");

void PTCutOff::Init() {

  static ClassDocumentation<PTCutOff> documentation
    ("There is no documentation for the PTCutOff class");


  static Parameter<PTCutOff,Energy> interfacepTmin
    ("pTmin",
     "The minimum pT if using a cut-off on the pT",
     &PTCutOff::pTmin_, GeV, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

}

void PTCutOff::doinit() {
  pT2min_ = sqr(pTmin_);
  SudakovCutOff::doinit();
}

const vector<Energy> & PTCutOff::virtualMasses(const IdList & ids) {
  static vector<Energy> output;
  output.clear();
  for(auto id : ids) 
      output.push_back(id->mass());
  return output;
}
