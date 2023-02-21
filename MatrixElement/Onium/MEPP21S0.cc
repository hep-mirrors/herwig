// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP21S0 class.
//

#include "MEPP21S0.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MEPP21S0::MEPP21S0() {}

IBPtr MEPP21S0::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP21S0::fullclone() const {
  return new_ptr(*this);
}

void MEPP21S0::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2);
}

void MEPP21S0::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP21S0,MEPP2OniumPowheg>
describeHerwigMEPP21S0("Herwig::MEPP21S0", "HwMEHadronOnium.so");

void MEPP21S0::Init() {

  static ClassDocumentation<MEPP21S0> documentation
    ("The MEPP21S0 class implements the matrix elements for the production of 1S0 quarkonium states.");

}

Energy2 MEPP21S0::leadingOrderME2() const {
  return 2.*sqr(Constants::pi)*O1_/9./sqrt(sHat());
}
 
void MEPP21S0::doinit() {
  // // get the non-perturbative ME
  O1_ = parameters()->singletMEProduction<0>(state(),principalQuantumNumber(),0,0);
  setParticleData(1);
  MEPP2OniumPowheg::doinit();
}
