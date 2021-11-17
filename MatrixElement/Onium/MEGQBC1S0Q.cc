// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGQBC1S0Q class.
//

#include "MEGQBC1S0Q.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

using namespace Herwig;

double MEGQBC1S0Q::me2() const {
  assert(false);
  return 0.0;
}

IBPtr MEGQBC1S0Q::clone() const {
  return new_ptr(*this);
}

IBPtr MEGQBC1S0Q::fullclone() const {
  return new_ptr(*this);
}

void MEGQBC1S0Q::doinit() {
  MEGQBCQBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<0>(bcbar,principleQuantumNumber(),0,0);
}

void MEGQBC1S0Q::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2);
}

void MEGQBC1S0Q::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGQBC1S0Q,MEGQBCQBase>
describeHerwigMEGQBC1S0Q("Herwig::MEGQBC1S0Q", "HwMEHadronOnium.so");

void MEGQBC1S0Q::Init() {

  static ClassDocumentation<MEGQBC1S0Q> documentation
    ("The MEGQBC1S0Q class implements the matrix element for g c -> B_c(1S0) b ");

}

