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
describeHerwigMEPP21S0("Herwig::MEPP21S0", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP21S0::Init() {

  static ClassDocumentation<MEPP21S0> documentation
    ("The MEPP21S0 class implements the matrix elements for the production of 1S0 quarkonium states.");

}

Energy2 MEPP21S0::leadingOrderME2() const {
  return 2.*sqr(Constants::pi)*O1_/9./sqrt(sHat());
}

double MEPP21S0::ggME(Energy2 s, Energy2 t, Energy2 u) const {
  Energy6 Q(s*t*u);
  Energy4 P(s*t+t*u+u*s);
  Energy2 M2 = sqr(mass_);
  return 32./3.*sqr(Constants::pi)*sqr(s)*Constants::pi*O1_/mass_/sqr(s)/Q/sqr(Q-M2*P)*sqr(P)*(pow<4,1>(M2)-2*sqr(M2)*P+sqr(P)+2.*M2*Q);
}

double MEPP21S0::qgME(Energy2 s, Energy2 t, Energy2 u) const {
  Energy2 M2 = sqr(mass_);
  return -64./3.*sqr(Constants::pi)*Constants::pi*O1_*(sqr(s)+sqr(t))/(9.*mass_*u*sqr(u-M2));
}

double MEPP21S0::qbargME(Energy2 s, Energy2 t, Energy2 u) const {
  Energy2 M2 = sqr(mass_);
  return -64./3.*sqr(Constants::pi)*Constants::pi*O1_*(sqr(s)+sqr(t))/(9.*mass_*u*sqr(u-M2));
}

void MEPP21S0::doinit() {
  // // get the non-perturbative ME
  O1_ = parameters()->singletMEProduction<0>(state(),principalQuantumNumber(),0,0);
  setParticleData(1);
  MEPP2OniumPowheg::doinit();
}
