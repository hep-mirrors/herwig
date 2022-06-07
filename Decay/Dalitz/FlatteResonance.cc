// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FlatteResonance class.
//

#include "FlatteResonance.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;
void FlatteResonance::persistentOutput(PersistentOStream & os) const {
  os << gpi_ << gK_;
}

void FlatteResonance::persistentInput(PersistentIStream & is, int) {
  is >> gpi_ >> gK_;
}

//The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<FlatteResonance,DalitzResonance>
describeHerwigFlatteResonance("Herwig::FlatteResonance", "HwDalitzDecay.so");

void FlatteResonance::Init() {

  static ClassDocumentation<FlatteResonance> documentation
    ("The FlatteResonance class implements the Flatte lineshape for Dalitz decays.");

}

void FlatteResonance::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << gpi_ << " " << gK_; 
}

Complex FlatteResonance::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  assert(type==ResonanceType::Flattef0);
  static const Complex ii = Complex(0.,1.);
  Energy mpi = CurrentGenerator::current().getParticleData(111)->mass();
  Energy mK  = CurrentGenerator::current().getParticleData(321)->mass();
  Energy Gamma_pi = gpi_*sqrt(0.25*sqr(mAB)-sqr(mpi));
  Energy2 arg = 0.25*sqr(mAB)-sqr(mK);
  complex<Energy> Gamma_K  = arg>=ZERO ? gK_*sqrt(arg) : gK_*ii*sqrt(-arg);
  return GeV2/(sqr(mass)-sqr(mAB)-ii*mass*(Gamma_pi+Gamma_K));
}
