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
  os << g1_ << g2_;
}

void FlatteResonance::persistentInput(PersistentIStream & is, int) {
  is >> g1_ >> g2_;
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
  output << " " << g1_ << " " << g2_; 
}

Complex FlatteResonance::BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const {
  static const Complex ii = Complex(0.,1.);
  Energy mpi = CurrentGenerator::current().getParticleData(111)->mass();
  Energy mK  = CurrentGenerator::current().getParticleData(321)->mass();
  if(type==ResonanceType::Flattef0) {
    Energy Gamma_pi = g1_*sqrt(0.25*sqr(mAB)-sqr(mpi));
    Energy2 arg = 0.25*sqr(mAB)-sqr(mK);
    complex<Energy> Gamma_K  = arg>=ZERO ? g2_*sqrt(arg) : g2_*ii*sqrt(-arg);
    return GeV2/(sqr(mass)-sqr(mAB)-ii*mass*(Gamma_pi+Gamma_K));
  }
  else if(type==ResonanceType::Flattea0) {
    Energy meta = CurrentGenerator::current().getParticleData(221)->mass();
    Energy2 q2=sqr(mAB);
    Energy Gamma_pi = mAB>meta+mpi ? g1_*0.5/mAB*sqrt((q2-sqr(mpi+meta))*(q2-sqr(mpi-meta))) : ZERO;
    Energy2 arg = 0.25*sqr(mAB)-sqr(mK);
    complex<Energy> Gamma_K  = arg>=ZERO ? g2_*sqrt(arg) : g2_*ii*sqrt(-arg);
    return GeV2/(sqr(mass)-sqr(mAB)-ii*mass*(Gamma_pi+Gamma_K));
  }
  else if(type==ResonanceType::FlatteKstar0) {
    Energy metaP = CurrentGenerator::current().getParticleData(331)->mass();
    Energy2 q2=sqr(mAB);
    Energy Gamma_pi = g1_*0.5/mAB*sqrt((q2-(sqr(mpi+mK)))*(q2-(sqr(mpi-mK))));
    Energy2 arg = 0.5/q2*(q2-(sqr(mK+metaP)))*(q2-(sqr(mK-metaP)));
    complex<Energy> Gamma_eta  = mAB>mK+metaP ? g2_*sqrt(arg) : g2_*ii*sqrt(abs(arg));
    return GeV2/(sqr(mass)-sqr(mAB)-ii*mass*(Gamma_pi+Gamma_eta));
  }
  else {
    assert(false);
  }
}
