// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MIPWA class.
//

#include "MIPWA.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void MIPWA::persistentOutput(PersistentOStream & os) const {
  os << ounit(energy_,GeV) << mag_ << phase_;
}

void MIPWA::persistentInput(PersistentIStream & is, int) {
  is >> iunit(energy_,GeV) >> mag_ >> phase_;
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<MIPWA,DalitzResonance>
  describeHerwigMIPWA("Herwig::MIPWA", "HwDalitzDecay.so");

void MIPWA::Init() {

  static ClassDocumentation<MIPWA> documentation
    ("The MIPWA class allows the use on experimental extractions from"
     " Model Independent Partial Wave Analyses. ");

}

Complex MIPWA::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  static Complex ii(0.,1.);
  if(!iMag_) {
    iMag_   = make_InterpolatorPtr(mag_  ,energy_,3);
    iPhase_ = make_InterpolatorPtr(phase_,energy_,3);
  }
  return (*iMag_)(mAB)*exp(ii*(*iPhase_)(mAB));
}

void MIPWA::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << energy_.size();
  for(unsigned int ix=0;ix<energy_.size();++ix)
    output << " " << energy_[ix]/GeV << " " << mag_[ix] << " " << phase_[ix];
}
