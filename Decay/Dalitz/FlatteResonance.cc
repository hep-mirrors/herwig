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

Complex FlatteResonance::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  assert(false);
}
