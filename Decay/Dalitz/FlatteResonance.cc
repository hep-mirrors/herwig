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
  os << ounit(g1_,GeV) << ounit(g2_,GeV);
}

void FlatteResonance::persistentInput(PersistentIStream & is, int) {
  is >> iunit(g1_,GeV) >> iunit(g2_,GeV);
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
  output << " " << g1_/GeV << " " << g2_/GeV; 
}

namespace {

double rho2(const Energy2 & q2, const Energy & m1, const Energy & m2) {
  return (1.-sqr(m1+m2)/q2)*(1.-sqr(m1-m2)/q2);
}

}

Complex FlatteResonance::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  static const Complex ii = Complex(0.,1.);
  Energy mpi = CurrentGenerator::current().getParticleData(111)->mass();
  Energy mK  = CurrentGenerator::current().getParticleData(321)->mass();
  Energy2 q2=sqr(mAB);
  if(type==ResonanceType::Flattef0) {
    complex<Energy2> MGamma = sqr(g1_)*sqrt(max(0.,rho2(q2,mpi,mpi)));
    double arg = rho2(q2,mK,mK);
    MGamma += mAB>2.*mK ? sqr(g2_)*sqrt(arg) : sqr(g2_)*ii*sqrt(abs(arg));
    return GeV2/(sqr(mass)-sqr(mAB)-ii*MGamma);
  }
  else if(type==ResonanceType::Flattea0) {
    Energy meta = CurrentGenerator::current().getParticleData(221)->mass();
    complex<Energy2> MGamma = sqr(g1_)*sqrt(max(0.,rho2(q2,meta,mpi)));
    double arg = rho2(q2,mK,mK);
    MGamma += mAB>2.*mK ? sqr(g2_)*sqrt(arg) : sqr(g2_)*ii*sqrt(abs(arg));
    return GeV2/(sqr(mass)-sqr(mAB)-ii*MGamma);
  }
  else if(type==ResonanceType::FlatteKstar0) {
    Energy metaP = CurrentGenerator::current().getParticleData(331)->mass();
    complex<Energy2> MGamma = sqr(g1_)*sqrt(max(0.,rho2(q2,mK,mpi)));
    double arg = rho2(q2,mK,metaP);
    MGamma += mAB>mK+metaP ? sqr(g2_)*sqrt(arg) : sqr(g2_)*ii*sqrt(abs(arg));
    return GeV2/(sqr(mass)-sqr(mAB)-ii*MGamma);
  }
  else {
    assert(false);
  }
}
