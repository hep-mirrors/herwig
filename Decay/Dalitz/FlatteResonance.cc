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
  os << ounit(g_,GeV);
}

void FlatteResonance::persistentInput(PersistentIStream & is, int) {
  is >> iunit(g_,GeV);
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
  for(const Energy & g : g_)
    output << " " << g/GeV; 
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
    assert(g_.size()==2);
    complex<Energy2> MGamma = sqr(g_[0])*sqrt(max(0.,rho2(q2,mpi,mpi)));
    double arg = rho2(q2,mK,mK);
    MGamma += mAB>2.*mK ? sqr(g_[1])*sqrt(arg) : sqr(g_[1])*ii*sqrt(abs(arg));
    return GeV2/(sqr(mass)-sqr(mAB)-ii*MGamma);
  }
  else if(type==ResonanceType::Flattea0) {
    assert(g_.size()==2 || g_.size()==3);
    Energy meta = CurrentGenerator::current().getParticleData(221)->mass();
    complex<Energy2> MGamma = sqr(g_[0])*sqrt(max(0.,rho2(q2,meta,mpi)));
    double arg = rho2(q2,mK,mK);
    MGamma += mAB>2.*mK ? sqr(g_[1])*sqrt(arg) : sqr(g_[1])*ii*sqrt(abs(arg));
    if(g_.size()==3 ) {
      Energy metap = CurrentGenerator::current().getParticleData(331)->mass();
      arg = rho2(q2,metap,mpi);
      MGamma += mAB>mpi+metap ? sqr(g_[2])*sqrt(arg) : sqr(g_[2])*ii*sqrt(abs(arg));
    }
    return GeV2/(sqr(mass)-sqr(mAB)-ii*MGamma);
  }
  else if(type==ResonanceType::FlatteKstar0) {
    assert(g_.size()==2);
    Energy metaP = CurrentGenerator::current().getParticleData(331)->mass();
    complex<Energy2> MGamma = sqr(g_[0])*sqrt(max(0.,rho2(q2,mK,mpi)));
    double arg = rho2(q2,mK,metaP);
    MGamma += mAB>mK+metaP ? sqr(g_[1])*sqrt(arg) : sqr(g_[1])*ii*sqrt(abs(arg));
    return GeV2/(sqr(mass)-sqr(mAB)-ii*MGamma);
  }
  else {
    assert(false);
  }
}
