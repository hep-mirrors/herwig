// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzSigma class.
//

#include "DalitzSigma.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void DalitzSigma::persistentOutput(PersistentOStream & os) const {
  os << ounit(a_,GeV2) << ounit(b1_,GeV) << ounit(b2_,1./GeV) << ounit(g4Pi_,GeV);
}

void DalitzSigma::persistentInput(PersistentIStream & is, int) {
  is >> iunit(a_,GeV2) >> iunit(b1_,GeV) >> iunit(b2_,1./GeV) >> iunit(g4Pi_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzSigma,DalitzResonance>
  describeHerwigDalitzSigma("Herwig::DalitzSigma", "HwDalitzDecay.so");

void DalitzSigma::Init() {

  static ClassDocumentation<DalitzSigma> documentation
    ("The DalitzSigma class implements the model of Bou and Zou for the sigma propagator");

}

void DalitzSigma::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << a_/GeV2 << " " << b1_/GeV << " " << b2_*GeV << " " << g4Pi_/GeV;
}

Complex DalitzSigma::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  static const Complex II(0.,1.);
  Energy2 s(sqr(mAB));
  Energy mpi = CurrentGenerator::current().getParticleData(111)->mass();
  Energy2 sA=0.5*sqr(mpi);
  // two pion width
  Energy gamma = (b1_+b2_*s)*exp(-(s-sqr(mass))/a_)*(s-sA)/(sqr(mass)-sA)*
    sqrt(1.-4.*sqr(mpi)/s)/sqrt(1.-4.*sqr(mpi/mass));
  if(mAB>4.*mpi) {
    gamma+= g4Pi_*rho4pi(s,mpi)/rho4pi(sqr(mass),mpi);
  }
  // return the propagator
  return GeV2/(sqr(mass)-s-II*mass*gamma);
}
