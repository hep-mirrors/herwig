// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PiPiI2 class.
//

#include "PiPiI2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void PiPiI2::persistentOutput(PersistentOStream & os) const {
  os << ounit(a_,1./GeV) << ounit(b_,1./GeV2) << ounit(c_,1./GeV2/GeV2) << ounit(d_,1./GeV2/GeV2/GeV2)
     << ounit(mmin_,GeV) << ounit(mmax_,GeV) << deltaEta_;
}

void PiPiI2::persistentInput(PersistentIStream & is, int) {
  is >> iunit(a_,1./GeV) >> iunit(b_,1./GeV2) >> iunit(c_,1./GeV2/GeV2) >> iunit(d_,1./GeV2/GeV2/GeV2)
     >> iunit(mmin_,GeV) >> iunit(mmax_,GeV) >> deltaEta_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PiPiI2,DalitzResonance>
describeHerwigPiPiI2("Herwig::PiPiI2", "HwDalitzDecay.so");

void PiPiI2::Init() {

  static ClassDocumentation<PiPiI2> documentation
    ("The PiPiI2 class provides an implementation of the I=2 s-eave for pipi.");

}

Complex PiPiI2::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  static Complex ii(0.,1.);
  Energy2 m2 = sqr(mAB);
  Energy mpi = CurrentGenerator::current().getParticleData(111)->mass();
  double delta = -a_*sqrt(0.25*m2-sqr(mpi))/(1.+m2*(b_+m2*(c_+d_*m2)));
  double eta = 1.;
  if(mAB>mmax_) {
    eta = 1. - deltaEta_;
  }
  else if(mAB>mmin_) {
    eta = 1. - 0.5*deltaEta_*(1.-cos(Constants::pi*(mAB-mmin_)/(mmax_-mmin_)));
  }
  return -0.5*ii*(eta*exp(2.*ii*delta)-1.);
}

void PiPiI2::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << a_*GeV << " " << b_*GeV2 << " " << c_*GeV2*GeV2 << " " << d_*GeV2*GeV2*GeV2
	 << " " << mmin_/GeV << " " << mmax_/GeV << " " << deltaEta_;
}
