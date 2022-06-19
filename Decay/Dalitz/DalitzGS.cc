// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzGS class.
//

#include "DalitzGS.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Decay/ResonanceHelpers.h"

using namespace Herwig;

DalitzGS::DalitzGS(long pid, ResonanceType::Type rtype, Energy m, Energy w,
		   unsigned int d1, unsigned int d2, unsigned int s,
		   double mag, double phi, InvEnergy rr)
    : DalitzResonance(pid,rtype,m,w,d1,d2,s,mag,phi,rr) {
  mpi_ = CurrentGenerator::current().getParticleData(211)->mass();
  hres_ = Resonance::Hhat(sqr(mass),mass,width,mpi_,mpi_);
  dh_   = Resonance::dHhatds(mass,width,mpi_,mpi_);
  h0_   = Resonance::H(ZERO,mass,width,mpi_,mpi_,dh_,hres_);
}

void DalitzGS::persistentOutput(PersistentOStream & os) const {
  os << ounit(mpi_,GeV) << dh_ << ounit(hres_,GeV2) << ounit(h0_,GeV2);
}

void DalitzGS::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mpi_,GeV) >> dh_ >> iunit(hres_,GeV2) >> iunit(h0_,GeV2);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzGS,DalitzResonance>
describeHerwigDalitzGS("Herwig::DalitzGS", "HwDalitzDalitz.so");

void DalitzGS::Init() {

  static ClassDocumentation<DalitzGS> documentation
    ("The DalitzGS class implements the Gounaris and Sakurai Phys. Rev. "
     "Lett. 21, 244 (1968) form for the propagator.");

}

Complex DalitzGS::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  return GeV2/sqr(mass)*Resonance::BreitWignerGS(sqr(mAB),mass,width,mpi_,mpi_,h0_,dh_,hres_);
}
