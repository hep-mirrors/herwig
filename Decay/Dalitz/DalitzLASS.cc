// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzLASS class.
//

#include "DalitzLASS.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void DalitzLASS::persistentOutput(PersistentOStream & os) const {
  os << opt_ << FNR_ << phiNR_ << FRes_ << phiRes_
     << ounit(aScat_,1./GeV) << ounit(rEff_,1./GeV);
}
void DalitzLASS::persistentInput(PersistentIStream & is, int) {
  is >> opt_ >> FNR_ >> phiNR_ >> FRes_ >> phiRes_
     >> iunit(aScat_,1./GeV) >> iunit(rEff_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzLASS,DalitzResonance>
describeHerwigDalitzLASS("Herwig::DalitzLASS", "HwDalitzDecay.so");

void DalitzLASS::Init() {

  static ClassDocumentation<DalitzLASS> documentation
    ("The DalitzLASS class implements the LASS parameterization of the Kpi s-wave.");

}

Complex DalitzLASS::BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const {
  static const Complex ii = Complex(0.,1.);
  // momenta for the resonance decay
  // off-shell
  Energy pAB=sqrt(0.25*sqr(sqr(mAB) -sqr(mA)-sqr(mB)) - sqr(mA*mB))/mAB;
  // on-shell
  Energy  pR=sqrt(0.25*sqr( mass*mass - sqr(mA) - sqr(mB)) - sqr(mA*mB))/mass;
  // non-resonant phase
  double NRphase = phiNR_+atan(1./(1./(aScat_*pAB)+0.5*rEff_*pAB));
  // resonant phase
  Energy Gamma  = width*(pAB/pR)*mass/mAB;
  double Rphase = atan(mass*Gamma/(sqr(mass)-sqr(mAB)));
  // return the result
  // BABar/BES and hopefully right form
  if (opt_==0) {
    return double(mAB/pAB)*(FNR_*sin(NRphase)*exp(ii*NRphase) +FRes_*sin(Rphase)*exp(ii*(Rphase+phiRes_+2.*NRphase)));
  }
  // BELLE form
  else if(opt_==1) {
    return (FNR_*sin(NRphase)*exp(ii*NRphase) +FRes_*sin(Rphase)*exp(ii*(Rphase+phiRes_+2.*NRphase)));
  }
  else
    assert(false);
}

void DalitzLASS::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << " " << opt_ << " " << FNR_ << " " << phiNR_ << " " << FRes_ << " " << phiRes_
	 << " " << aScat_*GeV << " " << rEff_*GeV;
}
