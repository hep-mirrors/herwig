// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzResonance class.
//

#include "DalitzResonance.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"

using namespace Herwig;

void DalitzResonance::persistentOutput(PersistentOStream & os) const {
  os << id << oenum(type) << ounit(mass,GeV) << ounit(width,GeV)
     << daughter1 << daughter2 << spectator
     << amp << ounit(R,1./GeV);
}

void DalitzResonance::persistentInput(PersistentIStream & is, int) {
  is >> id >> ienum(type) >> iunit(mass,GeV) >> iunit(width,GeV)
     >> daughter1 >> daughter2 >> spectator
     >> amp >> iunit(R,1./GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzResonance,Base>
  describeHerwigDalitzResonance("Herwig::DalitzResonance", "HwDalitzDecay.so");

void DalitzResonance::Init() {

  static ClassDocumentation<DalitzResonance> documentation
    ("The DalitzResonance class provides a container class for"
     " information on resonances in multi-body dalitz decays.");

}

Complex DalitzResonance::BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const {
  static const Complex ii = Complex(0.,1.);
  // non-resonant pieces
  if(abs(type)/10==10) return 1.;
  // momenta for the resonance decay
  // off-shell
  Energy pAB=sqrt(0.25*sqr(sqr(mAB) -sqr(mA)-sqr(mB)) - sqr(mA*mB))/mAB;
  if(type==ResonanceType::BABARf0) {
    double rho = 2.*pAB/mAB;
    return GeV2/(sqr(mass)-sqr(mAB)-ii*mass*width*rho);
  }
  else if (type==ResonanceType::Flattef0) {
    assert(false);
    // Energy mpi = getParticleData(111)->mass();
    // Energy mK  = getParticleData(321)->mass();
    // Energy Gamma_pi = f0gpi_*sqrt(0.25*sqr(mAB)-sqr(mpi));
    // Energy2 arg = 0.25*sqr(mAB)-sqr(mK);
    // complex<Energy> Gamma_K  = arg>=ZERO ? f0gK_*sqrt(arg) : f0gK_*ii*sqrt(-arg);
    // output *= GeV2/(sqr(mass)-sqr(mAB)-ii*mass*(Gamma_pi+Gamma_K));
    // return output;
  }
  else if (type==ResonanceType::Spin0Complex) {
    complex<Energy> sR(mass,width);
    return GeV2/(sqr(sR)-sqr(mAB));
  }
  //  on-shell
  Energy  pR=sqrt(0.25*sqr( mass*mass - sqr(mA) - sqr(mB)) - sqr(mA*mB))/mass;
  // Blatt-Weisskopf factors
  double fR=1;
  unsigned int power(1);
  if(type!=ResonanceType::Spin0 &&
     type!=ResonanceType::Spin0E691) {
    double r1A(R*pR),r1B(R*pAB);
    // Blatt-Weisskopf factors and spin piece
    switch (type) {
    case ResonanceType::Spin0Gauss:
      fR = exp(-(r1B-r1A)/12.);
      // cerr << "testing scalar B " <<   exp(+r1A/12.) << "\n";
      break;
    case ResonanceType::Spin1: case ResonanceType::Spin1E691 :
      fR=sqrt( (1. + sqr(r1A)) / (1. + sqr(r1B)) );
      // cerr << "testing vector B " << sqrt(1. + sqr(r1A))   << "\n";
      power=3;
      break;
    case ResonanceType::Spin2: case ResonanceType::Spin2E691:
      fR = sqrt( (9. + sqr(r1A)*(3.+sqr(r1A))) / (9. + sqr(r1B)*(3.+sqr(r1B))));
      // cerr << "testing tensor B " <<  sqrt( (9. + sqr(r1A)*(3.+sqr(r1A))))  << "\n";
      power=5;
      break;
    default :
      assert(false);
    }
  }
  // multiply by Breit-Wigner piece and return
  if (type/10 == 1 ) {
    return fR*sqrt(0.5*width/GeV/Constants::pi)*GeV/(mAB-mass-complex<Energy>(ZERO,0.5*width));
  }
  else {
    Energy gam = width*pow(pAB/pR,power)*(mass/mAB)*fR*fR;
    return fR*GeV2/(sqr(mass)-sqr(mAB)-mass*gam*ii);
  }
}

void DalitzResonance::dataBaseOutput(ofstream & output) {
  output << id << " " << oenum(type) << " "
	 << mass/GeV << " " << width/GeV << " "
	 << daughter1 << " " << daughter2 << " "
	 << spectator << " " 
	 << abs(amp) << " " << arg(amp) << " "
	 << R*GeV; 
}
