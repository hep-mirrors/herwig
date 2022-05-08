// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HeavyMesonWidthGenerator class.
//

#include "HeavyMesonWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/HeavyMeson/HQETStrongDecayer.h"

using namespace Herwig;

HeavyMesonWidthGenerator::HeavyMesonWidthGenerator() 
  : fPi_(130.2*MeV), g_(0.566), h_(0.544), hp_(0.413), k_(0.407), kp_(0.242), gtilde_(0.283),
    psiL_(0.), psiS_(0.041), Lambda_(1.*GeV), couplingsSet_(false)
{}

IBPtr HeavyMesonWidthGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr HeavyMesonWidthGenerator::fullclone() const {
  return new_ptr(*this);
}

void HeavyMesonWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << ounit(fPi_,MeV) << g_ << h_ << hp_ << k_ << kp_ << gtilde_
     << psiL_ << psiS_ << ounit(Lambda_,GeV) << couplingsSet_;
}

void HeavyMesonWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> iunit(fPi_,MeV) >> g_ >> h_ >> hp_ >> k_ >> kp_ >> gtilde_
     >> psiL_ >> psiS_ >> iunit(Lambda_,GeV) >> couplingsSet_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HeavyMesonWidthGenerator,GenericWidthGenerator>
describeHerwigHeavyMesonWidthGenerator("Herwig::HeavyMesonWidthGenerator", "HwHMDecay.so");

void HeavyMesonWidthGenerator::Init() {

  static ClassDocumentation<HeavyMesonWidthGenerator> documentation
    ("The HeavyMesonWidthGenerator class calculates the width for heavy meson decays");

  static Parameter<HeavyMesonWidthGenerator,Energy> interfacefPi
    ("fPi",
     "The pion decay constant",
     &HeavyMesonWidthGenerator::fPi_, MeV, 130.2*MeV, 100.0*MeV, 200.0*MeV,
     false, false, Interface::limited);

  static Parameter<HeavyMesonWidthGenerator,double> interfaceg
    ("g",
     "The coupling for 1S (0-,1-) decays",
     &HeavyMesonWidthGenerator::g_, 0.566, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HeavyMesonWidthGenerator,double> interfaceh
    ("h",
     "The coupling for 1P (0+,1+) decays",
     &HeavyMesonWidthGenerator::h_, 0.544, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HeavyMesonWidthGenerator,double> interfacehp
    ("hp",
     "The coupling for 1P (1+,2+) decays",
     &HeavyMesonWidthGenerator::hp_, 0.413, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HeavyMesonWidthGenerator,double> interfacek
    ("k",
     "The coupling for 1D (2-,3-) decays",
     &HeavyMesonWidthGenerator::k_, 0.407, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HeavyMesonWidthGenerator,double> interfacekp
    ("kp",
     "The coupling for 1D (1-,2-) decays",
     &HeavyMesonWidthGenerator::kp_, 0.242, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HeavyMesonWidthGenerator,double> interfacegtilde
    ("gtilde",
     "The coupling for 2S (0-,1-) decays",
     &HeavyMesonWidthGenerator::gtilde_, 0.283, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<HeavyMesonWidthGenerator,Energy> interfacefLambda
    ("Lambda",
     "Strong decays momentum scale",
     &HeavyMesonWidthGenerator::Lambda_, GeV, 1.*GeV, .1*GeV, 2.*GeV,
     false, false, Interface::limited);

  static Parameter<HeavyMesonWidthGenerator,double> interfacefpsiL
    ("psiL",
     "D_1 mixing angle for up and down heavy mesons",
     &HeavyMesonWidthGenerator::psiL_, 0., -M_PI/2., M_PI/2.,
     false, false, Interface::limited);

  static Parameter<HeavyMesonWidthGenerator,double> interfacefpsiS
    ("psiS",
     "D_1 mixing angle for strange heavy mesons",
     &HeavyMesonWidthGenerator::psiS_, 0.041, -M_PI/2., M_PI/2.,
     false, false, Interface::limited);
}
 
void HeavyMesonWidthGenerator::setupMode(tcDMPtr, tDecayIntegratorPtr decayer,
					 unsigned int) {
  // cast the decayer
  Ptr<HQETStrongDecayer>::tcptr strong = dynamic_ptr_cast<Ptr<HQETStrongDecayer>::tcptr >(decayer);
  if(strong) {
    if(!couplingsSet_) {
      fPi_    = strong->fPi_   ;
      g_      = strong->g_     ;
      h_      = strong->h_     ;
      hp_     = strong->hp_    ;
      k_      = strong->k_     ;
      kp_     = strong->kp_    ;
      gtilde_ = strong->gtilde_;
      psiL_   = strong->psiL_  ;
      psiS_   = strong->psiS_  ;
      Lambda_ = strong->Lambda_;
      couplingsSet_ = true;
    }
  }
}

void HeavyMesonWidthGenerator::dataBaseOutput(ofstream & output, bool header) {
  if(header) output << "update Width_Generators set parameters=\"";
  // info from the base class
  GenericWidthGenerator::dataBaseOutput(output,false);
  // info from this class
  output << "newdef " << name() << ":fPi    " << fPi_/MeV    << "\n";
  output << "newdef " << name() << ":g      " << g_          << "\n";
  output << "newdef " << name() << ":h      " << h_          << "\n";
  output << "newdef " << name() << ":hp     " << hp_         << "\n";
  output << "newdef " << name() << ":k      " << k_          << "\n";
  output << "newdef " << name() << ":kp     " << kp_         << "\n";
  output << "newdef " << name() << ":gtilde " << gtilde_     << "\n";
  output << "newdef " << name() << ":Lambda " << Lambda_/GeV << "\n";
  output << "newdef " << name() << ":psiL   " << psiL_       << "\n";
  output << "newdef " << name() << ":psiS   " << psiS_       << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << name() << "\";" 
  		    << endl;
}

Energy HeavyMesonWidthGenerator::partial2BodyWidth(int imode, Energy q,Energy m1,
						   Energy m2) const {
  if(q<m1+m2) return ZERO;
  // mode from the base class
  int mecode(MEcode(imode));
  if(mecode<=100) { 
    return GenericWidthGenerator::partial2BodyWidth(imode,q,m1,m2);
  }
  // calcluate the decay momentum
  Energy2 q2(q*q),m12(m1*m1),m22(m2*m2),
    pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
  Energy pcm(sqrt(pcm2)),gam(ZERO),msum(q+m1);
  double test=0.;
  if(mecode==101) {
    test = 4.*sqr(g_)*q/m1*sqr(pcm)/sqr(fPi_);
  }
  else if(mecode==102) {
    test = 4.*sqr(g_)*m1*sqr(pcm)/3./sqr(fPi_)/q;
  }
  else if(mecode==103) {
    test = 8.*sqr(g_)*m1*sqr(pcm)/3./sqr(fPi_)/q;
  }
  else if(mecode==104) {
    long id = abs(particle()->id());
    double psi = (id%100)/10!=3 ? psiL_ : psiS_;
    unsigned int itemp = abs(id)-abs(id)%1000;
    double mix1(cos(psi)),mix2(sin(psi));
    if(itemp==20000) {
      swap(mix1,mix2);
      mix1 *=-1.;
    }
    InvEnergy2 A = -2.*sqrt(2.*m1/3./q)*hp_*mix1/fPi_/Lambda_;
    double     B = -h_*mix2/1/fPi_*sqrt(m1/q)/q*(sqr(q)-sqr(m1)+sqr(m2));
    B += A*sqr(pcm);
    test = (3.*sqr(A*q*sqr(pcm)) +sqr(B)/3.*(3.*sqr(m1)+sqr(pcm))
	    -A*B*sqr(pcm)*(sqr(q)+sqr(m1)-sqr(m2)))/sqr(m1);
  }
  else if(mecode==105) {
    test = 32.*sqr(hp_)*m1*sqr(sqr(pcm))/15./sqr(fPi_)/sqr(Lambda_)/q;
  }
  else if(mecode==106) {
    test = 16.*sqr(hp_)*m1*sqr(sqr(pcm))/5./sqr(fPi_)/sqr(Lambda_)/q;
  }
  else if(mecode==107) {
    test = sqr(h_)/sqr(fPi_)*m1/pow<3,1>(q)*sqr(sqr(q)-sqr(m1)+sqr(m2));
  }
  else if(mecode==108) {
     test = 8.*m1*sqr(kp_*pcm*(sqr(q)-sqr(m1)+sqr(m2)))/9./sqr(fPi_*Lambda_*q)/q;
  }
  else if(mecode==109) {
    test = 4.*m1*sqr(kp_*pcm*(sqr(q)-sqr(m1)+sqr(m2)))/9./sqr(fPi_*Lambda_*q)/q;
  }
  else if(mecode==110) {
    test = 2.*sqr(kp_*pcm*(sqr(q)-sqr(m1)+sqr(m2)))/15./sqr(fPi_*Lambda_)/m2/pow<5,1>(q)*
      (sqr(q)*(sqr(q)+8.*sqr(m1)-2.*sqr(m2))+sqr(sqr(m1)-sqr(m2)));
  }
  else if(mecode==111) {
    test = 32.*sqr(k_)*pow<6,1>(pcm)/225./sqr(fPi_*sqr(Lambda_))/m1/pow<3,1>(q)*
      (sqr(q)*(16.*sqr(q)-2.*sqr(m1)+8.*sqr(m2))+sqr(sqr(m1)-sqr(m2)));
  }
  else if(mecode==112) {
    test = 32.*sqr(k_)*m1/35./q*pow<6,1>(pcm)/sqr(fPi_*sqr(Lambda_));
  }
  else if(mecode==113) {
    test = 128.*sqr(k_)*m1/q*pow<6,1>(pcm)/105./sqr(fPi_*sqr(Lambda_));
  }
  else if(mecode==114) {
    test = 4.*sqr(gtilde_)*q/m1*sqr(pcm)/sqr(fPi_);
  }
  else if(mecode==115) {
    test = 4.*sqr(gtilde_)*m1*sqr(pcm)/3./sqr(fPi_)/q;
  }
  else if(mecode==116) {
    test = 8.*sqr(gtilde_)*m1*sqr(pcm)/3./sqr(fPi_)/q;
  }
  else
    assert(false);
  return test*pcm/8./Constants::pi*MEcoupling(imode)*MEcoupling(imode);
}
