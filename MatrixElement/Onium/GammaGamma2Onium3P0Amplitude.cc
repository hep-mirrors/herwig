// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2Onium3P0Amplitude class.
//

#include "GammaGamma2Onium3P0Amplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/EventHandler.h"

using namespace Herwig;

IBPtr GammaGamma2Onium3P0Amplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2Onium3P0Amplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2Onium3P0Amplitude::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2*GeV2) << oenum(state_) << n_ << ounit(Lambda2_,GeV2) << mOpt_ << massGen_;
}

void GammaGamma2Onium3P0Amplitude::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2*GeV2) >> ienum(state_) >> n_ >> iunit(Lambda2_,GeV2) >> mOpt_ >> massGen_;
}

void GammaGamma2Onium3P0Amplitude::doinit() {
  GammaGammaAmplitude::doinit();
  // get the non-perturbative ME
  O1_ = params_->singletMEProduction<1>(state_,n_,1,0);
  // get the mass generator of the onium state
  unsigned int iq = 4+state_;
  long id = iq*110+10001 + (n_-1)*100000;
  tcPDPtr ps = getParticleData(id);
  if(!ps)
    throw Exception() << "No onium particle with id " << id << " in " << fullName();
  if(ps->massGenerator())
    massGen_=dynamic_ptr_cast<GenericMassGeneratorPtr>(ps->massGenerator());
  if(!massGen_) mOpt_=0;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GammaGamma2Onium3P0Amplitude,GammaGammaAmplitude>
describeHerwigGammaGamma2Onium3P0Amplitude("Herwig::GammaGamma2Onium3P0Amplitude",
					   "HwOniumParameters.so HwMEGammaGamma.so HwMEGammaGammaOnium.so");

void GammaGamma2Onium3P0Amplitude::Init() {

  static ClassDocumentation<GammaGamma2Onium3P0Amplitude> documentation
    ("The GammaGamma2Onium3P0Amplitude class implements the amplitude for gamma gamma -> 3P0");

  static Reference<GammaGamma2Onium3P0Amplitude,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &GammaGamma2Onium3P0Amplitude::params_, false, false, true, false, false);
  
  static Switch<GammaGamma2Onium3P0Amplitude,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &GammaGamma2Onium3P0Amplitude::state_, ccbar, false, false);
  static SwitchOption interfaceStateccbar
    (interfaceState,
     "ccbar",
     "Charmonium state",
     ccbar);
  static SwitchOption interfaceStatebbbar
    (interfaceState,
     "bbbar",
     "Bottomonium state",
     bbbar);
  
  static Parameter<GammaGamma2Onium3P0Amplitude,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &GammaGamma2Onium3P0Amplitude::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Parameter<GammaGamma2Onium3P0Amplitude,Energy2> interfaceLambda2
    ("Lambda2",
     "The value of Lambda^2 for the form-factor",
     &GammaGamma2Onium3P0Amplitude::Lambda2_, GeV2, sqr(3.0969*GeV), 0.0*GeV2, 200.0*GeV2,
     false, false, Interface::limited);
  
  static Switch<GammaGamma2Onium3P0Amplitude,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the generation of the onium mass",
     &GammaGamma2Onium3P0Amplitude::mOpt_, 0, false, false);
  static SwitchOption interfaceMassOptionOnShell
    (interfaceMassOption,
     "OnShell",
     "Generate the onium state on-shell",
     0);
  static SwitchOption interfaceMassOptionOffShell
    (interfaceMassOption,
     "OffShell",
     "Generate an off-shell onium state using the nass generator",
     1);

}

vector<DiagPtr> GammaGamma2Onium3P0Amplitude::getDiagrams(unsigned int iopt) const {
  // construct the meson PDG code from quark ids
  unsigned int iq = 4+state_;
  tcPDPtr ps = getParticleData(long(iq*110+10001 + (n_-1)*100000));
  // construct the diagrams  
  vector<DiagPtr> output;
  output.reserve(1);
  tcPDPtr g  = getParticleData(ParticleID::gamma );
  if(iopt==0) {
    output.push_back(new_ptr((Tree2toNDiagram(2), g, g, 1, ps, -1)));
  }
  else {
    cPDPair in = generator()->eventHandler()->incoming();
    if(in.first->charged() && in.second->charged())
      output.push_back(new_ptr((Tree2toNDiagram(4), in.first, g, g, in.second,
				1, in.first, 3, in.second, 2, ps, -1)));
  }
  return output;
}

ProductionMatrixElement GammaGamma2Onium3P0Amplitude::
helicityAmplitude(const vector<VectorWaveFunction> & v1,
		  const vector<VectorWaveFunction> & v2,
		  const Energy & M, double & output) const {
  Lorentz5Momentum pG1 = v1[0].momentum();
  Lorentz5Momentum pG2 = v2[0].momentum();
  ProductionMatrixElement me;
  if(v1.size()==4&&v2.size()==4) {
    me = ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
  				 PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  }
  else if(v1.size()==2&&v2.size()==2) {
    me = ProductionMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  }
  else
    assert(false);
  // calculate the matrix element
  output = 0;
  for(unsigned int ih1A=0;ih1A<v1.size()/2;++ih1A) {
    for(unsigned int ih1B=0;ih1B<2;++ih1B) {
      unsigned int ih1 = 2*ih1A+ih1B;
      complex<Energy> d1 = v1[ih1].wave()*pG2;
      for(unsigned int ih2A=0;ih2A<v1.size()/2;++ih2A) {
  	for(unsigned int ih2B=0;ih2B<2;++ih2B) {
  	  unsigned int ih2 = 2*ih2A+ih2B;
	  Complex amp = v1[ih1].wave()*v2[ih2].wave()  -2./sqr(M)*d1*(v2[ih2].wave()*pG1);
	  output += norm(amp);
  	  if(v1.size()==4) {
  	    me(ih1A,ih1B,ih2A,ih2B,0) = amp;
  	  }
  	  else {
  	    me(2*ih1B,2*ih2B,0) = amp;
  	  }
	}
      }
    }
  }
  return me;
}

Energy GammaGamma2Onium3P0Amplitude::generateW(double r, const tcPDVector & partons,Energy Wmax,Energy2 & jacW, Energy2 scale) {
  Energy Wmin = partons.back()->massMin();
  Wmax = min(Wmax,partons.back()->massMax());
  double wgt(0.);
  Energy output = massGen_->mass(wgt,*partons.back(),Wmin,Wmax,r);
  jacW = scale*wgt;
  return output;
}

double GammaGamma2Onium3P0Amplitude::
generateKinematics(const double * ,
		   const Energy2 & scale, 
		   vector<Lorentz5Momentum> & momenta,
		   const tcPDVector & ) {
  Energy M = sqrt(scale);
  double jac = scale*massGen_->BreitWignerWeight(M)/pow(Constants::twopi,3);
  momenta[0].setVect(Momentum3(ZERO,ZERO,ZERO));
  momenta[0].setE(M);
  momenta[0].rescaleMass();
  return jac;
}
