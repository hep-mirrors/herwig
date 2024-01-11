// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2Onium1S0Amplitude class.
//

#include "GammaGamma2Onium1S0Amplitude.h"
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
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/EventHandler.h"

using namespace Herwig;

IBPtr GammaGamma2Onium1S0Amplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2Onium1S0Amplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2Onium1S0Amplitude::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2) << oenum(state_) << n_ << ounit(Lambda2_,GeV2) << mOpt_ << massGen_;
}

void GammaGamma2Onium1S0Amplitude::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2) >> ienum(state_) >> n_ >> iunit(Lambda2_,GeV2) >> mOpt_ >> massGen_;
}

void GammaGamma2Onium1S0Amplitude::doinit() {
  GammaGammaAmplitude::doinit();
  // get the non-perturbative ME
  O1_ = params_->singletMEProduction<0>(state_,n_,0,0);
  // get the mass generator of the onium state
  unsigned int iq = 4+state_;
  long id = iq*110+1 + (n_-1)*100000;
  tcPDPtr ps = getParticleData(id);
  if(!ps)
    throw Exception() << "No onium particle with id " << id
		      << "in " << fullName();
  if(ps->massGenerator())
    massGen_=dynamic_ptr_cast<GenericMassGeneratorPtr>(ps->massGenerator());
  if(!massGen_) mOpt_=0;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GammaGamma2Onium1S0Amplitude,GammaGammaAmplitude>
describeHerwigGammaGamma2Onium1S0Amplitude("Herwig::GammaGamma2Onium1S0Amplitude",
					   "HwOniumParameters.so HwMEGammaGamma.so HwMEGammaGammaOnium.so");

void GammaGamma2Onium1S0Amplitude::Init() {

  static ClassDocumentation<GammaGamma2Onium1S0Amplitude> documentation
    ("The GammaGamma2Onium1S0Amplitude class implements the amplitude for gamma gamma -> 1S0");

  static Reference<GammaGamma2Onium1S0Amplitude,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &GammaGamma2Onium1S0Amplitude::params_, false, false, true, false, false);
  
  static Switch<GammaGamma2Onium1S0Amplitude,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &GammaGamma2Onium1S0Amplitude::state_, ccbar, false, false);
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
  
  static Parameter<GammaGamma2Onium1S0Amplitude,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &GammaGamma2Onium1S0Amplitude::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Parameter<GammaGamma2Onium1S0Amplitude,Energy2> interfaceLambda2
    ("Lambda2",
     "The value of Lambda^2 for the form-factor",
     &GammaGamma2Onium1S0Amplitude::Lambda2_, GeV2, sqr(3.0969*GeV), 0.0*GeV2, 200.0*GeV2,
     false, false, Interface::limited);
  
  static Switch<GammaGamma2Onium1S0Amplitude,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the generation of the onium mass",
     &GammaGamma2Onium1S0Amplitude::mOpt_, 0, false, false);
  static SwitchOption interfaceMassOptionOnShell
    (interfaceMassOption,
     "OnShell",
     "Generate the onium state on-shell",
     0);
  static SwitchOption interfaceMassOptionOffShell
    (interfaceMassOption,
     "OffShell",
     "Generate an off-shell onium state using the mass generator",
     1);

}

vector<DiagPtr> GammaGamma2Onium1S0Amplitude::getDiagrams(unsigned int iopt) const {
  // construct the meson PDG code from quark ids
  unsigned int iq = 4+state_;
  tcPDPtr ps = getParticleData(long(iq*110+1 + (n_-1)*100000));
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

ProductionMatrixElement GammaGamma2Onium1S0Amplitude::
helicityAmplitude(const vector<VectorWaveFunction> & v1,
		  const vector<VectorWaveFunction> & v2,
		  const Energy & M, double & output) const {
  vector<unsigned int> ihMax(4,0);
  ProductionMatrixElement me = bookME(ihMax,v1.size(),v2.size(),vector<PDT::Spin>(1,PDT::Spin0));
  // calculate the matrix element
  output = 0;
  Lorentz5Momentum pG1 = v1[0].momentum();
  Lorentz5Momentum pG2 = v2[0].momentum();
  for(unsigned int ih1A=0;ih1A<ihMax[0];++ih1A) {
    for(unsigned int ih1B=0;ih1B<ihMax[1];++ih1B) {
      auto vOff = Helicity::epsilon(v1[2*ih1A+ih1B].wave(),pG1,pG2);
      for(unsigned int ih2A=0;ih2A<ihMax[2];++ih2A) {
  	for(unsigned int ih2B=0;ih2B<ihMax[3];++ih2B) {
  	  Complex amp = (vOff*v2[2*ih2A+ih2B].wave())/sqr(M);
  	  output += norm(amp);
	  if(v1.size()==2 && v2.size()==2)  me(2*ih1B,2*ih2B,0) = amp;
	  else if(v1.size()==2)             me(2*ih1B,ih2A,ih2B,0) = amp;
	  else if(v2.size()==2)             me(ih1A,ih1B,2*ih2B,0) = amp;
	  else                              me(ih1A,ih1B,ih2A,ih2B,0) = amp;
  	}
      }
    }
  }
  return me;
}

Energy GammaGamma2Onium1S0Amplitude::generateW(double r, const tcPDVector & partons, Energy Wmin,
					       Energy Wmax,Energy2 & jacW, Energy2 scale) {
  Wmin = max(Wmin,partons.back()->massMin());
  Wmax = min(Wmax,partons.back()->massMax());
  double wgt(0.);
  Energy output = massGen_->mass(wgt,*partons.back(),Wmin,Wmax,r);
  jacW = scale*wgt;
  return output;
}

double GammaGamma2Onium1S0Amplitude::
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
