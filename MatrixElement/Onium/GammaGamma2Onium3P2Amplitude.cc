// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2Onium3P2Amplitude class.
//

#include "GammaGamma2Onium3P2Amplitude.h"
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
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Handlers/EventHandler.h"

using namespace Herwig;

IBPtr GammaGamma2Onium3P2Amplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2Onium3P2Amplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2Onium3P2Amplitude::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2*GeV2) << oenum(state_)
     << n_ << ounit(Lambda2_,GeV2) << mOpt_ << massGen_;
}

void GammaGamma2Onium3P2Amplitude::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2*GeV2) >> ienum(state_)
     >> n_ >> iunit(Lambda2_,GeV2) >> mOpt_ >> massGen_;
}

void GammaGamma2Onium3P2Amplitude::doinit() {
  GammaGammaAmplitude::doinit();
  // get the non-perturbative ME
  O1_ = params_->singletMEProduction<1>(state_,n_,1,2);
  // get the mass generator of the onium state
  unsigned int iq = 4+state_;
  long id = iq*110+5 + (n_-1)*100000;
  tcPDPtr ps = getParticleData(id);
  if(!ps)
    throw Exception() << "No onium particle with id " << id << " in " << fullName();
  if(ps->massGenerator())
    massGen_=dynamic_ptr_cast<GenericMassGeneratorPtr>(ps->massGenerator());
  if(!massGen_) mOpt_=0;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GammaGamma2Onium3P2Amplitude,GammaGammaAmplitude>
describeHerwigGammaGamma2Onium3P2Amplitude("Herwig::GammaGamma2Onium3P2Amplitude",
					   "HwOniumParameters.so HwMEGammaGamma.so HwMEGammaGammaOnium.so");

void GammaGamma2Onium3P2Amplitude::Init() {

  static ClassDocumentation<GammaGamma2Onium3P2Amplitude> documentation
    ("The GammaGamma2Onium3P2Amplitude class implements the amplitude for gamma gamma -> 3P2");


  static Reference<GammaGamma2Onium3P2Amplitude,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &GammaGamma2Onium3P2Amplitude::params_, false, false, true, false, false);
  
  static Switch<GammaGamma2Onium3P2Amplitude,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &GammaGamma2Onium3P2Amplitude::state_, ccbar, false, false);
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

  static Parameter<GammaGamma2Onium3P2Amplitude,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &GammaGamma2Onium3P2Amplitude::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Parameter<GammaGamma2Onium3P2Amplitude,Energy2> interfaceLambda2
    ("Lambda2",
     "The value of Lambda^2 for the form-factor",
     &GammaGamma2Onium3P2Amplitude::Lambda2_, GeV2, sqr(3.0969*GeV), 0.0*GeV2, 200.0*GeV2,
     false, false, Interface::limited);

  static Switch<GammaGamma2Onium3P2Amplitude,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the generation of the onium mass",
     &GammaGamma2Onium3P2Amplitude::mOpt_, 0, false, false);
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

vector<DiagPtr> GammaGamma2Onium3P2Amplitude::getDiagrams(unsigned int iopt) const {
  // construct the meson PDG code from quark ids
  unsigned int iq = 4+state_;
  tcPDPtr ps = getParticleData(long(iq*110+5 + (n_-1)*100000));
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

ProductionMatrixElement GammaGamma2Onium3P2Amplitude::
helicityAmplitude(const vector<VectorWaveFunction> & v1,
		  const vector<VectorWaveFunction> & v2,
		  const vector<TensorWaveFunction> & twave,
		  const Energy & M, double & output) const {
  vector<unsigned int> ihMax(4,0);
  ProductionMatrixElement me = bookME(ihMax,v1.size(),v2.size(),vector<PDT::Spin>(1,PDT::Spin2));
  Lorentz5Momentum pG1 = v1[0].momentum();
  Lorentz5Momentum pG2 = v2[0].momentum();
  Lorentz5Momentum pDiff = pG1-pG2;
  // calculate the amplitudes
  for(unsigned int ix=0;ix<5;++ix) {
    auto vPre  = twave[ix].wave().preDot (pDiff);
    auto vPost = twave[ix].wave().postDot(pDiff);
    complex<Energy2> v1v2=vPre*pDiff;
    for(unsigned int ih1A=0;ih1A<ihMax[0];++ih1A) {
      for(unsigned int ih1B=0;ih1B<ihMax[1];++ih1B) {
	unsigned int ih1 = 2*ih1A+ih1B;
	auto vEps1 = twave[ix].wave().preDot(v1[ih1].wave())+twave[ix].wave().postDot(v1[ih1].wave());
	complex<Energy> dPreEps1  = vPre*v1[ih1].wave();
	complex<Energy> dPostEps1 = vPost*v1[ih1].wave();
	complex<Energy> d1 = v1[ih1].wave()*pG2;
	for(unsigned int ih2A=0;ih2A<ihMax[2];++ih2A) {
	  for(unsigned int ih2B=0;ih2B<ihMax[3];++ih2B) {
	    unsigned int ih2 = 2*ih2A+ih2B;
	    complex<Energy> d2 = v2[ih2].wave()*pG1;
	    Complex d12 = v1[ih1].wave()*v2[ih2].wave();
	    Complex amp = ((dPreEps1+dPostEps1)*d2 -(vPre*v2[ih2].wave()+vPost*v2[ih2].wave())*d1-v1v2*d12)/sqr(M) + vEps1*v2[ih2].wave();
	    output += norm(amp);
	    if(v1.size()==2 && v2.size()==2)  me(2*ih1B,2*ih2B,ix) = amp;
	    else if(v1.size()==2)             me(2*ih1B,ih2A,ih2B,ix) = amp;
	    else if(v2.size()==2)             me(ih1A,ih1B,2*ih2B,ix) = amp;
	    else                              me(ih1A,ih1B,ih2A,ih2B,ix) = amp;
	  }
	}
      }
    }
  }
  return me;
}

double GammaGamma2Onium3P2Amplitude::me2(const vector<VectorWaveFunction> & v1,
					 const vector<VectorWaveFunction> & v2,
					 const Energy2 & t1, const Energy2 & t2,
					 const Energy2 & scale, 
					 const vector<Lorentz5Momentum> & momenta,
					 const cPDVector & partons, DVector &  ) const {
  double output(0.);
  Energy M  = momenta.back().mass();
  // wavefunction for the tensor meson
  vector<TensorWaveFunction> twave(5);
  for(unsigned int i=0;i<5;++i) {
    twave[i] = TensorWaveFunction(momenta.back(), partons.back(),i,outgoing);
  }
  // calculate the matrix element
  helicityAmplitude(v1,v2,twave,M,output);
  // coupling factors
  double eQ = state_==ccbar ? 2./3. : -1./3.;
  double alpha = generator()->standardModel()->alphaEM();
  return 128./5.*output*O1_/pow<3,1>(M)/scale*sqr(Constants::pi*alpha*sqr(eQ)/(1.-t1/Lambda2_)/(1.-t2/Lambda2_));
}

Energy GammaGamma2Onium3P2Amplitude::generateW(double r, const tcPDVector & partons, Energy Wmin,
					       Energy Wmax,Energy2 & jacW, Energy2 scale) {
  Wmin = max(Wmin,partons.back()->massMin());
  Wmax = min(Wmax,partons.back()->massMax());
  double wgt(0.);
  Energy output = massGen_->mass(wgt,*partons.back(),Wmin,Wmax,r);
  jacW = scale*wgt;
  return output;
}

ProductionMatrixElement GammaGamma2Onium3P2Amplitude::me(const vector<VectorWaveFunction> & v1,
							 const vector<VectorWaveFunction> & v2,
							 tParticleVector & particles) const {
  vector<TensorWaveFunction> t3;
  TensorWaveFunction(t3,particles[0],outgoing,true,false);
  double output(0);
  return helicityAmplitude(v1,v2,t3,particles[0]->mass(),output);
}


double GammaGamma2Onium3P2Amplitude::
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
