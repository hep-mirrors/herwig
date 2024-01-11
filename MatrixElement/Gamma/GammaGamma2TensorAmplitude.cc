// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2TensorAmplitude class.
//

#include "GammaGamma2TensorAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Handlers/EventHandler.h"

using namespace Herwig;

IBPtr GammaGamma2TensorAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2TensorAmplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2TensorAmplitude::persistentOutput(PersistentOStream & os) const {
  os << particle_ << ounit(FTT0_,1./GeV) << ounit(FTT2_,GeV) << ounit(LambdaP2_,GeV2) << mOpt_ << massGen_;
}

void GammaGamma2TensorAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> particle_ >> iunit(FTT0_,1./GeV) >> iunit(FTT2_,GeV) >> iunit(LambdaP2_,GeV2) >> mOpt_ >> massGen_;
}

void GammaGamma2TensorAmplitude::doinit() {
  GammaGammaAmplitude::doinit();
  if(particle_->iSpin()!=PDT::Spin2)
    throw Exception() << "Mustr have a spin-2 particle in GammaGamma2TensorAmplitude"
		      << Exception::runerror;
  if(particle_->massGenerator())
    massGen_=dynamic_ptr_cast<GenericMassGeneratorPtr>(particle_->massGenerator());
  if(!massGen_) mOpt_=0;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GammaGamma2TensorAmplitude,GammaGammaAmplitude>
describeHerwigGammaGamma2TensorAmplitude("Herwig::GammaGamma2TensorAmplitude",
					       "HwMEGammaGamma.so");

void GammaGamma2TensorAmplitude::Init() {

  static ClassDocumentation<GammaGamma2TensorAmplitude> documentation
    ("The GammaGamma2TensorAmplitude class implements"
     " the amplitude for gamma gamma -> pseudoscalar");
  
  static Reference<GammaGamma2TensorAmplitude,ParticleData> interfaceParticle
    ("Particle",
     "The particle produced by the amplitude",
     &GammaGamma2TensorAmplitude::particle_, false, false, true, false, false);

  static Parameter<GammaGamma2TensorAmplitude,InvEnergy> interfaceFTT0
    ("FTT0",
     "The form factor at zero momentum transfer",
     &GammaGamma2TensorAmplitude::FTT0_, 1./GeV, 8.89/MeV, 0./GeV, 100./GeV,
     false, false, Interface::limited);

  static Parameter<GammaGamma2TensorAmplitude,Energy> interfaceFTT2
    ("FTT2",
     "The form factor at zero momentum transfer",
     &GammaGamma2TensorAmplitude::FTT2_, GeV, 8.89*MeV, 0.*GeV, 100.*GeV,
     false, false, Interface::limited);

  static Parameter<GammaGamma2TensorAmplitude,Energy2> interfaceLambdaP2
    ("LambdaP2",
     "The square of the pole mass for the form factor",
     &GammaGamma2TensorAmplitude::LambdaP2_, GeV2, 0.6*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);
  
  static Switch<GammaGamma2TensorAmplitude,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the generation of the scalar mass",
     &GammaGamma2TensorAmplitude::mOpt_, 1, false, false);
  static SwitchOption interfaceMassOptionOnShell
    (interfaceMassOption,
     "OnShell",
     "Generate the scalar state on-shell",
     0);
  static SwitchOption interfaceMassOptionOffShell
    (interfaceMassOption,
     "OffShell",
     "Generate an off-shell scalar state using the mass generator",
     1);

}

vector<DiagPtr> GammaGamma2TensorAmplitude::getDiagrams(unsigned int iopt) const {
  vector<DiagPtr> output;
  output.reserve(3);
  tcPDPtr g  = getParticleData(ParticleID::gamma );

  if(iopt==0) {
    output.push_back(new_ptr((Tree2toNDiagram(2), g, g, 1, particle_, -1)));
  }
  else {
    cPDPair in = generator()->eventHandler()->incoming();
    if(in.first->charged() && in.second->charged())
      output.push_back(new_ptr((Tree2toNDiagram(4), in.first, g, g, in.second,
				1, in.first, 3, in.second, 2, particle_, -1)));
  }
  return output;
}

ProductionMatrixElement GammaGamma2TensorAmplitude::
helicityAmplitude(const vector<VectorWaveFunction> & v1,
		  const vector<VectorWaveFunction> & v2,
		  const vector<TensorWaveFunction> & ten,
		  const Energy2 & t1, const Energy2 & t2,
		  const Energy & M, double & output) const {
  vector<unsigned int> ihMax(4,0);
  ProductionMatrixElement me = bookME(ihMax,v1.size(),v2.size(),vector<PDT::Spin>(1,PDT::Spin2));
  // calculate the matrix element
  output = 0;
  Lorentz5Momentum pG1 = v1[0].momentum();
  Lorentz5Momentum pG2 = v2[0].momentum();
  Energy2 p1p2=pG1*pG2;
  Energy4 X = sqr(p1p2)-t1*t2;
  Lorentz5Momentum pDiff=pG1-pG2;
  for(unsigned int ih1A=0;ih1A<ihMax[0];++ih1A) {
    for(unsigned int ih1B=0;ih1B<ihMax[1];++ih1B) {
      complex<Energy> d1[2]={ v1[2*ih1A+ih1B].wave()*pG1, v1[2*ih1A+ih1B].wave()*pG2};
      LorentzPolarizationVector Q1 = (t1*pG2-p1p2*pG1)*d1[1]/X+v1[2*ih1A+ih1B].wave();
      for(unsigned int ih2A=0;ih2A<ihMax[2];++ih2A) {
     	for(unsigned int ih2B=0;ih2B<ihMax[2];++ih2B) {
	  complex<Energy> d2[2]={ v2[2*ih2A+ih2B].wave()*pG1, v2[2*ih2A+ih2B].wave()*pG2};
	  LorentzPolarizationVector Q2 = (t2*pG1-p1p2*pG2)*d2[0]/X+v2[2*ih2A+ih2B].wave();
	  Complex amp0 = v1[2*ih1A+ih1B].wave()*v2[2*ih2A+ih2B].wave()
	    + (-p1p2*(d1[1]*d2[0]+d1[0]*d2[1]) + t2*d1[0]*d2[0] + t1*d1[1]*d2[1])/X;
	  for(unsigned int ih3=0;ih3<5;++ih3) {
	    Complex amp = FTT0_*amp0/M*ten[ih3].wave().preDot(pDiff)*pDiff
	      +FTT2_/M*ten[ih3].wave().preDot(Q1)*Q2;
	    output += norm(amp);
	    if(v1.size()==2 && v2.size()==2)  me(2*ih1B,2*ih2B,ih3) = amp;
	    else if(v1.size()==2)             me(2*ih1B,ih2A,ih2B,ih3) = amp;
	    else if(v2.size()==2)             me(ih1A,ih1B,2*ih2B,ih3) = amp;
	    else                              me(ih1A,ih1B,ih2A,ih2B,ih3) = amp;
	  }
      	}
      }
    }
  }
  return me;
}

double GammaGamma2TensorAmplitude::me2(const vector<VectorWaveFunction> & v1,
				       const vector<VectorWaveFunction> & v2,
				       const Energy2 & t1, const Energy2 & t2,
				       const Energy2 & scale, 
				       const vector<Lorentz5Momentum> & momenta,
				       const cPDVector & partons, DVector &  ) const {
  Energy M  = momenta.back().mass();
  // wavefunction for the tensor meson
  vector<TensorWaveFunction> twave(5);
  for(unsigned int i=0;i<5;++i) {
    twave[i] = TensorWaveFunction(momenta.back(), partons.back(),i,outgoing);
  }
  // calculate the matrix element
  double output(0.);
  helicityAmplitude(v1,v2,twave,t1,t2,M,output);
  // form factor
  double form = 1./(1.-t1/LambdaP2_)/(1.-t2/LambdaP2_);
  // coupling factors
  double alpha = generator()->standardModel()->alphaEM();
  return 0.25*sqr(M*form)/scale*output*sqr(alpha*4.*Constants::pi);
}

Energy GammaGamma2TensorAmplitude::generateW(double r, const tcPDVector & partons, Energy Wmin,
					     Energy Wmax,Energy2 & jacW, Energy2 scale) {
  Wmin = max(Wmin,partons.back()->massMin());
  Wmax = min(Wmax,partons.back()->massMax());
  double wgt(0.);
  Energy output = massGen_->mass(wgt,*partons.back(),Wmin,Wmax,r);
  jacW = scale*wgt;
  return output;
}

double GammaGamma2TensorAmplitude::generateKinematics(const double * ,
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
