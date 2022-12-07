// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2ScalarAmplitude class.
//

#include "GammaGamma2ScalarAmplitude.h"
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

IBPtr GammaGamma2ScalarAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2ScalarAmplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2ScalarAmplitude::persistentOutput(PersistentOStream & os) const {
  os << particle_ << ounit(FTT_,GeV) << ounit(LambdaP2_,GeV2) << mOpt_ << massGen_;
}

void GammaGamma2ScalarAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> particle_ >> iunit(FTT_,GeV) >> iunit(LambdaP2_,GeV2) >> mOpt_ >> massGen_;
}

void GammaGamma2ScalarAmplitude::doinit() {
  GammaGammaAmplitude::doinit();
  if(particle_->iSpin()!=PDT::Spin0)
    throw Exception() << "Must have a spin-0 particle in GammaGamma2ScalarAmplitude"
		      << Exception::runerror;
  if(particle_->massGenerator())
    massGen_=dynamic_ptr_cast<GenericMassGeneratorPtr>(particle_->massGenerator());
  if(!massGen_) mOpt_=0;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GammaGamma2ScalarAmplitude,GammaGammaAmplitude>
describeHerwigGammaGamma2ScalarAmplitude("Herwig::GammaGamma2ScalarAmplitude",
					 "HwMEGammaGamma.so");

void GammaGamma2ScalarAmplitude::Init() {

  static ClassDocumentation<GammaGamma2ScalarAmplitude> documentation
    ("The GammaGamma2ScalarAmplitude class implements"
     " the amplitude for gamma gamma -> pseudoscalar");
  
  static Reference<GammaGamma2ScalarAmplitude,ParticleData> interfaceParticle
    ("Particle",
     "The particle produced by the amplitude",
     &GammaGamma2ScalarAmplitude::particle_, false, false, true, false, false);

  static Parameter<GammaGamma2ScalarAmplitude,Energy> interfaceFTT
    ("FTT",
     "The form factor at zero momentum transfer",
     &GammaGamma2ScalarAmplitude::FTT_, GeV, 8.89*MeV, 0.*GeV, 100.*GeV,
     false, false, Interface::limited);

  static Parameter<GammaGamma2ScalarAmplitude,Energy2> interfaceLambdaP2
    ("LambdaP2",
     "The square of the pole mass for the form factor",
     &GammaGamma2ScalarAmplitude::LambdaP2_, GeV2, 0.6*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);
  
  static Switch<GammaGamma2ScalarAmplitude,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the generation of the scalar mass",
     &GammaGamma2ScalarAmplitude::mOpt_, 1, false, false);
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

vector<DiagPtr> GammaGamma2ScalarAmplitude::getDiagrams(unsigned int iopt) const {
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

ProductionMatrixElement GammaGamma2ScalarAmplitude::
helicityAmplitude(const vector<VectorWaveFunction> & v1,
		  const vector<VectorWaveFunction> & v2,
		  const Energy2 & t1, const Energy2 & t2,
		  const Energy &, double & output) const {
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
  Lorentz5Momentum pG1 = v1[0].momentum();
  Lorentz5Momentum pG2 = v2[0].momentum();
  Energy2 p1p2=pG1*pG2;
  Energy4 X = sqr(p1p2)-t1*t2;
  for(unsigned int ih1A=0;ih1A<v1.size()/2;++ih1A) {
    for(unsigned int ih1B=0;ih1B<2;++ih1B) {
      complex<Energy> d1[2]={ v1[2*ih1A+ih1B].wave()*pG1, v1[2*ih1A+ih1B].wave()*pG2};
      for(unsigned int ih2A=0;ih2A<v1.size()/2;++ih2A) {
      	for(unsigned int ih2B=0;ih2B<2;++ih2B) {
	  complex<Energy> d2[2]={ v2[2*ih2A+ih2B].wave()*pG1, v2[2*ih2A+ih2B].wave()*pG2};
	  Complex amp = v1[2*ih1A+ih1B].wave()*v2[2*ih2A+ih2B].wave()
	    + (-p1p2*(d1[1]*d2[0]+d1[0]*d2[1]) + t2*d1[0]*d2[0] + t1*d1[1]*d2[1])/X;
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

double GammaGamma2ScalarAmplitude::me2(const vector<VectorWaveFunction> & v1,
					     const vector<VectorWaveFunction> & v2,
					     const Energy2 & t1, const Energy2 & t2,
					     const Energy2 & scale, 
					     const vector<Lorentz5Momentum> & momenta,
					     const cPDVector & , DVector &  ) const {
  Energy M  = momenta.back().mass();
  // calculate the matrix element
  double output(0.);
  helicityAmplitude(v1,v2,t1,t2,M,output);
  // form factor
  Energy form = FTT_/(1.-t1/LambdaP2_)/(1.-t2/LambdaP2_);
  // // coupling factors
  double alpha = generator()->standardModel()->alphaEM();
  return 0.25*sqr(form)/scale*output*sqr(alpha*4.*Constants::pi);
}

Energy GammaGamma2ScalarAmplitude::generateW(double r, const tcPDVector & partons,
					     Energy Wmax,Energy2 & jacW, Energy2 scale) {
  Energy Wmin = partons.back()->massMin();
  Wmax = min(Wmax,partons.back()->massMax());
  double wgt(0.);
  Energy output = massGen_->mass(wgt,*partons.back(),Wmin,Wmax,r);
  jacW = scale*wgt;
  return output;
}

double GammaGamma2ScalarAmplitude::generateKinematics(const double * ,
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
