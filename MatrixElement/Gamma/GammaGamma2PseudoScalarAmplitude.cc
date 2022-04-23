// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2PseudoScalarAmplitude class.
//

#include "GammaGamma2PseudoScalarAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

IBPtr GammaGamma2PseudoScalarAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2PseudoScalarAmplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2PseudoScalarAmplitude::persistentOutput(PersistentOStream & os) const {
  os << particle_ << ounit(F00_,1./GeV) << ounit(LambdaP2_,GeV2);
}

void GammaGamma2PseudoScalarAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> particle_ >> iunit(F00_,1./GeV) >> iunit(LambdaP2_,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GammaGamma2PseudoScalarAmplitude,GammaGammaAmplitude>
describeHerwigGammaGamma2PseudoScalarAmplitude("Herwig::GammaGamma2PseudoScalarAmplitude",
					       "HwMEGammaGamma.so");

void GammaGamma2PseudoScalarAmplitude::Init() {

  static ClassDocumentation<GammaGamma2PseudoScalarAmplitude> documentation
    ("The GammaGamma2PseudoScalarAmplitude class implements"
     " the amplitude for gamma gamma -> pseudoscalar");
  
  static Reference<GammaGamma2PseudoScalarAmplitude,ParticleData> interfaceParticle
    ("Particle",
     "The particle produced by the amplitude",
     &GammaGamma2PseudoScalarAmplitude::particle_, false, false, true, false, false);

  static Parameter<GammaGamma2PseudoScalarAmplitude,InvEnergy> interfaceF00
    ("F00",
     "The form factor at zero momentum transfer",
     &GammaGamma2PseudoScalarAmplitude::F00_, 1/GeV, 0.274/GeV, 0./GeV, 100./GeV,
     false, false, Interface::limited);

  static Parameter<GammaGamma2PseudoScalarAmplitude,Energy2> interfaceLambdaP2
    ("LambdaP2",
     "The square of the pole mass for the form factor",
     &GammaGamma2PseudoScalarAmplitude::LambdaP2_, GeV2, 0.6*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);

}

vector<DiagPtr> GammaGamma2PseudoScalarAmplitude::getDiagrams(unsigned int iopt) const {
  vector<DiagPtr> output;
  output.reserve(3);
  tcPDPtr g  = getParticleData(ParticleID::gamma );

  if(iopt==0) {
    output.push_back(new_ptr((Tree2toNDiagram(2), g, g, 1, particle_, -1)));
  }
  else {
    tcPDPtr ep = getParticleData(ParticleID::eplus );
    tcPDPtr em = getParticleData(ParticleID::eminus);
    output.push_back(new_ptr((Tree2toNDiagram(4), em, g, g, ep, 1, em, 3, ep, 2, particle_, -1)));
  }
  return output;
}

ProductionMatrixElement GammaGamma2PseudoScalarAmplitude::
helicityAmplitude(const vector<VectorWaveFunction> & v1,
		  const vector<VectorWaveFunction> & v2,
		  const Energy & M, double & output) const {
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
  for(unsigned int ih1A=0;ih1A<v1.size()/2;++ih1A) {
    for(unsigned int ih1B=0;ih1B<2;++ih1B) {
      auto vOff = Helicity::epsilon(v1[2*ih1A+ih1B].wave(),pG1,pG2);
      for(unsigned int ih2A=0;ih2A<v1.size()/2;++ih2A) {
  	for(unsigned int ih2B=0;ih2B<2;++ih2B) {
  	  Complex amp = (vOff*v2[2*ih2A+ih2B].wave())/sqr(M);
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

double GammaGamma2PseudoScalarAmplitude::me2(const vector<VectorWaveFunction> & v1,
					     const vector<VectorWaveFunction> & v2,
					     const Energy2 & t1, const Energy2 & t2,
					     const Energy2 & scale, 
					     const vector<Lorentz5Momentum> & momenta,
					     const cPDVector & , DVector &  ) const {
  Energy M  = momenta.back().mass();
  // calculate the matrix element
  double output(0.);
  helicityAmplitude(v1,v2,M,output);
  // form factor
  InvEnergy form = F00_/(1.-t1/LambdaP2_)/(1.-t2/LambdaP2_);
  // coupling factors
  double alpha = generator()->standardModel()->alphaEM();
  return 0.25*pow<4,1>(M)/scale*sqr(form)*output*sqr(alpha*4.*Constants::pi);
}
