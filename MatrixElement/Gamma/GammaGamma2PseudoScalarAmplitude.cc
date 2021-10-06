// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2PseudoScalarAmplitude class.
//

#include "GammaGamma2PseudoScalarAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

using namespace Herwig;

IBPtr GammaGamma2PseudoScalarAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2PseudoScalarAmplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2PseudoScalarAmplitude::persistentOutput(PersistentOStream & os) const {
  os << ounit(F00_,1./GeV) << ounit(LambdaP2_,GeV2);
}

void GammaGamma2PseudoScalarAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> iunit(F00_,1./GeV) >> iunit(LambdaP2_,GeV2);
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

  static ParVector<GammaGamma2PseudoScalarAmplitude,InvEnergy> interfaceF00
    ("F00",
     "The form factor at zero momentum transfer",
     &GammaGamma2PseudoScalarAmplitude::F00_, 1./GeV, 3, 1./GeV, 0./GeV, 100./GeV,
     false, false, Interface::limited);

  static ParVector<GammaGamma2PseudoScalarAmplitude,Energy2> interfaceLambdaP2
    ("LambdaP2",
     "The square of the pole mass for the form factor",
     &GammaGamma2PseudoScalarAmplitude::LambdaP2_, GeV2, 3, 1.0*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);

}

vector<DiagPtr> GammaGamma2PseudoScalarAmplitude::getDiagrams(unsigned int iopt) const {
  vector<DiagPtr> output;
  output.reserve(3);
  tcPDPtr g  = getParticleData(ParticleID::gamma );
  if(iopt==0) {
    for(long pid=111; pid<340; pid+=110) {
      tcPDPtr ps = getParticleData(pid);
      output.push_back(new_ptr((Tree2toNDiagram(2), g, g, 1, ps, -1)));
    }
  }
  else {
    tcPDPtr ep = getParticleData(ParticleID::eplus );
    tcPDPtr em = getParticleData(ParticleID::eminus);
    for(long pid=111; pid<340; pid+=110) {
      tcPDPtr ps = getParticleData(pid);
      output.push_back(new_ptr((Tree2toNDiagram(4), em, g, g, ep, 1, em, 3, ep, 2, ps, -1)));
    }
  }
  return output;
}

double GammaGamma2PseudoScalarAmplitude::me2(const vector<VectorWaveFunction> & v1,
					     const vector<VectorWaveFunction> & v2,
					     const Energy2 & t1, const Energy2 & t2,
					     const Energy2 & scale, 
					     const vector<Lorentz5Momentum> &,
					     const cPDVector & partons,
					     DVector &  ) const {
  // form factor
  int iloc = partons[0]->id()/100 -1;
  InvEnergy form = F00_[iloc]/(1.-t1/LambdaP2_[iloc])/(1.-t2/LambdaP2_[iloc]);
  // calculate the matrix element
  double output(0.);
  Lorentz5Momentum pG1 = v1[0].momentum();
  Lorentz5Momentum pG2 = v2[0].momentum();
  Energy rs = sqrt(scale);
  for(unsigned int ih1=0;ih1<v1.size();++ih1) {
    auto vOff = Helicity::epsilon(v1[ih1].wave(),pG1,pG2);
    for(unsigned int ih2=0;ih2<v2.size();++ih2) {
      Complex amp = form*(vOff*v2[ih2].wave())/rs;
      output += norm(amp);
    }
  }
  // coupling factors
  double alpha = generator()->standardModel()->alphaEM();
  return 0.25*output*sqr(alpha*4.*Constants::pi);
}
