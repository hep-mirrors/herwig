// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2PiPiAmplitude class.
//

#include "GammaGamma2PiPiAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

using namespace Herwig;

IBPtr GammaGamma2PiPiAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2PiPiAmplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2PiPiAmplitude::persistentOutput(PersistentOStream & os) const {
  os << mode_;
}

void GammaGamma2PiPiAmplitude::persistentInput(PersistentIStream & is , int) {
  is >> mode_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GammaGamma2PiPiAmplitude,GammaGammaAmplitude>
describeHerwigGammaGamma2PiPiAmplitude("Herwig::GammaGamma2PiPiAmplitude",
				       "HwMEGammaGamma.so");

void GammaGamma2PiPiAmplitude::Init() {

  static ClassDocumentation<GammaGamma2PiPiAmplitude> documentation
    ("The GammaGamma2PiPiAmplitude class implements the amplitude for gamma gamma -> pi+pi-");

  
  static Switch<GammaGamma2PiPiAmplitude,unsigned int> interfaceMode
    ("Mode",
     "Which particles to produce",
     &GammaGamma2PiPiAmplitude::mode_, 0, false, false);
  static SwitchOption interfaceModeAll
    (interfaceMode,
     "All",
     "Produce all pions and kaons",
     0);
  static SwitchOption interfaceModePiPi
    (interfaceMode,
     "PiPi",
     "Produce pi+pi-",
     1);
  static SwitchOption interfaceModeKK
    (interfaceMode,
     "KK",
     "Produce K+K-",
     2);

}

vector<DiagPtr> GammaGamma2PiPiAmplitude::getDiagrams(unsigned int iopt) const {
  vector<DiagPtr> output;
  output.reserve(2);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  vector<long> ids; ids.reserve(2);
  if(mode_==0 || mode_==1) ids.push_back(211);
  if(mode_==0 || mode_==2) ids.push_back(321);
  // gamma gamma process
  if(iopt==0) {
    // first t-channel
    for (long id : ids) {
      tcPDPtr pip = getParticleData(id);
      tcPDPtr pim = pip->CC();
      output.push_back(new_ptr((Tree2toNDiagram(3),gamma,pip,gamma,1,pip, 2,pim,-1)));
      // interchange
      output.push_back(new_ptr((Tree2toNDiagram(3),gamma,pip,gamma,2,pip, 1,pim,-2)));
    }
  }
  // e+ e- > e+ e- pi+pi-
  else {
    tcPDPtr ep = getParticleData(ParticleID::eplus );
    tcPDPtr em = getParticleData(ParticleID::eminus);
    for (long id : ids) {
      tcPDPtr pip = getParticleData(id);
      tcPDPtr pim = pip->CC();
      // first t-channel
      output.push_back(new_ptr((Tree2toNDiagram(5),em,gamma,pip,gamma,ep, 1, em, 4, ep, 2,pip, 3,pim,-1)));
      // interchange
      output.push_back(new_ptr((Tree2toNDiagram(5),em,gamma,pip,gamma,ep, 1, em, 4, ep, 3,pip, 2,pim,-2)));
    }
  }
  return output;
}

double GammaGamma2PiPiAmplitude::me2(const vector<VectorWaveFunction> & v1,
				     const vector<VectorWaveFunction> & v2,
				     const Energy2 &, const Energy2 &,
				     const Energy2 & , 
				     const vector<Lorentz5Momentum> & out,
				     const cPDVector &,
				     DVector & dweights ) const {
  // scale (external photons so scale in couplings is 0)
  Energy2 mt(0.*GeV2);
  double output(0.),sumdiag[2]={0.,0.};
  Complex diag[2];
  Energy2 off = -0.5*v1[0].momentum().m2(); 
  Energy2 d1  = off+v1[0].momentum()*out[0], d2 = off+v1[0].momentum()*out[1];
  for(unsigned int ihel1=0;ihel1<v1.size();++ihel1) {
    complex<Energy> a1 = out[0]*v1[ihel1].wave();
    complex<Energy> a2 = out[1]*v1[ihel1].wave();
    for(unsigned int ihel2=0;ihel2<v2.size();++ihel2) {
      complex<Energy> b1 = out[0]*v2[ihel2].wave();
      complex<Energy> b2 = out[1]*v2[ihel2].wave();
      // t-channel diagram
      diag[0] = a1*b2/d1;
      // u-channel diagram
      diag[1] = a2*b1/d2;
      sumdiag[0] += norm(diag[0]);
      sumdiag[1] += norm(diag[1]);
      Complex amp = v1[ihel1].wave()*v2[ihel2].wave()-diag[0]-diag[1];
      output += norm(amp);
    }
  }
  // diagrams
  dweights.push_back(sumdiag[0]);
  dweights.push_back(sumdiag[1]);
  double alpha = generator()->standardModel()->alphaEM();
  return output*sqr(alpha*4.*Constants::pi);
}

Energy GammaGamma2PiPiAmplitude::generateW(double r, const tcPDVector & partons,Energy Wmax,Energy2 & jacW) {
  // Energy Wmin = 2.*partons[0]->constituentMass();
  // Energy W = Wmin+r*(Wmax-Wmin);
  // jacW = 2.*W*(Wmax-Wmin);
  // return W;
  Energy Wmin = 2.*partons[0]->constituentMass();
  Energy W = Wmin*pow(Wmax/Wmin,r);
  jacW = 2.*sqr(W)*log(Wmax/Wmin);
  return W;
}
