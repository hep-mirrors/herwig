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
#include "ThePEG/Handlers/EventHandler.h"

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
    cPDPair in = generator()->eventHandler()->incoming();
    if(in.first->charged() && in.second->charged()) {
      for (long id : ids) {
	tcPDPtr pip = getParticleData(id);
	tcPDPtr pim = pip->CC();
	// first t-channel
	output.push_back(new_ptr((Tree2toNDiagram(5),in.first,gamma,pip,gamma,in.second,
				  1, in.first, 4, in.second, 2,pip, 3,pim,-1)));
	// interchange
	output.push_back(new_ptr((Tree2toNDiagram(5),in.first,gamma,pip,gamma,in.second,
				  1, in.first, 4, in.second, 3,pip, 2,pim,-2)));
      }
    }
  }
  return output;
}

ProductionMatrixElement
GammaGamma2PiPiAmplitude::helicityAmplitude(const vector<VectorWaveFunction> & v1,
					    const vector<VectorWaveFunction> & v2,
					    const vector<Lorentz5Momentum> & out,
					    double & output,
					    DVector & dweight) const {
  dweight = {0.,0.};
  vector<unsigned int> ihMax(4,0);
  ProductionMatrixElement me = bookME(ihMax,v1.size(),v2.size(),vector<PDT::Spin>(2,PDT::Spin0));
  // calculate the matrix element
  Energy2 off = -0.5*v1[0].momentum().m2(); 
  Energy2 d1  = off+v1[0].momentum()*out[0], d2 = off+v1[0].momentum()*out[1];
  Complex diag[2];
  for(unsigned int ih1A=0;ih1A<ihMax[0];++ih1A) {
    for(unsigned int ih1B=0;ih1B<ihMax[1];++ih1B) {
      unsigned int ih1 = 2*ih1A+ih1B;
      complex<Energy> a1 = out[0]*v1[ih1].wave();
      complex<Energy> a2 = out[1]*v1[ih1].wave();
      for(unsigned int ih2A=0;ih2A<ihMax[2];++ih2A) {
  	for(unsigned int ih2B=0;ih2B<ihMax[3];++ih2B) {
	  unsigned int ih2 = 2*ih2A+ih2B;
	  complex<Energy> b1 = out[0]*v2[ih2].wave();
	  complex<Energy> b2 = out[1]*v2[ih2].wave();
	  // t-channel diagram
	  diag[0] = a1*b2/d1;
	  // u-channel diagram
	  diag[1] = a2*b1/d2;
	  dweight[0] += norm(diag[0]);
	  dweight[1] += norm(diag[1]);
	  Complex amp = v1[ih1].wave()*v2[ih2].wave()-diag[0]-diag[1];
	  output += norm(amp);
	  if(v1.size()==2 && v2.size()==2)  me(2*ih1B,2*ih2B,0,0) = amp;
	  else if(v1.size()==2)             me(2*ih1B,ih2A,ih2B,0,0) = amp;
	  else if(v2.size()==2)             me(ih1A,ih1B,2*ih2B,0,0) = amp;
	  else                              me(ih1A,ih1B,ih2A,ih2B,0,0) = amp;
	}
      }
    }
  }
  return me;
}


double GammaGamma2PiPiAmplitude::me2(const vector<VectorWaveFunction> & v1,
				     const vector<VectorWaveFunction> & v2,
				     const Energy2 &, const Energy2 &,
				     const Energy2 & , 
				     const vector<Lorentz5Momentum> & out,
				     const cPDVector &,
				     DVector & dweights ) const {
  // matrix element
  double output(0.);
  helicityAmplitude(v1,v2,out,output,dweights);
  // prefactors
  double alpha = generator()->standardModel()->alphaEM();
  return output*sqr(alpha*4.*Constants::pi);
}

Energy GammaGamma2PiPiAmplitude::generateW(double r, const tcPDVector & partons, Energy Wmin, Energy Wmax,Energy2 & jacW, Energy2) {
  Wmin = max(Wmin,2.*partons[0]->constituentMass());
  Energy W = Wmin*pow(Wmax/Wmin,r);
  jacW = 2.*sqr(W)*log(Wmax/Wmin);
  return W;
}
