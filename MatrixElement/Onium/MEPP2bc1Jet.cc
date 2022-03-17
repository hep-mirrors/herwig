// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2bc1Jet class.
//

#include "MEPP2bc1Jet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

double MEPP2bc1Jet::me2() const {
  // return value
  Complex me2Sum(0.);
  vector<double> flows(2,0.);
  // colour matrix
  static const double cMatrix[2][2]= {{8.,-4.},{-4.,8.}};
  // diquark wavefunction
  VectorWaveFunction diquarkw(rescaledMomenta()[2],mePartonData()[2],outgoing);
  vector<VectorWaveFunction> v3;
  for(unsigned int ix=0;ix<3;++ix) {
    diquarkw.reset(ix);
    v3.push_back(diquarkw);
  }
  Energy2 sh=(rescaledMomenta()[0]+rescaledMomenta()[1]).m2();
  Energy2 th=(rescaledMomenta()[0]-rescaledMomenta()[2]).m2();
  Energy2 uh=(rescaledMomenta()[0]-rescaledMomenta()[3]).m2();
  Energy  M = rescaledMomenta()[2].mass();
  Energy2 M2=sqr(M);
  // g  initiated
  if(mePartonData()[0]->id()==ParticleID::g) {
    // gluon wavefunction
    VectorWaveFunction g1w(rescaledMomenta()[0],mePartonData()[0],incoming);
    vector<VectorWaveFunction> g1;
    for(unsigned int ix=0;ix<2;++ix) {
      // gauge test g1w.reset(10);
      g1w.reset(2*ix);
      g1.push_back(g1w);
    }
    SpinorBarWaveFunction q2w(rescaledMomenta()[1],mePartonData()[1],incoming);
    SpinorWaveFunction    q4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
    vector<SpinorBarWaveFunction> vbar2;
    vector<SpinorWaveFunction> v4;
    for(unsigned int ix=0;ix<2;++ix) {
      q2w.reset(ix);
      vbar2.push_back(q2w);
      q4w.reset(ix);
      v4.push_back(q4w);
    }
    // g q -> diquark qbar
    // mass ratio
    double a1  = rescaledMomenta()[1].mass()/M;
    double a12 = sqr(a1);
    double a2  = 1.-a1;
    double a22 = sqr(a2);
    if(mePartonData()[1]->id()>0) {
      // matrix element squared
      // loop over the helicities
      for(unsigned int ih2=0;ih2<2;++ih2) {
       	for(unsigned int ih4=0;ih4<2;++ih4) {
	  complex<Energy> dot6=v4[ih4].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	  complex<Energy2> dot12 = vec1*rescaledMomenta()[0];
      	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    complex<Energy> dot1 = rescaledMomenta()[0]*v3[ih3].wave();
	    complex<Energy> dot3 = rescaledMomenta()[1]*v3[ih3].wave();
	    complex<Energy> dot7 = vec1*v3[ih3].wave();
	    LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	    complex<Energy2> dot9 = vec2*rescaledMomenta()[0];
       	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      complex<Energy> dot2 = rescaledMomenta()[1]*g1[ih1].wave();
	      complex<Energy> dot4 = rescaledMomenta()[3]*g1[ih1].wave();
	      Complex dot5 = g1[ih1].wave()*v3[ih3].wave();
	      LorentzPolarizationVectorE vec3 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	      LorentzPolarizationVectorE vec4 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g1[ih1].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	      complex<Energy2> dot8 = vec4*rescaledMomenta()[0];
	      complex<Energy> dot10 = vec2*g1[ih1].wave();
	      complex<Energy2> dot11 = vec3*rescaledMomenta()[0];
	      complex<Energy> dot13 = vec1*g1[ih1].wave();
	      // now for the diagrams
	      Complex diag[5];
	      diag[0] = -(2.*dot2*dot7 - dot8)*M2/(a2*sqr(-(a12*M2) + sh));
	      diag[1] =  M*(4.*dot11*dot3 + 2.*dot1*(dot11 + 2.*dot2*dot6) - 2.*dot4*dot9 + 2.*a1*dot4*dot7*M + dot8*M + dot10*M2 - 2.*a1*dot10*M2 + 2.*a12*dot5*dot6*M2 - 2.*dot2*(dot9 + a1*dot7*M) + dot10*sh - 2.*dot5*dot6*sh - dot10*uh)/(a2*a1*(a12*M2 - sh)*(M2 - th));
	      diag[2] = -M*(4.*dot11*dot3 - 2.*dot4*dot9 - 2.*dot12*dot5*M - 2.*dot4*dot7*M + 2.*a1*dot4*dot7*M + dot8*M + dot10*M2 - 2.*a1*dot10*M2 + 2.*a12*dot5*dot6*M2 + 2.*dot1*(dot11 + 2.*dot2*dot6 + dot13*M) - 2.*dot2*(dot9 + -a2*dot7*M) + dot10*sh - 2.*dot5*dot6*sh - dot10*uh)/(a2*a1*(M2 - th)*(a22*M2 - uh));
	      diag[3] = -(2.*dot1*dot13 - 2.*dot12*dot5 - 2.*dot4*dot7 + dot8)*M2/(a1*sqr(-(a22*M2) + uh));
	      diag[4] = -M*(4.*dot11*dot3 + 2.*dot1*(dot11 + 2.*dot2*dot6) - 2.*dot4*dot9 + 2.*a1*dot4*dot7*M + dot10*M2 - 2.*a1*dot10*M2 + 2.*a12*dot5*dot6*M2 - 2.*dot2*(dot9 + -a2*dot7*M) + dot10*sh - 2.*dot5*dot6*sh - dot10*uh)/(a2*a1*(a12*M2 - sh)*(a22*M2 - uh));
	      Complex flow[2];
	      flow[0] = 4.*diag[0]/3.+diag[2]+diag[4]+diag[1]/3.;
	      flow[1] = 4.*diag[3]/3.+diag[1]-diag[4]+diag[2]/3.;
	      for(unsigned int ix=0;ix<2;++ix) {
		flows[ix]+=norm(flow[ix]);
		for(unsigned int iy=0;iy<2;++iy) {
		  me2Sum+=cMatrix[ix][iy]*flow[ix]*std::conj(flow[iy]);
		}
	      }
	    }
	  }
	}
      }
    }
    // g qbar -> antidiquark q
    else {
      me2Sum=0.;
      // matrix element squared
      // loop over the helicities
      for(unsigned int ih2=0;ih2<2;++ih2) {
       	for(unsigned int ih4=0;ih4<2;++ih4) {
	  complex<Energy> dot6=v4[ih4].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	  complex<Energy2> dot12 = vec1*rescaledMomenta()[0];
      	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    complex<Energy> dot1 = rescaledMomenta()[0]*v3[ih3].wave();
	    complex<Energy> dot3 = rescaledMomenta()[1]*v3[ih3].wave();
	    LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	    complex<Energy2> dot9 = vec2*rescaledMomenta()[0];
       	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      complex<Energy> dot2 = rescaledMomenta()[1]*g1[ih1].wave();
	      complex<Energy> dot4 = rescaledMomenta()[3]*g1[ih1].wave();
	      Complex dot5 = g1[ih1].wave()*v3[ih3].wave();
	      LorentzPolarizationVectorE vec3 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	      LorentzPolarizationVectorE vec4 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g1[ih1].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	      complex<Energy> dot7 = vec1*v3[ih3].wave();
	      complex<Energy2> dot8 = vec4*rescaledMomenta()[0];
	      complex<Energy> dot10 = vec2*g1[ih1].wave();
	      complex<Energy2> dot11 = vec3*rescaledMomenta()[0];
	      complex<Energy> dot13 = vec1*g1[ih1].wave();
	      // now for the diagrams
	      Complex diag[5];
	      diag[0] = -(2.*dot2*dot7 - dot8)*M2/(a2*sqr(-(a12*M2) + sh));
	      diag[1] = M*(4.*dot11*dot3 + 2.*dot1*(dot11 + 2.*dot2*dot6) - 2.*dot4*dot9 + 2.*a1*dot4*dot7*M + dot8*M + dot10*M2 - 2.*a1*dot10*M2 + 2.*a12*dot5*dot6*M2 - 2.*dot2*(dot9 + a1*dot7*M) + dot10*sh - 2.*dot5*dot6*sh - dot10*uh)/(a2*a1*(a12*M2 - sh)*(M2 - th));
	      diag[2] = (M*(4.*dot11*dot3 - 2.*dot4*dot9 - 2.*dot12*dot5*M - 2.*dot4*dot7*M + 2.*a1*dot4*dot7*M + dot8*M + dot10*M2 - 2.*a1*dot10*M2 + 2.*a12*dot5*dot6*M2 + 2.*dot1*(dot11 + 2.*dot2*dot6 + dot13*M) - 2.*dot2*(dot9 + -a2*dot7*M) + dot10*sh - 2.*dot5*dot6*sh - dot10*uh))/(-a2*a1*(M2 - th)*(a22*M2 - uh));
	      diag[3] = -(2.*dot1*dot13 - 2.*dot12*dot5 - 2.*dot4*dot7 + dot8)*M2/(a1*sqr(-(a22*M2) + uh));
	      diag[4] = -M*(4.*dot11*dot3 + 2.*dot1*(dot11 + 2.*dot2*dot6) - 2.*dot4*dot9 + 2.*a1*dot4*dot7*M + dot10*M2 - 2.*a1*dot10*M2 + 2.*a12*dot5*dot6*M2 - 2.*dot2*(dot9 + -a2*dot7*M) + dot10*sh - 2.*dot5*dot6*sh - dot10*uh)/(a2*a1*(a12*M2 - sh)*(a22*M2 - uh));
	      Complex flow[2];
	      flow[0] = 4.*diag[0]/3.+diag[2]+diag[4]+diag[1]/3.;
	      flow[1] = 4.*diag[3]/3.+diag[1]-diag[4]+diag[2]/3.;
	      for(unsigned int ix=0;ix<2;++ix) {
		flows[ix]+=norm(flow[ix]);
		for(unsigned int iy=0;iy<2;++iy) {
		  me2Sum+=cMatrix[ix][iy]*flow[ix]*std::conj(flow[iy]);
		}
	      }
	    }
	  }
	}
      }
    }
    // prefactors, spin and colour average
    me2Sum *= R02_*sqr(Constants::pi)*pow<3,1>(standardModel()->alphaS(scale())/M)/48.;
  }
  else {
    // gluon wavefunction
    VectorWaveFunction g4w(rescaledMomenta()[3],mePartonData()[3],incoming);
    vector<VectorWaveFunction> g4;
    for(unsigned int ix=0;ix<2;++ix) {
      // gauge test g4w.reset(10);
      g4w.reset(2*ix);
      g4.push_back(g4w);
    }
    SpinorBarWaveFunction q1w(rescaledMomenta()[0],mePartonData()[0],incoming);
    SpinorWaveFunction    q2w(rescaledMomenta()[1],mePartonData()[1],outgoing);
    vector<SpinorBarWaveFunction> vbar1;
    vector<SpinorWaveFunction> v2;
    for(unsigned int ix=0;ix<2;++ix) {
      q1w.reset(ix);
      vbar1.push_back(q1w);
      q2w.reset(ix);
      v2.push_back(q2w);
    }
    // mass ratio
    double a1 = rescaledMomenta()[0].mass()/M;
    double a12 = sqr(a1);
    double a2 = 1.-a1;
    double a22 = sqr(a2);
    me2Sum=0.;
    // q q -> diquark g
    if(mePartonData()[0]->id()>0) {
      for(unsigned int ih1=0;ih1<2;++ih1) {
	for(unsigned int ih2=0;ih2<2;++ih2) {
	  complex<Energy> dot6=v2[ih2].dimensionedWave().scalar(vbar1[ih1].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = v2[ih2].dimensionedWave().vectorCurrent(vbar1[ih1].dimensionedWave());
	  complex<Energy2> dot13 = vec1*rescaledMomenta()[3];
	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    complex<Energy> dot2 = rescaledMomenta()[0]*v3[ih3].wave();
	    complex<Energy> dot4 = rescaledMomenta()[1]*v3[ih3].wave();
	    LorentzPolarizationVectorE vec3 = v2[ih2].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(vbar1[ih1].dimensionedWave());
	    complex<Energy2> dot10 = vec3*rescaledMomenta()[3];
	    for(unsigned int ih4=0;ih4<2;++ih4) {
	      complex<Energy> dot1 = rescaledMomenta()[0]*g4[ih4].wave();
	      complex<Energy> dot3 = rescaledMomenta()[1]*g4[ih4].wave();
	      Complex dot5 = v3[ih3].wave()*g4[ih4].wave();
	      LorentzPolarizationVectorE vec2 = v2[ih2].dimensionedWave().slash(g4[ih4].wave()).vectorCurrent(vbar1[ih1].dimensionedWave());
	      LorentzPolarizationVectorE vec4 = v2[ih2].dimensionedWave().slash(g4[ih4].wave()).slash(v3[ih3].wave()).vectorCurrent(vbar1[ih1].dimensionedWave());
	      complex<Energy> dot7 = vec1*v3[ih3].wave();
	      complex<Energy2> dot8 = vec2*rescaledMomenta()[3];
	      complex<Energy> dot9 = vec2*v3[ih3].wave();
	      complex<Energy> dot11 = vec1*g4[ih4].wave();
	      complex<Energy2> dot12 = vec4*rescaledMomenta()[3];
	      Complex diag[5];
	      diag[0] = -M*(-2.*dot1*dot10 + 2.*dot10*dot3 + 4.*dot1*(dot2 + dot4)*dot6 + 2.*dot2*dot8 - 2.*dot4*dot8 - 2.*a2*dot1*dot7*M + 2.*a1*dot3*dot7*M - 2.*dot5*dot6*M2 + 4.*a1*dot5*dot6*M2 - 2.*a12*dot5*dot6*M2 + dot9*M2 - 2.*a1*dot9*M2 + 2.*dot5*dot6*th - dot9*th + dot9*uh)/(a2*a1*(a22*M2 - th)*(a12*M2 - uh));
	      diag[1] = -(((dot12 - 2*dot11*(dot2 + dot4) + 2*dot3*dot7)*M2)/(a1*sqr(-(a22*M2) + th)));
	      diag[2] = -M*(-2.*dot10*dot3 + 2.*(-dot2 + dot4)*dot8 + dot12*M - 2.*dot11*(dot2 + dot4)*M + 2.*dot3*dot7*M - 2.*a1*dot3*dot7*M + 2.*dot5*dot6*M2 - 4.*a1*dot5*dot6*M2 + 2.*a12*dot5*dot6*M2 - dot9*M2 + 2.*a1*dot9*M2 + 2.*dot1*(dot10 - 2.*(dot2 + dot4)*dot6 + a2*dot7*M) - 2.*dot5*dot6*th + dot9*th - dot9*uh)/(a2*a1*(M2 - sh)*(a22*M2 - th));
	      diag[3] = -(dot12 - 2.*(dot13*dot5 + dot1*dot7))*M2/(a2*sqr(-(a12*M2) + uh));
	      diag[4] = -M*(-2.*dot10*dot3 - 2.*dot2*dot8 + 2.*dot4*dot8 + dot12*M - 2.*dot13*dot5*M - 2.*a1*dot3*dot7*M + 2.*dot5*dot6*M2 - 4.*a1*dot5*dot6*M2 + 2.*a12*dot5*dot6*M2 - dot9*M2 + 2.*a1*dot9*M2 + 2.*dot1*(dot10 - 2.*(dot2 + dot4)*dot6 - a1*dot7*M) - 2.*dot5*dot6*th + dot9*th - dot9*uh)/(a2*a1*(M2 - sh)*(a12*M2 - uh));
	      Complex flow[2];
	      flow[0] =  diag[0]-diag[2]+4.*diag[3]/3.+diag[4]/3.;
	      flow[1] =-(diag[0]-diag[4]+4.*diag[1]/3.+diag[2]/3.);
	      for(unsigned int ix=0;ix<2;++ix) {
		flows[ix]+=norm(flow[ix]);
		for(unsigned int iy=0;iy<2;++iy) {
		  me2Sum+=cMatrix[ix][iy]*flow[ix]*std::conj(flow[iy]);
		}
	      }
	    }
	  }
	}
      }
    }
    // qbar qbar -> antidiquark g
    else {
      for(unsigned int ih1=0;ih1<2;++ih1) {
	for(unsigned int ih2=0;ih2<2;++ih2) {
	  complex<Energy> dot6=v2[ih2].dimensionedWave().scalar(vbar1[ih1].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = v2[ih2].dimensionedWave().vectorCurrent(vbar1[ih1].dimensionedWave());
	  complex<Energy2> dot13 = vec1*rescaledMomenta()[3];
	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    complex<Energy> dot2 = rescaledMomenta()[0]*v3[ih3].wave();
	    complex<Energy> dot4 = rescaledMomenta()[1]*v3[ih3].wave();
	    LorentzPolarizationVectorE vec3 = v2[ih2].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(vbar1[ih1].dimensionedWave());
	    complex<Energy2> dot10 = vec3*rescaledMomenta()[3];
	    for(unsigned int ih4=0;ih4<2;++ih4) {
	      complex<Energy> dot1 = rescaledMomenta()[0]*g4[ih4].wave();
	      complex<Energy> dot3 = rescaledMomenta()[1]*g4[ih4].wave();
	      Complex dot5 = v3[ih3].wave()*g4[ih4].wave();
	      LorentzPolarizationVectorE vec2 = v2[ih2].dimensionedWave().slash(g4[ih4].wave()).vectorCurrent(vbar1[ih1].dimensionedWave());
	      LorentzPolarizationVectorE vec4 = v2[ih2].dimensionedWave().slash(g4[ih4].wave()).slash(v3[ih3].wave()).vectorCurrent(vbar1[ih1].dimensionedWave());
	      complex<Energy> dot7 = vec1*v3[ih3].wave();
	      complex<Energy2> dot8 = vec2*rescaledMomenta()[3];
	      complex<Energy> dot9 = vec2*v3[ih3].wave();
	      complex<Energy> dot11 = vec1*g4[ih4].wave();
	      complex<Energy2> dot12 = vec4*rescaledMomenta()[3];
	      Complex diag[5];
	      diag[0] = -M*(-2.*dot10*dot3 - 2.*dot2*dot8 + 2.*dot4*dot8 - 2.*a1*dot3*dot7*M + 2.*dot5*dot6*M2 - 4.*a1*dot5*dot6*M2 + 2.*a12*dot5*dot6*M2 - dot9*M2 + 2.*a1*dot9*M2 + 2.*dot1*(dot10 - 2.*(dot2 + dot4)*dot6 + a2*dot7*M) - 2.*dot5*dot6*th + dot9*th - dot9*uh)/(a2*a1*(a22*M2 - th)*(a12*M2 - uh));
	      diag[1] = ((dot12 - 2.*dot11*(dot2 + dot4) + 2.*dot3*dot7)*M2)/(a1*sqr(-(a22*M2) + th));
	      diag[2] = -M*(-2.*dot10*dot3 + 2.*(-dot2 + dot4)*dot8 + dot12*M - 2.*dot11*(dot2 + dot4)*M + 2.*dot3*dot7*M - 2.*a1*dot3*dot7*M + 2.*dot5*dot6*M2 - 4.*a1*dot5*dot6*M2 + 2.*a12*dot5*dot6*M2 - dot9*M2 + 2.*a1*dot9*M2 + 2.*dot1*(dot10 - 2.*(dot2 + dot4)*dot6 + a2*dot7*M) - 2.*dot5*dot6*th + dot9*th - dot9*uh)/(a2*a1*(M2 - sh)*(a22*M2 - th));
	      diag[3] = -(-dot12 + 2.*dot13*dot5 + 2.*dot1*dot7)*M2/(a2*sqr(-(a12*M2) + uh));
	      diag[4] = -M*(2.*dot10*dot3 + 2.*dot2*dot8 - 2.*dot4*dot8 - dot12*M + 2.*dot13*dot5*M + 2.*a1*dot3*dot7*M - 2.*dot5*dot6*M2 + 4.*a1*dot5*dot6*M2 - 2.*a12*dot5*dot6*M2 + dot9*M2 - 2.*a1*dot9*M2 + dot1*(-2.*dot10 + 4.*(dot2 + dot4)*dot6 + 2.*a1*dot7*M) + 2.*dot5*dot6*th - dot9*th + dot9*uh)/(a2*a1*(M2 - sh)*(a12*M2 - uh));
	      Complex flow[2];
	      flow[0] =  diag[0]+diag[2]+4.*diag[3]/3.+diag[4]/3.;
	      flow[1] =-(diag[0]-diag[4]+4.*diag[1]/3.-diag[2]/3.);
	      for(unsigned int ix=0;ix<2;++ix) {
		flows[ix]+=norm(flow[ix]);
		for(unsigned int iy=0;iy<2;++iy) {
		  me2Sum+=cMatrix[ix][iy]*flow[ix]*std::conj(flow[iy]);
		}
	      }
	    }
	  }
	}
      }
    }
    // prefactors, spin and colour average
    me2Sum *= R02_*sqr(Constants::pi)*pow<3,1>(standardModel()->alphaS(scale())/M)/18.;
  }
  // return the answer
  meInfo(flows);
  return me2Sum.real();
}

IBPtr MEPP2bc1Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2bc1Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2bc1Jet::doinit() {
  setDiquark(5403);
  MEPP2DiquarkJet::doinit();
  R02_ = oniumParameters()->radialWaveFunctionSquared(bc,1);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2bc1Jet,MEPP2DiquarkJet>
describeHerwigMEPP2bc1Jet("Herwig::MEPP2bc1Jet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2bc1Jet::Init() {

  static ClassDocumentation<MEPP2bc1Jet> documentation
    ("The MEPP2bc1Jet class implements the 2->2 matrix elements for vector bc production");

}

void MEPP2bc1Jet::persistentOutput(PersistentOStream & os) const {
  os << ounit(R02_,GeV*GeV2);
}

void MEPP2bc1Jet::persistentInput(PersistentIStream & is, int) {
  is >> iunit(R02_,GeV*GeV2);
}
