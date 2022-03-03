// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2cc1Jet class.
//

#include "MEPP2cc1Jet.h"
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

double MEPP2cc1Jet::me2() const {
  // return value
  Complex me2Sum(0.);
  vector<double> flows(2,0.);
  // check the right diquark
  assert(abs(mePartonData()[2]->id())==4403);
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
      g1w.reset(2*ix);
      g1.push_back(g1w);
    }
    // g q -> diquark q bar
    if(mePartonData()[1]->id()>0) {
      SpinorWaveFunction      q2w(rescaledMomenta()[1],mePartonData()[1],incoming);
      SpinorBarWaveFunction   q4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
      vector<SpinorWaveFunction> q2;
      vector<SpinorBarWaveFunction> q4;
      for(unsigned int ix=0;ix<2;++ix) {
	q2w.reset(ix);
	q2.push_back(q2w);
	q4w.reset(ix);
	q4.push_back(q4w);
      }
      // matrix element squared
      // loop over the helicities
      for(unsigned int ih2=0;ih2<2;++ih2) {
	for(unsigned int ih4=0;ih4<2;++ih4) {
	  complex<Energy> amp1 = q2[ih2].dimensionedWave().scalar(q4[ih4].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = q2[ih2].dimensionedWave().vectorCurrent(q4[ih4].dimensionedWave());
	  complex<Energy2> amp5 = vec1*rescaledMomenta()[0];
	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    complex<Energy> amp4 = vec1*v3[ih3].wave();
	    complex<Energy>  dot2 = v3[ih3].wave()*rescaledMomenta()[0];
	    complex<Energy>  dot4 = v3[ih3].wave()*rescaledMomenta()[1];
	    LorentzPolarizationVectorE vec2 = q2[ih2].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(q4[ih4].dimensionedWave());
	    complex<Energy2> amp8 = vec2*rescaledMomenta()[0];
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      complex<Energy> amp2 = vec1*g1[ih1].wave();
	      complex<Energy> amp3 = vec2*g1[ih1].wave();
	      LorentzPolarizationVectorE vec3 = q2[ih2].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(q4[ih4].dimensionedWave());
	      LorentzPolarizationVectorE vec4 = q2[ih2].dimensionedWave().slash(v3[ih3].wave()).slash(g1[ih1].wave()).vectorCurrent(q4[ih4].dimensionedWave());
	      complex<Energy2> amp6 = vec3*rescaledMomenta()[0];
	      complex<Energy2> amp7 = vec4*rescaledMomenta()[0];
	      Complex dot1 = g1[ih1].wave()*v3[ih3].wave();
	      complex<Energy>  dot3 = g1[ih1].wave()*rescaledMomenta()[1];
	      complex<Energy>  dot5 = g1[ih1].wave()*rescaledMomenta()[3];
	      // now for the diagrams
	      Complex diag[5];
	      diag[0] = -64.*(amp7 - 2*amp5*dot1 + 2*amp2*dot2 + 2*amp4*dot3)*M2/sqr(M2-4*sh);
	      diag[1] =  16.*M*(-4*amp6*(dot2 + 2*dot4) - 8*amp1*dot2*dot5 + 4*amp8*(dot3 + dot5) + 2*amp7*M - 2*amp4*dot3*M + 2*amp4*dot5*M + amp1*dot1*M2 - 2*amp3*sh + 2*amp3*uh - 4*amp1*dot1*uh)/((M2 - th)*(M2 - 4*uh));
	      diag[2] =  16.*M*(4*amp6*(dot2 + 2*dot4) + 8*amp1*dot2*dot5 - 4*amp8*(dot3 + dot5) - 2*amp7*M + 4*amp5*dot1*M - 4*amp2*dot2*M - 2*amp4*dot3*M + 2*amp4*dot5*M - amp1*dot1*M2 + 2*amp3*sh - 2*amp3*uh + 4*amp1*dot1*uh)/((M2 - 4*sh)*(M2 - th));
	      diag[3] =  64.*(amp7 + 2*amp4*dot5)*M2/sqr(M2 - 4*uh);
	      diag[4] = -64.*M*(4*amp6*(dot2 + 2*dot4) + 8*amp1*dot2*dot5 - 4*amp8*(dot3 + dot5) + 2*amp4*dot3*M + 2*amp4*dot5*M - amp1*dot1*M2 + 2*amp3*sh - 2*amp3*uh + 4*amp1*dot1*uh)/((M2 - 4*sh)*(M2 - 4*uh));
	      Complex flow[2];
	      flow[0] = 4.*diag[0]/3.+diag[1]+diag[4]+diag[2]/3.;
	      flow[1] = 4.*diag[3]/3.+diag[2]-diag[4]+diag[1]/3.;
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
      // prefactors, spin and colour average
      me2Sum *= R02_*sqr(Constants::pi)*pow<3,1>(standardModel()->alphaS(scale())/M)/96.;
    }
    else
      assert(false);
  }
  else
    assert(false);
  // return the answer
  meInfo(flows);
  return me2Sum.real();
}

IBPtr MEPP2cc1Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2cc1Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2cc1Jet::doinit() {
  MEPP2DiquarkJet::doinit();
  R02_ = oniumParameters()->radialWaveFunctionSquared(cc,1);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2cc1Jet,MEPP2DiquarkJet>
describeHerwigMEPP2cc1Jet("Herwig::MEPP2cc1Jet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2cc1Jet::Init() {

  static ClassDocumentation<MEPP2cc1Jet> documentation
    ("The MEPP2cc1Jet class implements the 2->2 matrix elements for vector cc production");

}

void MEPP2cc1Jet::persistentOutput(PersistentOStream & os) const {
  os << ounit(R02_,GeV*GeV2);
}

void MEPP2cc1Jet::persistentInput(PersistentIStream & is, int) {
  is >> iunit(R02_,GeV*GeV2);
}
