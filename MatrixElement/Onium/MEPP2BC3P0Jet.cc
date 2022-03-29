// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2BC3P0Jet class.
//

#include "MEPP2BC3P0Jet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include <numeric>

using namespace Herwig;

double MEPP2BC3P0Jet::me2() const {
  // invariants
  Energy M  = rescaledMomenta()[2].mass();
  Energy2 M2 = sqr(M);
  Energy2 sh((rescaledMomenta()[0]+rescaledMomenta()[1]).m2());
  Energy2 th((rescaledMomenta()[0]-rescaledMomenta()[2]).m2());
  Energy2 uh((rescaledMomenta()[0]-rescaledMomenta()[3]).m2());
  double a1,a2;
  double a12,a22;
  // storage of the result
  vector<Complex> diag(5);
  DVector save(2,0.);
  double meSum(0.);
  // gluon initiated processes
  if(mePartonData()[0]->id()==ParticleID::g) {
    // gluon wavefunction
    VectorWaveFunction g1w(rescaledMomenta()[0],mePartonData()[0],incoming);
    vector<VectorWaveFunction> g1;
    for(unsigned int ix=0;ix<2;++ix) {
      //g1w.reset(10);
      g1w.reset(2*ix);
      g1.push_back(g1w);
    }
    // matrix element
    ProductionMatrixElement me(PDT::Spin1,PDT::Spin1Half,PDT::Spin0,PDT::Spin1Half);
    // g c -> B_c b
    if(mePartonData()[1]->id()>0) {
      a1 = rescaledMomenta()[1].mass()/M;
      a2 = rescaledMomenta()[3].mass()/M;
      a12=sqr(a1);
      a22=sqr(a2);
      SpinorWaveFunction      q2w(rescaledMomenta()[1],mePartonData()[1],incoming);
      SpinorBarWaveFunction   q4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
      vector<SpinorWaveFunction> u2;
      vector<SpinorBarWaveFunction> ubar4;
      for(unsigned int ix=0;ix<2;++ix) {
        q2w.reset(ix);
        u2.push_back(q2w);
        q4w.reset(ix);
        ubar4.push_back(q4w);
      }
      for(unsigned int ih2=0;ih2<2;++ih2) {
	for(unsigned int ih4=0;ih4<2;++ih4) {
	  complex<Energy> dot3=u2[ih2].dimensionedWave().scalar(ubar4[ih4].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = u2[ih2].dimensionedWave().vectorCurrent(ubar4[ih4].dimensionedWave());
	  complex<Energy2> dot6 = vec1*rescaledMomenta()[0];
	  for(unsigned int ih1=0;ih1<2;++ih1) {
	    complex<Energy> dot1 = rescaledMomenta()[1]*g1[ih1].wave();
	    complex<Energy> dot2 = rescaledMomenta()[3]*g1[ih1].wave();
	    LorentzPolarizationVectorE vec2 = u2[ih2].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(ubar4[ih4].dimensionedWave());
	    complex<Energy> dot4 = vec1*g1[ih1].wave();
	    complex<Energy2> dot5 = vec2*rescaledMomenta()[0];
	    // diagrams
	    diag[0]=(32.*a12*(2.*dot1*dot3+dot5)*pow<4,1>(M))/pow<3,1>(-(a12*M2)+sh)+
	      (4.*(((-3.+2.*a1)*(1.+4.*a1)*(2.*dot1*dot3+dot5))/a2+4.*a1*dot4*M)*M2)/sqr(-(a12*M2)+sh)-(4.*(3.-4.*a1)*dot4*M)/(a2*(-(a12*M2)+sh));
	    diag[1]=(-2.*(1.-2.*a1)*dot4*M)/(a1*a2*(-M2+th))
	      +((-8.*(2.*dot1*dot3+dot5)*M2)/a2-(16.*pow<3,1>(M)*(-(a2*dot1*dot6)+2.*a12*dot1*dot3*M+a1*dot5*M-a1*dot2*(dot6+2.*a1*dot3*M)))/(-M2+th))/sqr(-(a12*M2)+sh)+
	      ((-2.*(1.-2.*a1)*dot4*M)/(a1*a2)+(16.*(-dot1+dot2)*dot6*pow<3,1>(M))/sqr(-M2+th)-
	       (2.*M*(2.*(a1*(-3.+4.*a1)*dot1-dot2+(5.-4.*a1)*a1*dot2)*dot6+
		      (2.*a1*(1.+2.*a1*(-5.+4.*a1))*(dot1-dot2)*dot3-sqr(1.-2.*a1)*dot5)*M+8.*a1*a2*dot4*M2))/(a1*a2*(-M2+th)))/(-(a12*M2)+sh);
	    diag[2]=(-2.*(-1.+2.*a1)*dot4*M)/(a1*a2*(-M2+th))
	      +((-8.*(2.*dot2*dot3+dot5)*M2)/a1+
		(16.*pow<3,1>(M)*((-(a2*dot1)-a1*dot2)*dot6-a2*(-2.*a2*(dot1-dot2)*dot3+dot5)*M))/(-M2+th))/sqr(-(a22*M2)+uh)+
	      ((-2.*(-1.+2.*a1)*dot4*M)/(a1*a2)+(16.*(-dot1+dot2)*dot6*pow<3,1>(M))/sqr(-M2+th)-
	       (2.*M*(2.*(a1*(-3.+4.*a1)*dot1-dot2+(5.-4.*a1)*a1*dot2)*dot6+(-dot5-2.*a2*((-1.-6.*a1+8.*a12)*(dot1-dot2)*dot3-2.*a1*dot5))*M+8.*a1*a2*dot4*M2))/(a1*a2*(-M2+th)))/(-(a22*M2)+uh);
	    diag[3]=(32.*a22*(2.*dot2*dot3+dot5)*pow<4,1>(M))/pow<3,1>(-(a22*M2)+uh)+(4.*((1.+2.*a1)*(-5.+4.*a1)*(2.*dot2*dot3+dot5)+4.*a1*a2*dot4*M)*M2)/(a1*sqr(-(a22*M2)+uh))
	      +(4.*(1.-4.*a1)*dot4*M)/(a1*(-(a22*M2)+uh));
	    diag[4]=(8.*(2.*dot2*dot3+dot5)*M2)/(a1*sqr(-(a22*M2)+uh))
	      -(2.*(dot4-2.*a1*dot4)*M)/(a1*a2*(-(a22*M2)+uh))
	      +((8.*(2.*dot1*dot3+dot5)*M2)/a2-(16.*(-(a2*dot1)-a1*dot2)*pow<3,1>(M)*(dot6+2.*a1*dot3*M))/(-(a22*M2)+uh))/sqr(-(a12*M2)+sh)+
	      ((-2.*(-1.+2.*a1)*dot4*M)/(a1*a2)+(16.*(-(a2*dot1)-a1*dot2)*pow<3,1>(M)*(dot6-2.*a2*dot3*M))/sqr(-(a22*M2)+uh)-
	       (4.*M*((a1*(-3.+4.*a1)*dot1-dot2+(5.-4.*a1)*a1*dot2)*dot6+
		      (-((-1.-6.*a1+8.*a12)*a2*dot1*dot3)+a1*(-1.+2.*(5.-4.*a1)*a1)*dot2*dot3+4.*a1*a2*dot5)*M+2.*a1*a2*dot4*M2))/(a1*a2*(-(a22*M2)+uh)))/(-(a12*M2)+sh);
 	    // diagram weights
 	    save[0]+=norm(diag[3]);
 	    save[1]+=norm(diag[1]);
	    Complex aSum = (-0.5*diag[0] + diag[1] + diag[2] - 0.5*diag[3])/3. + 1.5*(diag[0] + diag[3] + 2.*diag[4]);
	    meSum+=norm(aSum);
	    me(2*ih1,ih2,0,ih4)=aSum;
	  }
	}
      }
    }
    // g cbar -> B_c- bbar
    else {
      a1 = rescaledMomenta()[3].mass()/M;
      a2 = rescaledMomenta()[1].mass()/M;
      a12=sqr(a1);
      a22=sqr(a2);
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
      for(unsigned int ih2=0;ih2<2;++ih2) {
        for(unsigned int ih4=0;ih4<2;++ih4) {
	  complex<Energy> dot3=v4[ih4].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	  complex<Energy2> dot6 = vec1*rescaledMomenta()[0];
	  for(unsigned int ih1=0;ih1<2;++ih1) {
	    complex<Energy> dot1 = rescaledMomenta()[1]*g1[ih1].wave();
	    complex<Energy> dot2 = rescaledMomenta()[3]*g1[ih1].wave();
	    LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	    complex<Energy> dot4 = vec1*g1[ih1].wave();
	    complex<Energy2> dot5 = vec2*rescaledMomenta()[0];
	    // diagrams
	    diag[0]=(32.*a22*(2.*dot1*dot3-dot5)*pow<4,1>(M))/pow<3,1>(-(a22*M2)+sh)+(4.*((1.+2.*a1)*(-5.+4.*a1)*(2.*dot1*dot3-dot5)-4.*a1*a2*dot4*M)*M2)/(a1*sqr(-(a22*M2)+sh))
	      +(4.*(-1.+4.*a1)*dot4*M)/(a1*(-(a22*M2)+sh));
	    diag[1]=(-2.*(-1.+2.*a1)*dot4*M)/(a1*a2*(-M2+th))
	      +((8.*(-2.*dot2*dot3+dot5)*M2)/a2+(16.*pow<3,1>(M)*(dot2*dot6+2.*a12*(dot1-dot2)*dot3*M+a1*(dot1*dot6-dot2*dot6+dot5*M)))/(-M2+th))/sqr(-(a12*M2)+uh)+
	      ((-2.*(-1.+2.*a1)*dot4*M)/(a1*a2)
	       +(16.*(dot1-dot2)*dot6*pow<3,1>(M))/sqr(-M2+th)-
	       (2.*M*(-2.*(-((-1.+4.*a1)*a2*dot1)+(3.-4.*a1)*a1*dot2)*dot6+
		      (-2.*a1*(1.+2.*a1*(-5.+4.*a1))*(dot1-dot2)*dot3+dot5-4.*a1*a2*dot5)*M-8.*a1*a2*dot4*M2))/(a1*a2*(-M2+th)))/(-(a12*M2)+uh);
	    diag[2]=(-2.*(1.-2.*a1)*dot4*M)/(a1*a2*(-M2+th))
	      +((8.*(-2.*dot1*dot3+dot5)*M2)/a1+(16.*pow<3,1>(M)*(-((a1*(dot1-dot2)+dot2)*dot6)+a2*(-2.*a2*(dot1-dot2)*dot3+dot5)*M))/(-M2+th))/sqr(-(a22*M2)+sh)+
	      ((-2.*(dot4-2.*a1*dot4)*M)/(a1*a2)
	       +(16.*(dot1-dot2)*dot6*pow<3,1>(M))/sqr(-M2+th)-
	       (2.*M*(-2.*(-((-1.+4.*a1)*a2*dot1)+(3.-4.*a1)*a1*dot2)*dot6+
		      (dot5-2.*a2*(-((-1.-6.*a1+8.*a12)*(dot1-dot2)*dot3)+2.*a1*dot5))*M-8.*a1*a2*dot4*M2))/(a1*a2*(-M2+th)))/(-(a22*M2)+sh);
	    diag[3]=(32.*a12*(2.*dot2*dot3-dot5)*pow<4,1>(M))/pow<3,1>(-(a12*M2)+uh)+
	      (4.*(((-3.+2.*a1)*(1.+4.*a1)*(2.*dot2*dot3-dot5))/a2-4.*a1*dot4*M)*M2)/sqr(-(a12*M2)+uh)-(4.*(-3.+4.*a1)*dot4*M)/(a2*(-(a12*M2)+uh));
	    diag[4]=(-8.*(-2.*dot2*dot3+dot5)*M2)/(a2*sqr(-(a12*M2)+uh))
	      -(2.*(1.-2.*a1)*dot4*M)/(a1*a2*(-(a12*M2)+uh))
	      +((8.*(2.*dot1*dot3-dot5)*M2)/a1-(16.*(a1*(dot1-dot2)+dot2)*pow<3,1>(M)*(dot6-2.*a2*dot3*M))/(-(a12*M2)+uh))/sqr(-(a22*M2)+sh)+
	      ((-2.*(-1.+2.*a1)*dot4*M)/(a1*a2)
	       +(16.*(a1*(dot1-dot2)+dot2)*pow<3,1>(M)*(dot6+2.*a1*dot3*M))/sqr(-(a12*M2)+uh)-
	       (4.*M*(a1*(-3.+4.*a1)*dot2*dot6-a2*((-1.-6.*a1+8.*a12)*dot2*dot3+4.*a1*dot5)*M-
		      dot1*(-((-1.+4.*a1)*a2*dot6)+a1*(1.+2.*a1*(-5.+4.*a1))*dot3*M)-2.*a1*a2*dot4*M2))/(a1*a2*(-(a12*M2)+uh)))/(-(a22*M2)+sh);
	    // diagram weights
	    save[0]+=norm(diag[3]);
	    save[1]+=norm(diag[1]);
	    Complex aSum = (-0.5*diag[0] + diag[1] + diag[2] - 0.5*diag[3])/3. + 1.5*(diag[0] + diag[3] + 2.*diag[4]);
	    meSum+=norm(aSum);
	    me(2*ih1,ih2,0,ih4)=aSum;
	  }
	}
      }
    }
    // set matrix element
    setME(me);
    // final factors
    meSum *= 1./1296.;
  }
  // q qbar process
  else {
    a1 = rescaledMomenta()[0].mass()/M;
    a2 = rescaledMomenta()[1].mass()/M;
    a12=sqr(a1);
    a22=sqr(a2);
    // gluon wavefunction
    VectorWaveFunction g4w(rescaledMomenta()[3],mePartonData()[3],incoming);
    vector<VectorWaveFunction> g4;
    for(unsigned int ix=0;ix<2;++ix) {
      //g4w.reset(10);
      g4w.reset(2*ix);
      g4.push_back(g4w);
    }
    SpinorWaveFunction    q1w(rescaledMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction q2w(rescaledMomenta()[1],mePartonData()[1],incoming);
    vector<SpinorWaveFunction> u1;
    vector<SpinorBarWaveFunction> vbar2;
    for(unsigned int ix=0;ix<2;++ix) {
      q1w.reset(ix);
      u1.push_back(q1w);
      q2w.reset(ix);
      vbar2.push_back(q2w);
    }
    // matrix element
    ProductionMatrixElement me(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0,PDT::Spin1);
    for(unsigned int ih1=0;ih1<2;++ih1) { 
      for(unsigned int ih2=0;ih2<2;++ih2) {
     	for(unsigned int ih4=0;ih4<2;++ih4) {
	  complex<Energy> dot1 = rescaledMomenta()[0]*g4[ih4].wave();
	  complex<Energy> dot2 = rescaledMomenta()[1]*g4[ih4].wave();
	  complex<Energy> dot3=u1[ih1].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = u1[ih1].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	  LorentzPolarizationVectorE vec2 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	  complex<Energy2> dot4 = vec1*rescaledMomenta()[3];
	  complex<Energy>  dot5 = vec1*g4[ih4].wave();
	  complex<Energy2> dot6 = vec2*rescaledMomenta()[3];
	  // diagrams
	  diag[0]=(-8.*(2.*dot2*dot3+dot6)*M2)/(a1*sqr(-(a22*M2)+th))
	    -(2.*(dot5-2.*a1*dot5)*M)/(a1*a2*(-(a22*M2)+th))
	    +((-8.*(-2.*dot1*dot3+dot6)*M2)/a2-(16.*(-(a2*dot1)+a1*dot2)*pow<3,1>(M)*(-dot4+2.*a1*dot3*M))/(-(a22*M2)+th))/sqr(-(a12*M2)+uh)+
	    ((-2.*(-1.+2.*a1)*dot5*M)/(a1*a2)
	     -(16.*(-(a2*dot1)+a1*dot2)*pow<3,1>(M)*(dot4+2.*a2*dot3*M))/sqr(-(a22*M2)+th)-
	     (4.*M*(-((dot2+4.*a12*(dot1+dot2)-a1*(3.*dot1+5.*dot2))*dot4)+
		    (-((-1.-6.*a1+8.*a12)*a2*dot1*dot3)+a1*(1.+2.*a1*(-5.+4.*a1))*dot2*dot3-4.*a1*a2*dot6)*M+2.*a1*a2*dot5*M2))/(a1*a2*(-(a22*M2)+th)))/(-(a12*M2)+uh);
	  diag[1]=(-32.*a22*(2.*dot2*dot3+dot6)*pow<4,1>(M))/pow<3,1>(-(a22*M2)+th)+
	    (4.*(-((1.+2.*a1)*(-5.+4.*a1)*(2.*dot2*dot3+dot6))+4.*a1*a2*dot5*M)*M2)/(a1*sqr(-(a22*M2)+th))
	    +(4.*(1.-4.*a1)*dot5*M)/(a1*(-(a22*M2)+th));
	  diag[2]=(8.*(2.*dot2*dot3+dot6)*M2)/(a1*sqr(-(a22*M2)+th))
	    -(2.*(-1.+2.*a1)*dot5*M)/(a1*a2*(-(a22*M2)+th))
	    +(16.*(dot1+dot2)*dot4*pow<3,1>(M))/(sqr(-M2+sh)*(-(a22*M2)+th))
	    +((-2.*(-1.+2.*a1)*dot5*M)/(a1*a2)
	      +(16.*pow<3,1>(M)*((dot1-a1*dot1-a1*dot2)*dot4-a2*(-2.*a2*(dot1+dot2)*dot3-dot6)*M))/sqr(-(a22*M2)+th)-
	      (2.*M*(-2.*(dot2+4.*a12*(dot1+dot2)-a1*(3.*dot1+5.*dot2))*dot4+
		     (dot6-2.*a2*((-1.-6.*a1+8.*a12)*(dot1+dot2)*dot3+2.*a1*dot6))*M+8.*a1*a2*dot5*M2))/(a1*a2*(-(a22*M2)+th)))/(-M2+sh);
	  diag[3]=(32.*a12*(2.*dot1*dot3-dot6)*pow<4,1>(M))/pow<3,1>(-(a12*M2)+uh)+
	    (4.*(((-3.+2.*a1)*(1.+4.*a1)*(2.*dot1*dot3-dot6))/a2+4.*a1*dot5*M)*M2)/sqr(-(a12*M2)+uh)-(4.*(3.-4.*a1)*dot5*M)/(a2*(-(a12*M2)+uh));
	  diag[4]=(-8.*(2.*dot1*dot3-dot6)*M2)/(a2*sqr(-(a12*M2)+uh))
	    -(2.*(dot5-2.*a1*dot5)*M)/(a1*a2*(-(a12*M2)+uh))
	    +(16.*(dot1+dot2)*dot4*pow<3,1>(M))/(sqr(-M2+sh)*(-(a12*M2)+uh))
	    +((-2.*(dot5-2.*a1*dot5)*M)/(a1*a2)
	      +(16.*pow<3,1>(M)*((-(a2*dot1)+a1*dot2)*dot4+a1*(-2.*a1*(dot1+dot2)*dot3+dot6)*M))/sqr(-(a12*M2)+uh)-
	      (2.*M*(-2.*(dot2+4.*a12*(dot1+dot2)-a1*(3.*dot1+5.*dot2))*dot4+
		     (2.*a1*(1.+2.*a1*(-5.+4.*a1))*(dot1+dot2)*dot3+dot6-4.*a1*a2*dot6)*M+8.*a1*a2*dot5*M2))/(a1*a2*(-(a12*M2)+uh)))/(-M2+sh);
    	  // diagram weights
    	  save[0]+=norm(diag[2]);
    	  save[1]+=norm(diag[4]);
    	  Complex aSum =1.5*(2.*diag[0] + diag[1] + diag[3]) + (-0.5*diag[1] + diag[2] - 0.5*diag[3] + diag[4])/3.;
    	  meSum+=norm(aSum);
    	  me(ih1,ih2,0,2*ih4)=aSum;
    	}
      }
    }
    // set matrix element
    setME(me);
    meSum *= 1./486.;
  }
  // save the diagram weights
  meInfo(save);
  // final factors
  return meSum*O1_*pow<3,1>(Constants::pi*standardModel()->alphaS(scale())/M)/M2/sqr(a1*a2);
}

IBPtr MEPP2BC3P0Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2BC3P0Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2BC3P0Jet::doinit() {
  setBcState(10541);
  MEPP2BCJetBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<1>(bcbar,principleQuantumNumber(),1,0);

}

void MEPP2BC3P0Jet::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2*GeV2);
}

void MEPP2BC3P0Jet::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2*GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2BC3P0Jet,MEPP2BCJetBase>
describeHerwigMEPP2BC3P0Jet("Herwig::MEPP2BC3P0Jet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2BC3P0Jet::Init() {

  static ClassDocumentation<MEPP2BC3P0Jet> documentation
    ("The MEPP2BC3P0Jet class implements the matrix element for g q -> B_c(3P0) q and "
     "q bar -> Bc(3P0) g and charged conjugate");

}
