// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2BC3D1Jet class.
//

#include "MEPP2BC3D1Jet.h"
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

double MEPP2BC3D1Jet::me2() const {
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
  // B_C* wavefunction
  VectorWaveFunction Bcw(rescaledMomenta()[2],mePartonData()[2],outgoing);
  vector<VectorWaveFunction> v3;
  for(unsigned int ix=0;ix<3;++ix) {
    Bcw.reset(ix);
    v3.push_back(Bcw);
  }
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
    ProductionMatrixElement me(PDT::Spin1,PDT::Spin1Half,PDT::Spin1,PDT::Spin1Half);
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
	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    for(unsigned int ih1=0;ih1<2;++ih1)
	      {auto dot1 = rescaledMomenta()[0]*v3[ih3].wave();
		auto dot2 = rescaledMomenta()[1]*g1[ih1].wave();
		auto dot3 = rescaledMomenta()[1]*v3[ih3].wave();
		auto dot4 = rescaledMomenta()[3]*g1[ih1].wave();
		auto dot5 = g1[ih1].wave()*v3[ih3].wave();
		complex<Energy> dot6=u2[ih2].dimensionedWave().scalar(ubar4[ih4].dimensionedWave());
		auto vec1 = u2[ih2].dimensionedWave().vectorCurrent(ubar4[ih4].dimensionedWave());
		auto vec2 = u2[ih2].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(ubar4[ih4].dimensionedWave());
		auto vec3 = u2[ih2].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(ubar4[ih4].dimensionedWave());
		auto vec4 = u2[ih2].dimensionedWave().slash(v3[ih3].wave()).slash(g1[ih1].wave()).vectorCurrent(ubar4[ih4].dimensionedWave());
		auto dot7 = vec1*rescaledMomenta()[0];
		auto dot8 = vec1*g1[ih1].wave();
		auto dot9 = vec1*v3[ih3].wave();
		auto dot10 = vec2*g1[ih1].wave();
		auto dot11 = vec3*rescaledMomenta()[0];
		auto dot12 = vec4*rescaledMomenta()[0];
		auto dot13 = vec2*rescaledMomenta()[0];
	      // diagrams
		diag[0]=(256.*a12*(dot1+dot3)*(dot11+2.*dot2*dot6)*pow<5,1>(M))/(a2*pow<4,1>(-(a12*M2)+sh))
		  +(16.*pow<3,1>(M)*((-15.+4.*a1*(-7.+4.*a1))*(dot1+dot3)*(dot11+2.*dot2*dot6)-8.*a1*a2*(dot12-2.*dot5*dot7-(dot1+3.*dot3)*dot8+2.*dot2*dot9)*M))/(3.*a22*pow<3,1>(-(a12*M2)+sh))
		  +(8.*((15.-6.*a1)*dot12+2.*(-15.*dot5*dot7+6.*a1*dot5*dot7+2.*a1*dot1*dot8-15.*dot3*dot8+8.*a1*dot3*dot8+15.*dot2*dot9-6.*a1*dot2*dot9-5.*a2*(dot10-2.*dot5*dot6)*M))*M2)/(3.*a22*sqr(-(a12*M2)+sh));
		diag[1]=(-8.*((3.-2.*a1)*dot1*dot8+2.*a1*(5.*dot2*dot9-5.*dot4*dot9+4.*dot10*M-4.*dot5*dot6*M))*M2)/(3.*a12*a2*sqr(-M2+th))
		  +(8.*(-1.-2.*a1*a2)*(dot10-dot5*dot6)*M)/(3.*a12*a22*(-M2+th))
		  +((-64.*(dot1+dot3)*(dot11+2.*dot2*dot6)*pow<3,1>(M))/a22-
		    (128.*(dot1+dot3)*pow<4,1>(M)*((-(a2*dot2)-a1*dot4)*dot7+a1*(dot11+2.*a1*(dot2-dot4)*dot6)*M))/(a2*(-M2+th)))/pow<3,1>(-(a12*M2)+sh)+
		  ((4.*(-1.-2.*a1*a2)*(dot10-2.*dot5*dot6)*M)/(3.*a12*a22)
		   +(128.*dot1*(-dot2+dot4)*dot7*pow<4,1>(M))/(a1*pow<3,1>(-M2+th))
		   -(8.*M2*(2.*(a1*(-9.+8.*a1)*dot1*dot2-10.*a1*dot2*dot3+(-3.+8.*a1)*a2*dot1*dot4+10.*a1*dot3*dot4)*dot7+
			    ((-3.-4.*(-4.+a1)*a1)*dot1*dot11+2.*a1*(8.*dot11*dot3+dot13*(-9.*dot2+10.*a1*dot2+dot4-10.*a1*dot4))+
			     2.*a1*dot1*(dot2+2.*a1*(-13.+8.*a1)*dot2+(7.+2.*(13.-8.*a1)*a1)*dot4)*dot6)*M-
			    4.*a1*(dot12+(-5.+3.*a1)*dot5*dot7+2.*(-2.+3.*a1)*dot1*dot8+2.*a1*(dot2-dot4)*dot9)*M2))/(3.*a12*a2*sqr(-M2+th))
		   +(4.*M*(2.*(-1.-2.*a1*a2)*(2.*dot11*dot3-dot13*(dot2+dot4)+dot1*(dot11+2.*dot4*dot6))+
			   (dot12-8.*a1*a2*dot12+2.*(-1.-a1*(-7.+8.*a1)*a2)*dot5*dot7-
			    2.*(2.*(1.+a1+a12)*dot1+a1*(-11.+8.*a1)*dot3)*dot8+2.*a1*((-10.+a1*(-5.+6.*a1))*dot2+3.*(5.-2.*a1)*a1*dot4)*dot9)*M+
			   2.*a1*(2.*(-2.+a1*(-10.+7.*a1))*dot10+(19.+a1*(21.+2.*a1*(-23.+8.*a1)))*dot5*dot6)*M2))/(3.*a12*a22*(-M2+th)))/(-(a12*M2)+sh)
		  +((8.*(-5.*dot12+((-9.+2.*a1)*dot1+dot3+2.*a1*dot3)*dot8+2.*a1*((-6.+a1)*dot10+(7.-2.*a1)*dot5*dot6)*M)*M2)/(3.*a1*a22)
		    -(64.*pow<4,1>(M)*(a1*dot3*(dot2-dot4)*dot7-dot1*(dot2+(-3.+a1)*a1*dot2-(-2.+a1)*a1*dot4)*dot7+a1*a2*dot1*(dot11+2.*a1*(dot2-dot4)*dot6)*M))/(a1*a2*sqr(-M2+th))
		    +(8.*M2*(-2.*(-(a1*(dot1+dot3)*(15.*dot2-17.*dot4))+8.*a12*(dot1+dot3)*(dot2-dot4)+(-9.*dot1+dot3)*dot4)*dot7+
			     ((-3.-22.*a1+8.*a12)*dot1*dot11+(-11.+6.*a1*(-3.+2.*a1))*dot11*dot3+
			      2.*dot13*((5.+(2.-7.*a1)*a1)*dot2+a1*(7.+3.*a1)*dot4)-2.*a1*(9.+4.*a1*(-7.+4.*a1))*dot3*(dot2-dot4)*dot6-
			      2.*dot1*((12.+a1*(7.+2.*a1*(-19.+8.*a1)))*dot2+a1*(-5.+2.*(17.-8.*a1)*a1)*dot4)*dot6)*M-
			     24.*a12*a22*dot5*dot6*pow<3,1>(M)+4.*a1*a2*(3.*dot12+3.*(-3.+a1)*dot5*dot7-3.*(dot1+3.*dot3)*dot8+
									 2.*(dot2+2.*a1*dot2-2.*a1*dot4)*dot9)*M2))/(3.*a1*a22*(-M2+th)))/sqr(-(a12*M2)+sh);
		diag[2]=(8.*(dot1*(dot8+2.*a1*dot8)+2.*a2*(5.*dot2*dot9-5.*dot4*dot9+4.*dot10*M-4.*dot5*dot6*M))*M2)/(3.*a1*a22*sqr(-M2+th))
		  +(8.*(1.+2.*a1*a2)*(dot10-dot5*dot6)*M)/(3.*a12*a22*(-M2+th))
		  +((-64.*dot3*(dot11+2.*dot4*dot6)*pow<3,1>(M))/a12+(128.*dot3*pow<4,1>(M)*((-(a2*dot2)-a1*dot4)*dot7-a2*(dot11-2.*a2*(dot2-dot4)*dot6)*M))/(a1*(-M2+th)))/pow<3,1>(-(a22*M2)+uh)+
		  ((-8.*(5.*dot12+(-3.+2.*a1)*dot3*dot8-2.*(5.+a1)*a2*dot10*M-10.*dot5*(dot7-a2*dot6*M))*M2)/(3.*a12*a2)
		   +(64.*pow<4,1>(M)*((dot3*(-dot2+dot4)*dot7)/a1+dot1*(dot2*dot7+(a1*dot4*dot7)/a2+dot11*M-2.*a2*(dot2-dot4)*dot6*M)))/sqr(-M2+th)-
		   (8.*M2*(-20.*dot1*dot2*dot7+2.*dot3*((-10.+a1+8.*a12)*dot2+(7.+a1-8.*a12)*dot4)*dot7+
			   ((17.+6.*(1.-2.*a1)*a1)*dot11*dot3-2.*dot13*(-((-10.+3.*a1)*a2*dot2)+(12.-7.*a1)*a1*dot4)+2.*(3.+a1+4.*a12*(-5.+4.*a1))*dot3*(dot2-dot4)*dot6)*M-
			   4.*dot1*((-3.+a1)*a1*dot11-5.*a2*dot2*dot6+a1*(-1.+2.*a1)*dot4*dot6)*M+
			   24.*a12*a22*dot5*dot6*pow<3,1>(M)-4.*a1*a2*(3.*dot12+3.*a1*dot5*dot7-9.*dot3*dot8-4.*dot2*dot9+4.*a1*(dot2-dot4)*dot9+
								       6.*dot4*dot9)*M2))/(3.*a12*a2*(-M2+th)))/sqr(-(a22*M2)+uh)+
		  ((4.*(1.+2.*a1*a2)*dot10*M)/(3.*a12*a22)-(128.*dot1*(-dot2+dot4)*dot7*pow<4,1>(M))/(a2*pow<3,1>(-M2+th))
		   +(8.*M2*(2.*((10.+a1*(-15.+8.*a1))*dot1*dot2+(-9.+8.*a1)*a2*dot1*dot4-10.*a2*dot3*(-dot2+dot4))*dot7+
			    ((-7.-4.*(-2.+a1)*a1)*dot1*dot11-2.*a2*(8.*dot11*dot3+dot13*(-9.*dot2+10.*a1*dot2+dot4-10.*a1*dot4))-
			     2.*a2*dot1*(dot2+2.*a1*(-13.+8.*a1)*dot2+(7.+2.*(13.-8.*a1)*a1)*dot4)*dot6)*M+
			    4.*a2*(dot12+3.*a1*dot5*dot7+6.*a1*dot1*dot8-2.*a2*(dot2-dot4)*dot9)*M2))/(3.*a1*a22*sqr(-M2+th))
		   +(4.*M*(-4.*dot11*dot3-2.*(-1.-2.*a1*a2)*dot13*(dot2+dot4)+dot12*M+
			   2.*dot1*((-1.-2.*a1*a2)*(dot11+2.*dot4*dot6)+(-2.+a1)*(-5.+2.*a1)*dot8*M)-
			   2.*a2*(4.*a1*dot11*dot3+
				  (4.*a1*dot12+a1*(-9.+8.*a1)*dot5*dot7-8.*a1*dot3*dot8+9.*dot4*dot9+a1*((3.+6.*a1)*dot2+(7.-6.*a1)*dot4)*dot9-3.*(dot3*dot8+3.*dot2*dot9))*M+
				  (2.*(-5.+a1*(-4.+7.*a1))*dot10+(10.+a1*(-7.+2.*a1*(-15.+8.*a1)))*dot5*dot6)*M2)))/(3.*a12*a22*(-M2+th)))/(-(a22*M2)+uh);
		diag[3]=(256.*a22*dot3*(dot11+2.*dot4*dot6)*pow<5,1>(M))/(a1*pow<4,1>(-(a22*M2)+uh))
		  +(16.*pow<3,1>(M)*((-27.+4.*a1*(-1.+4.*a1))*dot3*(dot11+2.*dot4*dot6)-8.*a1*a2*(dot12-3.*dot3*dot8+2.*dot4*dot9)*M))/(3.*a12*pow<3,1>(-(a22*M2)+uh))
		  +(8.*((9.+6.*a1)*dot12-2.*(7.+8.*a1)*dot3*dot8+6.*(3.+2.*a1)*dot4*dot9+10.*a1*dot10*M)*M2)/(3.*a12*sqr(-(a22*M2)+uh));
		diag[4]=(64.*dot3*(dot11+2.*dot4*dot6)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+uh))
		  +(8.*(-5.*dot12+(7.+2.*a1)*dot3*dot8-10.*dot4*dot9-2.*a1*a2*dot10*M)*M2)/(3.*a12*a2*sqr(-(a22*M2)+uh))
		  +(4.*(-1.-2.*a1*a2)*dot10*M)/(3.*a12*a22*(-(a22*M2)+uh))
		  +((64.*(dot1+dot3)*(dot11+2.*dot2*dot6)*pow<3,1>(M))/a22-(128.*(dot1+dot3)*(-(a2*dot2)-a1*dot4)*pow<4,1>(M)*(dot7+2.*a1*dot6*M))/(a2*(-(a22*M2)+uh)))/pow<3,1>(-(a12*M2)+sh)+
		  ((-4.*(-1.-2.*a1*a2)*(dot10-2.*dot5*dot6)*M)/(3.*a12*a22)
		   +(128.*dot3*(-(a2*dot2)-a1*dot4)*pow<4,1>(M)*(dot7-2.*a2*dot6*M))/(a1*pow<3,1>(-(a22*M2)+uh))
		   -(16.*M2*(dot3*(-3.*dot4+a1*(dot2+8.*a1*dot2+dot4-8.*a1*dot4))*dot7+
			     (-(a1*(2.*(2.+a1)*dot1*dot11+2.*(-1.+4.*a1)*dot11*dot3+dot13*(-3.*a2*dot2+2.*dot4-7.*a1*dot4)))+
			      ((3.+a1+4.*a12*(-5.+4.*a1))*dot2*dot3-a1*(4.*(2.+a1)*dot1+3.*dot3-16.*a1*a2*dot3)*dot4)*dot6)*M+
			     12.*a12*a22*dot5*dot6*pow<3,1>(M)-2.*a1*a2*(3.*a1*dot5*dot7-3.*dot3*dot8-4.*dot2*dot9+4.*a1*(dot2-dot4)*dot9)*M2))/(3.*a12*a2*sqr(-(a22*M2)+uh))
		   +(8.*M*(-((-1.-2.*a1*a2)*(-2.*dot11*dot3+dot13*(dot2+dot4)))+(-1.-2.*a1*a2)*dot1*(dot11+2.*dot4*dot6)+a1*(1.+2.*a1)*dot1*dot8*M+
			   (3.*dot3*dot8-a1*a2*(5.*dot12+(-9.+8.*a1)*dot5*dot7-6.*dot3*dot8)-dot2*dot9+a1*((-2.-3.*a1+6.*a12)*dot2+(-10.+3.*(5.-2.*a1)*a1)*dot4)*dot9)*M-
			   a1*a2*(7.*(-1.+2.*a1)*dot10+(-7.+2.*a1*(-15.+8.*a1))*dot5*dot6)*M2))/(3.*a12*a22*(-(a22*M2)+uh)))/(-(a12*M2)+sh)
		  +((8.*(-5.*dot12+10.*dot5*dot7-dot1*dot8-2.*a1*dot1*dot8+9.*dot3*dot8-2.*a1*dot3*dot8-10.*dot2*dot9+2.*a1*a2*(dot10-2.*dot5*dot6)*M)*M2)/(3.*a1*a22)
		    +(64.*(-(a2*dot2)-a1*dot4)*pow<4,1>(M)*(-(dot3*dot7)+a1*(dot1+2.*dot3)*dot7-2.*a1*a2*(dot1+2.*dot3)*dot6*M))/(a1*a2*sqr(-(a22*M2)+uh))
		    +(16.*M2*(-((dot1+dot3)*((10.+a1*(-15.+8.*a1))*dot2+(-9.+8.*a1)*a2*dot4)*dot7)+
			      (-(a2*(5.*dot13*dot2-6.*dot11*dot3+a1*(6.*dot1*dot11-7.*dot13*dot2+8.*dot11*dot3+3.*dot13*dot4)))+
			       (a2*dot2*(dot1+2.*a1*(-13.+8.*a1)*dot1+3.*dot3-16.*a1*a2*dot3)+a1*((-3.+2.*a1)*(-5.+8.*a1)*dot1+(9.+4.*a1*(-7.+4.*a1))*dot3)*dot4)*dot6)*M-
			      12.*a12*a22*dot5*dot6*pow<3,1>(M)+2.*a1*a2*(-3.*a2*dot5*dot7-3.*(dot1+dot3)*dot8-4.*a2*dot2*dot9-4.*a1*dot4*dot9)*M2))/(3.*a1*a22*(-(a22*M2)+uh)))/sqr(-(a12*M2)+sh);
	      // diagram weights
	      save[0]+=norm(diag[3]);
	      save[1]+=norm(diag[1]);
	      Complex aSum =  (-0.5*diag[0] + diag[1] + diag[2] - 0.5*diag[3])/3. + 1.5*(diag[0] + diag[3] + 2.*diag[4]);
	      meSum+=norm(aSum);
	      me(2*ih1,ih2,ih3,ih4)=aSum;
	    }
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
	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      complex<Energy> dot1 = rescaledMomenta()[0]*v3[ih3].wave();
	      complex<Energy> dot2 = rescaledMomenta()[1]*g1[ih1].wave();
	      complex<Energy> dot3 = rescaledMomenta()[1]*v3[ih3].wave();
	      complex<Energy> dot4 = rescaledMomenta()[3]*g1[ih1].wave();
	      Complex dot5 = g1[ih1].wave()*v3[ih3].wave();
	      complex<Energy> dot6=v4[ih4].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	      LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	      LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	      LorentzPolarizationVectorE vec3 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	      LorentzPolarizationVectorE vec4 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g1[ih1].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	      complex<Energy> dot7 = vec1*g1[ih1].wave();
	      complex<Energy> dot8 = vec1*v3[ih3].wave();
	      complex<Energy> dot9 = vec2*g1[ih1].wave();
	      complex<Energy2> dot10 = vec3*rescaledMomenta()[0];
	      complex<Energy2> dot11 = vec4*rescaledMomenta()[0];
	      complex<Energy2> dot12 = vec1*rescaledMomenta()[0];
	      complex<Energy2> dot13 = vec2*rescaledMomenta()[0];
	      // diagrams
	      diag[0]=(256.*a22*(dot1+dot3)*(dot10-2.*dot2*dot6)*pow<5,1>(M))/(a1*pow<4,1>(-(a22*M2)+sh))
		+(16.*pow<3,1>(M)*((-27.+4.*a1*(-1.+4.*a1))*(dot1+dot3)*(dot10-2.*dot2*dot6)+8.*a1*a2*(dot11+3.*(dot1+dot3)*dot7-2.*dot2*dot8)*M))/(3.*a12*pow<3,1>(-(a22*M2)+sh))
		-(8.*((9.+6.*a1)*dot11+2.*((7.+8.*a1)*dot1*dot7+(7.+8.*a1)*dot3*dot7-3.*(3.+2.*a1)*dot2*dot8+5.*a1*dot9*M))*M2)/(3.*a12*sqr(-(a22*M2)+sh));
	      diag[1]=(-8.*((-3.+2.*a1)*dot1*dot7+2.*a1*(-5.*dot2*dot8+5.*dot4*dot8+4.*dot5*dot6*M-4.*dot9*M))*M2)/(3.*a12*a2*sqr(-M2+th))
		+(8.*(-1.-2.*a1*a2)*(dot5*dot6-dot9)*M)/(3.*a12*a22*(-M2+th))
		+((-64.*dot3*(dot10-2.*dot4*dot6)*pow<3,1>(M))/a22-(128.*dot3*pow<4,1>(M)*(dot12*(a1*(dot2-dot4)+dot4)+a1*(dot10+2.*a1*(dot2-dot4)*dot6)*M))/(a2*(-M2+th)))/pow<3,1>(-(a12*M2)+uh)+
		((8.*(5.*dot11+(10.*dot1+dot3+2.*a1*dot3)*dot7+2.*a1*((-7.+2.*a1)*dot5*dot6-(-6.+a1)*dot9)*M)*M2)/(3.*a1*a22)
		 -(64.*pow<4,1>(M)*(a1*dot12*dot2*(-(a2*dot1)+dot3)-dot12*(a22*dot1+a1*dot3)*dot4-a1*a2*dot1*(dot10+2.*a1*(dot2-dot4)*dot6)*M))/(a1*a2*sqr(-M2+th))
		 +(8.*M2*(2.*dot12*(10.*dot1*dot2+dot3*(dot2+(17.-8.*a1)*a1*dot2+a1*(-15.+8.*a1)*dot4))+
			  (-11.*dot10*dot3+2.*(7.*a1*dot13*dot2+3.*a12*dot13*dot2-9.*a1*dot10*dot3+6.*a12*dot10*dot3+5.*dot13*dot4+2.*a1*dot13*dot4-
					       7.*a12*dot13*dot4-a1*(9.+4.*a1*(-7.+4.*a1))*dot3*(dot2-dot4)*dot6)+
			   4.*dot1*((-2.+a1+a12)*dot10-a1*(2.+3.*a1)*dot2*dot6-(6.+5.*a1)*a2*dot4*dot6))*M+
			  24.*a12*a22*dot5*dot6*pow<3,1>(M)-4.*a1*a2*(3.*dot11+3.*(-3.+a1)*dot12*dot5+6.*dot1*dot7+9.*dot3*dot7+4.*a1*(dot2-dot4)*dot8-2.*dot4*dot8)*M2))/(3.*a1*a22*(-M2+th)))/sqr(-(a12*M2)+uh)+
		((4.*(-1.-2.*a1*a2)*(2.*dot5*dot6-dot9)*M)/(3.*a12*a22)
		 +(128.*dot1*dot12*(dot2-dot4)*pow<4,1>(M))/(a1*pow<3,1>(-M2+th))
		 -(8.*M2*(20.*a1*dot12*dot3*(-dot2+dot4)+2.*dot1*dot12*((-3.+a1-8.*a12)*dot2+a1*(1.+8.*a1)*dot4)+
			  ((3.+4.*a12)*dot1*dot10+2.*a1*(8.*dot10*dot3+dot13*(dot2-10.*a1*dot2-9.*dot4+10.*a1*dot4))+2.*a1*dot1*((7.+2.*(13.-8.*a1)*a1)*dot2+dot4+2.*a1*(-13.+8.*a1)*dot4)*dot6)*M+
			  4.*a1*(dot11+(-5.+3.*a1)*dot12*dot5+2.*(-2.+3.*a1)*dot1*dot7+2.*a1*(dot2-dot4)*dot8)*M2))/(3.*a12*a2*sqr(-M2+th))
		 +(8.*(-1.-2.*a1*a2)*(2.*dot10*dot3-dot13*(dot2+dot4)+dot1*(dot10+2.*dot2*dot6))*M-
		   8.*a1*((19.+a1*(21.+2.*a1*(-23.+8.*a1)))*dot5*dot6+2.*(-2.+a1*(-10.+7.*a1))*dot9)*pow<3,1>(M)-
		   4.*(dot11-8.*a1*a2*dot11-2.*(dot12*dot5+2.*dot1*dot7)+
		       2.*a1*(-((-7.+8.*a1)*a2*dot12*dot5)+(-13.+6.*a1)*dot1*dot7+(-11.+8.*a1)*dot3*dot7+10.*dot4*dot8+a1*(-15.*dot2+6.*a1*dot2+5.*dot4-6.*a1*dot4)*dot8))*M2)/(3.*a12*a22*(-M2+th)))/(-(a12*M2)+uh);
	      diag[2]=(-8.*(dot1*(dot7+2.*a1*dot7)+2.*a2*(5.*dot2*dot8-5.*dot4*dot8-4.*dot5*dot6*M+4.*dot9*M))*M2)/(3.*a1*a22*sqr(-M2+th))
		+(8.*(1.+2.*a1*a2)*(dot5*dot6-dot9)*M)/(3.*a12*a22*(-M2+th))
		+((-64.*(dot1+dot3)*(dot10-2.*dot2*dot6)*pow<3,1>(M))/a12+
		  (128.*(dot1+dot3)*pow<4,1>(M)*(dot12*(a1*(dot2-dot4)+dot4)-a2*(dot10-2.*a2*(dot2-dot4)*dot6)*M))/(a1*(-M2+th)))/pow<3,1>(-(a22*M2)+sh)+
		((-8.*(-5.*dot11+10.*dot12*dot5+(-3.+2.*a1)*(dot1+dot3)*dot7+2.*a2*(-5.*dot5*dot6+(5.+a1)*dot9)*M)*M2)/(3.*a12*a2)
		 +(64.*pow<4,1>(M)*(-(a2*dot12*dot3*(dot2-dot4))+dot1*dot12*((a12-a2)*dot2+dot4-a12*dot4)-a1*a2*dot1*(dot10-2.*a2*(dot2-dot4)*dot6)*M))/(a1*a2*sqr(-M2+th))
		 -(8.*M2*(2.*dot12*(-((7.+8.*a1)*a2*dot2*(dot1+dot3))+10.*dot3*dot4-a1*(1.+8.*a1)*(dot1+dot3)*dot4)+
			  ((17.-2.*a1*(3.+4.*a1))*dot1*dot10-24.*a1*dot13*dot2+14.*a12*dot13*dot2+17.*dot10*dot3+6.*a1*dot10*dot3-12.*a12*dot10*dot3-
			   20.*dot13*dot4+26.*a1*dot13*dot4-6.*a12*dot13*dot4+2.*(3.+a1+4.*a12*(-5.+4.*a1))*dot3*(dot2-dot4)*dot6+
			   2.*dot1*((3.+a1*(3.+8.*a1*(-3.+2.*a1)))*dot2+(7.+4.*a1*(-1.+4.*a1))*a2*dot4)*dot6)*M-
			  24.*a12*a22*dot5*dot6*pow<3,1>(M)+4.*a1*a2*(3.*dot11+3.*a1*dot12*dot5+9.*(dot1+dot3)*dot7-6.*dot2*dot8+4.*(a1*(dot2-dot4)+dot4)*dot8)*M2))/(3.*a12*a2*(-M2+th)))/sqr(-(a22*M2)+sh)+
		((4.*(-1.-2.*a1*a2)*dot9*M)/(3.*a12*a22)
		 -(128.*dot1*dot12*(dot2-dot4)*pow<4,1>(M))/(a2*pow<3,1>(-M2+th))
		 +(8.*M2*(dot1*(2.*dot12*(dot2+(7.-8.*a1)*a1*dot2+a1*(-5.+8.*a1)*dot4)+
				M*((-9.+4.*a1*(2.+a1))*dot10-2.*a2*((7.+2.*(13.-8.*a1)*a1)*dot2*dot6+(1.+2.*a1*(-13.+8.*a1))*dot4*dot6+12.*a1*dot7*M)))-
			  2.*a2*(10.*dot12*dot3*(-dot2+dot4)+(8.*dot10*dot3+dot13*(dot2-10.*a1*dot2+(-9.+10.*a1)*dot4))*M+2.*(dot11+3.*a1*dot12*dot5-2.*a2*(dot2-dot4)*dot8)*M2)))/(3.*a1*a22*sqr(-M2+th))
		 +(4.*M*(-4.*dot10*dot3-2.*(-1.-2.*a1*a2)*dot13*(dot2+dot4)+2.*(-1.-2.*a1*a2)*dot1*(dot10+2.*dot2*dot6)-dot11*M-2.*(7.+2.*a1*(-7.+5.*a1))*dot1*dot7*M-
			 2.*a2*(4.*a1*dot10*dot3+(-3.*dot3*dot7-a1*(4.*dot11+(-9.+8.*a1)*dot12*dot5+8.*dot3*dot7)+
						  9.*dot2*dot8+(-9.*dot4+a1*(7.*dot2-6.*a1*dot2+3.*dot4+6.*a1*dot4))*dot8)*M+
				((-10.+a1*(7.+2.*(15.-8.*a1)*a1))*dot5*dot6+2.*(5.+(4.-7.*a1)*a1)*dot9)*M2)))/(3.*a12*a22*(-M2+th)))/(-(a22*M2)+sh);
	      diag[3]=(256.*a12*dot3*(dot10-2.*dot4*dot6)*pow<5,1>(M))/(a2*pow<4,1>(-(a12*M2)+uh))
		+(16.*pow<3,1>(M)*((-15.+4.*a1*(-7.+4.*a1))*dot3*(dot10-2.*dot4*dot6)+8.*a1*a2*(dot11-2.*dot12*dot5+2.*dot1*dot7+3.*dot3*dot7-2.*dot4*dot8)*M))/(3.*a22*pow<3,1>(-(a12*M2)+uh))
		+(8.*(2.*a1*(3.*dot11-6.*dot12*dot5+6.*dot1*dot7+8.*dot3*dot7-6.*dot4*dot8+10.*dot5*dot6*M-5.*dot9*M)+
		      5.*(-3.*dot11+6.*dot12*dot5-6.*(dot1+dot3)*dot7+6.*dot4*dot8-4.*dot5*dot6*M+2.*dot9*M))*M2)/(3.*a22*sqr(-(a12*M2)+uh));
	      diag[4]=(64.*dot3*(dot10-2.*dot4*dot6)*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+uh))
		+(8.*(5.*dot11-10.*dot12*dot5+10.*dot1*dot7+(9.-2.*a1)*dot3*dot7-10.*dot4*dot8+2.*a1*a2*(2.*dot5*dot6-dot9)*M)*M2)/(3.*a1*a22*sqr(-(a12*M2)+uh))
		-(4.*(-1.-2.*a1*a2)*(2.*dot5*dot6-dot9)*M)/(3.*a12*a22*(-(a12*M2)+uh))
		+((64.*(dot1+dot3)*(dot10-2.*dot2*dot6)*pow<3,1>(M))/a12+(128.*(dot1+dot3)*(a1*(dot2-dot4)+dot4)*pow<4,1>(M)*(dot12-2.*a2*dot6*M))/(a1*(-(a12*M2)+uh)))/pow<3,1>(-(a22*M2)+sh)+
		((-8.*(-5.*dot11-(7.+2.*a1)*(dot1+dot3)*dot7+10.*dot2*dot8-2.*a1*a2*dot9*M)*M2)/(3.*a12*a2)
		 +(64.*(a1*(dot2-dot4)+dot4)*pow<4,1>(M)*(-(a2*dot1*dot12)+(-1.+2.*a1)*dot12*dot3-2.*a1*a2*(dot1+2.*dot3)*dot6*M))/(a1*a2*sqr(-(a12*M2)+uh))
		 -(16.*M2*(dot12*(dot1+dot3)*((3.+a1*(-1.+8.*a1))*dot2-a1*(1.+8.*a1)*dot4)+
			   (6.*a1*a2*dot1*dot10+2.*(1.-4.*a1)*a1*dot10*dot3+a1*dot13*((-2.+7.*a1)*dot2+3.*a2*dot4)+
			    a1*dot2*((-5.+4.*a1*(-5.+4.*a1))*dot1+3.*dot3-16.*a1*a2*dot3)*dot6-(3.+a1+4.*a12*(-5.+4.*a1))*(dot1+dot3)*dot4*dot6)*M-
			   12.*a12*a22*dot5*dot6*pow<3,1>(M)+2.*a1*a2*(3.*a1*dot12*dot5+3.*(dot1+dot3)*dot7+4.*(a1*(dot2-dot4)+dot4)*dot8)*M2))
		 /(3.*a12*a2*(-(a12*M2)+uh)))/sqr(-(a22*M2)+sh)+
		((4.*(1.+2.*a1*a2)*dot9*M)/(3.*a12*a22)
		 -(128.*dot3*(a1*(dot2-dot4)+dot4)*pow<4,1>(M)*(dot12+2.*a1*dot6*M))/(a2*pow<3,1>(-(a12*M2)+uh))
		 +(16.*M2*(dot12*dot3*((-9.+8.*a1)*a2*dot2+(10.+a1*(-15.+8.*a1))*dot4)+
			   (-(a2*(2.*(-3.+a1)*dot1*dot10-6.*dot10*dot3+5.*dot13*dot4+a1*(3.*dot13*dot2+8.*dot10*dot3-7.*dot13*dot4)))-
			    a1*dot2*(-6.*a2*dot1+(9.+4.*a1*(-7.+4.*a1))*dot3)*dot6-a2*(2.*(1.+5.*a1)*dot1+3.*dot3-16.*a1*a2*dot3)*dot4*dot6)*M+
			   12.*a12*a22*dot5*dot6*pow<3,1>(M)-2.*a1*a2*(-3.*a2*dot12*dot5+3.*dot3*dot7+4.*(a1*dot2+dot4-a1*dot4)*dot8)*M2))/(3.*a1*a22*sqr(-(a12*M2)+uh))
		 +(8.*M*(-((-1.-2.*a1*a2)*(-2.*dot10*dot3+dot13*(dot2+dot4)))+(-1.-2.*a1*a2)*dot1*(dot10+2.*dot2*dot6)+
			 (3.+8.*a1)*a2*dot1*dot7*M-(-3.*dot3*dot7-a1*a2*(5.*dot11+(-9.+8.*a1)*dot12*dot5+6.*dot3*dot7)+
						    (6.*pow(a1,3)*(dot2-dot4)+dot4+3.*a12*(-5.*dot2+dot4)+2.*a1*(5.*dot2+dot4))*dot8)*M+
			 a1*a2*((-7.+2.*a1*(-15.+8.*a1))*dot5*dot6+7.*(-1.+2.*a1)*dot9)*M2))/(3.*a12*a22*(-(a12*M2)+uh)))/(-(a22*M2)+sh);
	      // diagram weights
	      save[0]+=norm(diag[3]);
	      save[1]+=norm(diag[1]);
	      Complex aSum = (-0.5*diag[0] + diag[1] + diag[2] - 0.5*diag[3])/3. + 1.5*(diag[0] + diag[3] + 2.*diag[4]);
	      meSum+=norm(aSum);
	      me(2*ih1,ih2,ih3,ih4)=aSum;
	    }
	  }
	}
      }
    }
    // set matrix element
    setME(me);
    // final factors
    meSum *= 1./34560.;
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
    ProductionMatrixElement me(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1,PDT::Spin1);
    for(unsigned int ih1=0;ih1<2;++ih1) { 
      for(unsigned int ih2=0;ih2<2;++ih2) {
	complex<Energy> dot6=u1[ih1].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	LorentzPolarizationVectorE vec1 = u1[ih1].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	complex<Energy2> dot7 = vec1*rescaledMomenta()[3];
	for(unsigned int ih3=0;ih3<3;++ih3) {
	  complex<Energy> dot2 = rescaledMomenta()[0]*v3[ih3].wave();
	  complex<Energy> dot4 = rescaledMomenta()[1]*v3[ih3].wave();
	  LorentzPolarizationVectorE vec3 = u1[ih1].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	  complex<Energy> dot8 = vec1*v3[ih3].wave();
	  complex<Energy2> dot12 = vec3*rescaledMomenta()[3];
	  for(unsigned int ih4=0;ih4<2;++ih4) {
	    complex<Energy> dot1 = rescaledMomenta()[0]*g4[ih4].wave();
	    complex<Energy> dot3 = rescaledMomenta()[1]*g4[ih4].wave();
	    Complex dot5 = v3[ih3].wave()*g4[ih4].wave();
	    LorentzPolarizationVectorE vec2 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	    LorentzPolarizationVectorE vec4 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).slash(v3[ih3].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	    complex<Energy> dot9 = vec1*g4[ih4].wave();
	    complex<Energy2> dot10 = vec2*rescaledMomenta()[3];
	    complex<Energy> dot11 = vec2*v3[ih3].wave();
	    complex<Energy2> dot13 = vec4*rescaledMomenta()[3];
	    // diagrams
	    diag[0]=(-64.*dot2*(dot10+2.*dot3*dot6)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+th))
	      -(8.*(5.*dot13-10.*dot5*dot7-10.*dot3*dot8-7.*dot2*dot9-2.*a1*dot2*dot9-2.*a1*a2*(dot11-2.*dot5*dot6)*M)*M2)/(3.*a12*a2*sqr(-(a22*M2)+th))
	      -(4.*(-1.-2.*a1*a2)*(dot11-2.*dot5*dot6)*M)/(3.*a12*a22*(-(a22*M2)+th))
	      +((64.*dot4*(dot10-2.*dot1*dot6)*pow<3,1>(M))/a22+(128.*(-(a2*dot1)+a1*dot3)*dot4*pow<4,1>(M)*(-dot7+2.*a1*dot6*M))/(a2*(-(a22*M2)+th)))/pow<3,1>(-(a12*M2)+uh)+
	      ((8.*(-5.*dot13-10.*dot1*dot8+(10.*dot2+dot4+2.*a1*dot4)*dot9-2.*a1*a2*dot11*M)*M2)/(3.*a1*a22)
	       +(64.*(-(a2*dot1)+a1*dot3)*pow<4,1>(M)*((dot2-a1*dot2+a1*dot4)*dot7-2.*a1*a2*(dot2-dot4)*dot6*M))/(a1*a2*sqr(-(a22*M2)+th))
	       +(16.*M2*(-(((10.+a1*(-15.+8.*a1))*dot1-(-9.+8.*a1)*a2*dot3)*dot4*dot7)+
			 (-(a2*((-5.+7.*a1)*dot1*dot12+3.*a1*dot12*dot3+dot10*(6.*dot2-2.*a1*dot2+6.*a1*dot4)))+
			  (a1*dot3*(6.*a2*dot2+(-3.+2.*a1)*(-5.+8.*a1)*dot4)-a2*dot1*(dot4-2.*(dot2+5.*a1*dot2+(13.-8.*a1)*a1*dot4)))*dot6)*M-
			 12.*a12*a22*dot5*dot6*pow<3,1>(M)-2.*a1*a2*(-3.*a2*dot5*dot7+4.*dot1*dot8-4.*a1*(dot1+dot3)*dot8-3.*dot4*dot9)*M2))/(3.*a1*a22*(-(a22*M2)+th)))/sqr(-(a12*M2)+uh)+
	      ((4.*(-1.-2.*a1*a2)*dot11*M)/(3.*a12*a22)
	       -(128.*dot2*(-(a2*dot1)+a1*dot3)*pow<4,1>(M)*(dot7+2.*a2*dot6*M))/(a1*pow<3,1>(-(a22*M2)+th))
	       -(16.*M2*(-(dot2*(a1*(dot1-dot3)+3.*dot3+8.*a12*(dot1+dot3))*dot7)+
			 (a1*(-3.*a2*dot1*dot12-6.*a2*dot10*dot2+(-2.+7.*a1)*dot12*dot3-2.*(2.+a1)*dot10*dot4)+
			  ((3.+a1+4.*a12*(-5.+4.*a1))*dot1*dot2+a1*dot3*((-5.+4.*a1*(-5.+4.*a1))*dot2-4.*(2.+a1)*dot4))*dot6)*M+
			 12.*a12*a22*dot5*dot6*pow<3,1>(M)+2.*a1*a2*(3.*a1*dot5*dot7+4.*dot1*dot8-4.*a1*(dot1+dot3)*dot8+3.*dot2*dot9)*M2))/(3.*a12*a2*sqr(-(a22*M2)+th))
	       +(8.*M*((-1.-2.*a1*a2)*(dot1*dot12-dot12*dot3+dot10*(-dot2+dot4)+2.*dot3*(dot2+dot4)*dot6)-
		       (dot1*dot8-3.*dot2*dot9+a1*(5.*a2*dot13-(1.+8.*a1)*a2*dot5*dot7+2.*dot1*dot8-10.*dot3*dot8+
						   3.*a1*(dot1-2.*a1*dot1+(5.-2.*a1)*dot3)*dot8-5.*dot2*dot9+(dot4+2.*a1*(4.*dot2+dot4))*dot9))*M-
		       a1*a2*((7.-14.*a1)*dot11+(-21.+2.*a1*(-1.+8.*a1))*dot5*dot6)*M2))/(3.*a12*a22*(-(a22*M2)+th)))/(-(a12*M2)+uh);
	    diag[1]=(-256.*a22*dot2*(dot10+2.*dot3*dot6)*pow<5,1>(M))/(a1*pow<4,1>(-(a22*M2)+th))
	      +(16.*pow<3,1>(M)*(-((-27.+4.*a1*(-1.+4.*a1))*dot2*(dot10+2.*dot3*dot6))-8.*a1*a2*(dot13-2.*dot5*dot7-2.*dot3*dot8-3.*dot2*dot9)*M))/(3.*a12*pow<3,1>(-(a22*M2)+th))
	      +(8.*(3.*(3.+2.*a1)*(dot13-2.*(dot5*dot7+dot3*dot8))-2.*(7.+8.*a1)*dot2*dot9-10.*a1*(dot11-2.*dot5*dot6)*M)*M2)/(3.*a12*sqr(-(a22*M2)+th));
	    diag[2]=(64.*dot2*(dot10+2.*dot3*dot6)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+th))
	      -(8.*(5.*dot13+(-3.+2.*a1)*dot2*dot9+2.*a2*((5.+a1)*dot11-(5.+2.*a1)*dot5*dot6)*M)*M2)/(3.*a12*a2*sqr(-(a22*M2)+th))
	      +(4.*(-1.-2.*a1*a2)*(dot11-2.*dot5*dot6)*M)/(3.*a12*a22*(-(a22*M2)+th))
	      +(128.*(dot1+dot3)*(dot2+dot4)*dot7*pow<4,1>(M))/(a2*pow<3,1>(-M2+sh)*(-(a22*M2)+th))
	      +((8.*(-1.-2.*a1*a2)*(dot11-dot5*dot6)*M)/(3.*a12*a22)
	    	+(128.*dot2*pow<4,1>(M)*((dot1-a1*dot1-a1*dot3)*dot7+a2*(dot10+2.*a2*(dot1+dot3)*dot6)*M))/(a1*pow<3,1>(-(a22*M2)+th))
	    	-(8.*M2*(-2.*(a1*(1.+8.*a1)*dot1*dot2-(7.+8.*a1)*a2*dot2*dot3+10.*dot1*dot4)*dot7+
			 ((-17.+6.*a1+8.*a12)*dot10*dot2-4.*(-3.+a1)*a1*dot10*dot4-
			  2.*a2*dot1*((-10.+3.*a1)*dot12+(7.+4.*a1*(-1.+4.*a1))*dot2*dot6+10.*dot4*dot6)+
			  2.*dot3*(a1*(-12.+7.*a1)*dot12+(3.+a1*(3.+8.*a1*(-3.+2.*a1)))*dot2*dot6+2.*(1.-2.*a1)*a1*dot4*dot6))*M+
			 24.*a12*a22*dot5*dot6*pow<3,1>(M)+4.*a1*a2*(-3.*dot13+3.*(2.+a1)*dot5*dot7+4.*dot1*dot8+6.*dot3*dot8-4.*a1*(dot1+dot3)*dot8+9.*dot2*dot9)*M2))
		/(3.*a12*a2*sqr(-(a22*M2)+th))
	    	+(4.*M*(2.*(-1.-2.*a1*a2)*(dot1*dot12-dot12*dot3+dot10*(-dot2+dot4)+2.*dot3*(dot2+dot4)*dot6)+
			(dot13-8.*a1*a2*dot13-2.*(1.+a1-9.*a12+8.*pow(a1,3))*dot5*dot7-
			 2.*a2*(3.*(-3.+a1+2.*a12)*dot1+(-9.+a1*(-7.+6.*a1))*dot3)*dot8-2.*((7.+2.*a1*(-7.+5.*a1))*dot2+(-2.+a1)*(-5.+2.*a1)*dot4)*dot9)*M-
			2.*a2*(2.*(5.+(4.-7.*a1)*a1)*dot11+(-10.+a1*(-23.+2.*a1*(-1.+8.*a1)))*dot5*dot6)*M2))/(3.*a12*a22*(-(a22*M2)+th)))/(-M2+sh)
	      +((8.*(10.*a2*(dot1+dot3)*dot8-(1.+2.*a1)*(dot2+dot4)*dot9-8.*a2*(dot11-dot5*dot6)*M)*M2)/(3.*a1*a22)
	    	+(64.*pow<4,1>(M)*(-((((a12-a2)*dot2*dot3+a12*dot3*dot4-a2*dot1*(dot2+a1*dot2+a1*dot4))*dot7)/(a1*a2))+
				   (dot2+dot4)*(dot10+2.*a2*(dot1+dot3)*dot6)*M))/sqr(-(a22*M2)+th)+
		(8.*M2*(2.*(-(dot2*dot3)+10.*dot1*dot4+9.*dot3*dot4+8.*a12*(dot1+dot3)*(dot2+dot4)-a1*(5.*dot1*(dot2+3.*dot4)+dot3*(7.*dot2+17.*dot4)))*dot7-
			((-9.+4.*a1*(2.+a1))*dot10*dot2+(7.+4.*(-2.+a1)*a1)*dot10*dot4-2.*a2*dot3*((-1.+10.*a1)*dot12+(-7.+2.*a1*(-13.+8.*a1))*(dot2+dot4)*dot6)-
			 2.*a2*dot1*((-9.+10.*a1)*dot12+(1.+2.*a1*(-13.+8.*a1))*(dot2+dot4)*dot6))*M-
			4.*a2*(-dot13+(2.+3.*a1)*dot5*dot7+2.*a2*(dot1+dot3)*dot8+6.*a1*(dot2+dot4)*dot9)*M2))/(3.*a1*a22*(-(a22*M2)+th)))/sqr(-M2+sh);
	    diag[3]=(256.*a12*dot4*(dot10-2.*dot1*dot6)*pow<5,1>(M))/(a2*pow<4,1>(-(a12*M2)+uh))
	      +(16.*pow<3,1>(M)*((-15.+4.*a1*(-7.+4.*a1))*dot4*(dot10-2.*dot1*dot6)-8.*a1*a2*(dot13+2.*dot1*dot8+(-2.*dot2+dot4)*dot9)*M))/(3.*a22*pow<3,1>(-(a12*M2)+uh))
	      -(8.*(3.*(-5.+2.*a1)*dot13+6.*(-5.+2.*a1)*dot1*dot8+6.*(5.-2.*a1)*dot2*dot9+4.*a1*dot4*dot9-10.*a2*dot11*M)*M2)/(3.*a22*sqr(-(a12*M2)+uh));
	    diag[4]=(-64.*dot4*(dot10-2.*dot1*dot6)*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+uh))
	      -(8.*(5.*dot13-10.*dot5*dot7-10.*dot2*dot9-9.*dot4*dot9+2.*a1*dot4*dot9+2.*a1*((-6.+a1)*dot11+5.*dot5*dot6)*M)*M2)/(3.*a1*a22*sqr(-(a12*M2)+uh))
	      +(4.*(1.+2.*a1*a2)*dot11*M)/(3.*a12*a22*(-(a12*M2)+uh))
	      -(128.*(dot1+dot3)*(dot2+dot4)*dot7*pow<4,1>(M))/(a1*pow<3,1>(-M2+sh)*(-(a12*M2)+uh))
	      +((-8.*(-1.-2.*a1*a2)*(dot11-dot5*dot6)*M)/(3.*a12*a22)
	    	-(128.*dot4*pow<4,1>(M)*((-(a2*dot1)+a1*dot3)*dot7+a1*(dot10-2.*a1*(dot1+dot3)*dot6)*M))/(a2*pow<3,1>(-(a12*M2)+uh))
	    	+(8.*M2*(-2.*(10.*dot2*dot3+a1*(-15.+8.*a1)*dot1*dot4-(-9.+8.*a1)*a2*dot3*dot4)*dot7+
			 (-4.*(-2.+a1+a12)*dot10*dot2+(-3.-22.*a1+8.*a12)*dot10*dot4+
			  2.*a1*dot3*((7.+3.*a1)*dot12-2.*(2.+3.*a1)*dot2*dot6+(5.+2.*a1*(-17.+8.*a1))*dot4*dot6)+
			  2.*dot1*(-((5.+7.*a1)*a2*dot12)-2.*(-6.+a1+5.*a12)*dot2*dot6+(12.+a1*(7.+2.*a1*(-19.+8.*a1)))*dot4*dot6))*M-
			 24.*a12*a22*dot5*dot6*pow<3,1>(M)+
			 4.*a1*a2*(3.*dot13+3.*a2*dot5*dot7+2.*dot1*dot8+4.*a1*(dot1+dot3)*dot8+3.*(-2.*dot2+dot4)*dot9)*M2))/(3.*a1*a22*sqr(-(a12*M2)+uh))
	    	+(4.*M*(2.*(-1.-2.*a1*a2)*(dot1*dot12-dot12*dot3+dot10*(-dot2+dot4)+2.*dot3*(dot2+dot4)*dot6)+
			(dot13-8.*a1*a2*dot13+2.*a1*(dot5*dot7+(7.-8.*a1)*a1*dot5*dot7+(-10.+a1*(-5.+6.*a1))*dot1*dot8+3.*a1*(-5.+2.*a1)*dot3*dot8)+
			 2.*((2.+(13.-6.*a1)*a1)*dot2+2.*(1.+a1+a12)*dot4)*dot9)*M+
			2.*a1*((4.+2.*(10.-7.*a1)*a1)*dot11+(11.+a1*(-19.+2.*a1*(-9.+8.*a1)))*dot5*dot6)*M2))/(3.*a12*a22*(-(a12*M2)+uh)))/(-M2+sh)
	      +((-8.*(-3.*(dot2+dot4)*dot9+2.*a1*(5.*(dot1+dot3)*dot8+dot2*dot9+dot4*dot9-4.*dot11*M+4.*dot5*dot6*M))*M2)/(3.*a12*a2)
	    	+(64.*pow<4,1>(M)*(((-(a2*dot2*(-(a2*dot1)+a1*dot3))+(dot1+(-3.+a1)*a1*dot1+(-2.+a1)*a1*dot3)*dot4)*dot7)/(a1*a2)
	    			   +(dot2+dot4)*(-dot10+2.*a1*(dot1+dot3)*dot6)*M))/sqr(-(a12*M2)+uh)-
		(8.*M2*(2.*(a1*dot1*(dot2-9.*dot4)+3.*dot3*(dot2+dot4)+8.*a12*(dot1+dot3)*(dot2+dot4)-a1*dot3*(dot2+11.*dot4))*dot7-
			(3.*dot10*(dot2+dot4)-2.*a1*(9.*dot1*dot12+dot12*dot3+8.*dot10*dot4)+2.*a1*(dot1-7.*dot3)*(dot2+dot4)*dot6+
			 32.*pow(a1,3)*(dot1+dot3)*(dot2+dot4)*dot6+4.*a12*(dot10*(dot2+dot4)+(dot1+dot3)*(5.*dot12-13.*(dot2+dot4)*dot6)))*M-
			4.*a1*(dot13+3.*a2*dot5*dot7+2.*a1*(dot1+dot3)*dot8-2.*(-2.+3.*a1)*(dot2+dot4)*dot9)*M2))/(3.*a12*a2*(-(a12*M2)+uh)))/sqr(-M2+sh);
	    // diagram weights
	    save[0]+=norm(diag[2]);
	    save[1]+=norm(diag[4]);
	    Complex aSum =1.5*(2.*diag[0] + diag[1] + diag[3]) + (-0.5*diag[1] + diag[2] - 0.5*diag[3] + diag[4])/3.;
	    meSum+=norm(aSum);
	    me(ih1,ih2,ih3,2*ih4)=aSum;
	  }
	}
      }
    }
    // final factors
    meSum *=  1./12960.;
  }
  // save the diagram weights
  meInfo(save);
  // final factors
  return O1_*pow<3,1>(Constants::pi*standardModel()->alphaS(scale())/M)*meSum/sqr(a1*a2*M2);
}

IBPtr MEPP2BC3D1Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2BC3D1Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2BC3D1Jet::doinit() {
  setBcState(543);
  MEPP2BCJetBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<2>(bcbar,principleQuantumNumber(),1,1);

}

void MEPP2BC3D1Jet::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2*GeV2*GeV2);
}

void MEPP2BC3D1Jet::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2*GeV2*GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2BC3D1Jet,MEPP2BCJetBase>
describeHerwigMEPP2BC3D1Jet("Herwig::MEPP2BC3D1Jet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2BC3D1Jet::Init() {

  static ClassDocumentation<MEPP2BC3D1Jet> documentation
    ("The MEPP2BC3D1Jet class implements the matrix element for g c -> B_c(3D1) b q and "
     "q bar -> Bc(3D1) g and charged conjugate");

}

