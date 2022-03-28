// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2BC3D3Jet class.
//

#include "MEPP2BC3D3Jet.h"
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
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include <numeric>

using namespace Herwig;

double MEPP2BC3D3Jet::me2() const {
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
  // B_C2 wavefunction
  Rank3TensorWaveFunction Bcw(rescaledMomenta()[2],mePartonData()[2],outgoing);
  vector<Rank3TensorWaveFunction> ten3;
  for(unsigned int ix=0;ix<7;++ix) {
    Bcw.reset(ix);
    ten3.push_back(Bcw);
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
    ProductionMatrixElement me(PDT::Spin1,PDT::Spin1Half,PDT::Spin3,PDT::Spin1Half);
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
	  for(unsigned int ih3=0;ih3<7;++ih3) {
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      auto dot1 = rescaledMomenta()[1]*g1[ih1].wave();
	      auto dot2 = rescaledMomenta()[3]*g1[ih1].wave();
	      complex<Energy> dot3=u2[ih2].dimensionedWave().scalar(ubar4[ih4].dimensionedWave());
	      auto vec1 = u2[ih2].dimensionedWave().vectorCurrent(ubar4[ih4].dimensionedWave());
	      auto vec2 = u2[ih2].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(ubar4[ih4].dimensionedWave());
	      auto vec3 = u2[ih2].dimensionedWave().slash(rescaledMomenta()[0]).vectorCurrent(ubar4[ih4].dimensionedWave());
	      auto vec4 = u2[ih2].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[0]).vectorCurrent(ubar4[ih4].dimensionedWave());
	      auto vec5 = ten3[ih3].wave().dot(rescaledMomenta()[0],0).preDot(rescaledMomenta()[0]);
	      auto vec6 = ten3[ih3].wave().dot(rescaledMomenta()[0],0).preDot(rescaledMomenta()[1]);
	      auto vec7 = ten3[ih3].wave().dot(rescaledMomenta()[0],0).preDot(g1[ih1].wave());
	      auto vec8 = ten3[ih3].wave().dot(rescaledMomenta()[1],0).preDot(rescaledMomenta()[1]);
	      auto vec9 = ten3[ih3].wave().dot(rescaledMomenta()[1],0).preDot(g1[ih1].wave());
	      auto dot4 = vec1*vec5;
	      auto dot5 = vec1*vec6;
	      auto dot6 = vec1*vec8;
	      auto dot7 = vec4*vec5;
	      auto dot8 = vec4*vec6;
	      auto dot9 = vec4*vec8;
	      auto dot10 = vec1*vec7;
	      auto dot11 = vec1*vec9;
	      auto dot12 = vec2*vec5;
	      auto dot13 = vec2*vec6;
	      auto dot14 = vec2*vec8;
	      auto dot15 = vec3*vec5;
	      auto dot16 = vec3*vec6;
	      auto dot17 = vec3*vec8;
	      auto dot18 = vec5*rescaledMomenta()[0];
	      auto dot19 = vec2*rescaledMomenta()[0];
	      auto dot20 = vec5*rescaledMomenta()[1];
	      auto dot21 = vec5*g1[ih1].wave();
	      auto dot22 = vec6*rescaledMomenta()[1];
	      auto dot23 = vec6*g1[ih1].wave();
	      auto dot24 = vec8*rescaledMomenta()[1];
	      auto dot25 = vec8*g1[ih1].wave();
	      auto dot26 = vec1*g1[ih1].wave();
	      auto dot27 = vec1*rescaledMomenta()[0];
	      // diagrams
	      diag[0]=(64.*a1*(2.*dot1*(dot4+2.*dot5+dot6)+dot7+2.*dot8+dot9)*pow<4,1>(M))/(a22*pow<4,1>(-(a12*M2)+sh));
	      diag[1]=(-64.*(dot12-dot21*dot3)*pow<3,1>(M))/(a12*pow<3,1>(-M2+th))
		+((32.*pow<3,1>(M)*(2.*(dot15*(dot1+dot2)+dot19*(dot18+2.*dot20)-2.*dot1*dot18*dot3)-(2.*a1*(dot1-dot2)*dot4+dot7)*M))/(a12*pow<3,1>(-M2+th))
		  +(32.*pow<3,1>(M)*(-((1.+a1)*dot12)+a1*(-2.*dot13+2.*(dot21+dot23)*dot3-a2*dot10*M)))/(a12*a2*sqr(-M2+th)))/(-(a12*M2)+sh)
		+((-32.*(dot12+2.*dot13+dot14)*pow<3,1>(M))/a22+
		  (32.*pow<3,1>(M)*(2.*(dot15+2.*dot16+dot17)*(dot1+dot2)+2.*dot19*(dot18+4.*dot20+5.*dot22+2.*dot24)-
				    4.*dot1*(dot18+2.*dot20+dot22)*dot3-(2.*a1*(dot1-dot2)*(dot4+2.*dot5+dot6)+dot7+2.*dot8+dot9)*M))/(a22*(-M2+th)))/pow<3,1>(-(a12*M2)+sh)+
		((-32.*pow<3,1>(M)*(-2.*((dot15+dot16)*(dot1+dot2)+dot19*(dot18+3.*dot20+2.*dot22))+
				    4.*dot1*(dot18+dot20)*dot3+(2.*a1*(dot1-dot2)*(dot4+dot5)+dot7+dot8)*M))/(a1*a2*sqr(-M2+th))
		 +(-32.*(dot12+a1*dot12+dot13+3.*a1*dot13+2.*a1*dot14-2.*a1*(dot21+2.*dot23+dot25)*dot3)*pow<3,1>(M)-
		   32.*a1*a2*(dot10+dot11)*pow<4,1>(M))/(a1*a22*(-M2+th)))/sqr(-(a12*M2)+sh);
	      diag[2]=(64.*(dot12-dot21*dot3)*pow<3,1>(M))/(a22*pow<3,1>(-M2+th))
		+((32.*(dot14-2.*dot25*dot3)*pow<3,1>(M))/a12+
		  (32.*pow<3,1>(M)*(2.*dot17*(dot1+dot2)+2.*dot19*(dot22+2.*dot24)-4.*dot1*dot22*dot3+
				    (2.*dot22*dot26-2.*dot25*dot27+2.*a2*(dot1-dot2)*dot6-dot9)*M))/(a12*(-M2+th)))/pow<3,1>(-(a22*M2)+uh)+
		((-32.*pow<3,1>(M)*(2.*dot16*(dot1+dot2)+2.*dot19*(dot20+2.*dot22)-4.*dot1*dot20*dot3+(2.*dot20*dot26-2.*dot23*dot27+2.*a2*(dot1-dot2)*dot5-dot8)*M))/(a1*a2*sqr(-M2+th))
		 +(32.*pow<3,1>(M)*(2.*dot14-2.*dot25*dot3+a12*dot11*M-a1*(dot13+2.*dot14-2.*(dot23+dot25)*dot3+dot11*M)))/(a12*a2*(-M2+th)))/sqr(-(a22*M2)+uh)+
		((32.*pow<3,1>(M)*(2.*dot15*(dot1+dot2)+2.*dot19*(dot18+2.*dot20)-4.*dot1*dot18*dot3+(2.*dot18*dot26-2.*dot21*dot27+2.*a2*(dot1-dot2)*dot4-dot7)*M))/(a22*pow<3,1>(-M2+th))
		 -(32.*pow<3,1>(M)*(2.*dot13-2.*dot23*dot3+a12*dot10*M-a1*(dot12+2.*dot13-2.*(dot21+dot23)*dot3+dot10*M)))/(a1*a22*sqr(-M2+th)))/(-(a22*M2)+uh);
	      diag[3]=(-64.*a2*(2.*dot22*dot26-2.*dot25*dot27-2.*dot2*dot6-dot9)*pow<4,1>(M))/(a12*pow<4,1>(-(a22*M2)+uh));
	      diag[4]=(-32.*(dot14-2.*dot25*dot3)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+uh))
		+((64.*pow<3,1>(M)*(dot1*dot17+dot17*dot2+dot19*dot22+2.*dot19*dot24-2.*dot1*dot22*dot3+(dot1-a1*dot1+a1*dot2)*dot6*M))/(a12*pow<3,1>(-(a22*M2)+uh))
		  +(32.*pow<3,1>(M)*(dot14+a1*(-dot13-2.*dot14+2.*(dot23+dot25)*dot3-a2*dot11*M)))/(a12*a2*sqr(-(a22*M2)+uh)))/(-(a12*M2)+sh)
		+((-64.*pow<3,1>(M)*(-((dot16+dot17)*(dot1+dot2))-dot19*(dot20+3.*dot22+2.*dot24)+2.*dot1*(dot20+dot22)*dot3+(-(a2*dot1)-a1*dot2)*(dot5+dot6)*M))/(a1*a2*sqr(-(a22*M2)+uh))
		  +(32.*pow<3,1>(M)*(dot13+dot14+a1*(-dot12-3.*dot13-2.*dot14+2.*(dot21+2.*dot23+dot25)*dot3-a2*(dot10+dot11)*M)))/(a1*a22*(-(a22*M2)+uh)))/sqr(-(a12*M2)+sh)+
		((32.*(dot12+2.*dot13+dot14)*pow<3,1>(M))/a22+
		 (64.*pow<3,1>(M)*(dot19*(dot18+4.*dot20+5.*dot22+2.*dot24)+dot2*(dot15+2.*dot16+dot17+a1*(dot4+2.*dot5+dot6)*M)+
				   dot1*(dot15+2.*dot16+dot17-2.*(dot18+2.*dot20+dot22)*dot3+a2*(dot4+2.*dot5+dot6)*M)))/(a22*(-(a22*M2)+uh)))/pow<3,1>(-(a12*M2)+sh);
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
	  for(unsigned int ih3=0;ih3<7;++ih3) {
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      auto dot1 = rescaledMomenta()[1]*g1[ih1].wave();
	      auto dot2 = rescaledMomenta()[3]*g1[ih1].wave();
	      complex<Energy> dot3=v4[ih4].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	      auto vec1 = v4[ih4].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	      auto vec2 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	      auto vec3 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[0]).vectorCurrent(vbar2[ih2].dimensionedWave());
	      auto vec4 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[0]).vectorCurrent(vbar2[ih2].dimensionedWave());
	      auto vec5 = ten3[ih3].wave().dot(rescaledMomenta()[0],0).preDot(rescaledMomenta()[0]);
	      auto vec6 = ten3[ih3].wave().dot(rescaledMomenta()[0],0).preDot(rescaledMomenta()[1]);
	      auto vec7 = ten3[ih3].wave().dot(rescaledMomenta()[0],0).preDot(g1[ih1].wave());
	      auto vec8 = ten3[ih3].wave().dot(rescaledMomenta()[1],0).preDot(rescaledMomenta()[1]);
	      auto vec9 = ten3[ih3].wave().dot(rescaledMomenta()[1],0).preDot(g1[ih1].wave());
	      auto dot4 = vec1*vec5;
	      auto dot5 = vec1*vec6;
	      auto dot6 = vec1*vec8;
	      auto dot7 = vec4*vec5;
	      auto dot8 = vec4*vec6;
	      auto dot9 = vec4*vec8;
	      auto dot10 = vec1*g1[ih1].wave();
	      auto dot11 = vec5*rescaledMomenta()[0];
	      auto dot12 = vec5*rescaledMomenta()[1];
	      auto dot13 = vec1*rescaledMomenta()[0];
	      auto dot14 = vec5*g1[ih1].wave();
	      auto dot15 = vec6*rescaledMomenta()[1];
	      auto dot16 = vec6*g1[ih1].wave();
	      auto dot17 = vec8*g1[ih1].wave();
	      auto dot18 = vec1*vec7;
	      auto dot19 = vec1*vec9;
	      auto dot20 = vec2*vec5;
	      auto dot21 = vec2*vec6;
	      auto dot22 = vec2*vec8;
	      auto dot23 = vec3*vec5;
	      auto dot24 = vec3*vec6;
	      auto dot25 = vec3*vec8;
	      auto dot26 = vec2*rescaledMomenta()[0];
	      auto dot27 = vec8*rescaledMomenta()[1];
	      // diagrams
	      diag[0]=(64.*a2*(2.*dot10*(dot11+2.*dot12+dot15)-2.*dot13*(dot14+2.*dot16+dot17)+2.*dot1*(dot4+2.*dot5+dot6)-dot7-2.*dot8-dot9)*pow<4,1>(M))/(a12*pow<4,1>(-(a22*M2)+sh));
	      diag[1]=(64.*(dot20-dot14*dot3)*pow<3,1>(M))/(a12*pow<3,1>(-M2+th))
		+((32.*dot22*pow<3,1>(M))/a22+
		  (32.*pow<3,1>(M)*(2.*(dot1+dot2)*dot25+2.*dot26*(dot15+2.*dot27)-4.*dot15*dot2*dot3+2.*a1*(dot1-dot2)*dot6*M+dot9*M))/(a22*(-M2+th)))/pow<3,1>(-(a12*M2)+uh)+
		((-32.*pow<3,1>(M)*(2.*(dot1+dot2)*dot24+2.*(dot12+2.*dot15)*dot26-4.*dot12*dot2*dot3+2.*a1*(dot1-dot2)*dot5*M+dot8*M))/(a1*a2*sqr(-M2+th))
		 +(32.*pow<3,1>(M)*(-(a2*dot21)+a1*(2.*dot22-2.*dot17*dot3-a2*dot19*M)))/(a1*a22*(-M2+th)))/sqr(-(a12*M2)+uh)+
		((32.*pow<3,1>(M)*(2.*(dot1+dot2)*dot23+2.*(dot11+2.*dot12)*dot26-4.*dot11*dot2*dot3+2.*a1*(dot1-dot2)*dot4*M+dot7*M))/(a12*pow<3,1>(-M2+th))
		 -(32.*pow<3,1>(M)*(-(a2*dot20)+a1*(2.*dot21-2.*dot16*dot3-a2*dot18*M)))/(a12*a2*sqr(-M2+th)))/(-(a12*M2)+uh);
	      diag[2]=(-64.*(dot20-dot14*dot3)*pow<3,1>(M))/(a22*pow<3,1>(-M2+th))
		+((32.*pow<3,1>(M)*(2.*(dot11+2.*dot12)*dot26+(-2.*dot10*dot11+2.*dot13*dot14+dot7)*M+
				    2.*dot2*(dot23-2.*dot11*dot3+dot4*M-a1*dot4*M)+2.*dot1*(dot23-a2*dot4*M)))/(a22*pow<3,1>(-M2+th))
		  +(32.*pow<3,1>(M)*((-2.+a1)*dot20-2.*a2*dot21+2.*(dot14+dot16-a1*dot16)*dot3-a1*a2*dot18*M))/(a1*a22*sqr(-M2+th)))/(-(a22*M2)+sh)
		+((-32.*pow<3,1>(M)*(-2.*(dot1+dot2)*(dot23+dot24)-2.*(dot11+3.*dot12+2.*dot15)*dot26+4.*(dot11+dot12)*dot2*dot3+
				     (2.*dot10*(dot11+dot12)-2.*dot13*(dot14+dot16)+2.*a2*(dot1-dot2)*(dot4+dot5)-dot7-dot8)*M))/(a1*a2*sqr(-M2+th))
		  -(32.*pow<3,1>(M)*(-((-2.+a1)*dot20)+(4.-3.*a1)*dot21+2.*a2*dot22-
				     2.*(dot14+2.*dot16-a1*dot16+dot17-a1*dot17)*dot3+a1*a2*(dot18+dot19)*M))/(a12*a2*(-M2+th)))/sqr(-(a22*M2)+sh)+
		((-32.*(dot20+2.*dot21+dot22-2.*(dot14+2.*dot16+dot17)*dot3)*pow<3,1>(M))/a12+
		 (32.*pow<3,1>(M)*(2.*(dot1+dot2)*(dot23+2.*dot24+dot25)+2.*dot26*(dot11+4.*dot12+5.*dot15+2.*dot27)-
				   4.*(dot11+2.*dot12+dot15)*dot2*dot3+(-2.*dot10*(dot11+2.*dot12+dot15)+2.*dot13*(dot14+2.*dot16+dot17)-
									2.*a2*(dot1-dot2)*(dot4+2.*dot5+dot6)+dot7+2.*dot8+dot9)*M))/(a12*(-M2+th)))/pow<3,1>(-(a22*M2)+sh);
	      diag[3]=(64.*a1*(2.*dot2*dot6-dot9)*pow<4,1>(M))/(a22*pow<4,1>(-(a12*M2)+uh));
	      diag[4]=(-32.*dot22*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+uh))
		+((64.*pow<3,1>(M)*(dot26*(dot15+2.*dot27)+dot2*(dot25-2.*dot15*dot3+dot6*M-a1*dot6*M)+dot1*(dot25+a1*dot6*M)))/(a22*pow<3,1>(-(a12*M2)+uh))
		  +(32.*pow<3,1>(M)*(-(a2*dot21)-dot22+2.*a1*dot22-2.*a1*dot17*dot3-a1*a2*dot19*M))/(a1*a22*sqr(-(a12*M2)+uh)))/(-(a22*M2)+sh)
		+((-64.*pow<3,1>(M)*(-((dot1+dot2)*(dot24+dot25))-dot26*(dot12+3.*dot15+2.*dot27)+2.*(dot12+dot15)*dot2*dot3-
				     (a1*(dot1-dot2)+dot2)*(dot5+dot6)*M))/(a1*a2*sqr(-(a12*M2)+uh))
		  -(32.*pow<3,1>(M)*(dot20-a1*dot20+2.*dot21+dot22+a1*(-3.*dot21-2.*dot22+2.*(dot16+dot17)*dot3+a2*(dot18+dot19)*M)))/(a12*a2*(-(a12*M2)+uh)))/sqr(-(a22*M2)+sh)+
		((32.*(dot20+2.*dot21+dot22-2.*(dot14+2.*dot16+dot17)*dot3)*pow<3,1>(M))/a12+
		 (64.*pow<3,1>(M)*(dot26*(dot11+4.*dot12+5.*dot15+2.*dot27)+dot1*(dot23+2.*dot24+dot25+a1*(dot4+2.*dot5+dot6)*M)+
				   dot2*(dot23+2.*dot24+dot25-2.*(dot11+2.*dot12+dot15)*dot3+a2*(dot4+2.*dot5+dot6)*M)))/(a12*(-(a12*M2)+uh)))/pow<3,1>(-(a22*M2)+sh);
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
    meSum *= 1./12096.;
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
    ProductionMatrixElement me(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin3,PDT::Spin1);
    for(unsigned int ih1=0;ih1<2;++ih1) { 
      for(unsigned int ih2=0;ih2<2;++ih2) {
	for(unsigned int ih3=0;ih3<7;++ih3) {
	  for(unsigned int ih4=0;ih4<2;++ih4) {
	    auto dot1 = rescaledMomenta()[0]*g4[ih4].wave();
	    auto dot2 = rescaledMomenta()[1]*g4[ih4].wave();
	    complex<Energy> dot3=u1[ih1].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	    auto vec1 = u1[ih1].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	    auto vec2 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	    auto vec3 = u1[ih1].dimensionedWave().slash(rescaledMomenta()[3]).vectorCurrent(vbar2[ih2].dimensionedWave());
	    auto vec4 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).slash(rescaledMomenta()[3]).vectorCurrent(vbar2[ih2].dimensionedWave());
	    auto vec5 = ten3[ih3].wave().dot(rescaledMomenta()[0],0).preDot(rescaledMomenta()[0]);
	    auto vec6 = ten3[ih3].wave().dot(rescaledMomenta()[0],0).preDot(rescaledMomenta()[1]);
	    auto vec7 = ten3[ih3].wave().dot(rescaledMomenta()[0],0).preDot(g4[ih4].wave());
	    auto vec8 = ten3[ih3].wave().dot(rescaledMomenta()[1],0).preDot(rescaledMomenta()[1]);
	    auto vec9 = ten3[ih3].wave().dot(rescaledMomenta()[1],0).preDot(g4[ih4].wave());
	    auto dot4 = vec1*vec5;
	    auto dot5 = vec1*vec6;
	    auto dot6 = vec1*vec7;
	    auto dot7 = vec1*vec8;
	    auto dot8 = vec1*vec9;
	    auto dot9 = vec2*vec5;
	    auto dot10 = vec2*vec6;
	    auto dot11 = vec2*vec8;
	    auto dot12 = vec3*vec5;
	    auto dot13 = vec3*vec6;
	    auto dot14 = vec3*vec8;
	    auto dot15 = vec5*rescaledMomenta()[0];
	    auto dot16 = vec2*rescaledMomenta()[3];
	    auto dot17 = vec5*rescaledMomenta()[1];
	    auto dot18 = vec5*g4[ih4].wave();
	    auto dot19 = vec6*rescaledMomenta()[1];
	    auto dot20 = vec6*g4[ih4].wave();
	    auto dot21 = vec8*rescaledMomenta()[1];
	    auto dot22 = vec8*g4[ih4].wave();
	    auto dot23 = vec4*vec5;
	    auto dot24 = vec1*g4[ih4].wave();
	    auto dot25 = vec1*rescaledMomenta()[3];
	    auto dot26 = vec4*vec6;
	    auto dot27 = vec4*vec8;
	    // diagrams
	    diag[0]=(32.*(2.*dot18*dot3-dot9)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+th))
	      +((32.*dot11*pow<3,1>(M))/a22-
		(64.*pow<3,1>(M)*(-(dot14*dot2)+dot16*(dot19-dot21)+a1*dot2*dot7*M+dot1*(dot14-2.*(dot19+dot21)*dot3-a2*dot7*M)))/(a22*(-(a22*M2)+th)))/pow<3,1>(-(a12*M2)+uh)+
	      ((64.*pow<3,1>(M)*(dot16*(dot17-dot19)-dot13*dot2+a1*dot2*dot5*M+dot1*(dot13-2.*(dot17+dot19)*dot3-a2*dot5*M)))/(a1*a2*sqr(-(a22*M2)+th))
	       -(32.*pow<3,1>(M)*(dot10-a1*dot10+a1*(dot11-2.*dot22*dot3-a2*dot8*M)))/(a1*a22*(-(a22*M2)+th)))/sqr(-(a12*M2)+uh)+
	      ((-64.*pow<3,1>(M)*(dot15*dot16-dot16*dot17-dot12*dot2+a1*dot2*dot4*M+dot1*(dot12-2.*(dot15+dot17)*dot3-a2*dot4*M)))/(a12*pow<3,1>(-(a22*M2)+th))
	       +(32.*pow<3,1>(M)*(dot9+a1*(dot10-2.*dot20*dot3-dot9-a2*dot6*M)))/(a12*a2*sqr(-(a22*M2)+th)))/(-(a12*M2)+uh);
	    diag[1]=(-64.*a2*(dot23-2.*(dot15+dot17)*dot24+2.*dot18*dot25+2.*dot2*dot4)*pow<4,1>(M))/(a12*pow<4,1>(-(a22*M2)+th));
	    diag[2]=(32.*(-2.*dot18*dot3+dot9)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+th))
	      +((32.*pow<3,1>(M)*(2.*(-(dot1*dot12)-dot15*dot16+dot16*dot17+dot12*dot2+2.*dot1*(dot15+dot17)*dot3)+
				  (dot23-2.*(dot15+dot17)*dot24+2.*dot18*dot25+2.*a2*(dot1+dot2)*dot4)*M))/(a12*pow<3,1>(-(a22*M2)+th))
		+(32.*pow<3,1>(M)*(-2.*dot18*dot3+2.*dot9+a1*(dot10-2.*dot20*dot3-dot9-a2*dot6*M)))/(a12*a2*sqr(-(a22*M2)+th)))/(-M2+sh)
	      +((64.*(2.*dot10+dot11-(dot18+2.*dot20+dot22)*dot3+dot9)*pow<3,1>(M))/a22+
		(32.*pow<3,1>(M)*(-2.*(dot12+2.*dot13+dot14)*(dot1-dot2)+2.*dot16*(-dot15-dot17+dot19+dot21)+
				  4.*dot1*(dot15+3.*(dot17+dot19)+dot21)*dot3+
				  (dot23-2.*(dot15+3.*(dot17+dot19)+dot21)*dot24+2.*(dot18+2.*dot20+dot22)*dot25+
				   2.*dot26+dot27+2.*a2*(dot1+dot2)*(dot4+2.*dot5+dot7))*M))/(a22*(-(a22*M2)+th)))/pow<3,1>(-M2+sh)+
	      ((-32.*pow<3,1>(M)*(2.*(dot16*(dot15-dot19)+(dot12+dot13)*(dot1-dot2)-2.*dot1*(dot15+2.*dot17+dot19)*dot3)-
				  (dot23-2.*(dot15+2.*dot17+dot19)*dot24+2.*(dot18+dot20)*dot25+dot26+2.*a2*(dot1+dot2)*(dot4+dot5))*M))/(a1*a2*sqr(-(a22*M2)+th))
	       +(32.*pow<3,1>(M)*(2.*dot10-2.*(dot18+dot20)*dot3+2.*dot9+a1*(dot11-2.*(dot20+dot22)*dot3-dot9-a2*(dot6+dot8)*M)))/(a1*a22*(-(a22*M2)+th)))/sqr(-M2+sh);
	    diag[3]=(-64.*a1*(dot27-2.*dot1*dot7)*pow<4,1>(M))/(a22*pow<4,1>(-(a12*M2)+uh));
	    diag[4]=(-32.*dot11*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+uh))
	      +((-32.*pow<3,1>(M)*(-2.*dot14*dot2+2.*dot16*(dot19-dot21)-dot27*M+2.*a1*dot2*dot7*M+2.*dot1*(dot14-2.*(dot19+dot21)*dot3+a1*dot7*M)))/(a22*pow<3,1>(-(a12*M2)+uh))
		-(32.*pow<3,1>(M)*(dot10-a1*dot10+dot11+a1*(dot11-2.*dot22*dot3-a2*dot8*M)))/(a1*a22*sqr(-(a12*M2)+uh)))/(-M2+sh)
	      +((-64.*(2.*dot10+dot11-(dot18+2.*dot20+dot22)*dot3+dot9)*pow<3,1>(M))/a12+
		(32.*pow<3,1>(M)*(-2.*(dot12+2.*dot13+dot14)*(dot1-dot2)+2.*dot16*(-dot15-dot17+dot19+dot21)+
				  4.*dot1*(dot15+3.*(dot17+dot19)+dot21)*dot3+(dot23+2.*dot26+dot27-2.*a1*(dot1+dot2)*(dot4+2.*dot5+dot7))*M))/(a12*(-(a12*M2)+uh)))/pow<3,1>(-M2+sh)+
	      ((-32.*pow<3,1>(M)*(2.*((dot13+dot14)*(dot1-dot2)+dot16*(dot17-dot21)-2.*dot1*(dot17+2.*dot19+dot21)*dot3)-
				  (dot26+dot27-2.*a1*(dot1+dot2)*(dot5+dot7))*M))/(a1*a2*sqr(-(a12*M2)+uh))
	       -(32.*pow<3,1>(M)*(2.*dot10+(1.+a1)*dot11+dot9-a1*(2.*(dot20+dot22)*dot3+dot9)-a1*a2*(dot6+dot8)*M))/(a12*a2*(-(a12*M2)+uh)))/sqr(-M2+sh);
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
    meSum *= 1./4536.;
  }
  // save the diagram weights
  meInfo(save);
  // final factors
  return meSum*O1_*pow<3,1>(Constants::pi*standardModel()->alphaS(scale())/M)/sqr(M2*a1*a2);
}

IBPtr MEPP2BC3D3Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2BC3D3Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2BC3D3Jet::doinit() {
  setBcState(547);
  MEPP2BCJetBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<2>(bcbar,principleQuantumNumber(),1,3);


}

void MEPP2BC3D3Jet::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2*GeV2*GeV2);
}

void MEPP2BC3D3Jet::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2*GeV2*GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2BC3D3Jet,MEPP2BCJetBase>
describeHerwigMEPP2BC3D3Jet("Herwig::MEPP2BC3D3Jet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2BC3D3Jet::Init() {

  static ClassDocumentation<MEPP2BC3D3Jet> documentation
    ("The MEPP2BC3D3Jet class implements the matrix element for g c -> B_c(3D3) b q and "
     "q bar -> Bc(3D3) g and charged conjugate");

}

