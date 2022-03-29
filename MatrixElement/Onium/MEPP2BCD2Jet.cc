// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2BCD2Jet class.
//

#include "MEPP2BCD2Jet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
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
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/epsilon.h"
#include <numeric>

using namespace Herwig;

double MEPP2BCD2Jet::me2() const {
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
  // B_c1 wavefunction
  TensorWaveFunction Bcw(rescaledMomenta()[2],mePartonData()[2],outgoing);
  vector<TensorWaveFunction> ten3;
  for(unsigned int ix=0;ix<5;++ix) {
    Bcw.reset(ix);
    ten3.push_back(Bcw);
  }
  // stuff for mixing
  int itest = (abs(mePartonData()[2]->id())%100000)/10000;
  double mix1 = itest==1 ? sTheta_ :  cTheta_;
  double mix2 = itest==1 ? cTheta_ : -sTheta_;
  if(mePartonData()[2]->id()<0) mix2*=-1.;
  double rt(sqrt(6.));
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
    ProductionMatrixElement me(PDT::Spin1,PDT::Spin1Half,PDT::Spin2,PDT::Spin1Half);
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
	  for(unsigned int ih3=0;ih3<5;++ih3) {
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      auto dot1 = rescaledMomenta()[1]*g1[ih1].wave();
	      auto dot2 = rescaledMomenta()[3]*g1[ih1].wave();
	      complex<Energy> dot3=u2[ih2].dimensionedWave().pseudoScalar(ubar4[ih4].dimensionedWave());
	      auto vec1 = u2[ih2].dimensionedWave().generalCurrent(ubar4[ih4].dimensionedWave(),-1,1);
	      auto vec2 = u2[ih2].dimensionedWave().slash(g1[ih1].wave()).generalCurrent(ubar4[ih4].dimensionedWave(),1,-1);
	      auto vec3 = u2[ih2].dimensionedWave().slash(rescaledMomenta()[0]).generalCurrent(ubar4[ih4].dimensionedWave(),1,-1);
	      auto vec4 = u2[ih2].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[0]).generalCurrent(ubar4[ih4].dimensionedWave(),-1,1);
	      auto vec5 = ten3[ih3].wave().preDot(rescaledMomenta()[0]);
	      auto vec6 = ten3[ih3].wave().preDot(rescaledMomenta()[1]);
	      auto vec7 = ten3[ih3].wave().preDot(g1[ih1].wave());
	      auto vec8 = ten3[ih3].wave().preDot(rescaledMomenta()[0]+rescaledMomenta()[1]);
	      auto dot4 = vec1*vec8;
	      auto dot5 = vec4*vec8;
	      auto dot6 = vec1*g1[ih1].wave();
	      auto dot7 = vec5*rescaledMomenta()[0];
	      auto dot8 = vec5*rescaledMomenta()[1];
	      auto dot9 = vec6*rescaledMomenta()[1];
	      auto dot10 = vec2*rescaledMomenta()[0];
	      auto dot11 = vec8*rescaledMomenta()[0];
	      auto dot12 = vec8*rescaledMomenta()[1];
	      auto dot13 = vec1*vec5;
	      auto dot14 = vec1*vec6;
	      auto dot15 = vec1*vec7;
	      auto dot16 = vec2*vec5;
	      auto dot17 = vec2*vec6;
	      auto dot18 = vec3*vec5;
	      auto dot19 = vec3*vec6;
	      auto dot20 = vec4*vec5;
	      auto dot21 = vec4*vec6;
	      auto dot22 = vec1*rescaledMomenta()[0];
	      auto dot23 = vec5*g1[ih1].wave();
	      auto dot24 = vec6*g1[ih1].wave();
	      auto dot25 = vec3*vec7;
	      // diagrams 1D2
	      diag[0]=rt*mix1*((64.*a1*(dot11+dot12)*(dot10+2.*dot1*dot3)*pow<4,1>(M))/(a22*pow<4,1>(-(a12*M2)+sh))
			       +(16.*((2.-4.*a1)*dot1*dot4+dot5-2.*a1*dot5+4.*a1*dot6*(dot7+2.*dot8+dot9))*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+sh)));
	      diag[1]=rt*mix1*((32.*dot6*dot7*pow<3,1>(M))/(a12*pow<3,1>(-M2+th))
			       -(16.*(dot16-dot23*dot3)*M2)/(a12*a2*sqr(-M2+th))
			       +((32.*dot6*(dot7+2.*dot8+dot9)*pow<3,1>(M))/a22-
				 (32.*(dot7+2.*dot8+dot9)*pow<3,1>(M)*(2.*dot2*dot22+dot10*M+2.*a1*(dot1-dot2)*(dot22+dot3*M)))/(a22*(-M2+th)))/pow<3,1>(-(a12*M2)+sh)+
			       ((-32.*dot7*pow<3,1>(M)*(2.*dot2*dot22+dot10*M+2.*a1*(dot1-dot2)*(dot22+dot3*M)))/(a12*pow<3,1>(-M2+th))
				-(8.*M2*(-2.*dot18*(dot1+dot2)+4.*dot1*dot3*dot7-2.*dot10*(dot7+2.*dot8)+
					 (dot20-4.*dot6*dot7-2.*a1*((-1.+2.*a1)*dot13*(dot1-dot2)+dot20-2.*dot22*dot23+2.*a1*dot22*dot23+2.*dot6*dot8))*M+
					 4.*a1*a2*dot23*dot3*M2))/(a12*a2*sqr(-M2+th))
				-(8.*((1.+a1)*dot16+a1*(2.*dot17-2.*(dot23+dot24)*dot3-(-1.+2.*a1)*a2*dot15*M))*M2)/(a12*a22*(-M2+th)))/(-(a12*M2)+sh)
			       +((-8.*(dot16+dot17)*M2)/(a1*a22)
				 -(32.*(dot7+dot8)*pow<3,1>(M)*(2.*dot2*dot22+dot10*M+2.*a1*(dot1-dot2)*(dot22+dot3*M)))/(a1*a2*sqr(-M2+th))
				 +(8.*M2*(2.*((dot18+dot19)*dot2+dot1*(dot18+dot19-2.*dot3*(dot7+dot8))+dot10*(dot7+3.*dot8+2.*dot9))+
					  (-dot20-dot21+4.*a12*((dot13+dot14)*(dot1-dot2)+dot22*(dot23+dot24))+4.*dot6*(dot7+dot8)+
					   2.*a1*(-(dot1*(dot13+dot14))+(dot13+dot14)*dot2+dot20+dot21-2.*dot22*(dot23+dot24)+2.*dot6*(dot8+dot9)))*M-
					  4.*a1*a2*(dot23+dot24)*dot3*M2))/(a1*a22*(-M2+th)))/sqr(-(a12*M2)+sh));
	      diag[2]=rt*mix1*((-32.*dot6*dot7*pow<3,1>(M))/(a22*pow<3,1>(-M2+th))
			       -(16.*(dot16-dot23*dot3)*M2)/(a1*a22*sqr(-M2+th))
			       +((-32.*dot6*dot9*pow<3,1>(M))/a12-
				 (32.*dot9*pow<3,1>(M)*(2.*(a1*(dot1-dot2)+dot2)*dot22+
							dot10*M-2.*a2*(dot1-dot2)*dot3*M))/(a12*(-M2+th)))/pow<3,1>(-(a22*M2)+uh)+
			       ((8.*(dot17-2.*dot24*dot3)*M2)/(a12*a2)
				+(32.*dot8*pow<3,1>(M)*(2.*(a1*(dot1-dot2)+dot2)*dot22+dot10*M-2.*a2*(dot1-dot2)*dot3*M))/(a1*a2*sqr(-M2+th))
				-(8.*M2*(-2.*dot19*(dot1+dot2)+4.*dot1*dot3*dot8-2.*dot10*(dot8+2.*dot9)+
					 (2.*(-1.+2.*a1)*a2*dot1*dot14-2.*(-1.+2.*a1)*a2*dot14*dot2+dot21+2.*dot22*dot24-
					  2.*dot6*dot8+4.*dot6*dot9-2.*a1*(dot21+2.*a1*dot22*dot24+2.*dot6*dot9))*M+4.*a1*a2*dot24*dot3*M2))/(a12*a2*(-M2+th)))/sqr(-(a22*M2)+uh)+
			       ((-32.*dot7*pow<3,1>(M)*(2.*(a1*(dot1-dot2)+dot2)*dot22+dot10*M-2.*a2*(dot1-dot2)*dot3*M))/(a22*pow<3,1>(-M2+th))
				+(8.*M2*(-2.*dot18*(dot1+dot2)+4.*dot1*dot3*dot7-2.*dot10*(dot7+2.*dot8)+
					 (2.*(-1.+2.*a1)*a2*dot1*dot13-2.*(-1.+2.*a1)*a2*dot13*dot2+dot20+2.*dot22*dot23-
					  2.*dot6*dot7+4.*dot6*dot8-2.*a1*(dot20+2.*a1*dot22*dot23+2.*dot6*dot8))*M+4.*a1*a2*dot23*dot3*M2))/(a1*a22*sqr(-M2+th))
				-(8.*(-2.*dot17+2.*dot24*dot3-3.*a12*dot15*M+2.*pow(a1,3)*dot15*M+a1*(dot16+2.*dot17-2.*(dot23+dot24)*dot3+dot15*M))*M2)/(a12*a22*(-M2+th)))/(-(a22*M2)+uh));
	      diag[3]=rt*mix1*((64.*a2*(dot10+2.*dot2*dot3)*dot9*pow<4,1>(M))/(a12*pow<4,1>(-(a22*M2)+uh))
			       -(16.*((-1.+2.*a1)*(2.*dot14*dot2+dot21+2.*dot22*dot24-2.*dot6*dot8)+4.*a2*dot6*dot9)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+uh)));
	      diag[4]=rt*mix1*((32.*dot6*dot9*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+uh))
			       -(8.*(dot17-2.*dot24*dot3)*M2)/(a12*a2*sqr(-(a22*M2)+uh))
			       +((-32.*dot6*(dot7+2.*dot8+dot9)*pow<3,1>(M))/a22-
				 (64.*(dot7+2.*dot8+dot9)*pow<3,1>(M)*(dot2*dot22-dot1*dot3*M+a1*(dot1-dot2)*(dot22+dot3*M)))/(a22*(-(a22*M2)+uh)))/pow<3,1>(-(a12*M2)+sh)+
			       ((-64.*dot9*pow<3,1>(M)*(dot2*dot22-dot1*dot3*M+a1*(dot1-dot2)*(dot22+dot3*M)))/(a12*pow<3,1>(-(a22*M2)+uh))
				-(16.*M2*(-(dot19*dot2)-dot10*(dot8+2.*dot9)+
					  a1*((-1.+2.*a1)*dot14*dot2+2.*a2*dot22*dot24)*M+
					  2.*dot6*(dot9-a1*(dot8+2.*dot9))*M-
					  dot1*(dot19-2.*dot3*dot8-(-1.+2.*a1)*a2*dot14*M)+2.*a1*a2*dot24*dot3*M2))/(a12*a2*sqr(-(a22*M2)+uh))
				-(8.*(-dot17+a1*(dot16+2.*dot17-2.*(dot23+dot24)*dot3-(-1.+2.*a1)*a2*dot15*M))*M2)/(a12*a22*(-(a22*M2)+uh)))/(-(a12*M2)+sh)
			       +((8.*(dot16+dot17)*M2)/(a1*a22)
				 -(64.*(dot8+dot9)*pow<3,1>(M)*(dot2*dot22-dot1*dot3*M+a1*(dot1-dot2)*(dot22+dot3*M)))/(a1*a2*sqr(-(a22*M2)+uh))
				 +(16.*M2*((dot18+dot19)*dot2+dot10*(dot7+3.*dot8+2.*dot9)+
					   (2.*a12*(-((dot13+dot14)*dot2)+dot22*(dot23+dot24))-2.*dot6*(dot8+dot9)+
					    a1*((dot13+dot14)*dot2-2.*dot22*(dot23+dot24)+2.*dot6*(dot7+3.*dot8+2.*dot9)))*M+
					   dot1*(dot18+dot19-2.*dot3*(dot7+dot8)-(-1.+2.*a1)*a2*(dot13+dot14)*M)-
					   2.*a1*a2*(dot23+dot24)*dot3*M2))/(a1*a22*(-(a22*M2)+uh)))/sqr(-(a12*M2)+sh));
	      // diagrams 3D2
	      diag[0]+=mix2*((128.*a1*pow<4,1>(M)*(-((-1.+2.*a1)*(dot7+2.*dot8+dot9)*(dot10+2.*dot1*dot3))+
						   2.*a1*a2*(2.*dot1*(dot13+dot14)+dot20+dot21)*M))/(a22*pow<4,1>(-(a12*M2)+sh))
			     -(16.*pow<3,1>(M)*(3.*(2.*dot1*(dot13+dot14)+dot20+dot21)-8.*a12*(dot16+dot17)*M+
						4.*a1*(-2.*dot6*(dot7+2.*dot8+dot9)+2.*dot1*(dot13+dot14)+
						       dot20+dot21+2.*(dot16+dot17)*M)))/(a22*pow<3,1>(-(a12*M2)+sh)));
	      diag[1]+=mix2*((64.*(dot7*dot6+2.*a1*(dot1-dot2)*dot13)*pow<3,1>(M))/(a12*pow<3,1>(-M2+th))
			     -(16.*(-(dot23*dot3)+dot16-2.*a1*a2*dot15*M)*M2)/(a12*a2*sqr(-M2+th))
			     +((-64.*pow<3,1>(M)*(dot7*(2.*dot22*dot2+dot10*M)+2.*a1*(dot1-dot2)*(2.*dot8*dot22-dot18*M+dot7*(dot22+dot3*M))))/(a12*pow<3,1>(-M2+th))
			       -(8.*M2*(-2.*dot10*(dot7+2.*dot8)-2.*dot18*(dot1+dot2)+4.*dot1*dot7*dot3+
					(-8.*(1.-a1*a2)*dot7*dot6+dot20+
					 2.*a1*(4.*a2*dot6*dot8+4.*dot22*dot24-9.*dot1*dot13+5.*dot2*dot13-3.*dot20+
						2.*a1*(-2.*dot22*dot24-2.*dot1*dot14+2.*dot2*(dot13+dot14)+dot20)))*M+
					4.*a1*(-(a2*dot25)+2.*a1*(dot23*dot3-dot16))*M2))/(a12*a2*sqr(-M2+th))
			       -(8.*(dot16+a1*(2.*(dot23+dot24)*dot3-3.*dot16-2.*dot17-5.*a2*dot15*M))*M2)/(a12*a22*(-M2+th)))/(-(a12*M2)+sh)
			     +((8.*(dot16+dot17)*M2)/(a1*a22)
			       -(64.*pow<3,1>(M)*(2.*dot22*(a1*(dot7+2.*dot8+dot9)*(dot1-dot2)+dot7*dot2)+
						  ((-1.+3.*a1)*dot10*dot8-dot1*dot18+
						   a1*(-(dot1*dot18)+a1*dot1*dot18-dot1*dot19+4.*dot18*dot2-a1*dot18*dot2+dot19*dot2+
						       2.*a2*dot8*(dot1-dot2)*dot3)+dot7*((1.+a1)*dot10+2.*(dot1-a12*dot1+(-2.+a1)*a1*dot2)*dot3))*M+
						  a1*a2*(4.*dot7*dot6-4.*dot22*dot23+2.*a1*(dot1-dot2)*dot13-3.*dot20)*M2))/(a1*a2*sqr(-M2+th))
			       +(8.*M2*(-2.*(dot10*(dot7+3.*dot8+2.*dot9)+(dot18+dot19)*(dot1+dot2))+4.*dot1*(dot7+dot8)
					*dot3+(8.*dot22*dot23+5.*dot20+dot21+
					       2.*a1*(-8.*dot6*(dot8+dot9)+4.*a1*dot6*(dot7+2.*dot8+dot9)-4.*dot22*dot23+(7.*dot1-3.*dot2)
						      *(dot13+dot14)-dot20+dot21-2.*a1*(2.*dot2*(dot13+dot14)+dot20+dot21)))*M-
					8.*a12*a22*dot15*pow<3,1>(M)+4.*a1*(a22*dot25+2.*dot24*dot3-2.*a12*(dot23+dot24)*dot3+dot16-dot17
									    +a1*(-2.*(dot23+2.*dot24)*dot3+3.*dot16+5.*dot17))*M2))/(a1*a22*(-M2+th)))/sqr(-(a12*M2)+sh)+
			     ((-64.*pow<3,1>(M)*(dot6*(dot7+2.*dot8+dot9)-2.*dot22*(dot23+dot24)-dot20-dot21+2.*a1*((dot23+dot24)*dot3-dot16-dot17)*M))/a22+
			      (64.*pow<3,1>(M)*(-2.*a1*(dot7+2.*dot8+dot9)*dot22*(dot1-dot2)+2.*(-dot7+dot9)*dot22*dot2-
						(dot10*(-2.*dot8+8.*a1*dot8-3.*dot9+6.*a1*dot9)-2.*dot1*dot18+
						 6.*a1*(dot18+dot19)*dot2+2.*a1*(dot1*dot9-(4.*dot8+dot9)*dot2)*dot3-2.*dot1*(dot19-2.*dot8*dot3)+
						 2.*a12*(dot1-dot2)*(dot18+dot19-2.*(2.*dot8+dot9)*dot3)+
						 dot7*(dot10+2.*a1*dot10+4.*dot1*dot3-2.*a1*(dot1+2.*a1*dot1+(3.-2.*a1)*dot2)*dot3))*M-
						2.*a1*a2*(4.*dot6*(dot7+dot8)-4.*dot22*(dot23+dot24)+2.*a1*(dot1-dot2)*(dot13+dot14)-
							  3.*(dot20+dot21))*M2))/(a22*(-M2+th)))/pow<3,1>(-(a12*M2)+sh));
	      diag[2]+=mix2*((64.*(dot7*dot6+2.*a2*(dot1-dot2)*dot13)*pow<3,1>(M))/(a22*pow<3,1>(-M2+th))
			     +(16.*((-3.+2.*a1)*(dot23*dot3-dot16)+2.*a1*a2*dot15*M)*M2)/(a1*a22*sqr(-M2+th))
			     +((64.*(-(dot10*dot9)-dot19*dot2+2.*dot8*(dot2*dot3-a2*dot6*M)+
				     M*(-(dot2*dot14)+dot21-a1*(dot22*dot24+dot21)-a2*(dot24*dot3-dot17)*M))*M2)/a12+
			       (64.*pow<3,1>(M)*(2.*dot1*((-2.+a1)*dot8-a2*dot9)*dot22-2.*a1*(dot8+dot9)*dot22*dot2+
						 (dot10*dot9+2.*dot1*dot19-2.*a1*(dot10*(dot8+dot9)+dot19*dot2)+
						  2.*a1*(dot8+3.*dot9)*(-dot1+dot2)*dot3-2.*(dot1*dot8-dot1*dot9+dot9*dot2)*dot3-
						  2.*a12*(dot1-dot2)*(dot19-2.*(dot8+dot9)*dot3))*M-
						 2.*a1*a2*(-2.*a2*(dot1-dot2)*dot14-dot21)*M2))/(a12*(-M2+th)))/pow<3,1>(-(a22*M2)+uh)+
			     ((-8.*(-2.*a2*dot25+4.*a2*dot23*dot3-2.*dot24*dot3+4.*a1*dot24*dot3+5.*dot17-6.*a1*dot17-2.*a22*dot15*M)*M2)/(a12*a2)
			      -(64.*pow<3,1>(M)*(dot1*((-2.+a1)*dot7+dot8+2.*a2*dot9)*dot22-(dot8+a1*(dot7-2.*dot9)+2.*dot9)*dot22*dot2+
						 (dot1*dot18+dot19*dot2-dot1*(dot19+dot7*dot3)-a12*(dot1-dot2)*(dot18-2.*(dot7+dot8)*dot3)+
						  dot8*(dot10+2.*(dot1-dot2)*dot3)+
						  a1*(-(dot10*(dot7+dot8))+dot1*dot19-(dot18+dot19)*dot2-(dot7+4.*dot8)*(dot1-dot2)*dot3))*M-
						 a1*a2*(-2.*a2*(dot1-dot2)*dot13-dot20)*M2))/(a1*a2*sqr(-M2+th))
			      -(8.*M2*(2.*(-(dot10*(dot8+2.*a1*dot8-2.*dot9+4.*a1*dot9))+dot2*(dot19-2.*a1*(dot18+dot19)+4.*a1*dot7*dot3-8.*a2*dot8*dot3)+
					   dot1*(-2.*a2*dot18+dot19-2.*a1*dot19+4.*dot7*dot3+6.*dot8*dot3-4.*a1*(dot7+dot8)*dot3))+
				       ((22.-16.*a1)*dot6*dot8-8.*(-1.+a12)*dot6*dot9-6.*dot22*dot24-4.*dot1*dot13+6.*dot1*dot14-6.*dot2*dot14+
					4.*a12*(dot22*dot23-dot1*dot13+dot2*dot13+2.*dot1*dot14-dot21)-3.*dot21+
					a1*(4.*dot22*(-2.*dot23+dot24)+8.*dot1*dot13-8.*dot2*dot13-14.*dot1*dot14+6.*dot2*dot14+6.*dot21))*M+
				       8.*a12*a22*dot15*pow<3,1>(M)+4.*a2*(2.*dot24*dot3+a12*(dot25-2.*(dot23+dot24)*dot3)-2.*dot17+
									   a1*(-(dot23*dot3)+2.*dot24*dot3+dot16+dot17))*M2))/(a12*a2*(-M2+th)))/sqr(-(a22*M2)+uh)+
			     ((64.*pow<3,1>(M)*(-2.*dot1*((-2.+a1)*dot7-2.*a2*dot8)*dot22-2.*a2*(dot7+2.*dot8)*dot22*dot2+
						(-2.*a2*dot18*(dot1-dot2)+dot7*(dot10+2.*a2*(dot1-dot2)*dot3))*M))/(a22*pow<3,1>(-M2+th))
			      +(8.*M2*(-2.*(-3.+2.*a1)*(dot7*dot10+2.*dot10*dot8+dot18*(dot1+dot2))-4.*dot7*(dot1+2.*a2*dot2)*dot3-
				       (2.*(-7.+4.*a1)*dot7*dot6+6.*dot22*dot23-2.*(dot1-dot2)*(3.*dot13-4.*dot14)+5.*dot20+
					4.*a12*(2.*(dot6*dot8+dot22*(dot23+dot24)+dot1*dot14-dot2*(dot13+dot14))+dot20)-
					2.*a1*(4.*dot6*dot8+4.*dot22*(2.*dot23+dot24)-(3.*dot1+dot2)*dot13+8.*(dot1-dot2)*dot14+5.*dot20))*M-
				       4.*a2*(-(dot23*dot3)+a1*(dot25-2.*dot23*dot3)+dot16)*M2))/(a1*a22*sqr(-M2+th))
			      +(8.*(2.*(-2.+a1)*dot23*dot3-2.*(-1.+2.*a1)*a2*dot24*dot3+2.*(dot16-dot17+dot15*M)+
				    a1*(dot16-2.*a1*dot16+6.*dot17-4.*a1*dot17+dot15*M-3.*a1*dot15*M))*M2)/(a12*a22*(-M2+th)))/(-(a22*M2)+uh));
	      diag[3]+=mix2*((-128.*a2*pow<4,1>(M)*((-1.+2.*a1)*dot9*(dot10+2.*dot2*dot3)+
						    2.*a1*a2*(2.*dot6*dot8-2.*dot22*dot24-2.*dot2*dot14-dot21)*M))/(a12*pow<4,1>(-(a22*M2)+uh))
			     +(16.*pow<3,1>(M)*(8.*a2*dot6*dot9+
						(7.-4.*a1)*(2.*dot6*dot8-2.*dot22*dot24-2.*dot2*dot14-dot21)+
						8.*a1*a2*(2.*dot24*dot3-dot17)*M))/(a12*pow<3,1>(-(a22*M2)+uh)));
	      diag[4]+=mix2*((-64.*(dot6*(2.*dot8+dot9)-2.*dot22*dot24-2.*dot2*dot14-dot21)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+uh))
			     -(8.*(-2.*dot24*dot3+dot17)*M2)/(a12*a2*sqr(-(a22*M2)+uh))
			     +((128.*pow<3,1>(M)*(dot9*dot22*(a1*dot1-(1.+a1)*dot2)-(-(a2*dot1)-a1*dot2)*(dot9*dot3+a1*(dot19-2.*(dot8+dot9)*dot3))*M-
						  2.*a1*a2*(-(a2*dot1)-a1*dot2)*dot14*M2))/(a12*pow<3,1>(-(a22*M2)+uh))
			       -(16.*M2*(-(dot10*(dot8+2.*dot9))-dot1*dot19-dot19*dot2+2.*dot1*dot8*dot3-
					 (-4.*dot6*dot9+dot1*dot14-2.*a12*(-2.*dot6*(dot8+dot9)+2.*dot22*dot24+2.*dot1*dot14+dot21)+
					  a1*(-4.*dot6*(dot7+3.*dot8+dot9)+4.*dot22*dot23+8.*dot22*dot24+4.*dot2*dot13+3.*dot1*dot14+
					      3.*dot2*dot14+2.*dot20+4.*dot21))*M+4.*a12*a22*dot15*pow<3,1>(M)-
					 2.*a1*a2*(-(a1*dot25)-2.*dot24*dot3+2.*a1*(dot23+dot24)*dot3+dot17)*M2))/(a12*a2*sqr(-(a22*M2)+uh))
			       -(8.*(-dot17+a1*(2.*(dot23+dot24)*dot3-dot16-5.*a2*dot15*M))*M2)/(a12*a22*(-(a22*M2)+uh)))/(-(a12*M2)+sh)
			     +((64.*(-(dot6*(dot7+2.*dot8+dot9))+2.*dot1*(dot13+dot14)+dot20+dot21)*pow<3,1>(M))/a22+
			       (128.*pow<3,1>(M)*(-((dot7+2.*dot8+dot9)*dot22*((-2.+a1)*dot1+dot2-a1*dot2))-
						  (-(a2*dot1)-a1*dot2)*(-(a2*(dot18+dot19))-(-1.+2.*a1)*(dot7+2.*dot8+dot9)*dot3)*M-
						  2.*a1*a2*(-(a2*dot1)-a1*dot2)*(dot13+dot14)*M2))/(a22*(-(a22*M2)+uh)))/pow<3,1>(-(a12*M2)+sh)+
			     ((-8.*(dot16+dot17)*M2)/(a1*a22)
			      -(64.*pow<3,1>(M)*(-2.*(dot8+dot9)*dot22*(dot1-dot2)+
						 (-(a2*dot1)-a1*dot2)*(-dot19+2.*(dot8+dot9)*dot3+a1*(dot18+2.*dot19-2.*(dot7+3.*dot8+2.*dot9)*dot3))*M+
						 2.*a1*a2*(-(a2*dot1)-a1*dot2)*(dot13+2.*dot14)*M2))/(a1*a2*sqr(-(a22*M2)+uh))
			      +(16.*M2*(-(dot10*(dot7+3.*dot8+2.*dot9))-(dot18+dot19)*(dot1+dot2)+2.*dot1*(dot7+dot8)*dot3+
					(-4.*dot6*(dot8+dot9)-dot1*dot13+3.*dot1*dot14+
					 a1*(-4.*dot6*(2.*dot7+3.*dot8+dot9)+dot1*(dot13-3.*dot14)+
					     5.*dot2*(dot13+dot14)+2.*dot20)+2.*dot21+
					 a12*(4.*dot6*(dot7+2.*dot8+dot9)-2.*(2.*dot2*(dot13+dot14)+dot20+dot21)))*M-
					4.*a12*a22*dot15*pow<3,1>(M)-2.*a1*a2*(-(a2*dot25)+2.*a2*(dot23+dot24)*dot3-dot16-dot17)*M2))
			      /(a1*a22*(-(a22*M2)+uh)))/sqr(-(a12*M2)+sh));
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
	  for(unsigned int ih3=0;ih3<5;++ih3) {
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      auto dot1 = rescaledMomenta()[1]*g1[ih1].wave();
	      auto dot2 = rescaledMomenta()[3]*g1[ih1].wave();
	      complex<Energy> dot3=v4[ih4].dimensionedWave().pseudoScalar(vbar2[ih2].dimensionedWave());
	      auto vec1 = v4[ih4].dimensionedWave().generalCurrent(vbar2[ih2].dimensionedWave(),-1,1);
	      auto vec2 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),1,-1);
	      auto vec3 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[0]).generalCurrent(vbar2[ih2].dimensionedWave(),1,-1);
	      auto vec4 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[0]).generalCurrent(vbar2[ih2].dimensionedWave(),-1,1);
	      auto vec5 = ten3[ih3].wave().preDot(rescaledMomenta()[0]);
	      auto vec6 = ten3[ih3].wave().preDot(rescaledMomenta()[1]);
	      auto vec7 = ten3[ih3].wave().preDot(g1[ih1].wave());
	      auto dot4 = vec1*vec5;
	      auto dot5 = vec1*vec6;
	      auto dot6 = vec4*vec5;
	      auto dot7 = vec4*vec6;
	      auto dot8 = vec5*rescaledMomenta()[0];
	      auto dot9 = vec1*g1[ih1].wave();
	      auto dot10 = vec2*rescaledMomenta()[0];
	      auto dot11 = vec5*rescaledMomenta()[1];
	      auto dot12 = vec1*rescaledMomenta()[0];
	      auto dot13 = vec5*g1[ih1].wave();
	      auto dot14 = vec6*rescaledMomenta()[1];
	      auto dot15 = vec6*g1[ih1].wave();
	      auto dot16 = vec1*vec7;
	      auto dot17 = vec2*vec5;
	      auto dot18 = vec2*vec6;
	      auto dot19 = vec3*vec5;
	      auto dot20 = vec3*vec6;
	      auto dot21 = vec3*vec7;
	      // diagrams 1D2
	      diag[0]=rt*mix1*((-64.*a2*(dot10-2.*dot1*dot3)*(2.*dot11+dot14+dot8)*pow<4,1>(M))/(a12*pow<4,1>(-(a22*M2)+sh))
			       +(16.*(-((-1.+2.*a1)*(2.*dot12*(dot13+dot15)-2.*dot1*(dot4+dot5)+dot6+dot7))+
				      2.*((3.-2.*a1)*dot11+2.*a2*dot14+dot8)*dot9)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+sh)));
	      diag[1]=rt*mix1*((-32.*dot8*dot9*pow<3,1>(M))/(a12*pow<3,1>(-M2+th))
			       +(16.*(dot17-dot13*dot3)*M2)/(a12*a2*sqr(-M2+th))
			       +((-32.*dot14*dot9*pow<3,1>(M))/a22+
				 (32.*dot14*pow<3,1>(M)*(dot10*M-2.*a1*dot2*(dot12+dot3*M)+2.*dot1*(-(a2*dot12)+a1*dot3*M)))/(a22*(-M2+th)))/pow<3,1>(-(a12*M2)+uh)+
			       ((-8.*dot18*M2)/(a1*a22)
				-(32.*dot11*pow<3,1>(M)*(dot10*M-2.*a1*dot2*(dot12+dot3*M)+2.*dot1*(-(a2*dot12)+a1*dot3*M)))/(a1*a2*sqr(-M2+th))
				+(8.*M2*(-2.*dot10*(dot11+2.*dot14)-2.*(dot1+dot2)*dot20+
					 4.*dot11*dot2*dot3+
					 (-dot7+4.*dot11*dot9+2.*a1*(-2.*a2*dot12*dot15+(-1.+2.*a1)*(dot1-dot2)*dot5+dot7-2.*(dot11+dot14)*dot9))*M-
					 4.*a1*a2*dot15*dot3*M2))/(a1*a22*(-M2+th)))/sqr(-(a12*M2)+uh)+
			       ((32.*dot8*pow<3,1>(M)*(dot10*M-2.*a1*dot2*(dot12+dot3*M)+2.*dot1*(-(a2*dot12)+a1*dot3*M)))/(a12*pow<3,1>(-M2+th))
				-(8.*M2*(-2.*dot19*(dot1+dot2)+4.*dot2*dot3*dot8-2.*dot10*(2.*dot11+dot8)-4.*a1*a2*dot12*dot13*M+
					 (-1.+2.*a1)*(2.*a1*(dot1-dot2)*dot4+dot6)*M+4.*(dot8-a1*(dot11+dot8))*dot9*M-4.*a1*a2*dot13*dot3*M2))/(a12*a2*sqr(-M2+th))
				+(8.*(dot17-a1*dot17+a1*(-2.*dot18+2.*dot15*dot3-(-1.+2.*a1)*a2*dot16*M))*M2)/(a12*a22*(-M2+th)))/(-(a12*M2)+uh));
	      diag[2]=rt*mix1*((32.*dot8*dot9*pow<3,1>(M))/(a22*pow<3,1>(-M2+th))
			       +(16.*(dot17-dot13*dot3)*M2)/(a1*a22*sqr(-M2+th))
			       +((32.*(2.*dot11+dot14+dot8)*dot9*pow<3,1>(M))/a12+
				 (32.*(2.*dot11+dot14+dot8)*pow<3,1>(M)*((dot10+2.*dot2*dot3)*M-
									 2.*a2*dot1*(dot12+dot3*M)-2.*a1*dot2*(dot12+dot3*M)))/(a12*(-M2+th)))/pow<3,1>(-(a22*M2)+sh)+
			       ((32.*dot8*pow<3,1>(M)*((dot10+2.*dot2*dot3)*M-2.*a2*dot1*(dot12+dot3*M)-2.*a1*dot2*(dot12+dot3*M)))/(a22*pow<3,1>(-M2+th))
				+(8.*M2*(-2.*dot19*(dot1+dot2)+4.*dot2*dot3*dot8-2.*dot10*(2.*dot11+dot8)+
					 (2.*(-1.+2.*a12)*dot12*dot13+(-1.+2.*a1)*(-2.*a2*(dot1-dot2)*dot4+dot6)+
					  2.*(2.*dot11-2.*a1*dot11+3.*dot8-2.*a1*dot8)*dot9)*M-4.*a1*a2*dot13*dot3*M2))/(a1*a22*sqr(-M2+th))
				+(8.*(-((-2.+a1)*dot17)+2.*a2*dot18-2.*(dot13+dot15-a1*dot15)*dot3-a1*(-1.+2.*a1)*a2*dot16*M)*M2)/(a12*a22*(-M2+th)))/(-(a22*M2)+sh)
			       +((8.*(dot17+dot18-2.*(dot13+dot15)*dot3)*M2)/(a12*a2)
				 +(32.*(dot11+dot8)*pow<3,1>(M)*((dot10+2.*dot2*dot3)*M-2.*a2*dot1*(dot12+dot3*M)-2.*a1*dot2*(dot12+dot3*M)))/(a1*a2*sqr(-M2+th))
				 -(8.*M2*(2.*(dot1+dot2)*(dot19+dot20)-4.*dot2*dot3*(dot11+dot8)+
					  2.*dot10*(3.*dot11+2.*dot14+dot8)-2.*(-1.+2.*a12)*dot12*(dot13+dot15)*M
					  -(-1.+2.*a1)*(-2.*a2*(dot1-dot2)*(dot4+dot5)+dot6+dot7)*M+
					  2.*((-5.+4.*a1)*dot11-2.*dot14-3.*dot8+2.*a1*(dot14+dot8))*dot9*M+
					  4.*a1*a2*(dot13+dot15)*dot3*M2))/(a12*a2*(-M2+th)))/sqr(-(a22*M2)+sh));
	      diag[3]=rt*mix1*((-64.*a1*dot14*(dot10-2.*dot2*dot3)*pow<4,1>(M))/(a22*pow<4,1>(-(a12*M2)+uh))
			       +(16.*(2.*(-1.+2.*a1)*dot2*dot5+dot7-2.*a1*dot7-4.*a1*dot14*dot9)*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+uh)));
	      diag[4]=rt*mix1*((32.*dot14*dot9*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+uh))
			       +(8.*dot18*M2)/(a1*a22*sqr(-(a12*M2)+uh))
			       +((-32.*(2.*dot11+dot14+dot8)*dot9*pow<3,1>(M))/a12+
				 (64.*(2.*dot11+dot14+dot8)*pow<3,1>(M)*(-(a2*dot1*dot12)+a1*dot1*dot3*M+dot2*dot3*M-a1*dot2*(dot12+dot3*M)))/(a12*(-(a12*M2)+uh)))/pow<3,1>(-(a22*M2)+sh)+
			       ((64.*dot14*pow<3,1>(M)*(-(a2*dot1*dot12)+a1*dot1*dot3*M+dot2*dot3*M-a1*dot2*(dot12+dot3*M)))/(a22*pow<3,1>(-(a12*M2)+uh))
				+(16.*M2*(-(dot10*(dot11+2.*dot14))-dot1*dot20-dot2*dot20+2.*dot11*dot2*dot3-2.*a1*a2*dot12*dot15*M+
					  (-1.+2.*a1)*(a1*(dot1-dot2)+dot2)*dot5*M+2.*(dot11-a1*dot11+dot14-2.*a1*dot14)*dot9*M-2.*a1*a2*dot15*dot3*M2))/(a1*a22*sqr(-(a12*M2)+uh))
				+(8.*(dot17-a1*dot17+dot18-2.*a1*dot18+2.*a1*dot15*dot3-a1*(-1.+2.*a1)*a2*dot16*M)*M2)/(a12*a22*(-(a12*M2)+uh)))/(-(a22*M2)+sh)
			       +((-8.*(dot17+dot18-2.*(dot13+dot15)*dot3)*M2)/(a12*a2)
				 +(64.*(dot11+dot14)*pow<3,1>(M)*(-(a2*dot1*dot12)+a1*dot1*dot3*M+dot2*dot3*M-a1*dot2*(dot12+dot3*M)))/(a1*a2*sqr(-(a12*M2)+uh))
				 -(16.*M2*(dot1*dot19+dot19*dot2+dot1*dot20+dot2*dot20-2.*dot11*dot2*dot3-2.*dot2*dot3*dot8+
					   dot10*(3.*dot11+2.*dot14+dot8)+2.*a1*a2*dot12*(dot13+dot15)*M-
					   (-1.+2.*a1)*(a1*(dot1-dot2)+dot2)*(dot4+dot5)*M+2.*((-2.+3.*a1)*dot11+(-1.+2.*a1)*dot14-a2*dot8)*dot9*M+
					   2.*a1*a2*(dot13+dot15)*dot3*M2))/(a12*a2*(-(a12*M2)+uh)))/sqr(-(a22*M2)+sh));
	      // diagrams 3D2
	      diag[0]+=mix2*((-128.*a2*pow<4,1>(M)*(-((-1.+2.*a1)*(dot8+2.*dot11+dot14)*(dot10-2.*dot1*dot3))+
						    2.*a1*a2*(2.*dot9*(dot8+dot11)-2.*dot12*(dot13+dot15)+2.*dot1*(dot4+dot5)-dot6-dot7)*M))/(a12*pow<4,1>(-(a22*M2)+sh))
			     +(16.*pow<3,1>(M)*(6.*dot8*dot9+dot9*((-2.+8.*a1)*dot11-8.*a2*dot14)+
						(-7.+4.*a1)*(2.*dot12*(dot13+dot15)-2.*dot1*(dot4+dot5)+dot6+dot7)+
						8.*a1*a2*(2.*(dot13+dot15)*dot3-dot17-dot18)*M))/(a12*pow<3,1>(-(a22*M2)+sh)));
	      diag[1]+=mix2*((-64.*(dot8*dot9+2.*a1*(dot1-dot2)*dot4)*pow<3,1>(M))/(a12*pow<3,1>(-M2+th))
			     +(16.*(-(dot13*dot3)+dot17-2.*a1*a2*dot16*M)*M2)/(a12*a2*sqr(-M2+th))
			     +((64.*pow<3,1>(M)*(dot9*dot14+2.*dot12*dot15+dot7+2.*a1*(-(dot15*dot3)+dot18)*M))/a22+
			       (64.*pow<3,1>(M)*(2.*dot1*dot12*(2.*dot11+dot14+a1*dot14)-2.*a1*dot12*dot14*dot2+
						 (dot10*(-4.*a2*dot11-3.*dot14+6.*a1*dot14)+6.*a1*dot1*dot20-
						  2.*a1*(2.*dot1*dot11-dot1*dot14+2.*dot11*dot2+dot14*dot2)*dot3-2.*dot2*(dot20-2.*dot11*dot3)-
						  2.*a12*(dot1-dot2)*(dot20+2.*dot14*dot3))*M-
						 2.*a1*a2*(4.*dot9*dot11-4.*dot12*dot15+2.*a1*(dot1-dot2)*dot5-3.*dot7)*M2))/(a22*(-M2+th)))/pow<3,1>(-(a12*M2)+uh)+
			     ((8.*dot18*M2)/(a1*a22)
			      -(64.*pow<3,1>(M)*(2.*dot12*(dot1*(dot8-a1*dot14)+a1*dot14*dot2)+
						 ((-1.+3.*a1)*dot10*dot11+3.*a1*dot1*dot19-a12*dot1*dot19-dot19*dot2+a12*dot19*dot2-a1*dot1*dot20+
						  a1*dot2*dot20+2.*a1*a2*dot11*(dot1-dot2)*dot3-2.*dot8*(dot10-a1*dot10+a1*dot1*dot3-a2*dot2*dot3))*M-
						 a1*a2*(4.*dot8*dot9-4.*dot13*dot12+2.*a1*(dot1-dot2)*dot4-3.*dot6)*M2))/(a1*a2*sqr(-M2+th))
			      +(8.*M2*(2.*dot10*(dot11+2.*dot14)+2.*(dot1+dot2)*dot20-4.*dot11*dot2*dot3+
				       (-8.*dot13*dot12-4.*dot6+dot7+
					2.*a1*(8.*dot9*dot11+4.*dot13*dot12-4.*(-2.+a1)*dot9*dot14+(3.+4.*a1)*dot1*dot5-
					       7.*dot2*dot5+2.*dot6+dot7-2.*a1*dot7))*M+
				       8.*a12*a22*dot16*pow<3,1>(M)-
				       4.*a1*(a22*dot21-2.*a2*dot13*dot3-2.*dot15*dot3+4.*a1*dot15*dot3+
					      2.*a12*dot15*dot3+2.*dot17-2.*a1*dot17+dot18-5.*a1*dot18)*M2))/(a1*a22*(-M2+th)))/sqr(-(a12*M2)+uh)+
			     ((64.*pow<3,1>(M)*(-2.*dot1*(dot8+a1*dot8+2.*a1*dot11)*dot12+2.*a1*(dot8+2.*dot11)*dot12*dot2+
						2.*a1*dot19*(-dot1+dot2)*M+dot8*(dot10+2.*a1*(dot1-dot2)*dot3)*M))/(a12*pow<3,1>(-M2+th))
			      -(8.*M2*(-2.*dot10*(dot8+2.*dot11)-2.*dot19*(dot1+dot2)+4.*dot8*dot2*dot3+
				       (8.*dot8*dot9-dot6-4.*a12*(2.*(dot9*dot11+dot12*(dot13+dot15)+dot1*dot5-dot2*(dot4+dot5))+dot6)+
					2.*a1*(4.*dot9*dot11+4.*dot12*(dot13+dot15)+5.*dot1*dot4-9.*dot2*dot4+3.*dot6))*M+
				       4.*a1*(dot21-a1*dot21+2.*a1*(-(dot13*dot3)+dot17))*M2))/(a12*a2*sqr(-M2+th))
			      +(8.*(dot17+a1*(-2.*dot15*dot3-dot17+2.*dot18-5.*a2*dot16*M))*M2)/(a12*a22*(-M2+th)))/(-(a12*M2)+uh));
	      diag[2]+=mix2*((-64.*(dot8*dot9+2.*a2*(dot1-dot2)*dot4)*pow<3,1>(M))/(a22*pow<3,1>(-M2+th))
			     +(16.*((3.-2.*a1)*dot13*dot3-3.*dot17+2.*a1*(dot17-a2*dot16*M))*M2)/(a1*a22*sqr(-M2+th))
			     +((64.*pow<3,1>(M)*(2.*(a1*dot8-2.*a2*dot11)*dot12*dot2-2.*a2*dot19*dot2*M-
						 dot8*(dot10-2.*a2*dot2*dot3)*M+2.*a2*dot1*(2.*dot11*dot12+dot19*M+dot8*(dot12-dot3*M))))/(a22*pow<3,1>(-M2+th))
			       +(8.*M2*(-2.*(-3.+2.*a1)*(dot8*dot10+2.*dot10*dot11+dot19*(dot1+dot2))+
					4.*dot8*(-2.*a2*dot1-dot2)*dot3+
					(-2.*(7.+4.*(-2.+a1)*a1)*dot8*dot9+6.*dot13*dot12-2.*(dot1-dot2)*(7.*dot4+4.*dot5)+
					 2.*a1*(4.*dot9*dot11-4.*dot13*dot12+4.*dot12*dot15+7.*dot1*dot4-11.*dot2*dot4+
						8.*dot1*dot5-8.*dot2*dot5-5.*dot6)+
					 5.*dot6+4.*a12*(-2.*dot9*dot11-2.*dot12*dot15-2.*dot1*dot5+2.*dot2*(dot4+dot5)+dot6))*M+
					4.*a2*(-(dot13*dot3)+a1*(dot21-2.*dot13*dot3)+dot17)*M2))/(a1*a22*sqr(-M2+th))
			       +(8.*(2.*(3.+2.*(-2.+a1)*a1)*dot13*dot3-2.*(-1.+2.*a1)*a2*dot15*dot3-
				     2.*(2.*dot17+dot18+dot16*M)+a1*(5.*dot17-2.*a1*dot17+6.*dot18-4.*a1*dot18+(-1.+3.*a1)*dot16*M))*M2)/(a12*a22*(-M2+th)))/(-(a22*M2)+sh)
			     +((-8.*(2.*a2*dot21+2.*(-3.+4.*a1)*dot13*dot3-2.*dot15*dot3+4.*a1*dot15*dot3+5.*dot17-6.*a1*dot17+5.*dot18-6.*a1*dot18+2.*a22*dot16*M)*M2)/(a12*a2)
			       -(64.*pow<3,1>(M)*(dot1*dot12*(-(a2*dot8)+(-3.+4.*a1)*dot11-2.*a2*dot14)-dot12*(dot8+a1*dot8-3.*dot11+4.*a1*dot11-2.*a2*dot14)*dot2+
						  (-(a2*(-(dot10*dot11)+dot1*(dot19+a1*dot19+dot20)-dot2*((2.+a1)*dot19+dot20)-
							 2.*a2*dot11*(dot1-dot2)*dot3))+dot8*(dot10+(2.-3.*a1)*dot1*dot3-3.*a2*dot2*dot3))*M+
						  a1*a2*(-2.*a2*(dot1-dot2)*dot4-dot6)*M2))/(a1*a2*sqr(-M2+th))
			       -(8.*M2*(2.*dot10*(-5.*dot11+6.*a1*dot11-2.*dot14+4.*a1*dot14)+8.*a1*dot19*dot2+4.*a1*dot1*dot20-
					2.*dot2*dot20+4.*a1*dot2*dot20-2.*(3.*dot19*dot2+dot1*(dot19+dot20))+
					4.*dot11*(4.*a2*dot1+(-3.+2.*a1)*dot2)*dot3+
					(-3.*(2.*dot13*dot12-2.*dot1*dot4+dot6+dot7)+
					 2.*((3.-8.*a1*a2)*dot9*dot11+4.*(-1.+a12)*dot9*dot14-3.*dot12*dot15-5.*dot2*dot4+
					     3.*dot1*dot5-3.*dot2*dot5-2.*a12*(dot13*dot12-dot1*dot4+3.*dot2*dot4+2.*dot2*dot5+dot6+dot7)+
					     a1*(2.*dot12*(3.*dot13+dot15)-7.*dot1*dot4+11.*dot2*dot4-3.*dot1*dot5+7.*dot2*dot5+3.*(dot6+dot7))))*M-
					8.*a12*a22*dot16*pow<3,1>(M)+
					2.*dot8*((-3.+2.*a1)*dot10-2.*(2.*(-2.+a1)*dot1+dot2)*dot3+(7.+4.*(-2.+a1)*a1)*dot9*M)-
					4.*a2*(a12*(dot21+2.*dot15*dot3)-a1*(3.*dot13*dot3+2.*dot15*dot3+dot18)+
					       2.*(-((dot13+dot15)*dot3)+dot17+dot18))*M2))/(a12*a2*(-M2+th)))/sqr(-(a22*M2)+sh)+
			     ((64.*M2*(2.*dot10*dot11+dot10*dot14+dot1*(dot19+dot20)-2.*dot1*dot11*dot3+
				       (-2.*a2*dot9*dot11+dot1*(dot4+dot5)+dot6+dot7-a1*(dot12*(dot13+dot15)+dot6+dot7))*M+
				       dot8*(dot10-2.*dot1*dot3-2.*a2*dot9*M)-a2*((dot13+dot15)*dot3-dot17-dot18)*M2))/a12+
			      (64.*pow<3,1>(M)*(-2.*a1*dot12*(dot11+dot14)*(dot1-dot2)+2.*dot12*(dot8-dot14)*dot2-
						(dot10*(2.*dot11-2.*a1*dot11+dot14-2.*a1*dot14)-2.*a2*(a1*dot1-(1.+a1)*dot2)*(dot19+dot20)+
						 2.*((-1.+2.*a1)*dot1*((-2.+a1)*dot11-a2*dot14)+a2*(-3.*dot11+2.*a1*dot11-dot14+2.*a1*dot14)*dot2)*dot3+
						 dot8*(dot10+2.*(dot1-2.*a1*dot1-2.*a2*dot2)*dot3))*M-
						2.*a1*a2*(-2.*a2*(dot1-dot2)*(dot4+dot5)-dot6-dot7)*M2))/(a12*(-M2+th)))/pow<3,1>(-(a22*M2)+sh));
	      diag[3]+=mix2*((128.*a1*pow<4,1>(M)*((-1.+2.*a1)*dot14*(dot10-2.*dot2*dot3)-2.*a1*a2*(2.*dot2*dot5-dot7)*M))/(a22*pow<4,1>(-(a12*M2)+uh))
			     -(16.*pow<3,1>(M)*(-6.*dot2*dot5+3.*dot7+4.*a1*(2.*dot9*dot14-2.*dot2*dot5+dot7+2.*a2*dot18*M)))/(a22*pow<3,1>(-(a12*M2)+uh)));
	      diag[4]+=mix2*((64.*(dot9*dot14-2.*dot2*dot5+dot7)*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+uh))
			     -(8.*dot18*M2)/(a1*a22*sqr(-(a12*M2)+uh))
			     +((128.*pow<3,1>(M)*(dot12*dot14*(-(a2*dot1)-(-2.+a1)*dot2)-
						  (a1*(dot1-dot2)+dot2)*(-(a2*dot20)+(-1.+2.*a1)*dot14*dot3)*M-2.*a1*a2*(a1*(dot1-dot2)+dot2)*dot5*M2))/(a22*pow<3,1>(-(a12*M2)+uh))
			       +(16.*M2*(dot10*dot11+2.*dot10*dot14+(dot1+dot2)*dot20-2.*dot11*dot2*dot3+
					 (4.*dot9*(dot14+a2*(dot11+a1*dot14))+a1*(-5.+4.*a1)*dot1*dot5-a2*dot2*(4.*dot4+3.*dot5)+2.*(dot6+dot7)-2.*a1*(dot6+a1*dot7))*M+
					 4.*a12*a22*dot16*pow<3,1>(M)+2.*a1*a2*(-(a2*(dot21+2.*dot15*dot3))+dot18)*M2))/(a1*a22*sqr(-(a12*M2)+uh))
			       +(8.*(-2.*a1*dot15*dot3+dot17-a1*dot17+dot18-5.*a1*a2*dot16*M)*M2)/(a12*a22*(-(a12*M2)+uh)))/(-(a22*M2)+sh)
			     +((-64.*(dot8*dot9-dot9*dot14-2.*dot12*(dot13+dot15)+2.*dot1*(dot4+dot5)-dot6-dot7)*pow<3,1>(M))/a12+
			       (128.*pow<3,1>(M)*(-(dot12*(dot8+2.*dot11+dot14)*(dot1+a1*dot1-a1*dot2))-
						  (a1*(dot1-dot2)+dot2)*(-((dot8+2.*dot11+dot14)*dot3)+a1*(dot19+dot20+2.*(dot11+dot14)*dot3))*M-
						  2.*a1*a2*(a1*(dot1-dot2)+dot2)*(dot4+dot5)*M2))/(a12*(-(a12*M2)+uh)))/pow<3,1>(-(a22*M2)+sh)+
			     ((-8.*(-2.*(dot13+dot15)*dot3+dot17+dot18)*M2)/(a12*a2)
			      -(64.*pow<3,1>(M)*(2.*dot12*(dot11+dot14)*(dot1-dot2)+
						 (a1*(dot1-dot2)+dot2)*(-(a2*dot19)+(-1.+2.*a1)*dot20+2.*(-(a2*dot11)+(-1.+2.*a1)*dot14)*dot3)*M+
						 2.*a1*a2*(a1*(dot1-dot2)+dot2)*(dot4+2.*dot5)*M2))/(a1*a2*sqr(-(a12*M2)+uh))
			      -(16.*M2*(3.*dot10*dot11+2.*dot10*dot14+(dot1+dot2)*(dot19+dot20)-2.*dot11*dot2*dot3+
					(4.*dot9*(-dot14-a2*((2.+a1)*dot11+a1*dot14))+dot2*(dot4+dot5)+
					 2.*a12*(2.*dot12*(dot13+dot15)-2.*dot2*(dot4+dot5)+dot6+dot7)-
					 a1*(4.*dot12*(dot13+2.*dot15)+(dot1-3.*dot2)*dot4-3.*(dot1+dot2)*dot5+2.*(dot6+2.*dot7)))*M-
					4.*a12*a22*dot16*pow<3,1>(M)+dot8*(dot10-2.*dot2*dot3-4.*a2*dot9*M)-
					2.*a1*a2*(-2.*(dot13+dot15)*dot3+a1*(dot21+2.*dot15*dot3)+dot17+dot18)*M2))/(a12*a2*(-(a12*M2)+uh)))/sqr(-(a22*M2)+sh));
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
    meSum *= 1./51840.;
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
    ProductionMatrixElement me(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin2,PDT::Spin1);
    for(unsigned int ih1=0;ih1<2;++ih1) { 
      for(unsigned int ih2=0;ih2<2;++ih2) {
    	for(unsigned int ih3=0;ih3<5;++ih3) {
    	  for(unsigned int ih4=0;ih4<2;++ih4) {
	    auto dot1 = rescaledMomenta()[0]*g4[ih4].wave();
	    auto dot2 = rescaledMomenta()[1]*g4[ih4].wave();
	    complex<Energy> dot3=u1[ih1].dimensionedWave().pseudoScalar(vbar2[ih2].dimensionedWave());
	    auto vec1 = u1[ih1].dimensionedWave().generalCurrent(vbar2[ih2].dimensionedWave(),-1,1);
	    auto vec2 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),1,-1);
	    auto vec3 = u1[ih1].dimensionedWave().slash(rescaledMomenta()[3]).generalCurrent(vbar2[ih2].dimensionedWave(),1,-1);
	    auto vec4 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).slash(rescaledMomenta()[3]).generalCurrent(vbar2[ih2].dimensionedWave(),-1,1);
	    auto vec5 = ten3[ih3].wave().preDot(rescaledMomenta()[0]);
	    auto vec6 = ten3[ih3].wave().preDot(rescaledMomenta()[1]);
	    auto vec7 = ten3[ih3].wave().preDot(g4[ih4].wave());
	    auto dot4 = vec1*vec5;
	    auto dot5 = vec1*vec6;
	    auto dot6 = vec1*vec7;
	    auto dot7 = vec2*vec5;
	    auto dot8 = vec2*vec6;
	    auto dot9 = vec3*vec5;
	    auto dot10 = vec3*vec6;
	    auto dot11 = vec5*rescaledMomenta()[0];
	    auto dot12 = vec1*rescaledMomenta()[3];
	    auto dot13 = vec1*g4[ih4].wave();
	    auto dot14 = vec2*rescaledMomenta()[3];
	    auto dot15 = vec5*rescaledMomenta()[1];
	    auto dot16 = vec5*g4[ih4].wave();
	    auto dot17 = vec6*rescaledMomenta()[1];
	    auto dot18 = vec6*g4[ih4].wave();
	    auto dot19 = vec4*vec5;
	    auto dot20 = vec4*vec6;
	    auto dot21 = vec3*vec7;
    	    // diagrams 1D2
	    diag[0]=rt*mix1*((32.*dot11*dot13*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+th))
			     -(8.*(-2.*dot16*dot3+dot7)*M2)/(a12*a2*sqr(-(a22*M2)+th))
			     +((-32.*dot13*dot17*pow<3,1>(M))/a22+
			       (64.*dot17*pow<3,1>(M)*(-(dot12*dot2)+dot1*dot3*M+a1*(dot1+dot2)*(dot12-dot3*M)))/(a22*(-(a22*M2)+th)))/pow<3,1>(-(a12*M2)+uh)+
			     ((-8.*dot8*M2)/(a1*a22)
			      -(64.*dot15*pow<3,1>(M)*(-(dot12*dot2)+dot1*dot3*M+a1*(dot1+dot2)*(dot12-dot3*M)))/(a1*a2*sqr(-(a22*M2)+th))
			      +(16.*M2*(dot14*dot15-dot14*dot17-dot10*dot2+2.*dot13*(dot15-a1*dot15+a1*dot17)*M+a1*(-2.*a2*dot12*dot18+(1.-2.*a1)*dot2*dot5)*M+
					dot1*(dot10-2.*(dot15+dot17)*dot3+(-1.+(3.-2.*a1)*a1)*dot5*M)+2.*a1*a2*dot18*dot3*M2))/(a1*a22*(-(a22*M2)+th)))/sqr(-(a12*M2)+uh)+
			     ((64.*dot11*pow<3,1>(M)*(-(dot12*dot2)+dot1*dot3*M+a1*(dot1+dot2)*(dot12-dot3*M)))/(a12*pow<3,1>(-(a22*M2)+th))
			      -(16.*M2*(-(dot14*dot15)-2.*dot1*dot15*dot3+dot1*dot9-dot2*dot9+
					2.*a1*(dot13*dot15-a2*dot12*dot16)*M-(-1.+2.*a1)*(-(a2*dot1)+a1*dot2)*dot4*M+
					dot11*(dot14-2.*dot1*dot3+2.*a2*dot13*M)+2.*a1*a2*dot16*dot3*M2))/(a12*a2*sqr(-(a22*M2)+th))
			      -(8.*(-dot7+a1*(2.*dot18*dot3+dot7-dot8-(-1.+2.*a1)*a2*dot6*M))*M2)/(a12*a22*(-(a22*M2)+th)))/(-(a12*M2)+uh));
	    diag[1]=rt*mix1*((-64.*a2*dot11*(dot14+2.*dot2*dot3)*pow<4,1>(M))/(a12*pow<4,1>(-(a22*M2)+th))
			     +(16.*(-2.*dot11*dot13+(1.-2.*a1)*(2.*dot13*dot15-2.*dot12*dot16-dot19-2.*dot2*dot4))*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+th)));
	    diag[2]=rt*mix1*((-32.*dot11*dot13*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+th))
			     -(8.*(2.*dot16*dot3-dot7)*M2)/(a12*a2*sqr(-(a22*M2)+th))
			     +((-32.*dot13*(dot11+2.*dot15+dot17)*pow<3,1>(M))/a22+
			       (32.*(dot11+2.*dot15+dot17)*pow<3,1>(M)*(-2.*dot12*dot2+2.*a1*dot12*(dot1+dot2)+dot14*M+2.*a2*(dot1+dot2)*dot3*M))/(a22*(-(a22*M2)+th)))/pow<3,1>(-M2+sh)+
			     ((32.*dot11*pow<3,1>(M)*(-2.*dot12*dot2+2.*a1*dot12*(dot1+dot2)+dot14*M+2.*a2*(dot1+dot2)*dot3*M))/(a12*pow<3,1>(-(a22*M2)+th))
			      -(8.*M2*(-2.*(dot14*dot15+2.*dot1*dot15*dot3-dot1*dot9+dot2*dot9)+
				       (2.*dot13*dot15+2.*(-1.+2.*a12)*dot12*dot16+(-1.+2.*a1)*(dot19+2.*a2*(dot1+dot2)*dot4))*M+
				       2.*dot11*(dot14-2.*dot1*dot3+(3.-2.*a1)*dot13*M)+4.*a1*a2*dot16*dot3*M2))/(a12*a2*sqr(-(a22*M2)+th))
			      -(8.*(2.*dot16*dot3-2.*dot7+a1*(2.*dot18*dot3+dot7-dot8-(-1.+2.*a1)*a2*dot6*M))*M2)/(a12*a22*(-(a22*M2)+th)))/(-M2+sh)
			     +((-16.*((dot16+dot18)*dot3-dot7-dot8)*M2)/(a1*a22)
			       +(32.*(dot11+dot15)*pow<3,1>(M)*(-2.*dot12*dot2+2.*a1*dot12*(dot1+dot2)+dot14*M+2.*a2*(dot1+dot2)*dot3*M))/(a1*a2*sqr(-(a22*M2)+th))
			       +(8.*M2*(2.*dot14*(-dot11+dot17)+2.*dot2*(dot10+dot9)+
					(2.*(-3.+2.*a1)*dot11*dot13+2.*dot12*dot16+dot19+dot20-
					 2.*(dot13*(-2.*(-2.+a1)*dot15+dot17)-dot12*dot18+2.*a12*dot12*(dot16+dot18)+a1*(dot19+dot20)))*M-
					2.*(-1.+2.*a1)*a2*dot2*(dot4+dot5)*M+2.*dot1*(-dot10+2.*(dot11+2.*dot15+dot17)*dot3-dot9-(-1.+2.*a1)*a2*(dot4+dot5)*M)-
					4.*a1*a2*(dot16+dot18)*dot3*M2))/(a1*a22*(-(a22*M2)+th)))/sqr(-M2+sh));
	    diag[3]=rt*mix1*((-64.*a1*dot17*(dot14-2.*dot1*dot3)*pow<4,1>(M))/(a22*pow<4,1>(-(a12*M2)+uh))
			     +(16.*(dot20-2.*dot1*dot5+a1*(4.*dot13*dot17-2.*dot20+4.*dot1*dot5))*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+uh)));
	    diag[4]=rt*mix1*((32.*dot13*dot17*pow<3,1>(M))/(a22*pow<3,1>(-(a12*M2)+uh))
			     +(8.*dot8*M2)/(a1*a22*sqr(-(a12*M2)+uh))
			     +((32.*dot13*(dot11+2.*dot15+dot17)*pow<3,1>(M))/a12+
			       (32.*(dot11+2.*dot15+dot17)*pow<3,1>(M)*(-2.*dot12*dot2+dot14*M+2.*a1*(dot1+dot2)*(dot12-dot3*M)))/(a12*(-(a12*M2)+uh)))/pow<3,1>(-M2+sh)+
			     ((32.*dot17*pow<3,1>(M)*(-2.*dot12*dot2+dot14*M+2.*a1*(dot1+dot2)*(dot12-dot3*M)))/(a22*pow<3,1>(-(a12*M2)+uh))
			      +(8.*M2*(2.*dot14*(dot15-dot17)-2.*dot10*dot2+2.*dot1*(dot10-2.*(dot15+dot17)*dot3)+
				       4.*dot13*(dot15-a1*dot15+dot17)*M-dot20*M+
				       2.*a1*(-2.*a2*dot12*dot18+dot20-(-1.+2.*a1)*(dot1+dot2)*dot5)*M+4.*a1*a2*dot18*dot3*M2))/(a1*a22*sqr(-(a12*M2)+uh))
			      -(8.*(-dot7-dot8+a1*(2.*dot18*dot3+dot7-dot8-(-1.+2.*a1)*a2*dot6*M))*M2)/(a12*a22*(-(a12*M2)+uh)))/(-M2+sh)
			     +((-16.*((dot16+dot18)*dot3-dot7-dot8)*M2)/(a12*a2)
			       +(32.*(dot15+dot17)*pow<3,1>(M)*(-2.*dot12*dot2+dot14*M+2.*a1*(dot1+dot2)*(dot12-dot3*M)))/(a1*a2*sqr(-(a12*M2)+uh))
			       -(8.*M2*(2.*dot14*(-dot11+dot17)+2.*dot2*(dot10+dot9)-
					2.*dot1*(dot10-2.*(dot11+2.*dot15+dot17)*dot3+dot9)+
					(-4.*a2*dot11*dot13+4.*(-2.+a1)*dot13*dot15-4.*dot13*dot17+dot19+dot20-
					 2.*a1*(-2.*a2*dot12*(dot16+dot18)+dot19+dot20)+2.*a1*(-1.+2.*a1)*(dot1+dot2)*dot4+
					 2.*a1*(-1.+2.*a1)*(dot1+dot2)*dot5)*M-4.*a1*a2*(dot16+dot18)*dot3*M2))/(a12*a2*(-(a12*M2)+uh)))/sqr(-M2+sh));
    	    // 3D2 diagrams
	    diag[0]+=mix2*((-64.*(dot19-dot11*dot13-2.*dot13*dot15+2.*dot12*dot16+2.*dot2*dot4)*pow<3,1>(M))/(a12*pow<3,1>(-(a22*M2)+th))
			   -(8.*(-2.*dot16*dot3+dot7)*M2)/(a12*a2*sqr(-(a22*M2)+th))
			   +((64.*(dot20-dot13*dot17-2.*dot1*dot5)*pow<3,1>(M))/a22+
			     (128.*pow<3,1>(M)*(dot12*((-2.+a1)*dot1-a2*dot2)*dot17-
						(-(a2*dot1)+a1*dot2)*(-(a2*dot10)+(1.-2.*a1)*dot17*dot3)*M+
						2.*a1*a2*(-(a2*dot1)+a1*dot2)*dot5*M2))/(a22*(-(a22*M2)+th)))/pow<3,1>(-(a12*M2)+uh)+
			   ((8.*dot8*M2)/(a1*a22)
			    -(64.*pow<3,1>(M)*(-2.*dot12*dot15*(dot1+dot2)+
					       (-(a2*dot1)+a1*dot2)*(-2.*dot15*dot3+a1*(dot10+2.*dot15*dot3-2.*dot17*dot3-dot9)+dot9)*M+
					       2.*a1*a2*(-(a2*dot1)+a1*dot2)*(dot4-dot5)*M2))/(a1*a2*sqr(-(a22*M2)+th))
			    +(16.*(dot10*dot2+dot14*(-dot15+dot17)-
				   2.*a2*(dot19-a1*dot20)*M+4.*dot13*(dot15-a1*dot15+(-2.+a1)*a1*dot17)*M+
				   dot1*(-dot10+2.*(dot15+dot17)*dot3+a2*(4.*dot4+dot5)*M)+
				   a1*M*((5.-4.*a1)*dot2*dot5-2.*a2*M*(dot21-a1*dot21+dot8-2.*a2*(dot18*dot3-a1*dot6*M))))*M2)/(a1*a22*(-(a22*M2)+th)))/sqr(-(a12*M2)+uh)+
			   ((128.*pow<3,1>(M)*(-(dot11*dot12*(dot2+a1*(dot1+dot2)))-
					       (-(a2*dot1)+a1*dot2)*(dot11*dot3+2.*a1*dot15*dot3-a1*dot9)*M-
					       2.*a1*a2*(-(a2*dot1)+a1*dot2)*dot4*M2))/(a12*pow<3,1>(-(a22*M2)+th))
			    -(16.*M2*(dot11*dot14-dot14*dot15-2.*dot1*dot11*dot3-2.*dot1*dot15*dot3+dot1*dot9-dot2*dot9-2.*a1*(-(a2*dot19)+dot20)*M-
				      (-4.*a2*dot11*dot13+dot1*dot4+a1*(-4.*dot13*(-(a2*dot15)+dot17)+4.*dot12*(-(a2*dot16)+dot18)+
									3.*dot1*dot4-4.*a1*dot1*dot4+dot2*dot4+4.*dot2*dot5))*M+4.*a12*a22*dot6*pow<3,1>(M)-
				      2.*a1*a2*(-2.*dot16*dot3+a1*(dot21-2.*dot18*dot3)+dot7)*M2))/(a12*a2*sqr(-(a22*M2)+th))
			    -(8.*(-dot7+a1*(-2.*dot18*dot3+dot7+dot8-5.*a2*dot6*M))*M2)/(a12*a22*(-(a22*M2)+th)))/(-(a12*M2)+uh));
	    diag[1]+=mix2*((-128.*a2*pow<4,1>(M)*(-((-1.+2.*a1)*dot11*(dot14+2.*dot2*dot3))+2.*a1*a2*(dot19-2.*dot13*(dot11+dot15)+2.*dot12*dot16+2.*dot2*dot4)*M))/(a12*pow<4,1>(-(a22*M2)+th))
			   +(16.*pow<3,1>(M)*((7.-4.*a1)*dot19-6.*dot11*dot13+2.*(-7.+4.*a1)*(dot13*dot15-dot12*dot16-dot2*dot4)+8.*a1*a2*(2.*dot16*dot3-dot7)*M))/(a12*pow<3,1>(-(a22*M2)+th)));
	    diag[2]+=mix2*((64.*M2*(dot11*dot14+2.*dot11*dot2*dot3+2.*dot15*dot2*dot3-dot2*dot9+
				    (-(a2*dot19)+2.*a2*dot13*(dot11+dot15)+a1*dot12*dot16+dot2*dot4)*M-a2*(dot16*dot3-dot7)*M2))/(a12*pow<3,1>(-(a22*M2)+th))
			   -(8.*(2.*a2*dot21+2.*(-3.+4.*a1)*dot16*dot3-4.*dot18*dot3+4.*a1*dot18*dot3+5.*dot7-6.*a1*dot7-2.*a22*dot6*M)*M2)/(a12*a2*sqr(-(a22*M2)+th))
			   +((64.*(dot13*(dot11+2.*dot15+dot17)-2.*a2*(dot1+dot2)*(dot4+dot5))*pow<3,1>(M))/a22+
			     (64.*pow<3,1>(M)*(-2.*dot12*(a1*(dot1+dot2)*(dot11-dot17)+dot2*(-dot11+dot17)+2.*dot1*(dot15+dot17))+
					       (-(dot14*(dot11+2.*dot15+dot17))-2.*a2*(dot1+dot2)*(dot10-(dot11+2.*dot15+dot17)*dot3+dot9))*M))/(a22*(-(a22*M2)+th)))/pow<3,1>(-M2+sh)+
			   ((64.*pow<3,1>(M)*(2.*dot12*(-(dot1*dot11)+(-2.+a1)*dot1*dot15+a1*dot15*dot2)-
					      (-2.*dot1*dot15*dot3+dot11*(dot14-4.*a2*dot1*dot3+2.*(-1.+2.*a1)*dot2*dot3)+
					       2.*a1*dot15*(dot14+(-1.+2.*a1)*(dot1+dot2)*dot3)+2.*a2*(dot1+a1*dot1+a1*dot2)*dot9)*M-
					      2.*a1*a2*(dot19-2.*a2*(dot1+dot2)*dot4)*M2))/(a12*pow<3,1>(-(a22*M2)+th))
			    -(8.*M2*(2.*(-3.*dot11*dot14-8.*dot11*dot2*dot3-dot15*(dot14+8.*dot2*dot3)+
					 2.*a1*(dot11*dot14-dot14*dot15+dot10*dot2+2.*dot2*(dot11-dot17)*dot3)+
					 dot2*dot9+dot1*(-2.*a2*dot10-2.*(dot11+(-1.+2.*a1)*dot15-2.*a2*dot17)*dot3+(-3.+4.*a1)*dot9))+
				     ((3.-6.*a1+4.*a12)*dot19-22.*dot13*dot15+6.*dot2*dot4+2.*(-7.*dot11*dot13+3.*dot12*dot16+5.*dot1*dot4)+
				      4.*dot1*dot5+4.*a12*(-2.*dot11*dot13+dot12*(dot16+dot18)+(3.*dot1+dot2)*dot4+(dot1+dot2)*dot5)-
				      2.*a1*(-8.*dot13*(dot11+dot15)+6.*dot12*dot16+4.*dot12*dot18+11.*dot1*dot4+7.*dot2*dot4+
					     4.*(dot1+dot2)*dot5))*M+
				     8.*a12*a22*dot6*pow<3,1>(M)-
				     4.*a2*(-2.*dot16*dot3+2.*dot7+a1*(-((3.*dot16+dot18)*dot3)+a1*(dot21-2.*dot18*dot3)+dot8))*M2))/(a12*a2*sqr(-(a22*M2)+th))
			    +(8.*(2.*(3.+2.*(-2.+a1)*a1)*dot16*dot3-2.*(-2.+a1)*dot18*dot3-2.*(2.*dot7+dot8-dot6*M)+
				  a1*((5.-2.*a1)*dot7-dot8+2.*a1*dot8+dot6*M-3.*a1*dot6*M))*M2)/(a12*a22*(-(a22*M2)+th)))/(-M2+sh)
			   +((16.*(-((-3.+2.*a1)*((dot16+dot18)*dot3-dot7-dot8))+2.*a1*a2*dot6*M)*M2)/(a1*a22)
			     -(64.*pow<3,1>(M)*(dot12*(dot1*(dot11+a1*dot11+5.*dot15-2.*a1*dot15+2.*dot17-a1*dot17)+
						       dot2*(-(a2*dot11)+dot15-a1*(2.*dot15+dot17)))+
						(dot14*(dot11+dot15)+
						 a1*(dot14*(dot15+dot17)+dot2*(dot10+(3.*dot11+2.*dot15-dot17)*dot3))+
						 dot2*(-2.*(dot11+dot15)*dot3+dot9)-a12*dot2*(dot10-2.*(dot15+dot17)*dot3+dot9)+
						 a2*dot1*((1.+a1)*dot10-(3.*dot11+2.*(2.+a1)*dot15+dot17+2.*a1*dot17)*dot3+(2.+a1)*dot9))*M+
						a1*a2*(dot19+dot20-2.*a2*(dot1+dot2)*(dot4+dot5))*M2))/(a1*a2*sqr(-(a22*M2)+th))
			     +(8.*M2*(6.*dot1*dot10-4.*a1*dot1*dot10+6.*dot11*dot14-4.*a1*dot11*dot14-6.*dot10*dot2+4.*a1*dot10*dot2-
				      6.*dot14*dot17+4.*a1*dot14*dot17-4.*(dot1-2.*a2*dot2)*(dot11+2.*dot15+dot17)*dot3-
				      2.*(-3.+2.*a1)*(dot1-dot2)*dot9-5.*dot19*M-
				      (5.*dot20-14.*dot13*(dot11+2.*dot15+dot17)+6.*dot12*(dot16+dot18)+
				       2.*(dot1+dot2)*(7.*dot4+3.*dot5)+
				       4.*a12*(dot19+dot20-2.*dot13*(dot11+dot15)+2.*dot12*dot18+2.*dot1*dot4-2.*dot2*dot5)-
				       2.*a1*(5.*dot19+5.*dot20+4.*dot12*dot16-4.*dot13*(2.*dot11+3.*dot15+dot17)+8.*dot12*dot18+
					      11.*dot1*dot4+7.*dot2*dot4+3.*dot1*dot5-dot2*dot5))*M+
				      4.*a2*(-((dot16+dot18)*dot3)+a1*(dot21-2.*(dot16+dot18)*dot3)+dot7+dot8)*M2))/(a1*a22*(-(a22*M2)+th)))/sqr(-M2+sh));
	    diag[3]+=mix2*((128.*a1*pow<4,1>(M)*((-1.+2.*a1)*dot17*(dot14-2.*dot1*dot3)+2.*a1*a2*(dot20-2.*dot1*dot5)*M))/(a22*pow<4,1>(-(a12*M2)+uh))
			   +(16.*pow<3,1>(M)*(-((3.+4.*a1)*dot20)+6.*dot1*dot5+8.*a1*(dot13*dot17+dot1*dot5+dot8*M-a1*dot8*M)))/(a22*pow<3,1>(-(a12*M2)+uh)));
	    diag[4]+=mix2*((64.*pow<3,1>(M)*(dot20-dot13*dot17+2.*dot12*dot18+2.*a1*dot18*dot3*M-2.*a1*dot8*M))/(a22*pow<3,1>(-(a12*M2)+uh))
			   -(8.*dot8*M2)/(a1*a22*sqr(-(a12*M2)+uh))
			   +((64.*(dot13*(dot11+2.*dot15+dot17)-2.*a1*(dot1+dot2)*(dot4+dot5))*pow<3,1>(M))/a12-
			     (64.*pow<3,1>(M)*((dot11+2.*dot15+dot17)*(2.*dot12*dot2-dot14*M)+
					       2.*a1*(dot1+dot2)*(dot12*(dot11-dot17)+
								  (-dot10+(dot11+2.*dot15+dot17)*dot3-dot9)*M)))/(a12*(-(a12*M2)+uh)))/pow<3,1>(-M2+sh)+
			   ((64.*pow<3,1>(M)*(2.*dot12*(-2.*dot15*dot2-dot2*dot17+a1*(dot1+dot2)*dot17)+
					      (dot14*(4.*a2*dot15+dot17+2.*a1*dot17)-
					       2.*a1*dot2*((-3.+a1)*dot10+2.*dot15*dot3+(3.-2.*a1)*dot17*dot3)+
					       2.*dot1*(dot10-a12*dot10-2.*a2*dot15*dot3+(-2.+a1+2.*a12)*dot17*dot3))*M+
					      2.*a1*a2*(3.*dot20-4.*dot13*(dot15+dot17)+4.*dot12*dot18+2.*a1*(dot1+dot2)*dot5)*M2))/(a22*pow<3,1>(-(a12*M2)+uh))
			    +(8.*M2*(2.*(-(dot1*dot10)+dot10*dot2+dot14*(-dot15+dot17)+2.*dot1*(dot15+dot17)*dot3)-
				     (-4.*a2*dot19-5.*dot20-8.*dot12*dot16-8.*dot12*dot18+
				      4.*a12*(dot20-2.*dot13*dot17+2.*dot2*dot5)+2.*a1*(dot20-8.*dot13*dot15+4.*dot12*(dot16+dot18)+(7.*dot1+3.*dot2)*dot5))*M-
				     8.*a12*a22*dot6*pow<3,1>(M)-4.*a1*(a22*dot21-2.*a2*dot16*dot3+2.*dot7-2.*a1*((1.+a1)*dot18*dot3+dot7)+dot8+3.*a1*dot8)*M2))/(a1*a22*sqr(-(a12*M2)+uh))
			    +(8.*(dot7+dot8-a1*(-2.*dot18*dot3+dot7+3.*dot8-5.*a2*dot6*M))*M2)/(a12*a22*(-(a12*M2)+uh)))/(-M2+sh)
			   +((-16.*((dot16+dot18)*dot3-dot7-dot8-2.*a1*a2*dot6*M)*M2)/(a12*a2)
			     -(64.*pow<3,1>(M)*(2.*dot12*(dot11*dot2+2.*dot15*dot2+(dot2-a1*(dot1+dot2))*dot17)+
						(-2.*dot11*dot14-dot14*(3.*dot15+dot17)+
						 a1*(dot14*(dot15-dot17)+2.*dot11*(dot14+dot2*dot3)+
						     dot2*(-4.*dot10+6.*dot15*dot3+4.*dot17*dot3-3.*dot9))+
						 a12*dot2*(dot10-2.*(dot15+dot17)*dot3+dot9)+
						 dot1*((-1.-a1*a2)*dot10+2.*a2*(dot11+(2.+a1)*dot15+dot17+a1*dot17)*dot3+(-1.+a12)*dot9))*M-
						a1*a2*(3.*dot19+3.*dot20-4.*dot13*(dot11+2.*dot15+dot17)+
						       4.*dot12*(dot16+dot18)+2.*a1*(dot1+dot2)*(dot4+dot5))*M2))/(a1*a2*sqr(-(a12*M2)+uh))
			     -(8.*M2*(2.*(-(dot11*dot14)+dot14*dot17+dot2*(dot10+dot9)-dot1*(dot10-2.*(dot11+2.*dot15+dot17)*dot3+dot9))+
				      ((1.-6.*a1+4.*a12)*dot19+dot20-8.*dot13*(dot11+2.*dot15)-8.*dot13*dot17+
				       4.*a12*(dot20+2.*dot12*dot16-2.*dot13*(dot15+dot17)-2.*dot1*dot4+2.*dot2*dot5)+
				       2.*a1*(-3.*dot20-4.*dot12*dot16+4.*dot13*(dot15+dot17)+
					      (9.*dot1+5.*dot2)*(dot4+dot5)))*M+
				      4.*a1*(dot21-a1*dot21+2.*a1*(-((dot16+dot18)*dot3)+dot7+dot8))*M2))
			     /(a12*a2*(-(a12*M2)+uh)))/sqr(-M2+sh));
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
    meSum *=  1./19440.;
  }
  // save the diagram weights
  meInfo(save);
  // final factors
  return meSum*O1_*pow<3,1>(Constants::pi*standardModel()->alphaS(scale())/M)/sqr(M2*a1*a2);
}

IBPtr MEPP2BCD2Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2BCD2Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2BCD2Jet::doinit() {
  setBcState(10545+mode_*10000);
  MEPP2BCJetBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<2>(bcbar,principleQuantumNumber(),0,2);
  double theta = oniumParameters()->singletTripletMixing(principleQuantumNumber(),2);
  sTheta_ = sin(theta);
  cTheta_ = cos(theta);
}

void MEPP2BCD2Jet::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2*GeV2*GeV2) << mode_ << theta_ << sTheta_ << cTheta_;
}

void MEPP2BCD2Jet::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2*GeV2*GeV2) >> mode_ >> theta_ >> sTheta_ >> cTheta_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2BCD2Jet,MEPP2BCJetBase>
describeHerwigMEPP2BCD2Jet("Herwig::MEPP2BCD2Jet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2BCD2Jet::Init() {

  static ClassDocumentation<MEPP2BCD2Jet> documentation
    ("The MEPP2BCD2Jet class implements the matrix element for g c -> B_c(D2) b q and "
     "q bar -> Bc(D2) g and charged conjugate");

  static Switch<MEPP2BCD2Jet,unsigned int> interfaceMode
    ("Mode",
     "Which b_c state to produce",
     &MEPP2BCD2Jet::mode_, 0, false, false);
  static SwitchOption interfaceModeB_c1
    (interfaceMode,
     "B_c2",
     "Produce the lighter B_c2 state",
     0);
  static SwitchOption interfaceModeB_c1Prime
    (interfaceMode,
     "B_c2Prime",
     "Produce the heavier B_c2 state",
     1);

}

