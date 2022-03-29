// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2BCP1Jet class.
//

#include "MEPP2BCP1Jet.h"
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
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/epsilon.h"
#include <numeric>

using namespace Herwig;

double MEPP2BCP1Jet::me2() const {
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
  VectorWaveFunction Bcw(rescaledMomenta()[2],mePartonData()[2],outgoing);
  vector<VectorWaveFunction> v3;
  for(unsigned int ix=0;ix<3;++ix) {
    Bcw.reset(ix);
    v3.push_back(Bcw);
  }
  // stuff for mixing
  int itest = (abs(mePartonData()[2]->id())%100000)/10000;
  double mix1 = itest==1 ? sTheta_ :  cTheta_;
  double mix2 = itest==1 ? cTheta_ : -sTheta_;
  if(mePartonData()[2]->id()<0) mix2*=-1.;
  double rt(sqrt(2.));
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
	  complex<Energy> dot6=u2[ih2].dimensionedWave().pseudoScalar(ubar4[ih4].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = u2[ih2].dimensionedWave().generalCurrent(ubar4[ih4].dimensionedWave(),-1,1);
	  complex<Energy2> dot7 = vec1*rescaledMomenta()[0];
	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    complex<Energy> dot1 = rescaledMomenta()[0]*v3[ih3].wave();
	    complex<Energy> dot3 = rescaledMomenta()[1]*v3[ih3].wave();
	    complex<Energy> dot9 = vec1*v3[ih3].wave();
	    LorentzPolarizationVectorE vec3 = u2[ih2].dimensionedWave().slash(v3[ih3].wave()).generalCurrent(ubar4[ih4].dimensionedWave(),1,-1);
	    complex<Energy2> dot12 = vec3*rescaledMomenta()[0];
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      complex<Energy> dot2 = rescaledMomenta()[1]*g1[ih1].wave();
	      complex<Energy> dot4 = rescaledMomenta()[3]*g1[ih1].wave();
	      Complex dot5 = g1[ih1].wave()*v3[ih3].wave();
	      LorentzPolarizationVectorE vec2 = u2[ih2].dimensionedWave().slash(g1[ih1].wave()).generalCurrent(ubar4[ih4].dimensionedWave(),1,-1);
	      LorentzPolarizationVectorE vec4 = u2[ih2].dimensionedWave().slash(v3[ih3].wave()).slash(g1[ih1].wave()).generalCurrent(ubar4[ih4].dimensionedWave(),-1,1);
	      complex<Energy> dot8 = vec1*g1[ih1].wave();
	      complex<Energy2> dot10 = vec2*rescaledMomenta()[0];
	      complex<Energy2> dot11 = vec4*rescaledMomenta()[0];
	      complex<Energy> dot13 = vec3*g1[ih1].wave();
	      // diagrams 1P1
	      diag[0]=rt*mix1*((-16.*a1*(dot1+dot3)*(dot10+2.*dot2*dot6)*pow<3,1>(M))/(a2*pow<3,1>(-(a12*M2)+sh))
			       -(4.*(dot11-2.*a1*dot11+2.*((-1.+2.*a1)*dot5*dot7+dot1*dot8+2.*a1*dot3*dot8+dot2*dot9-2.*a1*dot2*dot9))*M2)/(a2*sqr(-(a12*M2)+sh)));
	      diag[1]=rt*mix1*((-8.*dot1*dot8*M2)/(a1*sqr(-M2+th))
			       -(4.*(dot13-dot5*dot6)*M)/(a1*a2*(-M2+th))
			       +((-8.*(dot1+dot3)*dot8*M2)/a2+
				 (8.*(dot1+dot3)*(2.*dot4*dot7+dot10*M+2.*a1*(dot2-dot4)*(dot7+dot6*M))*M2)/(a2*(-M2+th)))/sqr(-(a12*M2)+sh)+
			       ((-2.*(dot13-2.*dot5*dot6)*M)/(a1*a2)
				+(8.*dot1*(2.*dot4*dot7+dot10*M+2.*a1*(dot2-dot4)*(dot7+dot6*M))*M2)/(a1*sqr(-M2+th))
				-(2.*M*(4.*dot10*dot3-2.*dot12*(dot2+dot4)+2.*dot1*(dot10+2.*dot4*dot6)+(-1.+2.*a1)*dot11*M+
					2.*((1.+2.*(-2.+a1)*a1)*dot5*dot7+dot1*dot8+2.*a1*(dot1+dot3)*dot8+a1*(-1.+2.*a1)*(dot2-dot4)*dot9)*M-
					4.*a1*a2*dot5*dot6*M2))/(a1*a2*(-M2+th)))/(-(a12*M2)+sh));
	      diag[2]=rt*mix1*((-8.*dot1*dot8*M2)/(a2*sqr(-M2+th))
			       +(4.*(dot13-dot5*dot6)*M)/(a1*a2*(-M2+th))
			       +((8.*dot3*dot8*M2)/a1+
				 (8.*dot3*(2.*(a1*(dot2-dot4)+dot4)*dot7+(dot10-2.*a2*(dot2-dot4)*dot6)*M)*M2)/(a1*(-M2+th)))/sqr(-(a22*M2)+uh)+
			       ((2.*dot13*M)/(a1-a12)
				-(8.*dot1*(2.*(a1*(dot2-dot4)+dot4)*dot7+(dot10-2.*a2*(dot2-dot4)*dot6)*M)*M2)/(a2*sqr(-M2+th))
				-(2.*M*(4.*dot10*dot3-2.*dot12*(dot2+dot4)+(-1.+2.*a1)*dot11*M-
					4.*a2*(a1*dot5*dot7+dot3*dot8)*M-2.*a2*M*((-1.+2.*a1)*(dot2-dot4)*dot9+2.*a1*dot5*dot6*M)+
					2.*dot1*(dot10+2.*dot4*dot6+2.*a1*dot8*M)))/(a1*a2*(-M2+th)))/(-(a22*M2)+uh));
	      diag[3]=rt*mix1*((-16.*a2*dot3*(dot10+2.*dot4*dot6)*pow<3,1>(M))/(a1*pow<3,1>(-(a22*M2)+uh))
			       -(4.*(-4.*a2*dot3*dot8-(-1.+2.*a1)*(dot11+2.*dot4*dot9))*M2)/(a1*sqr(-(a22*M2)+uh)));
	      diag[4]=rt*mix1*((-8.*dot3*dot8*M2)/(a1*sqr(-(a22*M2)+uh))
			       -(2.*dot13*M)/(a1*a2*(-(a22*M2)+uh))
			       +((8.*(dot1+dot3)*dot8*M2)/a2+
				 (16.*(dot1+dot3)*(dot4*dot7-dot2*dot6*M+a1*(dot2-dot4)*(dot7+dot6*M))*M2)/(a2*(-(a22*M2)+uh)))/sqr(-(a12*M2)+sh)+
			       ((2.*(dot13-2.*dot5*dot6)*M)/(a1*a2)
				+(16.*dot3*(dot4*dot7-dot2*dot6*M+a1*(dot2-dot4)*(dot7+dot6*M))*M2)/(a1*sqr(-(a22*M2)+uh))
				-(4.*M*(2.*dot10*dot3-dot12*(dot2+dot4)-2.*a1*a2*dot5*dot7*M+
					(-1.+2.*a1)*(2.*dot3*dot8-a2*dot2*dot9-a1*dot4*dot9)*M+dot1*(dot10+2.*dot4*dot6+2.*a1*dot8*M)-
					2.*a1*a2*dot5*dot6*M2))/(a1*a2*(-(a22*M2)+uh)))/(-(a12*M2)+sh));
	      // diagrams 3P1
	      diag[0]+=mix2*((-16.*a1*pow<3,1>(M)*((-1.+2.*a1)*(dot1+dot3)*(dot10+2.*dot2*dot6)-
						  2.*a1*a2*(dot11-2.*dot5*dot7+2.*dot1*dot8+2.*dot2*dot9)*M))/(a2*pow<3,1>(-(a12*M2)+sh))
			     -(8.*((1.+a1)*dot11+2.*(-(dot5*dot7)-a1*dot5*dot7+dot1*dot8-a1*dot3*dot8+
						      dot2*dot9+a1*dot2*dot9-a1*a2*(dot13-2.*dot5*dot6)*M))*M2)/(a2*sqr(-(a12*M2)+sh)));
	      diag[1]+=mix2*((8.*(dot1*dot8+2.*a1*(dot2-dot4)*dot9)*M2)/(a1*sqr(-M2+th))
			     +((-8.*(dot1*(2.*dot4*dot7+dot10*M)+2.*a1*(dot2-dot4)*(2.*dot3*dot7+dot12*M+dot1*(dot7-dot6*M)))*M2)/(a1*sqr(-M2+th))
			       -(8.*(-(dot1*dot8)+a1*(-(a2*dot11)+dot5*dot7+3.*dot3*dot8-3.*dot2*dot9+dot4*dot9-dot5*dot6*M+
							a1*(-2.*dot3*dot8+(dot2+dot4)*dot9+2.*dot13*M-dot5*(dot7+dot6*M))))*M2)/(a1*a2*(-M2+th)))/(-(a12*M2)+sh)
			     +((-8.*(-dot11-dot1*dot8+dot3*dot8+2.*a1*(dot13-dot5*dot6)*M)*M2)/a2-
			       (8.*M2*(2.*(a1*(dot1+dot3)*(dot2-dot4)+(dot1-dot3)*dot4)*dot7+
				       (2.*dot12*(dot2-a12*dot2+(-3.+a1)*a1*dot4)+dot1*(dot10+2.*a1*dot10-2.*a1*(dot2-3.*dot4)*dot6)+
					(-1.+2.*a1)*dot3*(3.*dot10+2.*a1*(-dot2+dot4)*dot6))*M-
				       2.*a1*a2*(3.*dot11-2.*dot5*dot7+2.*dot1*dot8+2.*a1*(-dot2+dot4)*dot9)*M2))/(a2*(-M2+th)))/sqr(-(a12*M2)+sh));
	      diag[2]+=mix2*((8.*(-((dot1*dot8)/a2)+2.*(-dot2+dot4)*dot9)*M2)/sqr(-M2+th)+
			     (8.*(dot13-dot5*dot6)*M)/(a1*(-M2+th))
			     +((8.*M*(-(dot10*dot3)+dot12*dot4+M*(dot11-a1*dot11-2.*dot5*dot7+a1*dot5*dot7-dot4*dot9-a2*(dot13-dot5*dot6)*M)))/a1+
			       (8.*M2*(2.*((-2.+a1)*dot1*dot2-a2*dot2*dot3-a1*(dot1+dot3)*dot4)*dot7+
				       (-2.*dot12*dot2+dot10*dot3+2.*dot2*(dot1+dot3)*dot6-2.*dot3*dot4*dot6+
					2.*a12*(dot2-dot4)*(dot12+2.*dot3*dot6)-
					2.*a1*(dot10*(dot1+dot3)-dot12*dot4+3.*dot3*(dot2-dot4)*dot6+dot1*(dot2+dot4)*dot6))*M+
				       2.*a1*a2*(dot11-2.*dot5*dot7+2.*dot1*dot8+2.*a2*(dot2-dot4)*dot9)*M2))/(a1*(-M2+th)))/sqr(-(a22*M2)+uh)+
			     ((8.*(dot13-dot5*dot6)*M)/a1+
			      (8.*(-2.*(dot2-dot4)*(2.*dot3*dot7+dot12*M)-(dot1*(-2.*((-2.+a1)*dot2+dot4-a1*dot4)*dot7+dot10*M-2.*a2*(dot2-dot4)*dot6*M))/a2)*M2)/sqr(-M2+th)+
			      (8.*M*(-2.*dot10*dot3+dot12*(dot2+dot4)+dot1*(-dot10-2.*dot2*dot6+(-2.*a1-1/a2)*dot8*M)+
				     M*(dot11-a1*dot11-(1.+2.*a1)*dot3*dot8-dot2*dot9+(dot4+a1*(dot2+dot4))*dot9-dot13*M-a2*dot5*(dot7-dot6*M))))/(a1*(-M2+th)))/(-(a22*M2)+uh));
	      diag[3]+=mix2*((-16.*a2*pow<3,1>(M)*((-1.+2.*a1)*dot3*(dot10+2.*dot4*dot6)-2.*a1*a2*(dot11+2.*dot4*dot9)*M))/(a1*pow<3,1>(-(a22*M2)+uh))
			     -(8.*(-2.*a2*dot3*dot8-(-2.+a1)*(dot11+2.*dot4*dot9)-2.*a1*a2*dot13*M)*M2)/(a1*sqr(-(a22*M2)+uh)));
	      diag[4]+=mix2*((8.*(dot11-dot3*dot8+2.*dot4*dot9)*M2)/(a1*sqr(-(a22*M2)+uh))
			     +((16.*M2*(dot3*(a1*dot2-(1.+a1)*dot4)*dot7+(-(a2*dot2)-a1*dot4)*(a1*dot12+(-1.+2.*a1)*dot3*dot6)*M-
					2.*a1*a2*(-(a2*dot2)-a1*dot4)*dot9*M2))/(a1*sqr(-(a22*M2)+uh))
			       -(8.*(dot3*dot8+a12*(dot11-dot5*dot7-2.*dot3*dot8+(dot2+dot4)*dot9-dot13*M+dot5*dot6*M)-
				     a1*(dot11-dot5*dot7-(dot1+2.*dot3)*dot8+(dot2+dot4)*dot9-dot13*M+dot5*dot6*M))*M2)/(a1*a2*(-(a22*M2)+uh)))/(-(a12*M2)+sh)
			     +((8.*(dot11-2.*dot5*dot7+dot1*dot8-dot3*dot8+2.*dot2*dot9)*M2)/a2-
			       (16.*M2*((dot1+dot3)*((-2.+a1)*dot2+dot4-a1*dot4)*dot7-
					(-(a2*dot2)-a1*dot4)*(-(a2*dot12)+(dot1+(-1.+2.*a1)*dot3)*dot6)*M+
					2.*a1*a2*(-(a2*dot2)-a1*dot4)*dot9*M2))/(a2*(-(a22*M2)+uh)))/sqr(-(a12*M2)+sh));
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
	  complex<Energy> dot6=v4[ih4].dimensionedWave().pseudoScalar(vbar2[ih2].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().generalCurrent(vbar2[ih2].dimensionedWave(),-1,1);
	  complex<Energy2> dot11 = vec1*rescaledMomenta()[0];
	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    complex<Energy> dot1 = rescaledMomenta()[0]*v3[ih3].wave();
	    complex<Energy> dot3 = rescaledMomenta()[1]*v3[ih3].wave();
	    complex<Energy> dot8 = vec1*v3[ih3].wave();
	    LorentzPolarizationVectorE vec3 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),1,-1);
	    complex<Energy2> dot12 = vec3*rescaledMomenta()[0];
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      complex<Energy> dot2 = rescaledMomenta()[1]*g1[ih1].wave();
	      complex<Energy> dot4 = rescaledMomenta()[3]*g1[ih1].wave();
	      Complex dot5 = g1[ih1].wave()*v3[ih3].wave();
	      LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),1,-1);
	      LorentzPolarizationVectorE vec4 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g1[ih1].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),-1,1);
	      complex<Energy> dot7 = vec1*g1[ih1].wave();
	      complex<Energy2> dot9 = vec2*rescaledMomenta()[0];
	      complex<Energy2> dot10 = vec4*rescaledMomenta()[0];
	      complex<Energy> dot13 = vec3*g1[ih1].wave();
	      // diagrams 1P1
	      diag[0]=rt*mix1*((16.*a2*(dot1+dot3)*(2.*dot2*dot6-dot9)*pow<3,1>(M))/(a1*pow<3,1>(-(a22*M2)+sh))
			       -(4.*((-1.+2.*a1)*dot10-4.*a2*(dot1+dot3)*dot7+2.*(1.-2.*a1)*dot2*dot8)*M2)/(a1*sqr(-(a22*M2)+sh)));
	      diag[1]=rt*mix1*((8.*dot1*dot7*M2)/(a1*sqr(-M2+th))
			       +(4.*(dot13-dot5*dot6)*M)/(a1*a2*(-M2+th))
			       +((-8.*dot3*dot7*M2)/a2+
				 (8.*dot3*(-2.*a2*dot11*dot2-2.*a1*dot11*dot4+(2.*a1*(dot2-dot4)*dot6+dot9)*M)*M2)/(a2*(-M2+th)))/sqr(-(a12*M2)+uh)+
			       ((2.*(dot13-2.*dot5*dot6)*M)/(a1*a2)
				-(8.*dot1*(-2.*a2*dot11*dot2-2.*a1*dot11*dot4+(2.*a1*(dot2-dot4)*dot6+dot9)*M)*M2)/(a1*sqr(-M2+th))
				-(2.*M*(-2.*dot12*(dot2+dot4)+4.*dot3*dot9+2.*dot1*(2.*dot2*dot6+dot9)+
					(dot10-2.*a1*dot10-2.*(1.+2.*(-2.+a1)*a1)*dot11*dot5-2.*dot1*dot7+4.*a1*dot3*dot7-2.*a1*(-1.+2.*a1)*(dot2-dot4)*dot8)*M+
					4.*a1*a2*dot5*dot6*M2))/(a1*a2*(-M2+th)))/(-(a12*M2)+uh));
	      diag[2]=rt*mix1*((8.*dot1*dot7*M2)/(a2*sqr(-M2+th))
			       -(4.*(dot13-dot5*dot6)*M)/(a1*a2*(-M2+th))
			       +((8.*(dot1+dot3)*dot7*M2)/a1+
				 (8.*(dot1+dot3)*(-2.*a2*dot11*dot2-2.*a1*dot11*dot4+(-2.*a2*(dot2-dot4)*dot6+dot9)*M)*M2)/(a1*(-M2+th)))/sqr(-(a22*M2)+sh)+
			       ((-2.*dot13*M)/(a1*a2)
				+(8.*dot1*(-2.*a2*dot11*dot2-2.*a1*dot11*dot4+(-2.*a2*(dot2-dot4)*dot6+dot9)*M)*M2)/(a2*sqr(-M2+th))
				-(2.*M*(-2.*dot12*(dot2+dot4)+4.*dot3*dot9+2.*dot1*(2.*dot2*dot6+dot9)+
					(dot10-2.*a1*dot10+4.*a1*a2*dot11*dot5-4.*(dot1+dot3-a1*dot3)*dot7+
					 2.*(-1.+2.*a1)*a2*(dot2-dot4)*dot8)*M+4.*a1*a2*dot5*dot6*M2))/(a1*a2*(-M2+th)))/(-(a22*M2)+sh));
	      diag[3]=rt*mix1*((-16.*a1*dot3*(-2.*dot4*dot6+dot9)*pow<3,1>(M))/(a2*pow<3,1>(-(a12*M2)+uh))
			       -(4.*((-1.+2.*a1)*dot10+2.*dot11*(dot5-2.*a1*dot5)-
				     2.*dot1*dot7+4.*a1*(dot1+dot3)*dot7+2.*dot4*dot8-4.*a1*dot4*dot8)*M2)/(a2*sqr(-(a12*M2)+uh)));
	      diag[4]=rt*mix1*((8.*dot3*dot7*M2)/(a2*sqr(-(a12*M2)+uh))
			       -(2.*(dot13-2.*dot5*dot6)*M)/(a1*a2*(-(a12*M2)+uh))
			       +((-8.*(dot1+dot3)*dot7*M2)/a1+
				 (16.*(dot1+dot3)*(-(a2*dot11*dot2)-a1*dot11*dot4+(a1*(dot2-dot4)+dot4)*dot6*M)*M2)/(a1*(-(a12*M2)+uh)))/sqr(-(a22*M2)+sh)+
			       ((2.*dot13*M)/(a1-a12)
				+(16.*dot3*(-(a2*dot11*dot2)-a1*dot11*dot4+(a1*(dot2-dot4)+dot4)*dot6*M)*M2)/(a2*sqr(-(a12*M2)+uh))
				-(4.*M*(-(dot12*(dot2+dot4))+2.*dot3*dot9+2.*a1*a2*dot11*dot5*M+
					(-1.+2.*a1)*(2.*dot3*dot7-(a1*dot2+dot4-a1*dot4)*dot8)*M+
					dot1*(2.*dot2*dot6+dot9-2.*a2*dot7*M)+2.*a1*a2*dot5*dot6*M2))/(a1*a2*(-(a12*M2)+uh)))/(-(a22*M2)+sh));
	      // diagrams 3P1
	      diag[0]+=mix2*((-16.*a2*pow<3,1>(M)*((-1.+2.*a1)*(dot1+dot3)*(dot9-2.*dot2*dot6)+2.*a1*a2*(dot10-2.*dot2*dot8)*M))/(a1*pow<3,1>(-(a22*M2)+sh))
			     -(8.*(-2.*dot10+a1*dot10-2.*a2*dot1*dot7+4.*dot2*dot8-2.*a1*dot2*dot8-2.*a2*(dot3*dot7-a1*dot13*M))*M2)/(a1*sqr(-(a22*M2)+sh)));
	      diag[1]+=mix2*((-8.*(dot1*dot7+2.*a1*(dot2-dot4)*dot8)*M2)/(a1*sqr(-M2+th))
			     +((-8.*(dot10+2.*dot1*dot7+dot3*dot7+2.*a1*dot5*dot6*M-2.*a1*dot13*M)*M2)/a2-
			       (8.*M2*(4.*dot1*dot11*dot2+2.*dot11*dot3*(dot2+a1*dot2-a1*dot4)+
				       (-4.*a2*dot1*dot9-3.*dot9*dot3+2.*dot12*dot4-4.*a1*dot1*((-2.+a1)*dot2+dot4-a1*dot4)*dot6+
					2.*a12*(dot2-dot4)*(dot12-2.*dot3*dot6)+a1*(-6.*dot12*dot2+6.*dot9*dot3+2.*dot3*(dot2-dot4)*dot6))*M-
				       2.*a1*a2*(-3.*dot10+2.*dot11*dot5-2.*dot1*dot7+2.*a1*(dot2-dot4)*dot8)*M2))/(a2*(-M2+th)))/sqr(-(a12*M2)+uh)+
			     ((8.*(-2.*dot1*dot11*(dot2+a1*dot2-a1*dot4)+dot1*(dot9+2.*a1*(-dot2+dot4)*dot6)*M+2.*a1*(dot2-dot4)*(-2.*dot11*dot3+dot12*M))*M2)/(a1*sqr(-M2+th))
			      -(8.*(dot1*dot7+a1*(dot10-dot11*dot5+3.*(dot1+dot3)*dot7+(dot2-3.*dot4)*dot8+dot5*dot6*M)+
				    a12*(-dot10+dot11*dot5-2.*(dot1+dot3)*dot7+(dot2+dot4)*dot8+dot5*dot6*M-2.*dot13*M))*M2)/(a1*a2*(-M2+th)))/(-(a12*M2)+uh));
	      diag[2]+=mix2*((8.*((dot1*dot7)/a2+2.*(dot2-dot4)*dot8)*M2)/sqr(-M2+th)+(8.*(dot5*dot6-dot13)*M)/(a1*(-M2+th))
			     +((8.*(dot5*dot6-dot13)*M)/a1+
			       (8.*(2.*dot1*dot11*(-(a2*dot2)-a1*dot4)+dot1*(dot9-2.*a2*(dot2-dot4)*dot6)*M-2.*a2*(dot2-dot4)*(2.*dot11*dot3-dot12*M))*M2)/(a2*sqr(-M2+th))
			       +(8.*M*(-((dot1*(dot9+2.*dot4*dot6))/a1)+(dot1*dot7*M)/a2+
				       (-2.*dot9*dot3+dot12*(dot2+dot4)-a2*(dot10-dot11*dot5)*M-(1.+2.*a1)*dot3*dot7*M+
					M*((1.+a1)*dot2*dot8-a2*dot4*dot8+(-(a2*dot5*dot6)+dot13)*M))/a1))/(-M2+th))/(-(a22*M2)+sh)
			     +((8.*M*(dot12*dot2-dot9*(dot1+dot3)+M*(-(a2*dot10)-(-2.+a1)*dot11*dot5-dot2*dot8-a2*(dot5*dot6-dot13)*M)))/a1+
			       (8.*M2*(2.*a1*dot11*dot2*dot3-2.*dot11*(dot1-a2*dot3)*dot4+
				       (dot9*(dot3-2.*a1*dot3)+2.*a12*dot12*dot4-2.*dot12*(-(a1*a2*dot2)+dot4)-
					2.*(-1.+2.*a1)*a2*dot3*(dot2-dot4)*dot6+dot1*(dot9+2.*(dot2+2.*(-2.+a1)*a1*dot2+2.*a1*a2*dot4)*dot6))*M-
				       2.*a1*a2*(dot10-2.*dot11*dot5+2.*dot1*dot7+2.*a2*(dot2-dot4)*dot8)*M2))/(a1*(-M2+th)))/sqr(-(a22*M2)+sh));
	      diag[3]+=mix2*((-16.*a1*pow<3,1>(M)*((-1.+2.*a1)*dot3*(dot9-2.*dot4*dot6)+2.*a1*a2*(dot10-2.*dot11*dot5+2.*dot1*dot7-2.*dot4*dot8)*M))/(a2*pow<3,1>(-(a12*M2)+uh))
			     -(8.*(-((1.+a1)*dot10)+2.*(1.+a1)*dot11*dot5-2.*(dot1+a1*dot1+a1*dot3)*dot7+2.*(1.+a1)*dot4*dot8-2.*a1*a2*(2.*dot5*dot6-dot13)*M)*M2)/(a2*sqr(-(a12*M2)+uh)));
	      diag[4]+=mix2*((-8.*(dot10-2.*dot11*dot5+2.*dot1*dot7+dot3*dot7-2.*dot4*dot8)*M2)/(a2*sqr(-(a12*M2)+uh))
			     +((-16.*(dot11*dot3*(-(a2*dot2)-(-2.+a1)*dot4)+
				      (a1*(dot2-dot4)+dot4)*M*(-(a2*dot12)+2.*dot1*dot6+dot3*dot6+2.*a12*dot8*M-2.*a1*((dot1+dot3)*dot6+dot8*M)))*M2)/(a2*sqr(-(a12*M2)+uh))
			       -(8.*((dot1+dot3)*dot7+a1*(dot10-dot11*dot5+(dot1+2.*dot3)*dot7-(dot2+dot4)*dot8+dot5*dot6*M-dot13*M)+
				     a12*(-dot10+dot11*dot5-2.*(dot1+dot3)*dot7+(dot2+dot4)*dot8-dot5*dot6*M+dot13*M))*M2)/(a1*a2*(-(a12*M2)+uh)))/(-(a22*M2)+sh)
			     +((-8.*(dot10+(dot1+dot3)*dot7-2.*dot2*dot8)*M2)/a1+
			       (16.*M2*(dot11*(dot1+dot3)*(dot2+a1*dot2-a1*dot4)-
					(a1*(dot2-dot4)+dot4)*((dot1+dot3)*dot6+a1*(dot12-2.*(dot1+dot3)*dot6))*M+
					2.*a1*a2*(a1*(dot2-dot4)+dot4)*dot8*M2))/(a1*(-(a12*M2)+uh)))/sqr(-(a22*M2)+sh));
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
    meSum *= 1./2592.;
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
    	complex<Energy> dot6=u1[ih1].dimensionedWave().pseudoScalar(vbar2[ih2].dimensionedWave());
    	LorentzPolarizationVectorE vec1 = u1[ih1].dimensionedWave().generalCurrent(vbar2[ih2].dimensionedWave(),-1,1);
    	complex<Energy2> dot7 = vec1*rescaledMomenta()[3];
    	for(unsigned int ih3=0;ih3<3;++ih3) {
    	  complex<Energy> dot2 = rescaledMomenta()[0]*v3[ih3].wave();
    	  complex<Energy> dot4 = rescaledMomenta()[1]*v3[ih3].wave();
    	  complex<Energy> dot8 = vec1*v3[ih3].wave();
    	  LorentzPolarizationVectorE vec3 = u1[ih1].dimensionedWave().slash(v3[ih3].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),1,-1);
    	  complex<Energy2> dot12 = vec3*rescaledMomenta()[3];
    	  for(unsigned int ih4=0;ih4<2;++ih4) {
    	    complex<Energy> dot1 = rescaledMomenta()[0]*g4[ih4].wave();
    	    complex<Energy> dot3 = rescaledMomenta()[1]*g4[ih4].wave();
    	    Complex         dot5 = v3[ih3].wave()*g4[ih4].wave();
    	    LorentzPolarizationVectorE vec2 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),1,-1);
    	    LorentzPolarizationVectorE vec4 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).slash(v3[ih3].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),-1,1);
    	    complex<Energy> dot9 = vec1*g4[ih4].wave();
    	    complex<Energy2> dot10 = vec2*rescaledMomenta()[3];
    	    complex<Energy> dot11 = vec2*v3[ih3].wave();
    	    complex<Energy2> dot13 = vec4*rescaledMomenta()[3];
    	    // diagrams 1P1
	    diag[0]=rt*mix1*((-8.*dot2*dot9*M2)/(a1*sqr(-(a22*M2)+th))
			     +(2.*(dot11-2.*dot5*dot6)*M)/(a1*a2*(-(a22*M2)+th))
			     +((-8.*dot4*dot9*M2)/a2-
			       (16.*dot4*(dot3*dot7-dot1*dot6*M-a1*(dot1+dot3)*(dot7-dot6*M))*M2)/(a2*(-(a22*M2)+th)))/sqr(-(a12*M2)+uh)+
			     ((-2.*dot11*M)/(a1*a2)
			      +(16.*dot2*(dot3*dot7-dot1*dot6*M-a1*(dot1+dot3)*(dot7-dot6*M))*M2)/(a1*sqr(-(a22*M2)+th))
			      -(4.*M*(dot1*dot12-dot12*dot3+dot10*(-dot2+dot4)+2.*dot2*dot3*dot6+2.*dot3*dot4*dot6+2.*a1*a2*dot5*dot7*M+
				      (-1.+2.*a1)*(-(a2*dot1)+a1*dot3)*dot8*M+
				      2.*(-(a2*dot2)-a1*dot4)*dot9*M-2.*a1*a2*dot5*dot6*M2))/(a1*a2*(-(a22*M2)+th)))/(-(a12*M2)+uh));
	    diag[1]=rt*mix1*((16.*a2*dot2*(dot10+2.*dot3*dot6)*pow<3,1>(M))/(a1*pow<3,1>(-(a22*M2)+th))
			     +(4.*((-1.+2.*a1)*(dot13-2.*(dot5*dot7+dot3*dot8))+4.*a2*dot2*dot9)*M2)/(a1*sqr(-(a22*M2)+th)));
	    diag[2]=rt*mix1*((8.*dot2*dot9*M2)/(a1*sqr(-(a22*M2)+th))
			     -(2.*(dot11-2.*dot5*dot6)*M)/(a1*a2*(-(a22*M2)+th))
			     +((8.*(dot2+dot4)*dot9*M2)/a2-
			       (8.*(dot2+dot4)*(-2.*dot3*dot7+2.*a1*(dot1+dot3)*dot7+dot10*M+2.*a2*(dot1+dot3)*dot6*M)*M2)/(a2*(-(a22*M2)+th)))/sqr(-M2+sh)+
			     ((-4.*(dot11-dot5*dot6)*M)/(a1*a2)
			      +(8.*dot2*(2.*dot3*dot7-2.*a1*(dot1+dot3)*dot7-dot10*M-2.*a2*(dot1+dot3)*dot6*M)*M2)/(a1*sqr(-(a22*M2)+th))
			      -(2.*M*(2.*(dot1*dot12-dot12*dot3+dot10*(-dot2+dot4)+2.*dot3*(dot2+dot4)*dot6)+
				      ((-1.+2.*a1)*dot13+2.*dot5*(dot7-2.*a12*dot7)-2.*(-1.+2.*a1)*a2*(dot1+dot3)*dot8-4.*(dot2+a1*dot4)*dot9)*M-
				      4.*a1*a2*dot5*dot6*M2))/(a1*a2*(-(a22*M2)+th)))/(-M2+sh));
	    diag[3]=rt*mix1*((-16.*a1*dot4*(dot10-2.*dot1*dot6)*pow<3,1>(M))/(a2*pow<3,1>(-(a12*M2)+uh))
			     -(4.*(-((-1.+2.*a1)*(dot13+2.*dot1*dot8))-2.*(dot2-2.*a1*dot2+dot4)*dot9)*M2)/(a2*sqr(-(a12*M2)+uh)));
	    diag[4]=rt*mix1*((8.*dot4*dot9*M2)/(a2*sqr(-(a12*M2)+uh))
			     +(2.*dot11*M)/((a1-a12)*(-(a12*M2)+uh))
			     +((8.*(dot2+dot4)*dot9*M2)/a1+(8.*(dot2+dot4)*(-2.*dot3*dot7+dot10*M+2.*a1*(dot1+dot3)*(dot7-dot6*M))*M2)/(a1*(-(a12*M2)+uh)))/sqr(-M2+sh)+
			     ((4.*(dot11-dot5*dot6)*M)/(a1*a2)
			      -(8.*dot4*(2.*dot3*dot7-dot10*M-2.*a1*(dot1+dot3)*(dot7-dot6*M))*M2)/(a2*sqr(-(a12*M2)+uh))
			      -(2.*M*(2.*(dot1*dot12-dot12*dot3+dot10*(-dot2+dot4)+2.*dot3*(dot2+dot4)*dot6)+
				      ((-1.+2.*a1)*dot13+4.*a1*a2*dot5*dot7+2.*a1*(-1.+2.*a1)*(dot1+dot3)*dot8-2.*(dot2+dot4+2.*a1*dot4)*dot9)*M-
				      4.*a1*a2*dot5*dot6*M2))/(a1*a2*(-(a12*M2)+uh)))/(-M2+sh));
    	    // 3P1 diagrams
	    diag[0]+=mix2*((8.*(dot13-2.*dot5*dot7-2.*dot3*dot8-dot2*dot9)*M2)/(a1*sqr(-(a22*M2)+th))
			   +((8.*(dot13+2.*dot1*dot8-(2.*dot2+dot4)*dot9)*M2)/a2-
			     (16.*(dot3*(-(a2*dot4*dot7)+a1*dot4*dot6*M-a1*a2*M*(dot12-2.*dot2*dot6-2.*a1*dot8*M))+
				   dot1*((-2.+a1)*dot4*dot7-a2*dot4*dot6*M+a22*M*(dot12-2.*dot2*dot6-2.*a1*dot8*M)))*M2)/(a2*(-(a22*M2)+th)))/sqr(-(a12*M2)+uh)+
			   ((16.*M2*(-(dot2*(dot3+a1*(dot1+dot3))*dot7)-(-(a2*dot1)+a1*dot3)*(dot2*dot6+a1*(dot12-2.*dot2*dot6))*M-2.*a1*a2*(-(a2*dot1)+a1*dot3)*dot8*M2))/(a1*sqr(-(a22*M2)+th))
			    -(8.*(dot2*dot9+a1*(-(a2*dot13)-dot1*dot8+a1*dot1*dot8+dot3*dot8-a1*dot3*dot8+dot2*dot9-
						2.*a1*dot2*dot9-dot4*dot9-a2*dot11*M+a2*dot5*(dot7+dot6*M)))*M2)/(a1*a2*(-(a22*M2)+th)))/(-(a12*M2)+uh));
	    diag[1]+=mix2*((-16.*a2*pow<3,1>(M)*(-((-1.+2.*a1)*dot2*(dot10+2.*dot3*dot6))-2.*a1*a2*(dot13-2.*(dot5*dot7+dot3*dot8))*M))/(a1*pow<3,1>(-(a22*M2)+th))
			   +(8.*((-2.+a1)*(dot13-2.*(dot5*dot7+dot3*dot8))+2.*a2*dot2*dot9-2.*a1*a2*(dot11-2.*dot5*dot6)*M)*M2)/(a1*sqr(-(a22*M2)+th)));
	    diag[2]+=mix2*((8.*M*(dot10*dot2+dot12*dot3+M*(dot13-a1*dot13+a1*dot5*dot7+dot3*dot8+a2*(dot11-dot5*dot6)*M)))/(a1*sqr(-(a22*M2)+th))
			   -(8.*(dot11-dot5*dot6)*M)/(a1*(-(a22*M2)+th))
			   +((8.*(-2.*a2*(dot1+dot3)*dot8+(dot2+dot4)*dot9)*M2)/a2-
			     (8.*(2.*(-(dot2*dot3)+a1*(dot1+dot3)*(dot2-dot4)+(2.*dot1+dot3)*dot4)*dot7+
				  (dot10*(dot2+dot4)-2.*a2*(dot1+dot3)*(dot12-(dot2+dot4)*dot6))*M)*M2)/(a2*(-(a22*M2)+th)))/sqr(-M2+sh)+
			   ((-8.*(dot11-dot5*dot6)*M)/a1+
			    (8.*M2*(2.*(-(dot1*dot2)+(-2.+a1)*dot1*dot4+a1*dot3*dot4)*dot7-
				    (dot10*(dot2+2.*a1*dot4)-2.*dot2*dot3*dot6+2.*a1*dot3*(-(a2*dot12)+(4.*dot2-2.*a1*dot2+dot4)*dot6)-
				     2.*a2*dot1*((1.+a1)*dot12-(2.*a1*dot2+dot4)*dot6))*M-
				    2.*a1*a2*(-dot13-2.*a2*(dot1+dot3)*dot8+2.*(dot2+dot4)*dot9)*M2))/(a1*sqr(-(a22*M2)+th))
			    -(8.*M*(-(a2*(dot12*dot3+dot10*(dot2-dot4)))-
				    (a22*(dot13-dot5*dot7)+(-1.+a12)*dot3*dot8+(dot4+a1*(dot2+2.*a2*dot4))*dot9)*M-
				    a2*dot1*(-dot12+2.*(dot2+dot4)*dot6-a2*dot8*M)-a2*(dot11-(1.+a1)*dot5*dot6)*M2))/(a1*a2*(-(a22*M2)+th)))/(-M2+sh));
	    diag[3]+=mix2*((-16.*a1*pow<3,1>(M)*((-1.+2.*a1)*dot4*(dot10-2.*dot1*dot6)-2.*a1*a2*(dot13+2.*dot1*dot8-2.*(dot2+dot4)*dot9)*M))/(a2*pow<3,1>(-(a12*M2)+uh))
			   -(8.*((1.+a1)*dot13+2.*(1.+a1)*dot1*dot8-2.*(dot2+a1*dot2+dot4)*dot9+2.*a1*a2*dot11*M)*M2)/(a2*sqr(-(a12*M2)+uh)));
	    diag[4]+=mix2*((-8.*(-dot13+2.*dot2*dot9+dot4*dot9-2.*a1*dot11*M+2.*dot5*(dot7+a1*dot6*M))*M2)/(a2*sqr(-(a12*M2)+uh))
			   +((-8.*M2*(2.*(-dot3+a1*(dot1+dot3))*dot4*dot7+
				      (dot10*dot4+2.*((-1.+a12)*dot1*dot12+(-3.+a1)*a1*dot12*dot3+a1*dot1*dot4*dot6+a1*dot4*(dot10+3.*dot3*dot6)))*M-
				      4.*dot2*dot3*(dot7+(-2.+a1)*a1*dot6*M)+4.*a2*dot2*M*(dot10+a1*dot1*dot6+a1*dot9*M)+
				      2.*a1*a2*(-3.*dot13+4.*dot5*dot7+2.*a1*(dot1+dot3)*dot8+2.*dot4*dot9)*M2))/(a2*sqr(-(a12*M2)+uh))
			     -(8.*((dot2+dot4)*dot9-a1*(dot13-dot5*dot7+3.*dot1*dot8+dot3*dot8-3.*dot2*dot9+dot5*dot6*M)+
				   a12*(dot13-dot5*dot7+dot1*dot8-dot3*dot8-2.*dot2*dot9-2.*dot11*M+3.*dot5*dot6*M))*M2)/(a1*a2*(-(a12*M2)+uh)))/(-M2+sh)
			   +((8.*(2.*a1*(dot1+dot3)*dot8-(dot2+dot4)*dot9)*M2)/a1+
			     (8.*((dot2+dot4)*(2.*dot3*dot7-dot10*M)+
				  2.*a1*(dot1+dot3)*((dot2-dot4)*dot7+(dot12-(dot2+dot4)*dot6)*M))*M2)/(a1*(-(a12*M2)+uh)))/sqr(-M2+sh));
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
    meSum *=  1./972.;
  }
  // save the diagram weights
  meInfo(save);
  // final factors
  return meSum/M2*O1_*pow<3,1>(Constants::pi*standardModel()->alphaS(scale())/M)/sqr(a1*a2);
}

IBPtr MEPP2BCP1Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2BCP1Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2BCP1Jet::doinit() {
  setBcState(10543+mode_*10000);
  MEPP2BCJetBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<1>(bcbar,principleQuantumNumber(),0,1);
  double theta = oniumParameters()->singletTripletMixing(principleQuantumNumber(),1);
  sTheta_ = sin(theta);
  cTheta_ = cos(theta);
}

void MEPP2BCP1Jet::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2*GeV2) << mode_ << theta_ << sTheta_ << cTheta_;
}

void MEPP2BCP1Jet::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2*GeV2) >> mode_ >> theta_ >> sTheta_ >> cTheta_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2BCP1Jet,MEPP2BCJetBase>
describeHerwigMEPP2BCP1Jet("Herwig::MEPP2BCP1Jet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2BCP1Jet::Init() {

  static ClassDocumentation<MEPP2BCP1Jet> documentation
    ("The MEPP2BCP1Jet class implements the matrix element for g c -> B_c(P1) b q and "
     "q bar -> Bc(P1) g and charged conjugate");

  static Switch<MEPP2BCP1Jet,unsigned int> interfaceMode
    ("Mode",
     "Which b_c state to produce",
     &MEPP2BCP1Jet::mode_, 0, false, false);
  static SwitchOption interfaceModeB_c1
    (interfaceMode,
     "B_c1",
     "Produce the lighter B_c1 state",
     0);
  static SwitchOption interfaceModeB_c1Prime
    (interfaceMode,
     "B_c1Prime",
     "Produce the heavier B_c1 state",
     1);

}

