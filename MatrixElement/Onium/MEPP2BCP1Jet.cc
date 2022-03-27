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
  static bool first=true;
  if(first) {
    cerr << "testing mixing " << mix1 << " " << mix2 << "\n";
    first=false;
  }
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
	      diag[0]=rt*mix1*(4.*M2*(a1*M*(4.*(dot1+dot3)*(dot10+2.*dot2*dot6)+a1*(-1.+2.*a1)*(dot11-2.*dot5*dot7)*M-
					    2.*a1*(dot1*dot8+2.*a1*dot3*dot8+dot2*dot9-2.*a1*dot2*dot9)*M)+
				      (dot11-2.*a1*dot11+2.*((-1.+2.*a1)*dot5*dot7+dot1*dot8+
							     2.*a1*dot3*dot8+dot2*dot9-2.*a1*dot2*dot9))*sh))/(a2*pow<3,1>(a12*M2-sh));
	      diag[1]=rt*mix1*(-2.*M*(4.*pow(a1,5)*dot1*dot8*pow<5,1>(M)-
				      2.*pow(a1,4)*pow<3,1>(M)*(4.*dot1*(-3.*dot4*dot7-2.*dot4*dot6*M+3.*dot2*(dot7+dot6*M)+3.*dot8*M2)-
								(2.*dot5*dot7+2.*dot2*dot9-2.*dot4*dot9+dot5*dot6*M)*(M2-th))+
				      2.*pow(a1,3)*pow<3,1>(M)*((dot11-4.*dot5*dot7-2.*dot3*dot8-dot2*dot9+dot4*dot9-dot13*M)*(M2-th)+
								2.*dot1*(-(dot10*M)+6.*dot2*(dot7+dot6*M)-4.*dot4*(3.*dot7+dot6*M)+dot8*(8.*M2+sh-3.*uh)))+
				      a12*M*((M2-th)*(-2.*dot12*dot4*M+dot13*pow<3,1>(M)-2.*dot5*dot6*pow<3,1>(M)+8.*dot2*dot3*(dot7+dot6*M)-
						      8.*dot3*dot4*(dot7+dot6*M)+4.*dot3*M*(dot10+2.*dot8*M)-dot11*M2+2.*dot5*dot7*M2-
						      4.*dot5*dot7*sh+4.*dot4*dot9*sh+dot13*M*sh-2.*dot5*dot6*M*sh-2.*dot2*(dot12*M+2.*dot9*sh)-dot13*M*uh+2.*dot5*dot6*M*uh)-
					     2.*dot1*(dot10*M*(-3.*M2+th)+3.*dot8*M2*(3.*M2-sh-3.*uh)+4.*dot2*(dot7+dot6*M)*(M2-2.*sh-uh)+
						      2.*dot4*dot6*M*(-M2+sh+uh)+4.*dot4*dot7*(-4.*M2+2.*sh+uh)))+
				      sh*(-2.*dot1*(4.*dot4*dot7*M+dot10*(3.*M2-th)-dot8*M*(M2+sh-uh)+2.*dot4*dot6*(-M2+sh+uh))+
					  (M2-th)*(-4.*dot10*dot3+2.*dot12*(dot2+dot4)+dot11*M-dot13*(M2+sh-uh)-2.*dot5*(M*(dot7-dot6*M)+dot6*uh)))+
				      2.*a1*M*(2.*dot1*(dot8*pow<4,1>(M)-2.*dot2*dot7*sh-2.*dot2*dot6*M*sh-2.*dot8*M2*sh-dot8*sqr(sh)+dot10*M*(M2+sh-th)-
							2.*dot8*M2*uh+dot8*sh*uh+dot8*sqr(uh)+2.*dot4*dot7*(-M2+3.*sh+uh))+
					       (M2-th)*((-dot11+4.*dot5*dot7+dot2*dot9-dot4*dot9+dot13*M)*sh+
							2.*dot3*(2.*dot4*dot7+M*(dot10-dot8*M)+dot8*uh)))))/(a1*a2*sqr(-(a12*M2)+sh)*sqr(M2-th));
	      diag[2]=rt*mix1*(-2.*M*(-4.*pow(a1,5)*dot1*dot8*pow<5,1>(M)+
				      2.*pow(a1,4)*pow<3,1>(M)*(4.*dot1*(dot4*dot7-dot2*(dot7+dot6*M)+dot8*M2)+
								(2.*dot5*dot7+2.*dot2*dot9-2.*dot4*dot9+dot5*dot6*M)*(M2-th))+
				      4.*dot3*(M2-th)*(M*(2.*dot10*M-2.*dot2*dot6*M+2.*dot4*(dot7+dot6*M)-dot8*sh)-dot10*uh)+
				      (M2-uh)*((M2-th)*(-2.*dot12*(dot2+dot4)-dot11*M+2.*dot2*dot9*M-2.*dot4*dot9*M+dot13*M2-
							2.*dot5*dot6*M2+dot13*sh-dot13*uh+2.*dot5*dot6*uh)+2.*dot1*(dot10*(M2-th)+2.*dot4*dot6*(-M2+sh+uh)))+
				      2.*pow(a1,3)*pow<3,1>(M)*((dot11-6.*dot5*dot7-2.*dot3*dot8-7.*dot2*dot9+
								 7.*dot4*dot9-(dot13+2.*dot5*dot6)*M)*(M2-th)+
								2.*dot1*(-6.*dot4*dot7-dot10*M+dot2*(4.*dot7+6.*dot6*M)+dot8*(-M2+sh+uh)))+
				      a12*M*((M2-th)*(8.*dot3*dot4*(dot7+dot6*M)+4.*dot3*M*(dot10+dot8*M)+
						      M*(-2.*dot12*dot4-5.*dot11*M+M*(12.*dot5*dot7-18.*dot4*dot9+5.*dot13*M)+dot13*sh)-
						      (4.*dot5*dot7-4.*dot4*dot9+dot13*M)*uh-2.*dot2*(4.*dot3*(dot7+dot6*M)+M*(dot12-9.*dot9*M)+2.*dot9*uh))-
					     2.*dot1*(4.*dot2*(dot7+3.*dot6*M)*M2+M*(4.*dot8*M*sh+dot10*(-5.*M2+th))-
						      4.*dot2*(dot7+dot6*M)*uh-2.*dot4*(dot7*(6.*M2-2.*uh)+dot6*M*(-M2+sh+uh))))+
				      2.*a1*M*(-((M2-th)*(M*(-2.*dot12*dot4+M*(-2.*dot11+2.*dot5*dot7-5.*dot4*dot9+
									       2.*dot13*M-2.*dot5*dot6*M)+dot13*sh)+
							  dot3*(6.*dot10*M+8.*dot4*(dot7+dot6*M)-2.*dot8*sh)+
							  (dot11-2.*dot5*dot7+3.*dot4*dot9-2.*dot13*M+2.*dot5*dot6*M)*uh))+
					       dot2*(4.*dot1*dot6*pow<3,1>(M)-5.*dot9*pow<4,1>(M)+4.*dot3*(dot7+2.*dot6*M)*(M2-th)+
						     5.*dot9*M2*th+2.*dot12*(pow<3,1>(M)-M*th)-4.*dot1*dot6*M*uh+3.*dot9*M2*uh-3.*dot9*th*uh)+
					       2.*dot1*(dot8*sh*(M2-uh)+dot10*M*(-2.*M2+th+uh)-
							2.*dot4*(dot7*(M2-uh)+dot6*M*(-M2+sh+uh))))))/(a1*a2*sqr(M2-th)*sqr(-(a22*M2)+uh));
	      diag[3]=rt*mix1*(-4.*M2*(-4.*a2*dot10*dot3*M-(-1.+2.*a1)*(dot11+2.*dot4*dot9)*(a22*M2-uh)-
				       4.*a2*dot3*(M*(2.*dot4*dot6+a22*dot8*M)-dot8*uh)))/(a1*pow<3,1>(a22*M2-uh));
	      diag[4]=rt*mix1*(2.*M*(-((a1*M*(-4.*(dot1+dot3)*dot8+a1*(dot13-2.*dot5*dot6)*M)-(dot13-2.*dot5*dot6)*sh)/a2)-
				     (4.*dot3*M*(a12*M2-sh)*(2.*dot4*dot7-2.*dot2*dot6*M+2.*a1*(dot2-dot4)*(dot7+dot6*M)+a12*dot8*M2-dot8*sh))/sqr(-(a22*M2)+uh)+
				     (2.*pow(a1,3)*pow<3,1>(M)*(2.*dot1*dot8+4.*dot3*dot8-3.*dot2*dot9+dot4*dot9-2.*dot5*(dot7+dot6*M))+
				      pow(a1,4)*pow<3,1>(M)*(4.*dot2*dot9-4.*dot4*dot9-dot13*M+4.*dot5*(dot7+dot6*M))-
				      sh*(4.*dot10*dot3-2.*dot12*(dot2+dot4)+2.*dot1*(dot10+2.*dot4*dot6)-4.*dot3*dot8*M+2.*dot2*dot9*M+dot13*sh)+
				      2.*a1*M*(4.*(dot1+dot3)*(dot4*dot7-dot2*dot6*M)+(-2.*dot1*dot8-4.*dot3*dot8+3.*dot2*dot9-dot4*dot9+2.*dot5*(dot7+dot6*M))*sh)+
				      2.*a12*M*(4.*dot1*(dot2-dot4)*dot7-4.*dot3*dot4*dot7+2.*dot10*dot3*M-dot12*dot4*M-4.*dot3*dot4*dot6*M+
						dot1*(dot10+4.*dot2*dot6-2.*dot4*dot6)*M-2.*dot3*dot8*M2-2.*dot5*dot7*sh+2.*dot4*dot9*sh+dot13*M*sh-
						2.*dot5*dot6*M*sh+dot2*(-(dot12*M)+4.*dot3*(dot7+dot6*M)+
									dot9*M2-2.*dot9*sh)))/(-(pow(a2,3)*M2)+uh-a1*uh)))/(a1*sqr(-(a12*M2)+sh));
	      // diagrams 3P1
	      diag[0]+= mix2*8.*M2*(-2.*a1*(dot1+dot3)*(dot10+2.*dot2*dot6)*M-2.*pow(a1,4)*(dot13-2.*dot5*dot6)*pow<3,1>(M)+
				    pow(a1,3)*(3.*dot11-6.*dot5*dot7+8.*dot1*dot8+2.*dot3*dot8+6.*dot2*dot9+2.*dot13*M-4.*dot5*dot6*M)*M2+
				    (dot11-2.*dot5*dot7+2.*dot1*dot8+2.*dot2*dot9)*sh+a1*(dot11-2.*(dot5*dot7+dot3*dot8-dot2*dot9+dot13*M-2.*dot5*dot6*M))*sh+
				    a12*M*(4.*dot10*dot3+8.*dot2*dot3*dot6-5.*dot11*M+10.*dot5*dot7*M-10.*dot2*dot9*M+2.*dot1*(2.*dot10+4.*dot2*dot6-5.*dot8*M)+
					   2.*dot13*sh-4.*dot5*dot6*sh))/(a2*pow<3,1>(a12*M2-sh));
	      diag[1]+=-mix2*8.*M2*(4.*pow(a1,6)*(-2.*dot2+dot4)*dot9*pow<4,1>(M)+pow(a1,5)*(5.*dot1*dot8+22.*dot2*dot9-14.*dot4*dot9)*pow<4,1>(M)+
				    dot1*sh*(2.*dot4*dot7+M*(dot10-dot8*M)+dot8*uh)+
				    pow(a1,4)*M2*(4.*dot3*dot4*dot7-2.*dot12*dot4*M-2.*dot13*pow<3,1>(M)+3.*dot5*dot6*pow<3,1>(M)+dot11*M2-
						  dot5*dot7*M2-2.*dot3*dot8*M2+13.*dot4*dot9*M2+
						  dot1*(-6.*dot4*dot7+14.*dot4*dot6*M+6.*dot2*(dot7-dot6*M)-7.*dot8*M2)+
						  3.*dot4*dot9*sh-dot11*th+dot5*dot7*th+2.*dot3*dot8*th+2.*dot13*M*th-3.*dot5*dot6*M*th-
						  3.*dot4*dot9*uh+dot2*(8.*dot3*dot7-2.*dot12*M-19.*dot9*M2+3.*dot9*sh+5.*dot9*uh))+
				    a12*(-8.*dot12*dot2*pow<3,1>(M)-6.*dot10*dot3*pow<3,1>(M)+6.*dot12*dot4*pow<3,1>(M)-2.*dot2*dot3*dot6*pow<3,1>(M)+
					 2.*dot3*dot4*dot6*pow<3,1>(M)+4.*dot11*pow<4,1>(M)-4.*dot5*dot7*pow<4,1>(M)-2.*dot13*pow<5,1>(M)+
					 2.*dot5*dot6*pow<5,1>(M)+2.*dot2*dot3*dot7*M2+18.*dot3*dot4*dot7*M2-6.*dot2*dot3*dot7*sh+
					 2.*dot12*dot4*M*sh-dot5*dot6*pow<3,1>(M)*sh-dot11*M2*sh+dot5*dot7*M2*sh+2.*dot3*dot8*M2*sh+
					 7.*dot2*dot9*M2*sh-dot4*dot9*M2*sh+dot2*dot9*sqr(sh)-3.*dot4*dot9*sqr(sh)-6.*dot3*dot4*dot7*th+
					 2.*dot12*dot2*M*th+6.*dot10*dot3*M*th-6.*dot12*dot4*M*th+2.*dot2*dot3*dot6*M*th-2.*dot3*dot4*dot6*M*th+
					 2.*dot13*pow<3,1>(M)*th-2.*dot5*dot6*pow<3,1>(M)*th-4.*dot11*M2*th+4.*dot5*dot7*M2*th+dot11*sh*th-
					 dot5*dot7*sh*th-2.*dot3*dot8*sh*th+dot5*dot6*M*sh*th-
					 (4.*dot3*dot4*dot7+dot4*dot9*sh+dot2*(2.*dot3*dot7-2.*dot12*M+dot9*sh)+2.*(dot13-dot5*dot6)*M*(-M2+th))*uh+
					 dot1*(dot10*M*(-3.*M2+2.*th)+2.*dot4*dot6*M*(3.*M2-4.*sh-3.*uh)+2.*dot2*(dot7-dot6*M)*(M2-2.*sh-uh)+
					       dot8*M2*(M2-sh-uh)+2.*dot4*dot7*(-4.*M2+2.*sh+uh)))+
				    a1*(2.*dot12*dot2*pow<3,1>(M)+dot11*pow<4,1>(M)-2.*dot12*dot4*M*sh+dot5*dot6*pow<3,1>(M)*sh-dot5*dot7*M2*sh-
					3.*dot2*dot9*M2*sh+dot4*dot9*M2*sh+dot2*dot9*sqr(sh)+dot4*dot9*sqr(sh)-dot11*M2*th+dot5*dot7*sh*th-
					dot5*dot6*M*sh*th-(2.*dot12*dot2*M+(-3.*dot2+dot4)*dot9*sh+dot11*(M2-th))*uh+
					dot1*(2.*dot2*dot7*sh-M*(dot8*pow<3,1>(M)+2.*dot2*dot6*sh-4.*dot8*M*sh+dot10*(M2+sh-th))+
					      2.*dot4*(dot6*M*sh+dot7*(M2-3.*sh-uh))+2.*dot8*(M2-sh)*uh-dot8*sqr(uh))+
					dot3*(4.*dot2*dot7*sh+(3.*dot10*M+dot8*(M2-3.*sh-th))*(M2-th)+2.*dot4*dot7*(-3.*M2+th+2.*uh)))+
				    pow(a1,3)*M*(-12.*dot3*dot4*dot7*M-5.*dot11*pow<3,1>(M)+5.*dot5*dot7*pow<3,1>(M)+3.*dot3*dot8*pow<3,1>(M)-
						 5.*dot4*dot9*pow<3,1>(M)+4.*dot13*pow<4,1>(M)-5.*dot5*dot6*pow<4,1>(M)-4.*dot3*dot4*dot6*M2+
						 dot4*dot9*M*sh+dot1*M*(2.*dot4*(6.*dot7-7.*dot6*M)-6.*dot2*(dot7-dot6*M)+M*(dot10+2.*dot8*M)-2.*dot8*sh)+
						 2.*dot12*dot4*th+4.*dot3*dot4*dot6*th+5.*dot11*M*th-5.*dot5*dot7*M*th-3.*dot3*dot8*M*th-4.*dot13*M2*th+
						 5.*dot5*dot6*M2*th+5.*dot4*dot9*M*uh-
						 dot2*(8.*dot3*dot7*M+2.*dot12*(-4.*M2+th)+4.*dot3*dot6*(-M2+th)+dot9*M*(-7.*M2+11.*sh+7.*uh))))/
		(a1*a2*sqr(-(a12*M2)+sh)*sqr(M2-th));
	      diag[2]+= mix2*8.*M*(-(dot1*(pow(a1,5)*dot8*pow<5,1>(M)-2.*pow(a1,4)*pow<3,1>(M)*(dot4*dot7+dot4*dot6*M+3.*dot8*M2-dot8*th)+
					   pow(a1,3)*M2*(2.*dot4*M*(3.*dot7+dot6*M)+dot10*(-2.*M2+th)+2.*dot8*M*(8.*M2-5.*th-uh))+
					   (dot10+dot8*M)*(M2-th)*(M2-uh)+
					   a1*(-6.*dot10*pow<4,1>(M)-dot10*th*uh+dot10*M2*(5.*th+2.*uh)-2.*dot4*M*(dot6*M*sh+dot7*(-th+uh))+
					       dot8*M*(5.*pow<4,1>(M)-4.*M2*(th+uh)+uh*(2.*th+uh)))+
					   a12*M*(dot10*(7.*pow<3,1>(M)-5.*M*th)+2.*dot4*(dot6*M*sh+dot7*(-2.*M2-th+uh))+
						  dot8*(-17.*pow<4,1>(M)-2.*th*uh+M2*(13.*th+6.*uh)))+
					   2.*dot2*((-2.+a1)*dot7*M*((1.+(-2.+a1)*a12)*M2-a2*th-a1*uh)-
						    a2*dot6*((-2.+a12)*a2*pow<4,1>(M)-th*uh+M2*(-((-2.+a1)*a2*th)+uh+a1*uh)))))-
				   a2*(-2.*pow(a1,5)*dot4*dot9*pow<5,1>(M)+pow(a1,4)*pow<4,1>(M)*(8.*dot4*dot9*M+(dot13+dot5*dot6)*(M2-th))+
				       pow(a1,3)*pow<3,1>(M)*(-2.*dot12*dot4*M-3.*dot5*dot6*pow<3,1>(M)-dot11*M2+dot5*dot7*M2-17.*dot4*dot9*M2+dot11*th-
							      dot5*dot7*th+5.*dot4*dot9*th+3.*dot5*dot6*M*th-2.*dot3*(2.*dot4*dot7+dot8*(-M2+th))+4.*dot4*dot9*uh)+
				       a12*M2*(-(dot11*pow<3,1>(M))+dot5*dot7*pow<3,1>(M)+17.*dot4*dot9*pow<3,1>(M)+2.*dot5*dot6*pow<4,1>(M)-
					       2.*dot13*M2*sh-dot5*dot6*M2*sh-dot3*(4.*dot4*dot6+3.*dot8*M)*(M2-th)+dot11*M*th-dot5*dot7*M*th-
					       9.*dot4*dot9*M*th-dot13*M2*th-2.*dot5*dot6*M2*th+2.*dot13*sh*th+dot5*dot6*sh*th+dot13*sqr(th)+
					       dot12*dot4*(-M2+th)-8.*dot4*dot9*M*uh-dot5*dot6*M2*uh+dot5*dot6*th*uh)+
				       (M2-th)*(-(dot5*dot7*pow<3,1>(M))-dot12*dot4*sh-dot11*M*sh+2.*dot5*dot7*M*sh-
						2.*dot5*dot6*M2*sh+dot13*sqr(sh)-dot4*dot9*M*th+dot13*sh*th+dot10*dot3*(2.*M2+sh-uh)+
						dot5*dot7*M*uh+dot4*dot9*M*uh+dot5*dot6*sh*uh+dot3*M*(-2.*dot4*dot6*M+dot8*M2-dot8*uh))-
				       a1*M*(2.*dot5*dot7*pow<4,1>(M)+5.*dot4*dot9*pow<4,1>(M)+dot13*pow<5,1>(M)-2.*dot12*dot4*M*sh-
					     3.*dot5*dot6*pow<3,1>(M)*sh+dot5*dot7*M2*sh-dot11*(2.*M2+sh)*(M2-th)-2.*dot13*pow<3,1>(M)*th-
					     2.*dot5*dot7*M2*th-3.*dot4*dot9*M2*th-dot5*dot7*sh*th+3.*dot5*dot6*M*sh*th+dot13*M*sqr(th)-
					     5.*dot4*dot9*M2*uh+dot4*dot9*th*uh+2.*dot4*dot9*sqr(uh)+
					     2.*dot3*(dot4*dot7*(3.*M2-2.*sh-3.*th)+3.*dot4*dot6*M*(-M2+th)+(M2-th)*(2.*dot10*M+dot8*uh)))+
				       dot2*(2.*dot3*M*(-((-1.+2.*a1)*a2*dot6*M*(M2-th))+dot7*((-1.+3.*a1+2.*pow(a1,3))*M2-2.*a1*sh+th-3.*a1*th))-
					     dot12*((3.+a12*(-5.+2.*a1))*pow<4,1>(M)+th*uh+M2*(-((3.+a1)*a2*th)-(1.+2.*a1)*uh))+
					     dot9*M*((1.+a1*(5.+2.*(-2.+a1)*a1))*a22*pow<4,1>(M)+uh*(th-a1*th+2.*a1*uh)+a2*M2*(-uh-a2*(th+3.*a1*th+4.*a1*uh))))))
		/(a1*a2*sqr(M2-th)*sqr(-(a22*M2)+uh));
	      diag[3]+=-mix2*8.*M2*(-2.*(-1.+2.*a1)*a2*dot10*dot3*M+a22*((2.+3.*a1)*(dot11+2.*dot4*dot9)-2.*a1*a2*dot13*M)*M2+
				    ((-2.+a1)*(dot11+2.*dot4*dot9)+2.*a1*a2*dot13*M)*uh-
				    2.*a2*dot3*(M*(2.*(-1.+2.*a1)*dot4*dot6+a22*dot8*M)-dot8*uh))/(a1*pow<3,1>(a22*M2-uh));
	      diag[4]+=-mix2*8.*M2*(pow(a1,5)*pow<4,1>(M)*(dot11-2.*dot3*dot8-9.*dot2*dot9+3.*dot4*dot9+3.*dot13*M-3.*dot5*(dot7+dot6*M))+
				    pow(a1,6)*pow<4,1>(M)*(-dot11+2.*dot3*dot8+(dot2+dot4)*dot9-dot13*M+dot5*(dot7+dot6*M))-
				    sh*((dot11+2.*dot4*dot9)*sh+dot3*(-2.*dot4*dot7+2.*dot2*dot6*M-dot8*M2+dot8*th))+
				    pow(a1,4)*M2*(4.*dot12*(-dot2+dot4)*M-4.*dot2*dot3*dot6*M+4.*dot3*dot4*dot6*M-3.*dot13*pow<3,1>(M)+
						  3.*dot5*dot6*pow<3,1>(M)-2.*dot11*M2+3.*dot5*dot7*M2+2.*dot3*dot8*M2+23.*dot2*dot9*M2-
						  7.*dot4*dot9*M2+2.*dot1*(-(dot2*dot7)+dot4*dot7+dot2*dot6*M-dot4*dot6*M+dot8*M2)+
						  2.*dot11*sh-4.*dot3*dot8*sh+3.*dot2*dot9*sh-5.*dot4*dot9*sh+dot13*M*sh-dot5*dot6*M*sh+
						  dot11*th-dot5*dot7*th-2.*dot3*dot8*th+(-2.*dot5*dot7-5.*dot2*dot9+3.*dot4*dot9+dot13*M-dot5*dot6*M)*uh)+
				    pow(a1,3)*M*(2.*dot2*dot3*dot7*M+2.*dot3*dot4*dot7*M-dot11*pow<3,1>(M)-dot5*dot7*pow<3,1>(M)-25.*dot2*dot9*pow<3,1>(M)+
						 3.*dot4*dot9*pow<3,1>(M)+dot13*pow<4,1>(M)-dot5*dot6*pow<4,1>(M)+10.*dot12*dot2*M2-6.*dot12*dot4*M2+
						 10.*dot2*dot3*dot6*M2-6.*dot3*dot4*dot6*M2+2.*dot12*dot2*sh-2.*dot12*dot4*sh+4.*dot2*dot3*dot6*sh-
						 4.*dot3*dot4*dot6*sh-2.*dot11*M*sh+4.*dot3*dot8*M*sh-9.*dot2*dot9*M*sh+7.*dot4*dot9*M*sh-
						 3.*dot13*M2*sh+3.*dot5*dot6*M2*sh+dot11*M*th+3.*dot5*dot7*M*th+
						 (-2.*(dot2-dot4)*(dot12+2.*dot3*dot6)+(2.*dot5*dot7+13.*dot2*dot9-3.*dot4*dot9)*M+(-dot13+dot5*dot6)*M2)*uh+
						 dot1*M*(8.*dot2*dot7-6.*dot4*dot7-6.*dot2*dot6*M+4.*dot4*dot6*M-5.*dot8*M2-dot8*sh+dot8*uh))+
				    a1*(2.*dot12*dot2*pow<3,1>(M)-dot11*pow<4,1>(M)-2.*dot2*dot9*pow<4,1>(M)+2.*dot12*dot2*M*sh-
					dot13*pow<3,1>(M)*sh+dot5*dot6*pow<3,1>(M)*sh+dot11*M2*sh-3.*dot2*dot9*M2*sh+dot4*dot9*M2*sh+
					dot11*sqr(sh)+2.*dot4*dot9*sqr(sh)+2.*dot11*M2*th+dot5*dot7*M2*th-dot11*sh*th-dot11*sqr(th)-
					dot1*(2.*dot4*dot7+dot2*(-4.*dot7+2.*dot6*M)+dot8*(M2+sh-uh))*(M2-uh)-
					(2.*dot12*dot2*M-4.*dot2*dot9*M2+dot2*dot9*sh+dot4*dot9*sh-dot13*M*sh+dot5*(dot7*M2+dot6*M*sh+dot7*th))*uh+
					(dot5*dot7-2.*dot2*dot9)*sqr(uh)+
					dot3*(dot8*(-2.*sqr(sh)+sqr(M2-th))-2.*dot4*(dot6*M*sh+dot7*(-M2+sh+th))+
					      2.*dot2*(dot6*pow<3,1>(M)+4.*dot6*M*sh+dot7*th-(dot7+dot6*M)*uh)))+
				    a12*(2.*dot12*dot4*pow<3,1>(M)+2.*dot3*dot4*dot6*pow<3,1>(M)-dot3*dot8*pow<4,1>(M)-4.*dot3*dot4*dot7*M2+
					 2.*dot12*dot4*M*sh+6.*dot3*dot4*dot6*M*sh+3.*dot13*pow<3,1>(M)*sh-3.*dot5*dot6*pow<3,1>(M)*sh+
					 3.*dot11*M2*sh-2.*dot3*dot8*M2*sh-3.*dot4*dot9*M2*sh-dot11*sqr(sh)+2.*dot3*dot8*sqr(sh)+
					 2.*dot3*dot4*dot7*th-3.*dot5*dot7*M2*th+dot3*dot8*M2*th-dot11*sh*th+2.*dot3*dot8*sh*th-
					 2.*dot12*dot4*M*uh-2.*dot3*dot4*dot6*M*uh+dot5*dot7*M2*uh+dot4*dot9*sh*uh-dot13*M*sh*uh+
					 dot5*dot6*M*sh*uh+dot5*dot7*th*uh+dot5*dot7*sqr(uh)+
					 2.*dot1*(M2*(-5.*dot2*dot7+3.*dot4*dot7+3.*dot2*dot6*M-dot4*dot6*M+2.*dot8*M2+dot8*sh)+
						  (-(dot4*dot7)+dot4*dot6*M+dot2*(dot7-dot6*M)-2.*dot8*M2)*uh)+
					 dot2*(-2.*dot3*(dot7*(M2+th)+dot6*M*(4.*M2+5.*sh-3.*uh))-4.*dot12*M*(2.*M2+sh-uh)+
					       dot9*(12.*pow<4,1>(M)+3.*M2*(3.*sh-4.*uh)+sh*uh))))/(a1*a2*sqr(-(a12*M2)+sh)*sqr(-(a22*M2)+uh));
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
	      // diagrams
	      diag[0]=rt*mix1*(-4.*M2*(-(a2*M*(4.*dot3*(-2.*dot2*dot6+dot9)-
					       a2*((-1.+2.*a1)*dot10-4.*a2*dot3*dot7+2.*(1.-2.*a1)*dot2*dot8)*M))+
				       (dot10-2.*a1*dot10+4.*a2*dot3*dot7+2.*(-1.+2.*a1)*dot2*dot8)*sh-
				       4.*a2*dot1*(M*(-2.*dot2*dot6+dot9+a22*dot7*M)-dot7*sh)))/(a1*pow<3,1>(a22*M2-sh));
	      diag[1]=rt*mix1*(-2.*M*(4.*pow(a1,5)*dot1*dot7*pow<5,1>(M)-
				      2.*pow(a1,4)*pow<3,1>(M)*(4.*dot1*dot11*(-dot2+dot4)+2.*dot1*M*(2.*dot4*dot6+dot7*M)+
								(2.*dot11*dot5+2.*dot2*dot8-2.*dot4*dot8-dot13*M+dot5*dot6*M)*(M2-th))+
				      2.*a1*M*(2.*dot3*(M2-th)*(-2.*dot11*dot2+M*(dot9+dot7*M)-dot7*th)+
					       (dot1*dot11*(8.*dot2-4.*dot4)-2.*dot1*(2.*dot4*dot6+dot9)*M+
						(dot10-4.*dot11*dot5-2.*dot3*dot7-dot2*dot8+dot4*dot8)*(M2-th))*uh+2.*dot1*dot7*sqr(uh))+
				      2.*pow(a1,3)*pow<3,1>(M)*(-((dot10-4.*dot11*dot5-2.*dot3*dot7-dot2*dot8+dot4*dot8)*(M2-th))+
								dot1*(4.*dot11*(-2.*dot2+dot4)+4.*dot4*dot6*M+2.*dot9*M-4.*dot7*uh))-
				      uh*((M2-th)*(-2.*dot12*(dot2+dot4)+4.*dot3*dot9+dot10*M-2.*dot11*dot5*M+dot13*M2+
						   2.*dot5*dot6*M2-2.*dot5*dot6*sh-dot13*th-2.*dot13*uh)+
					  2.*dot1*(4.*dot11*dot2*M-dot7*pow<3,1>(M)-dot9*M2-dot9*th+dot7*M*th+2.*dot7*M*uh+2.*dot2*dot6*(-M2+sh+uh)))+
				      a12*M*((M2-th)*(M*(-2.*dot12*(dot2+dot4)+8.*dot2*dot3*dot6-8.*dot3*dot4*dot6+4.*dot3*dot9+M*(dot10+dot13*M)+
							 2.*dot5*dot6*M2-2.*dot5*dot6*sh-dot13*th)+2.*(2.*dot2*dot8-2.*dot4*dot8-2.*dot13*M+dot5*dot6*M)*uh+
						      dot11*(8.*dot2*dot3-8.*dot3*dot4-2.*dot5*M2+4.*dot5*uh))+
					     2.*dot1*(4.*dot11*(dot2*(M2-uh)+dot4*uh)+
						      M*(-(dot9*(M2+th))+4.*dot4*dot6*uh+2.*dot2*dot6*(-M2+sh+uh)+dot7*M*(-M2+th+4.*uh))))))
		/(a1*a2*sqr(M2-th)*sqr(-(a12*M2)+uh));
	      diag[2]=-rt*mix1*((-2.*M*(M2-th)*(-(dot10*pow<3,1>(M))+4.*a1*dot10*pow<3,1>(M)-5.*a12*dot10*pow<3,1>(M)+2.*pow(a1,3)*dot10*pow<3,1>(M)+
						8.*dot3*dot7*pow<3,1>(M)-16.*a1*dot3*dot7*pow<3,1>(M)+12.*a12*dot3*dot7*pow<3,1>(M)-
						4.*pow(a1,3)*dot3*dot7*pow<3,1>(M)+2.*dot2*dot8*pow<3,1>(M)-10.*a1*dot2*dot8*pow<3,1>(M)+
						18.*a12*dot2*dot8*pow<3,1>(M)-14.*pow(a1,3)*dot2*dot8*pow<3,1>(M)+4.*pow(a1,4)*dot2*dot8*pow<3,1>(M)-
						2.*dot4*dot8*pow<3,1>(M)+10.*a1*dot4*dot8*pow<3,1>(M)-18.*a12*dot4*dot8*pow<3,1>(M)+
						14.*pow(a1,3)*dot4*dot8*pow<3,1>(M)-4.*pow(a1,4)*dot4*dot8*pow<3,1>(M)+3.*dot13*pow<4,1>(M)-
						10.*a1*dot13*pow<4,1>(M)+13.*a12*dot13*pow<4,1>(M)-8.*pow(a1,3)*dot13*pow<4,1>(M)+
						2.*pow(a1,4)*dot13*pow<4,1>(M)-2.*dot5*dot6*pow<4,1>(M)+4.*a1*dot5*dot6*pow<4,1>(M)-
						4.*pow(a1,3)*dot5*dot6*pow<4,1>(M)+2.*pow(a1,4)*dot5*dot6*pow<4,1>(M)+8.*dot2*dot3*dot6*M2-
						16.*a1*dot2*dot3*dot6*M2+8.*a12*dot2*dot3*dot6*M2-8.*dot3*dot4*dot6*M2+16.*a1*dot3*dot4*dot6*M2-
						8.*a12*dot3*dot4*dot6*M2-8.*dot3*dot9*M2+12.*a1*dot3*dot9*M2-4.*a12*dot3*dot9*M2+
						2.*dot12*(dot2+dot4)*(a22*M2-sh)+4.*dot3*dot9*sh+dot10*M*sh-2.*a1*dot10*M*sh-4.*dot3*dot7*M*sh+
						4.*a1*dot3*dot7*M*sh-2.*dot2*dot8*M*sh+6.*a1*dot2*dot8*M*sh-4.*a12*dot2*dot8*M*sh+2.*dot4*dot8*M*sh-
						6.*a1*dot4*dot8*M*sh+4.*a12*dot4*dot8*M*sh-5.*dot13*M2*sh+8.*a1*dot13*M2*sh-4.*a12*dot13*M2*sh+
						4.*dot5*dot6*M2*sh-4.*a1*dot5*dot6*M2*sh+2.*dot13*sqr(sh)-2.*dot5*dot6*sqr(sh)-
						4.*a2*dot11*M*(-2.*a2*dot2*dot3+a1*(-2.*dot3*dot4+a22*dot5*M2-dot5*sh))-
						a2*M*(4.*dot3*dot7+dot13*M-a1*dot13*M)*th+dot13*sh*th)+
				 4.*dot1*M*(4.*dot4*dot6*pow<4,1>(M)-4.*a1*dot4*dot6*pow<4,1>(M)-8.*a12*dot4*dot6*pow<4,1>(M)+12.*pow(a1,3)*dot4*dot6*pow<4,1>(M)-
					    4.*pow(a1,4)*dot4*dot6*pow<4,1>(M)+3.*dot9*pow<4,1>(M)-2.*a1*dot9*pow<4,1>(M)-3.*a12*dot9*pow<4,1>(M)+
					    2.*pow(a1,3)*dot9*pow<4,1>(M)-4.*dot7*pow<5,1>(M)+4.*a1*dot7*pow<5,1>(M)+6.*a12*dot7*pow<5,1>(M)-
					    12.*pow(a1,3)*dot7*pow<5,1>(M)+8.*pow(a1,4)*dot7*pow<5,1>(M)-2.*pow(a1,5)*dot7*pow<5,1>(M)+
					    2.*dot7*pow<3,1>(M)*sh+4.*a1*dot7*pow<3,1>(M)*sh-8.*a12*dot7*pow<3,1>(M)*sh+4.*pow(a1,3)*dot7*pow<3,1>(M)*sh-
					    4.*a1*dot4*dot6*M2*sh+4.*a12*dot4*dot6*M2*sh-dot9*M2*sh-2.*a1*dot9*M2*sh-2.*a1*dot7*M*sqr(sh)-
					    a2*(4.*a2*dot4*dot6-(-3.+a1)*(dot9-2.*dot7*M))*M2*th+(dot9-2.*dot7*M)*sh*th-2.*a2*dot7*M*sqr(th)+
					    4.*dot11*(-(a2*dot2)-a1*dot4)*M*((1.+(-2.+a1)*a12)*M2-a1*sh-a2*th)+
					    2.*dot2*dot6*(a22*(-1.-2.*a1*a2)*pow<4,1>(M)+sh*th+
							  M2*((-1.+2.*a1*a2)*sh+a22*th))))/(a1*a2*sqr(-(a22*M2)+sh)*sqr(M2-th)));
	      diag[3]=rt*mix1*(-4.*M2*(a1*M*(a1*(-1.+2.*a1)*(dot10-2.*dot11*dot5+2.*dot1*dot7-2.*dot4*dot8)*M+
					     4.*dot3*(2.*dot4*dot6-dot9+a12*dot7*M))+
				       (dot10-2.*a1*dot10+2.*(-1.+2.*a1)*dot11*dot5+2.*dot1*dot7-
					4.*a1*(dot1+dot3)*dot7-2.*dot4*dot8+4.*a1*dot4*dot8)*uh))/(a2*pow<3,1>(a12*M2-uh));
	      diag[4]=rt*mix1*(2.*M*(2.*pow(a1,6)*pow<5,1>(M)*(2.*dot11*dot5+2.*dot2*dot8-2.*dot4*dot8+dot5*dot6*M)-
				     2.*pow(a1,5)*pow<5,1>(M)*(6.*dot11*dot5+5.*dot2*dot8-7.*dot4*dot8+dot13*M+2.*dot5*dot6*M)+
				     pow(a1,4)*pow<3,1>(M)*(8.*dot1*dot11*(dot2-dot4)+2.*dot12*(dot2+dot4)*M-4.*dot3*dot9*M+5.*dot13*pow<3,1>(M)+
							    2.*dot1*M*(2.*dot2*dot6-4.*dot4*dot6-dot9+4.*dot7*M)+12.*dot11*dot5*M2+
							    8.*dot2*dot8*M2-18.*dot4*dot8*M2-4.*dot11*dot5*sh-4.*dot2*dot8*sh+4.*dot4*dot8*sh-
							    dot13*M*sh-4.*dot11*dot5*uh-4.*dot2*dot8*uh+4.*dot4*dot8*uh+dot13*M*uh-2.*dot5*dot6*M*uh)+
				     2.*a1*M*(4.*dot11*dot2*dot3*(M2-sh)+2.*dot3*(-2.*dot4*dot6*M+dot7*(M2-sh))*(M2-sh)+
					      2.*dot11*(2.*(dot1+dot3)*(2.*dot2-dot4)+dot5*(M2-sh))*uh+
					      (M*(2.*dot12*(dot2+dot4)+4.*dot2*dot3*dot6-8.*(dot1+dot3)*dot4*dot6-2.*(dot1+2.*dot3)*dot9+
						  (6.*dot1*dot7+8.*dot3*dot7+dot2*dot8-5.*dot4*dot8)*M+2.*(dot13-dot5*dot6)*M2)-
					       (2.*dot1*dot7+4.*dot3*dot7+dot2*dot8-3.*dot4*dot8+2.*dot13*M-2.*dot5*dot6*M)*sh)*uh+
					      (2.*(dot1+dot3)*dot7+dot13*M)*sqr(uh))+
				     2.*pow(a1,3)*pow<3,1>(M)*(4.*dot11*dot3*(dot2-dot4)-2.*dot12*dot2*M-2.*dot12*dot4*M+
							       4.*dot2*dot3*dot6*M-4.*dot3*dot4*dot6*M+4.*dot3*dot9*M-
							       2.*dot13*pow<3,1>(M)+2.*dot5*dot6*pow<3,1>(M)+4.*dot3*dot7*M2-
							       dot2*dot8*M2+5.*dot4*dot8*M2+dot2*dot8*sh-3.*dot4*dot8*sh+2.*dot13*M*sh-
							       2.*dot5*dot6*M*sh+5.*dot2*dot8*uh-7.*dot4*dot8*uh+2.*dot5*dot6*M*uh+
							       2.*dot11*dot5*(-M2+sh+3.*uh)+2.*dot1*(2.*dot11*(-2.*dot2+dot4)+4.*dot4*dot6*M+
												     dot9*M-3.*dot7*M2+dot7*sh-dot7*uh))-
				     uh*(8.*dot11*dot2*dot3*M+4.*dot3*dot7*pow<3,1>(M)-2.*dot4*dot8*pow<3,1>(M)+dot13*pow<4,1>(M)-
					 2.*dot5*dot6*pow<4,1>(M)+2.*dot12*dot2*M2+2.*dot12*dot4*M2-8.*dot3*dot4*dot6*M2-4.*dot3*dot9*M2-
					 2.*dot12*dot2*sh-2.*dot12*dot4*sh+4.*dot3*dot9*sh-4.*dot3*dot7*M*sh+2.*dot4*dot8*M*sh-2.*dot13*M2*sh+
					 4.*dot5*dot6*M2*sh+dot13*sqr(sh)-2.*dot5*dot6*sqr(sh)+4.*dot3*dot7*M*uh+dot13*M2*uh-dot13*sh*uh+
					 2.*dot1*(4.*dot11*dot2*M+2.*dot7*pow<3,1>(M)-4.*dot4*dot6*M2-dot9*M2+dot9*sh-
						  2.*dot7*M*sh+2.*dot2*dot6*(-M2+sh)+2.*dot7*M*uh))-
				     a12*M*(pow<3,1>(M)*(-2.*dot12*(dot2+dot4)+4.*dot3*(2.*dot2*dot6-4.*dot4*dot6+dot9)+
							 2.*(6.*dot3*dot7+dot4*dot8)*M-(dot13-2.*dot5*dot6)*M2)+
					    2.*M*(dot12*(dot2+dot4)-2.*dot3*(2.*dot2*dot6-2.*dot4*dot6+dot9)-(6.*dot3*dot7+dot4*dot8)*M+
						  (dot13-2.*dot5*dot6)*M2)*sh-(dot13-2.*dot5*dot6)*M*sqr(sh)+
					    8.*dot11*dot3*(2.*dot2*M2-dot4*M2-dot2*sh+dot4*sh)+4.*dot11*(2.*dot2*dot3-2.*dot3*dot4+3.*dot5*M2-dot5*sh)*uh+
					    2.*(M*(dot12*(dot2+dot4)-2.*dot3*(-2.*dot2*dot6+2.*dot4*dot6+dot9)+
						   (6.*dot3*dot7+4.*dot2*dot8-9.*dot4*dot8)*M+2.*dot13*M2)+2.*(-dot2+dot4)*dot8*sh)*uh+dot13*M*sqr(uh)+
					    2.*dot1*(pow<3,1>(M)*(2.*dot2*dot6+4.*dot4*dot6+dot9-2.*dot7*M)-4.*dot11*dot2*M2-
						     M*(2.*dot2*dot6+dot9-2.*dot7*M)*sh+4.*dot11*(dot2-dot4)*uh+
						     M*(2.*dot2*dot6-4.*dot4*dot6-dot9+2.*dot7*M)*uh))))/(a1*a2*sqr(-(a22*M2)+sh)*sqr(-(a12*M2)+uh));
	      // diagrams 3P1
	      diag[0]+= mix2*8.*M2*(2.*(-1.+2.*a1)*a2*dot9*dot3*M-
				    a2*M*(2.*dot2*(2.*(-1.+2.*a1)*dot3*dot6+(2.+a1-3.*a12)*dot8*M)-a2*M*((2.+3.*a1)*dot10-2.*a2*(-(dot3*dot7)+a1*dot13*M)))+
				    ((-2.+a1)*dot10-2.*a2*dot3*dot7-2.*(-2.+a1)*dot2*dot8+2.*a1*a2*dot13*M)*sh-
				    2.*a2*dot1*(M*(dot9-2.*a1*dot9+2.*(-1.+2.*a1)*dot2*dot6-a22*dot7*M)+dot7*sh))/(a1*pow<3,1>(a22*M2-sh));
	      diag[1]+=-mix2*8.*M2*(2.*pow(a1,6)*(-dot2+dot4)*dot8*pow<4,1>(M)-pow(a1,5)*(dot1*dot7+2.*(-dot2+dot4)*dot8)*pow<4,1>(M)+
				    pow(a1,4)*M2*(-2.*dot12*dot2*M-2.*dot12*dot4*M+dot5*dot6*pow<3,1>(M)-2.*dot13*pow<3,1>(M)-dot10*M2-
						  2.*dot3*dot7*M2-3.*dot2*dot8*M2+5.*dot4*dot8*M2+dot10*th+2.*dot3*dot7*th+3.*dot2*dot8*th-
						  5.*dot4*dot8*th-dot5*dot6*M*th+2.*dot13*M*th+dot11*(4.*dot2*dot3-4.*dot3*dot4+dot5*M2-dot5*th)+
						  dot1*(2.*dot11*(dot2-dot4)-M*(-2.*dot2*dot6+2.*dot4*dot6+dot7*M)+2.*dot7*th)+4.*(dot2-dot4)*dot8*uh)+
				    a12*(-2.*dot11*(dot2*dot3-dot3*dot4-2.*dot5*M2)*(M2-th)+
					 dot1*M*(-2.*dot11*dot2*M-(8.*dot2*dot6-4.*dot4*dot6+3.*dot7*M)*(M2-th)+dot9*(-3.*M2+4.*th))+
					 2.*M*(3.*dot12*dot2*(M2-th)+dot12*dot4*(-4.*M2+sh+th)+
					       (M2-th)*(-3.*dot9*dot3-dot2*dot3*dot6+dot3*dot4*dot6-2.*dot10*M+dot5*dot6*M2-dot13*M2-dot5*dot6*th+dot13*th))+
					 (2.*dot12*dot2*M+(dot10+2.*dot3*dot7-(dot2+dot4)*dot8-dot5*dot6*M+2.*dot13*M)*(M2-th))*uh-
					 2.*dot1*((dot2-dot4)*(dot11+dot6*M)+dot7*th)*uh+dot11*(-4.*dot2*dot3+4.*dot3*dot4+dot5*(-M2+th))*uh+
					 2.*(-dot2+dot4)*dot8*sqr(uh))+
				    dot1*uh*(2.*dot11*dot2-dot9*M+dot7*(-M2+th+uh))+
				    pow(a1,3)*M*(5.*dot10*pow<3,1>(M)+3.*dot3*dot7*pow<3,1>(M)+5.*dot2*dot8*pow<3,1>(M)-7.*dot4*dot8*pow<3,1>(M)+
						 dot5*dot6*pow<4,1>(M)+8.*dot12*dot4*M2+4.*dot2*dot3*dot6*M2-4.*dot3*dot4*dot6*M2+2.*dot12*dot2*th-
						 2.*dot12*dot4*th-4.*dot2*dot3*dot6*th+4.*dot3*dot4*dot6*th-5.*dot10*M*th-3.*dot3*dot7*M*th-
						 5.*dot2*dot8*M*th+7.*dot4*dot8*M*th-dot5*dot6*M2*th+
						 dot11*M*(-4.*dot2*dot3+2.*dot1*dot4+4.*dot3*dot4-5.*dot5*M2+5.*dot5*th)+
						 4.*(-dot2+dot4)*dot8*M*uh+dot1*(7.*dot7*pow<3,1>(M)-dot9*M2-2.*dot4*dot6*M2+
										 2.*dot2*dot6*(M2-2.*th)+4.*dot4*dot6*th-7.*dot7*M*th+2.*dot7*M*uh))+
				    a1*(3.*dot9*dot3*pow<3,1>(M)+2.*dot12*dot4*pow<3,1>(M)-dot10*pow<4,1>(M)+dot3*dot7*pow<4,1>(M)-
					2.*dot12*dot4*M*sh+dot10*M2*sh-3.*dot9*dot3*M*th+dot10*M2*th-2.*dot3*dot7*M2*th-dot10*sh*th+dot3*dot7*sqr(th)+
					2.*dot11*dot2*dot3*(-M2+th)-M*(2.*dot12*dot2+M*(3.*dot3*dot7+dot2*dot8-3.*dot4*dot8+dot5*dot6*M))*uh+
					(3.*dot3*dot7+dot2*dot8-3.*dot4*dot8+dot5*dot6*M)*th*uh+dot11*(4.*dot2*dot3-4.*dot3*dot4+dot5*M2-dot5*th)*uh+
					2.*(dot2-dot4)*dot8*sqr(uh)+
					dot1*(2.*(M2-th)*(-2.*dot11*dot2+2.*dot9*M+dot7*M2-dot7*th)+
					      (-2.*dot11*dot4+M*(dot9+2.*dot2*dot6-2.*dot4*dot6-3.*dot7*M)+3.*dot7*th)*uh-dot7*sqr(uh))))
		/(a1*a2*sqr(M2-th)*sqr(-(a12*M2)+uh));
	      diag[2]+= mix2*8.*M*(-(a1*M*(a22*M2-sh)*(dot1*(-2.*a2*dot11*dot2-2.*a1*dot11*dot4+M*(dot9-a2*(2.*dot2*dot6-2.*dot4*dot6+dot7*M-a1*dot7*M))+dot7*sh)-
						       2.*a2*(dot2-dot4)*(2.*dot11*dot3+M*(-dot12+a22*dot8*M)-dot8*sh)))+
				   (pow(a1,5)*(dot5*dot6-dot13)*pow<4,1>(M)+
				    pow(a1,4)*pow<3,1>(M)*(-dot10+dot11*dot5+2.*dot3*dot7-5.*dot2*dot8+3.*dot4*dot8-6.*dot5*dot6*M+5.*dot13*M)+
				    M*(2.*dot11*(dot1-dot3)*dot4+(-3.*dot9*dot3+dot12*(dot2+3.*dot4)+2.*dot3*(-dot2+dot4)*dot6-
								  2.*dot1*(dot9+(dot2+dot4)*dot6))*M+2.*(-(dot5*dot6)+dot13)*pow<3,1>(M)-
				       (dot10-dot11*dot5+dot3*dot7-dot2*dot8+dot4*dot8)*M2)+
				    (2.*dot9*dot3-dot12*(dot2+dot4)+dot1*(dot9+2.*dot4*dot6)+(dot10-dot11*dot5+dot3*dot7-dot2*dot8+dot4*dot8)*M+3.*(dot5*dot6-dot13)*M2)*sh+
				    (-(dot5*dot6)+dot13)*sqr(sh)+
				    pow(a1,3)*M2*(2.*dot9*dot3+dot12*(-3.*dot2+dot4)+4.*dot2*dot3*dot6-4.*dot3*dot4*dot6+6.*dot10*M-
						  8.*dot11*dot5*M-5.*dot3*dot7*M+14.*dot2*dot8*M-8.*dot4*dot8*M+
						  dot1*(dot9+4.*dot2*dot6-2.*dot4*dot6+5.*dot7*M)+14.*dot5*dot6*M2-
						  11.*dot13*M2-2.*dot5*dot6*sh+2.*dot13*sh)+
				    a1*(M*(-2.*dot11*(dot2*dot3+(dot1-2.*dot3)*dot4)+
					   (4.*dot1*dot9+9.*dot9*dot3-5.*dot12*(dot2+dot4)+8.*dot3*(dot2-dot4)*dot6+2.*dot1*(5.*dot2+dot4)*dot6)*M+
					   (9.*dot5*dot6-8.*dot13)*pow<3,1>(M)+(6.*dot10-8.*dot11*dot5+5.*dot1*dot7+dot3*dot7+2.*dot2*dot8)*M2)+
					(-2.*dot9*dot3+dot12*(dot2+dot4)+(-2.*dot10+2.*dot11*dot5+dot3*dot7)*M-dot1*(dot9+2.*dot4*dot6+dot7*M)+
					 M*(-2.*dot4*dot8-8.*dot5*dot6*M+7.*dot13*M))*sh+(dot5*dot6-dot13)*sqr(sh))+
				    a12*M*(M*(-8.*dot9*dot3+dot12*(7.*dot2+dot4)+10.*dot3*(-dot2+dot4)*dot6-10.*dot10*M+3.*dot3*dot7*M+
					      6.*(-2.*dot2+dot4)*dot8*M-dot1*(3.*dot9+12.*dot2*dot6-2.*dot4*dot6+10.*dot7*M)+
					      (-16.*dot5*dot6+13.*dot13)*M2)+(dot10-2.*dot3*dot7+dot2*dot8+dot4*dot8+7.*dot5*dot6*M-
									      6.*dot13*M)*sh+dot11*(2.*dot2*dot3-2.*dot3*dot4+14.*dot5*M2-dot5*sh)))*(M2-th)+
				   a1*(-(dot12*dot2)+dot9*(dot1+dot3)+M*(-((-2.+a1)*dot10)+(-3.+a1)*dot11*dot5+dot2*dot8+
									 (5.+(-4.+a1)*a1)*(dot5*dot6-dot13)*M)+(-(dot5*dot6)+dot13)*sh)*sqr(M2-th)-
				   (-(dot12*dot2)+dot9*(dot1+dot3)+M*(dot10-2.*dot11*dot5+dot2*dot8+2.*dot5*dot6*M-2.*dot13*M)+
				    (-(dot5*dot6)+dot13)*sh)*sqr(M2-th))/(a1*a2*sqr(-(a22*M2)+sh)*sqr(M2-th));
	      diag[3]+=-mix2*8.*M2*(a1*M*(2.*(1.-2.*a1)*dot9*dot3-2.*dot3*((2.-4.*a1)*dot4*dot6+a12*dot7*M)+
					  a1*M*((-5.+3.*a1)*(dot10-2.*dot11*dot5+2.*dot1*dot7-2.*dot4*dot8)-2.*a1*a2*(2.*dot5*dot6-dot13)*M))+
				    ((1.+a1)*dot10-2.*(1.+a1)*dot11*dot5+2.*(dot1+a1*dot1+a1*dot3)*dot7-
				     2.*(1.+a1)*dot4*dot8+2.*a1*a2*(2.*dot5*dot6-dot13)*M)*uh)/(a2*pow<3,1>(a12*M2-uh));
	      diag[4]+= mix2*8.*M2*(-(pow(a1,5)*pow<4,1>(M)*(3.*dot10-5.*dot11*dot5+6.*(dot1+dot3)*dot7+3.*dot2*dot8-9.*dot4*dot8+3.*dot5*dot6*M-3.*dot13*M))+
				    pow(a1,6)*pow<4,1>(M)*(dot10-dot11*dot5+2.*(dot1+dot3)*dot7-(dot2+dot4)*dot8+dot5*dot6*M-dot13*M)+
				    uh*((dot1+dot3)*(2.*dot11*dot2-2.*dot4*dot6*M+dot7*(M2-sh))-(dot10+(dot1+dot3)*dot7-2.*dot2*dot8)*uh)+
				    pow(a1,3)*M*(-7.*dot10*pow<3,1>(M)-11.*dot1*dot7*pow<3,1>(M)-6.*dot3*dot7*pow<3,1>(M)-
						 3.*dot2*dot8*pow<3,1>(M)+25.*dot4*dot8*pow<3,1>(M)-dot5*dot6*pow<4,1>(M)+dot13*pow<4,1>(M)+
						 6.*dot12*dot2*M2-10.*dot12*dot4*M2-10.*dot1*dot2*dot6*M2-6.*dot2*dot3*dot6*M2+16.*dot1*dot4*dot6*M2+
						 10.*dot3*dot4*dot6*M2-2.*dot12*dot2*sh+2.*dot12*dot4*sh+4.*dot1*dot2*dot6*sh+4.*dot2*dot3*dot6*sh-
						 4.*dot1*dot4*dot6*sh-4.*dot3*dot4*dot6*sh+3.*dot10*M*sh+5.*dot1*dot7*M*sh+4.*dot3*dot7*M*sh+
						 3.*dot2*dot8*M*sh-13.*dot4*dot8*M*sh+dot5*dot6*M2*sh-dot13*M2*sh+
						 dot11*M*(-6.*dot2*dot3+2.*(dot1+5.*dot3)*dot4+dot5*(13.*M2-5.*sh-3.*uh))+
						 (2.*(dot2-dot4)*(dot12-2.*(dot1+dot3)*dot6)+
						  (dot10+3.*dot1*dot7+4.*dot3*dot7-7.*dot2*dot8+9.*dot4*dot8)*M+3.*(dot5*dot6-dot13)*M2)*uh)-
				    a1*(2.*dot11*(M2-sh)*(dot2*dot3-2.*dot3*dot4-dot5*M2+dot5*sh)+
					(M2-sh)*(M*(2.*dot12*dot4-2.*(2.*dot1+dot3)*dot4*dot6+(dot10+2.*dot1*dot7+dot3*dot7-2.*dot4*dot8)*M)-
						 (dot10+2.*dot1*dot7+dot3*dot7-2.*dot4*dot8)*sh)+
					dot11*(2.*(dot1+dot3)*dot4+dot5*(M2-sh))*uh+
					(M*(2.*dot12*dot4+2.*(dot1+dot3)*(dot2-4.*dot4)*dot6+
					    (-dot10+dot1*dot7+dot2*dot8-3.*dot4*dot8)*M+
					    (-(dot5*dot6)+dot13)*M2)+
					 (dot10+(dot1+2.*dot3)*dot7-(dot2+dot4)*dot8+(dot5*dot6-dot13)*M)*sh)*uh-
					(dot10+(dot1+dot3)*dot7-2.*dot2*dot8)*sqr(uh))+
				    pow(a1,4)*M2*(4.*dot11*dot3*(dot2-dot4)-4.*dot12*dot2*M+4.*dot12*dot4*M+4.*dot2*dot3*dot6*M-
						  4.*dot3*dot4*dot6*M+3.*dot5*dot6*pow<3,1>(M)-3.*dot13*pow<3,1>(M)+6.*dot10*M2+
						  8.*dot3*dot7*M2+7.*dot2*dot8*M2-23.*dot4*dot8*M2-dot10*sh-2.*dot3*dot7*sh-
						  3.*dot2*dot8*sh+5.*dot4*dot8*sh-dot5*dot6*M*sh+dot13*M*sh-
						  (dot10+2.*dot3*dot7-5.*dot2*dot8+3.*dot4*dot8+dot5*dot6*M-dot13*M)*uh+
						  dot11*dot5*(-11.*M2+sh+uh)+2.*dot1*(dot11*(dot2-dot4)+3.*dot2*dot6*M-3.*dot4*dot6*M+5.*dot7*M2-dot7*sh-dot7*uh))-
				    a12*(2.*(dot2-4.*dot4)*(dot12-dot3*dot6)*pow<3,1>(M)+(-4.*dot10-3.*dot3*dot7+12.*dot4*dot8)*pow<4,1>(M)+
					 M*(-2.*dot12*(dot2-2.*dot4)+2.*dot3*(dot2-3.*dot4)*dot6+(4.*dot10+3.*dot3*dot7-12.*dot4*dot8)*M)*sh+
					 dot1*M*(2.*dot11*dot2*M+(-4.*dot2*dot6+14.*dot4*dot6-7.*dot7*M)*M2+(4.*dot2*dot6-8.*dot4*dot6+7.*dot7*M)*sh)+
					 M*(2.*dot12*(dot2-2.*dot4)+2.*dot3*(-3.*dot2+5.*dot4)*dot6+
					    (dot10+3.*dot3*dot7-3.*dot2*dot8+9.*dot4*dot8)*M+3.*(dot5*dot6-dot13)*M2)*uh+
					 (-dot10-2.*dot3*dot7+(dot2+dot4)*dot8+(-(dot5*dot6)+dot13)*M)*sh*uh+
					 dot1*(2.*dot11*(dot2-dot4)-6.*dot2*dot6*M+10.*dot4*dot6*M+dot7*M2-2.*dot7*sh)*uh+
					 dot11*(2.*dot3*dot4*(5.*M2-sh-uh)+2.*dot2*dot3*(-2.*M2+sh+uh)+
						dot5*(8.*pow<4,1>(M)+sh*uh-M2*(8.*sh+3.*uh)))))/(a1*a2*sqr(-(a22*M2)+sh)*sqr(-(a12*M2)+uh));
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
    	    diag[0]= rt*mix1*(-2.*M*(-2.*pow(a1,5)*pow<5,1>(M)*(-6.*dot5*dot7+7.*dot1*dot8+5.*dot3*dot8-dot11*M+4.*dot5*dot6*M)+
    	    			     2.*pow(a1,6)*pow<5,1>(M)*(2.*(dot1+dot3)*dot8+dot5*(-2.*dot7+dot6*M))+
    	    			     2.*a1*M*(2.*dot4*(-2.*dot3*dot7+2.*dot1*dot6*M+dot9*(M2-th))*(M2-th)+
    	    				      (M*(-2.*dot12*dot3+2.*dot10*dot4-2.*dot5*dot7*M+2.*M*(dot4*dot9-dot11*M+dot5*dot6*M)+
    	    					  dot3*(4.*dot4*dot6+dot8*M))+
    	    				       (-(dot3*dot8)-2.*dot4*dot9+2.*dot11*M+2.*dot5*(dot7-dot6*M))*th+
    	    				       dot1*(2.*dot12*M+4.*dot2*(dot7-2.*dot6*M)+5.*dot8*M2-3.*dot8*th)+
    	    				       2.*dot2*(4.*dot3*dot7-M*(dot10+3.*dot9*M)+dot9*th))*uh-
    	    				      (2.*dot2*dot9+(dot11-2.*dot5*dot6)*M)*sqr(uh))-
    	    			     2.*pow(a1,3)*pow<3,1>(M)*(8.*dot2*dot3*dot7+12.*dot3*dot4*dot7-2.*dot10*dot2*M-2.*dot12*dot3*M+
    	    						       2.*dot10*dot4*M-4.*dot3*dot4*dot6*M-2.*dot11*pow<3,1>(M)+2.*dot5*dot6*pow<3,1>(M)-
    	    						       2.*dot5*dot7*M2+dot3*dot8*M2-6.*dot2*dot9*M2-10.*dot4*dot9*M2+2.*dot5*dot7*th-
    	    						       dot3*dot8*th+2.*dot2*dot9*th+2.*dot4*dot9*th+2.*dot11*M*th-2.*dot5*dot6*M*th+
    	    						       (6.*dot5*dot7-5.*dot3*dot8-2.*(dot2+dot4)*dot9-2.*dot5*dot6*M)*uh+
    	    						       dot1*(8.*dot4*dot7+2.*dot12*M-12.*dot4*dot6*M+4.*dot2*(dot7-2.*dot6*M)+
    	    							     5.*dot8*M2-3.*dot8*th-7.*dot8*uh))-
    	    			     a12*M*(-2.*dot1*dot12*pow<3,1>(M)+2.*dot12*dot3*pow<3,1>(M)-2.*dot10*dot4*pow<3,1>(M)+24.*dot1*dot4*dot6*pow<3,1>(M)+
    	    				    4.*dot3*dot4*dot6*pow<3,1>(M)-2.*dot1*dot8*pow<4,1>(M)+16.*dot4*dot9*pow<4,1>(M)+
    	    				    dot11*pow<5,1>(M)-8.*dot1*dot4*dot7*M2-24.*dot3*dot4*dot7*M2+8.*dot1*dot4*dot7*th+
    	    				    8.*dot3*dot4*dot7*th+2.*dot1*dot12*M*th-2.*dot12*dot3*M*th+2.*dot10*dot4*M*th-8.*dot1*dot4*dot6*M*th-
    	    				    4.*dot3*dot4*dot6*M*th-2.*dot11*pow<3,1>(M)*th+2.*dot1*dot8*M2*th-16.*dot4*dot9*M2*th+dot11*M*sqr(th)+
    	    				    2.*M*(dot1*dot12-dot12*dot3+dot10*dot4+2.*dot3*dot4*dot6-6.*dot5*dot7*M+
    	    					  (9.*dot1*dot8+4.*dot3*dot8+4.*dot4*dot9)*M-2.*(dot11-2.*dot5*dot6)*M2)*uh+
    	    				    4.*(dot5*dot7-(dot1+dot3)*dot8)*th*uh-(dot11-2.*dot5*dot6)*M*sqr(uh)+
    	    				    2.*dot2*(M*(-2.*dot3*M*(2.*dot7+dot6*M)+(dot10+4.*dot1*dot6+2.*dot9*M)*M2+2.*dot3*dot6*th-
    	    						(dot10+2.*dot9*M)*th)+(4.*(dot1+dot3)*dot7-(dot10+2.*(2.*dot1+dot3)*dot6)*M-2.*dot9*M2)*uh))+
    	    			     pow(a1,4)*pow<3,1>(M)*(8.*dot2*dot3*dot7+8.*dot3*dot4*dot7-2.*dot10*dot2*M-2.*dot12*dot3*M+2.*dot10*dot4*M-
    	    						    4.*dot2*dot3*dot6*M-4.*dot3*dot4*dot6*M-5.*dot11*pow<3,1>(M)+10.*dot5*dot6*pow<3,1>(M)-
    	    						    12.*dot5*dot7*M2+8.*dot3*dot8*M2-8.*dot2*dot9*M2-8.*dot4*dot9*M2+4.*dot5*dot7*th-
    	    						    4.*dot3*dot8*th+dot11*M*th-2.*dot5*dot6*M*th+4.*dot5*dot7*uh-4.*dot3*dot8*uh-dot11*M*uh+
    	    						    2.*dot1*(dot12*M+4.*dot2*(dot7-dot6*M)+4.*dot4*(dot7-dot6*M)+9.*dot8*M2-2.*dot8*(th+uh)))+
    	    			     uh*((M2-th)*(2.*dot12*dot3-2.*dot1*(dot12+dot8*M)+dot11*(M2-th+uh)-2.*(dot10*dot4+2.*dot3*dot4*dot6+dot5*dot6*uh))+
    	    				 dot2*(-4.*dot3*M*(2.*dot7+dot6*M)+2.*dot10*(M2-th)+4.*dot3*dot6*th+
    	    				       4.*M*(2.*dot1*dot6*M+dot9*(M2-th+uh))))))/(a1*a2*sqr(-(a22*M2)+th)*sqr(-(a12*M2)+uh));
    	    diag[1]=rt*mix1*(-4.*M2*(4.*a2*dot10*dot2*M-(-1.+2.*a1)*(dot13-2.*(dot5*dot7+dot3*dot8))*(a22*M2-th)-
    	    			     4.*a2*dot2*(-2.*dot3*dot6*M+a22*dot9*M2-dot9*th)))/(a1*pow<3,1>(a22*M2-th));
    	    diag[2]=rt*mix1*(-2.*M*(4.*pow(a1,5)*(-dot2+dot4)*dot9*pow<5,1>(M)+
    	    			    2.*pow(a1,4)*pow<3,1>(M)*(-4.*dot3*dot4*dot7+dot5*dot6*pow<3,1>(M)+4.*dot2*dot3*(dot7+dot6*M)-2.*dot5*dot7*M2+
    	    						      2.*dot3*dot8*M2-4.*dot4*dot9*M2+2.*dot5*dot7*sh-2.*dot3*dot8*sh-dot5*dot6*M*sh+
    	    						      2.*dot1*(-2.*dot4*dot7+2.*dot4*dot6*M+2.*dot2*(dot7-dot6*M)+dot8*M2-dot8*sh))+
    	    			    (M2-th)*(-(dot13*pow<3,1>(M))+2.*dot5*dot7*pow<3,1>(M)+2.*dot3*dot8*pow<3,1>(M)-dot11*pow<4,1>(M)-
    	    				     2.*dot12*dot3*M2+2.*dot10*dot4*M2-4.*dot3*dot4*dot6*M2+2.*dot1*(dot12+dot8*M)*(M2-sh)+2.*dot12*dot3*sh-
    	    				     2.*dot10*dot4*sh+dot13*M*sh-2.*dot5*dot7*M*sh-2.*dot3*dot8*M*sh+dot11*M2*sh+4.*dot3*dot4*dot6*th+
    	    				     dot11*M2*th-dot11*sh*th+4.*dot3*dot4*dot6*uh-(dot11-2.*dot5*dot6)*(M2-sh)*uh)-
    	    			    2.*dot2*(M2-sh)*(-4.*dot3*dot7*M+2.*dot3*dot6*(M2+th)+dot10*(M2+sh+uh)+2.*M*(2.*dot1*dot6*M+dot9*uh))+
    	    			    2.*pow(a1,3)*pow<3,1>(M)*(8.*dot1*dot4*dot7+12.*dot3*dot4*dot7-2.*dot10*dot4*M-12.*dot1*dot4*dot6*M+dot11*pow<3,1>(M)-
    	    						      4.*dot5*dot6*pow<3,1>(M)+dot13*M2+4.*dot5*dot7*M2-7.*dot1*dot8*M2-7.*dot3*dot8*M2+
    	    						      2.*dot4*dot9*M2-dot13*sh-4.*dot5*dot7*sh+7.*dot1*dot8*sh+7.*dot3*dot8*sh-dot11*M*sh+
    	    						      4.*dot5*dot6*M*sh-2.*dot4*dot9*th-2.*dot4*dot9*uh+
    	    						      2.*dot2*(dot10*M+2.*dot1*dot6*M-2.*dot3*(dot7+3.*dot6*M)+2.*dot9*(-M2+sh+uh)))-
    	    			    2.*a1*M*(-2.*dot12*dot3*pow<3,1>(M)+4.*dot10*dot4*pow<3,1>(M)-4.*dot3*dot4*dot6*pow<3,1>(M)-2.*dot13*pow<4,1>(M)+
    	    				     2.*dot5*dot7*pow<4,1>(M)+5.*dot3*dot8*pow<4,1>(M)-2.*dot11*pow<5,1>(M)+2.*dot5*dot6*pow<5,1>(M)-
    	    				     4.*dot3*dot4*dot7*M2+2.*dot12*dot3*M*sh-2.*dot10*dot4*M*sh+2.*dot11*pow<3,1>(M)*sh-2.*dot5*dot6*pow<3,1>(M)*sh+
    	    				     2.*dot13*M2*sh-2.*dot5*dot7*M2*sh-5.*dot3*dot8*M2*sh+4.*dot3*dot4*dot7*th-2.*dot10*dot4*M*th+
    	    				     4.*dot3*dot4*dot6*M*th+2.*dot11*pow<3,1>(M)*th-2.*dot5*dot6*pow<3,1>(M)*th+dot13*M2*th-3.*dot3*dot8*M2*th-
    	    				     dot13*sh*th+3.*dot3*dot8*sh*th-2.*dot11*M*sh*th+2.*dot5*dot6*M*sh*th+4.*dot2*dot3*(dot6*M*(sh-th)+dot7*(3.*M2-3.*sh-uh))+
    	    				     M*(4.*dot3*dot4*dot6+2.*dot4*dot9*M-(dot11-2.*dot5*dot6)*(M2-sh))*uh-2.*dot4*dot9*th*uh-
    	    				     2.*dot2*(2.*M2-2.*sh-uh)*(dot10*M+dot9*uh)+
    	    				     dot1*(4.*dot4*dot6*pow<3,1>(M)+5.*dot8*pow<4,1>(M)-5.*dot8*M2*sh+2.*dot12*(pow<3,1>(M)-M*sh)-
    	    					   4.*dot4*dot6*M*th-3.*dot8*M2*th+3.*dot8*sh*th+4.*dot2*((dot7-3.*dot6*M)*(M2-sh)+dot6*M*uh)))+
    	    			    a12*M*(-2.*dot12*dot3*pow<3,1>(M)+10.*dot10*dot4*pow<3,1>(M)-4.*dot3*dot4*dot6*pow<3,1>(M)-5.*dot13*pow<4,1>(M)-
    	    				   2.*dot5*dot7*pow<4,1>(M)+18.*dot3*dot8*pow<4,1>(M)-5.*dot11*pow<5,1>(M)+10.*dot5*dot6*pow<5,1>(M)-
    	    				   24.*dot3*dot4*dot7*M2+2.*dot12*dot3*M*sh-2.*dot10*dot4*M*sh+5.*dot11*pow<3,1>(M)*sh-10.*dot5*dot6*pow<3,1>(M)*sh+
    	    				   5.*dot13*M2*sh+2.*dot5*dot7*M2*sh-18.*dot3*dot8*M2*sh+8.*dot3*dot4*dot7*th+4.*dot3*dot4*dot6*M*th+
    	    				   dot11*pow<3,1>(M)*th-2.*dot5*dot6*pow<3,1>(M)*th+4.*dot5*dot7*M2*th-4.*dot3*dot8*M2*th-4.*dot5*dot7*sh*th+
    	    				   4.*dot3*dot8*sh*th-dot11*M*sh*th+2.*dot5*dot6*M*sh*th+
    	    				   2.*dot1*((M2-sh)*(M*(dot12+9.*dot8*M)-2.*dot8*th)-4.*dot4*(-3.*dot6*pow<3,1>(M)+dot7*M2-dot7*th+dot6*M*th)+
    	    					    4.*dot2*(dot7-dot6*M)*(2.*M2-2.*sh-uh))+M*(4.*dot3*dot4*dot6+8.*dot4*dot9*M-(dot11-2.*dot5*dot6)*(M2-sh))*uh+
    	    				   2.*dot2*(M*(dot10+2.*dot9*M)*(M2-sh)+2.*dot3*(4.*dot7*M2+dot6*M*(5.*M2+sh-2.*th)-2.*dot7*(2.*sh+uh))))))
	      /(a1*a2*sqr(M2-sh)*sqr(-(a22*M2)+th));
    	    diag[3]=rt*mix1*(-4.*M2*(a1*M*(-4.*dot4*(dot10-2.*dot1*dot6)+a1*(dot13-2.*a1*dot13-2.*((-1.+2.*a1)*dot1*dot8+(dot2-2.*a1*dot2+dot4)*dot9))*M)+
    	    			     ((-1.+2.*a1)*(dot13+2.*dot1*dot8)+2.*(dot2-2.*a1*dot2+dot4)*dot9)*uh))/(a2*pow<3,1>(a12*M2-uh));
    	    diag[4]=rt*mix1*(-2.*M*(4.*pow(a1,5)*(dot2-dot4)*dot9*pow<5,1>(M)+
    	    			    2.*pow(a1,4)*pow<3,1>(M)*(-4.*dot2*dot3*dot7-12.*dot3*dot4*dot7+8.*dot3*dot4*dot6*M+dot5*dot6*pow<3,1>(M)-
    	    						      2.*dot5*dot7*M2+2.*dot3*dot8*M2+12.*dot4*dot9*M2+2.*dot5*dot7*sh-2.*dot3*dot8*sh-
    	    						      dot5*dot6*M*sh+2.*dot1*(-2.*dot2*dot7-6.*dot4*dot7+2.*dot2*dot6*M+6.*dot4*dot6*M+dot8*M2-dot8*sh))-
    	    			    2.*pow(a1,3)*pow<3,1>(M)*(-24.*dot3*dot4*dot7+2.*dot10*dot4*M+8.*dot3*dot4*dot6*M-dot11*pow<3,1>(M)+2.*dot5*dot6*pow<3,1>(M)-
    	    						      dot13*M2-2.*dot5*dot7*M2+dot3*dot8*M2+16.*dot4*dot9*M2+dot13*sh+2.*dot5*dot7*sh-dot3*dot8*sh+dot11*M*sh-
    	    						      2.*dot5*dot6*M*sh+dot1*(-4.*dot2*dot7-12.*dot4*dot7+4.*dot2*dot6*M+12.*dot4*dot6*M+dot8*M2-dot8*sh)-
    	    						      6.*dot4*dot9*th+2.*dot4*dot9*uh+2.*dot2*(-4.*dot3*dot7+M*(dot10+dot9*M)+2.*dot9*uh))+
    	    			    2.*a1*M*(4.*dot3*dot4*dot7*(M2-th)-2.*dot4*dot9*sqr(M2-th)-dot3*(8.*dot2*dot7+12.*dot4*dot7-dot8*M2+dot8*sh)*uh+
    	    				     (-(dot11*pow<3,1>(M))+2.*dot5*dot6*pow<3,1>(M)-dot13*M2-2.*dot5*dot7*M2+2.*dot2*dot9*M2+4.*dot4*dot9*M2+dot13*sh+
    	    				      2.*dot5*dot7*sh+dot11*M*sh-2.*dot5*dot6*M*sh+dot1*(-4.*dot2*dot7-4.*dot4*dot7+4.*dot2*dot6*M+4.*dot4*dot6*M+dot8*M2-dot8*sh)-
    	    				      2.*dot4*dot9*th)*uh+2.*(dot2+dot4)*dot9*sqr(uh)+2.*dot10*M*(dot2*uh+dot4*(M2-sh+uh)))+
    	    			    uh*(8.*dot3*dot4*dot7*M+dot13*pow<3,1>(M)-2.*dot4*dot9*pow<3,1>(M)+dot11*pow<4,1>(M)-2.*dot1*dot12*M2-6.*dot10*dot4*M2+
    	    				2.*dot12*dot3*(M2-sh)+2.*dot1*dot12*sh+2.*dot10*dot4*sh-dot13*M*sh-dot11*M2*sh+2.*dot4*dot9*M*th-dot11*M2*th+dot11*sh*th+
    	    				4.*dot3*dot4*dot6*(M2-th-uh)-2.*dot4*dot9*M*uh+(dot11-2.*dot5*dot6)*(M2-sh)*uh-
    	    				2.*dot2*(-4.*dot3*dot7*M+dot10*(M2+sh)+dot9*M*(M2-th+uh)+2.*dot3*dot6*(-M2+th+uh)))+
    	    			    a12*M*(-2.*dot12*dot3*pow<3,1>(M)+6.*dot10*dot4*pow<3,1>(M)+4.*dot3*dot4*dot6*pow<3,1>(M)-
    	    				   dot13*pow<4,1>(M)+18.*dot4*dot9*pow<4,1>(M)-dot11*pow<5,1>(M)-32.*dot3*dot4*dot7*M2+
    	    				   2.*dot12*dot3*M*sh-2.*dot10*dot4*M*sh+dot11*pow<3,1>(M)*sh+dot13*M2*sh+8.*dot3*dot4*dot7*th-
    	    				   4.*dot3*dot4*dot6*M*th+dot11*pow<3,1>(M)*th-18.*dot4*dot9*M2*th-dot11*M*sh*th-
    	    				   8.*dot1*dot4*(dot7-dot6*M)*(M2-th-2.*uh)+8.*dot1*dot2*dot7*uh+4.*dot5*dot7*(M2-sh)*uh+
    	    				   M*(-(M*(6.*dot4*dot9+dot11*M))+dot11*sh)*uh+4.*dot3*(4.*dot4*dot7-dot4*dot6*M-dot8*M2+dot8*sh)*uh+
    	    				   2.*dot1*(M2-sh)*(dot12*M-2.*dot8*uh)+2.*dot2*M*(dot10*(M2+sh)-4.*dot1*dot6*uh+dot9*M*(M2-th+uh))+
    	    				   4.*dot2*dot3*(2.*dot7*(-M2+uh)+dot6*M*(-M2+th+uh)))))/(a1*a2*sqr(M2-sh)*sqr(-(a12*M2)+uh));
    	    // 3P1 diagrams
	    diag[0]+=mix2*8.*M2*(pow(a1,6)*pow<4,1>(M)*(dot13-dot1*dot8+dot3*dot8-2.*dot2*dot9-dot11*M+dot5*(dot7+dot6*M))-
				 pow(a1,5)*pow<4,1>(M)*(dot13-3.*(3.*dot1*dot8+dot3*dot8+dot11*M)+dot5*(dot7+3.*dot6*M))+
				 uh*(dot2*(-2.*dot3*dot7+2.*dot1*dot6*M-dot9*M2+dot9*sh)+(dot13-2.*(dot5*dot7+dot3*dot8))*uh)+
				 a1*(dot13*pow<4,1>(M)-2.*dot2*dot9*pow<4,1>(M)-dot4*dot9*pow<4,1>(M)-2.*dot3*dot4*dot7*M2-
				     2.*dot13*M2*sh+4.*dot2*dot9*M2*sh+dot13*sqr(sh)-2.*dot2*dot9*sqr(sh)+2.*dot3*dot4*dot7*th+
				     2.*dot4*dot9*M2*th-dot4*dot9*sqr(th)-
				     2.*dot1*(M2-th)*(dot4*(2.*dot7+dot6*M)-M*(dot12-2.*dot2*dot6+dot8*M)+dot8*th)+
				     (-2.*dot2*dot3*dot6*M+3.*dot2*dot9*(M2-sh)+dot13*(-M2+sh)+(dot3*dot8-dot4*dot9-dot11*M+dot5*(dot7+dot6*M))*(M2-th))*uh+
				     dot1*(2.*dot12*M-2.*dot2*(dot7+4.*dot6*M)+dot8*(3.*M2+th))*uh+(-dot13+2.*dot5*dot7+2.*dot3*dot8)*sqr(uh))+
				 pow(a1,4)*M2*(2.*dot3*dot4*dot7-4.*dot12*dot3*M+2.*dot3*dot4*dot6*M-3.*dot11*pow<3,1>(M)+3.*dot5*dot6*pow<3,1>(M)+
					       2.*dot13*M2+dot5*dot7*M2-7.*dot3*dot8*M2+2.*dot4*dot9*M2-dot13*sh-dot5*dot7*th+3.*dot3*dot8*th+dot11*M*th-
					       dot5*dot6*M*th-(2.*dot13+dot5*dot7+5.*dot3*dot8-dot11*M+dot5*dot6*M)*uh+
					       dot1*(-2.*dot2*dot7+2.*dot4*dot7-4.*dot12*M+6.*dot2*dot6*M+2.*dot4*dot6*M-23.*dot8*M2+5.*dot8*th-3.*dot8*uh)+
					       dot2*(-2.*dot3*dot7+6.*dot3*dot6*M-2.*dot9*M2+2.*dot9*sh+4.*dot9*uh))-
				 pow(a1,3)*M*(M*(-(dot11*pow<3,1>(M))+3.*dot2*dot9*M2+5.*dot4*dot9*M2-3.*dot2*dot9*sh+dot5*(dot7+dot6*M)*(M2-th)-
						 dot4*dot9*th+dot11*M*th+dot13*(-M2+sh-2.*uh)+dot4*dot9*uh+3.*dot11*M*uh+dot5*(dot7-3.*dot6*M)*uh)+
					      dot3*(-3.*dot8*pow<3,1>(M)+2.*dot4*M*(3.*dot7+2.*dot6*M)+10.*dot2*dot6*M2-4.*dot2*dot6*th+
						    3.*dot8*M*th+4.*dot2*dot6*uh-7.*dot8*M*uh-2.*dot12*(3.*M2-th+uh))+
					      dot1*(-2.*dot2*dot7*M-25.*dot8*pow<3,1>(M)+2.*dot4*M*(4.*dot7+3.*dot6*M)-10.*dot12*M2+
						    2.*dot12*th+13.*dot8*M*th-2.*dot12*uh-9.*dot8*M*uh+4.*dot2*dot6*(4.*M2-th+uh)))+
				 a12*(2.*M2*(dot3*(3.*dot4*dot7-dot12*M+dot4*dot6*M)+2.*dot4*dot9*M2+dot1*(5.*dot4*dot7-4.*dot12*M+3.*dot4*dot6*M-6.*dot8*M2))-
				      2.*(dot3*dot4*dot7-dot12*dot3*M+dot3*dot4*dot6*M+dot1*dot4*(dot7+dot6*M)-2.*dot1*M*(dot12+3.*dot8*M)+2.*dot4*dot9*M2)*th+
				      (-2.*dot12*dot3*M+3.*dot11*pow<3,1>(M)-3.*dot5*dot6*pow<3,1>(M)-3.*dot13*M2+dot5*dot7*M2-3.*dot3*dot8*M2+
				       2.*dot4*dot9*M2+dot13*sh+(dot5*dot7+dot3*dot8-dot11*M+dot5*dot6*M)*th-dot1*(4.*dot12*M+dot8*(9.*M2+th)))*uh+
				      dot13*sqr(uh)+dot2*(dot9*(pow<4,1>(M)-M2*(sh-2.*uh)-2.*uh*(sh+uh))+2.*dot3*(dot7*(M2+uh)+dot6*M*(2.*M2-2.*th+3.*uh))+
							  2.*dot1*(dot7*uh+dot6*M*(7.*M2-4.*th+5.*uh)))))/(a1*a2*sqr(-(a22*M2)+th)*sqr(-(a12*M2)+uh));
	    diag[1]+=mix2*8.*M2*(-2.*(-1.+2.*a1)*a2*dot10*dot2*M-2.*a1*dot11*pow<3,1>(M)+6.*a12*dot11*pow<3,1>(M)-6.*pow(a1,3)*dot11*pow<3,1>(M)+
				 2.*pow(a1,4)*dot11*pow<3,1>(M)+4.*a1*dot5*dot6*pow<3,1>(M)-12.*a12*dot5*dot6*pow<3,1>(M)+12.*pow(a1,3)*dot5*dot6*pow<3,1>(M)-
				 4.*pow(a1,4)*dot5*dot6*pow<3,1>(M)-2.*dot13*M2+a1*dot13*M2+4.*a12*dot13*M2-3.*pow(a1,3)*dot13*M2+4.*dot5*dot7*M2-
				 2.*a1*dot5*dot7*M2-8.*a12*dot5*dot7*M2+6.*pow(a1,3)*dot5*dot7*M2+4.*dot3*dot8*M2-2.*a1*dot3*dot8*M2-8.*a12*dot3*dot8*M2+
				 6.*pow(a1,3)*dot3*dot8*M2-(-2.+a1)*(dot13-2.*dot5*dot7-2.*dot3*dot8)*th+2.*a1*a2*(dot11-2.*dot5*dot6)*M*th-
				 2.*a2*dot2*(M*(2.*(-1.+2.*a1)*dot3*dot6-a22*dot9*M)+dot9*th))/(a1*pow<3,1>(a22*M2-th));
	    diag[2]+=mix2*8.*M*(2.*pow(a1,6)*dot3*dot8*pow<5,1>(M)+pow(a1,5)*pow<4,1>(M)*(-2.*dot3*dot8*M+M*((dot2+dot4)*dot9-dot11*M+dot5*dot6*M)+(dot11-dot5*dot6)*sh)+
				(M2-sh)*(-(dot5*dot7*pow<3,1>(M))+2.*dot3*dot8*pow<3,1>(M)+dot4*dot9*pow<3,1>(M)+dot10*dot2*M2-dot10*dot4*M2-
					 2.*dot2*dot3*dot6*M2-2.*dot3*dot8*M*sh+dot10*dot4*th+dot5*dot7*M*th-dot4*dot9*M*th+
					 (dot10*dot2+dot12*dot3+dot13*M-dot3*dot8*M+dot11*sh-dot5*dot6*sh)*uh+(dot11-dot5*dot6)*sqr(uh))+
				pow(a1,4)*pow<3,1>(M)*(2.*dot2*dot3*(dot7+dot6*M)+(-dot13+dot5*dot7-6.*dot4*dot9+dot11*M-2.*dot5*dot6*M)*M2+
						       (dot13-dot5*dot7+2.*dot4*dot9-dot11*M+2.*dot5*dot6*M)*sh+
						       dot3*(-2.*dot12*M+2.*dot4*(dot7+dot6*M)+7.*dot8*M2-7.*dot8*sh-4.*dot8*uh))+
				a12*M*(-7.*dot10*dot4*pow<3,1>(M)+3.*dot13*pow<4,1>(M)-3.*dot5*dot7*pow<4,1>(M)-17.*dot4*dot9*pow<4,1>(M)+dot11*pow<5,1>(M)-
				       2.*dot5*dot6*pow<5,1>(M)+5.*dot10*dot4*M*sh-3.*dot11*pow<3,1>(M)*sh+5.*dot5*dot6*pow<3,1>(M)*sh-3.*dot13*M2*sh+
				       3.*dot5*dot7*M2*sh+13.*dot4*dot9*M2*sh+2.*dot11*M*sqr(sh)-3.*dot5*dot6*M*sqr(sh)+
				       dot2*M*(dot10+8.*dot9*M)*(-M2+sh)+dot3*(M2-sh)*(4.*dot4*dot7+M*(dot12+14.*dot8*M)-dot8*sh)+6.*dot4*dot9*M2*th-
				       2.*dot4*dot9*sh*th+(dot13-dot5*dot7-2.*dot11*M+3.*dot5*dot6*M)*(M2-sh)*uh+
				       dot3*(2.*dot12*M-2.*dot4*(dot7+dot6*M)+3.*dot8*(-M2+sh))*uh+
				       2.*dot3*dot8*sqr(uh)-2.*dot2*dot3*(dot7*(-3.*M2+th)+dot6*M*(5.*M2-5.*sh+uh)))+
				pow(a1,3)*M2*(16.*dot4*dot9*pow<3,1>(M)+dot5*dot6*pow<4,1>(M)+2.*dot10*dot4*M2-dot10*dot4*sh-10.*dot4*dot9*M*sh+
					      dot11*M2*sh-2.*dot5*dot6*M2*sh-dot11*sqr(sh)+dot5*dot6*sqr(sh)-2.*dot4*dot9*M*th+
					      2.*(dot11-dot5*dot6)*(M2-sh)*uh+
					      dot3*(-2.*dot4*M*(dot7+dot6*M)+dot12*(M2+sh)+4.*dot8*M*(-4.*M2+4.*sh+uh))-
					      dot2*(dot10*sh+dot3*(6.*dot7*M-2.*dot6*M2+4.*dot6*sh)+dot9*M*(-5.*M2+5.*sh+2.*uh)))+
				a1*(6.*dot10*dot4*pow<4,1>(M)-2.*dot13*pow<5,1>(M)+3.*dot5*dot7*pow<5,1>(M)+5.*dot4*dot9*pow<5,1>(M)-
				    dot11*pow<6,1>(M)+2.*dot5*dot6*pow<6,1>(M)+2.*dot13*pow<3,1>(M)*sh-3.*dot5*dot7*pow<3,1>(M)*sh-
				    4.*dot4*dot9*pow<3,1>(M)*sh+2.*dot11*pow<4,1>(M)*sh-4.*dot5*dot6*pow<4,1>(M)*sh-5.*dot10*dot4*M2*sh-
				    dot11*M2*sqr(sh)+2.*dot5*dot6*M2*sqr(sh)-dot5*dot7*pow<3,1>(M)*th-4.*dot4*dot9*pow<3,1>(M)*th-
				    2.*dot10*dot4*M2*th+dot10*dot4*sh*th+dot5*dot7*M*sh*th+2.*dot4*dot9*M*sh*th+dot4*dot9*M*sqr(th)-
				    (M2-sh)*(2.*dot13*M+dot11*sh-dot5*(M*(dot7-dot6*M)+dot6*sh))*uh-(dot11-dot5*dot6)*(M2-sh)*sqr(uh)+
				    dot3*(-(M*(M2-sh)*(4.*dot4*dot7+7.*dot8*M2-3.*dot8*sh))+M*(2.*dot4*(dot7+dot6*M)+M*(-3.*dot12+4.*dot8*M))*uh+
					  (dot12-4.*dot8*M)*sh*uh-2.*dot8*M*sqr(uh))+
				    dot2*(-2.*(dot10-2.*dot9*M)*M2*(M2-sh)+dot10*sh*uh+
					  dot9*M*(-M2+sh)*uh+dot9*M*sqr(uh)+2.*dot3*M*(dot7*(-M2+th)+dot6*M*(4.*M2-4.*sh+uh))))+
				dot1*(2.*dot4*(-(a2*dot6*(pow(a1,3)*pow<4,1>(M)+a12*(pow<4,1>(M)-M2*sh)+a1*M2*(2.*M2-2.*sh-uh)-(M2-sh)*(sh+uh)))+
					       (-2.+a1)*dot7*M*((-1.+2.*a1+pow(a1,3))*M2+sh-a1*(2.*sh+uh)))-
				      2.*dot2*(-(a2*dot6*((1.+a1-3.*a12+pow(a1,3))*pow<4,1>(M)+(1.+a1)*M2*(-(a2*sh)-th)+sh*th))+
					       dot7*M*((a12+pow(a1,4)-a2)*M2+sh-a1*(sh+a1*sh+a1*uh)))-
				      a2*(dot12*((3.+a12*(-5.+2.*a1))*pow<4,1>(M)+sh*th+M2*(-((3.+a1)*a2*sh)-(1.+2.*a1)*th))+
					  dot8*M*((-1.+a1*(7.+a1*(-9.+9.*a1+2.*pow(a1,3))))*pow<4,1>(M)+
						  (sh+uh)*((-1.+3.*a1)*sh+2.*a1*uh)+M2*((2.+a1*(-10.+9.*a1*a2))*sh+uh-5.*a1*uh-
											4.*pow(a1,3)*uh)))))/(a1*a2*sqr(M2-sh)*sqr(-(a22*M2)+th));
	    diag[3]+=mix2*8.*M2*(a1*M*(2.*(-1.+2.*a1)*dot4*(dot10-2.*dot1*dot6)+
				       a1*((-5.+3.*a1)*(dot13+2.*dot1*dot8)+2.*((5.-3.*a1)*dot2+(5.-4.*a1)*dot4)*dot9)*M-
				       2.*a12*a2*dot11*M2)+((1.+a1)*dot13+2.*(1.+a1)*dot1*dot8-2.*(dot2+a1*dot2+dot4)*dot9+
							    2.*a1*a2*dot11*M)*uh)/(a2*pow<3,1>(a12*M2-uh));
	    diag[4]+=mix2*8.*M2*(4.*pow(a1,6)*(2.*dot1+dot3)*dot8*pow<4,1>(M)+
				 pow(a1,5)*(-22.*dot1*dot8-14.*dot3*dot8+(dot2+5.*dot4)*dot9)*pow<4,1>(M)+
				 pow(a1,4)*M2*(-6.*dot3*dot4*dot7+2.*dot12*dot3*M-14.*dot3*dot4*dot6*M-2.*dot11*pow<3,1>(M)+
					       dot5*dot6*pow<3,1>(M)-dot13*M2+dot5*dot7*M2+13.*dot3*dot8*M2-7.*dot4*dot9*M2+
					       dot13*sh-dot5*dot7*sh+2.*dot11*M*sh-dot5*dot6*M*sh+
					       dot2*(2.*dot3*dot7-2.*dot3*dot6*M+dot9*M2-2.*dot9*sh)-
					       3.*dot3*dot8*th+3.*dot3*dot8*uh-dot1*(-2.*dot2*dot7+6.*dot4*dot7+2.*dot12*M+
										     2.*dot2*dot6*M+6.*dot4*dot6*M-19.*dot8*M2+5.*dot8*th+3.*dot8*uh))-
				 a1*(-2.*dot1*dot12*pow<3,1>(M)-dot10*dot4*pow<3,1>(M)+dot13*pow<4,1>(M)-2.*dot5*dot7*pow<4,1>(M)+dot4*dot9*pow<4,1>(M)-
				     2.*dot3*dot4*dot7*M2+dot10*dot4*M*sh-dot13*M2*sh+2.*dot5*dot7*M2*sh+2.*dot3*dot4*dot7*th+2.*dot1*dot12*M*th-
				     dot13*M2*th+2.*dot5*dot7*M2*th-2.*dot4*dot9*M2*th+dot13*sh*th-2.*dot5*dot7*sh*th+dot4*dot9*sqr(th)+
				     (-(dot10*dot4*M)+dot5*dot6*pow<3,1>(M)+2.*dot1*dot4*(dot7+dot6*M)+dot5*dot7*M2-4.*dot4*dot9*M2-dot5*dot7*sh-
				      dot5*dot6*M*sh+2.*dot4*dot9*th+3.*dot1*dot8*(-M2+th)+dot3*(2.*dot4*(3.*dot7+dot6*M)-M*(2.*dot12+dot8*M)+dot8*th))*uh+
				     (dot1-dot3)*dot8*sqr(uh)+dot2*(2.*(2.*dot3*dot7-2.*dot10*M+dot9*(M2-sh))*(M2-sh)+
								    (-2.*dot1*(dot7-dot6*M)-M*(dot10-2.*dot3*dot6+3.*dot9*M)+3.*dot9*sh)*uh-dot9*sqr(uh)))-
				 uh*(dot4*(-2.*dot3*dot7+M*(dot10+dot9*M)-dot9*th)+dot2*(-2.*dot3*dot7+dot10*M+dot9*(-M2+sh+uh)))-
				 a12*(6.*dot12*dot3*pow<3,1>(M)-3.*dot10*dot4*pow<3,1>(M)+6.*dot3*dot4*dot6*pow<3,1>(M)+
				      4.*dot13*pow<4,1>(M)-4.*dot5*dot7*pow<4,1>(M)-dot4*dot9*pow<4,1>(M)+2.*dot11*pow<5,1>(M)-
				      2.*dot5*dot6*pow<5,1>(M)+8.*dot3*dot4*dot7*M2-6.*dot12*dot3*M*sh+2.*dot10*dot4*M*sh-
				      2.*dot11*pow<3,1>(M)*sh+2.*dot5*dot6*pow<3,1>(M)*sh-4.*dot13*M2*sh+4.*dot5*dot7*M2*sh+
				      dot2*M*(2.*dot3*M*(dot7-4.*dot6*M)+dot10*(3.*M2-4.*sh)-(4.*dot1*dot6+3.*dot9*M)*(M2-sh)+8.*dot3*dot6*sh)-
				      2.*dot3*dot4*dot7*th-6.*dot3*dot4*dot6*M*th-2.*dot11*pow<3,1>(M)*th+2.*dot5*dot6*pow<3,1>(M)*th+
				      dot4*dot9*M2*th+2.*dot11*M*sh*th-2.*dot5*dot6*M*sh*th-2.*dot2*(-((dot1+dot3)*(dot7-dot6*M))+dot9*sh)*uh+
				      ((-dot13+dot5*dot7+dot4*dot9-dot5*dot6*M)*M2+(dot13-dot5*dot7+dot5*dot6*M)*sh+
				       dot3*(2.*dot12*M-4.*dot4*(dot7+2.*dot6*M)+dot8*(M2+th)))*uh+
				      3.*dot3*dot8*sqr(uh)+dot1*(8.*dot12*pow<3,1>(M)-2.*dot12*M*(sh+th)+
								 2.*dot4*(dot7+dot6*M)*(M2-th-2.*uh)+dot8*uh*(7.*M2-th+uh)))+
				 pow(a1,3)*M*(-(M*(dot10*(dot2+dot4)*M-4.*dot11*pow<3,1>(M)+3.*dot5*dot6*pow<3,1>(M)+
						   5.*dot5*dot7*M2+7.*dot2*dot9*M2-2.*dot4*dot9*M2-5.*dot5*dot7*sh-7.*dot2*dot9*sh+
						   4.*dot11*M*sh-3.*dot5*dot6*M*sh+5.*dot13*(-M2+sh)+2.*(dot2+dot4)*dot9*uh))+
					      dot1*(-7.*dot8*pow<3,1>(M)-2.*dot2*M*(dot7+dot6*M)+6.*dot4*M*(dot7+dot6*M)+
						    8.*dot12*M2-2.*dot12*sh+4.*dot2*dot6*sh+7.*dot8*M*th+11.*dot8*M*uh)+
					      dot3*(2.*dot4*M*(6.*dot7+7.*dot6*M)-2.*dot2*dot6*(M2-2.*sh)-2.*dot12*sh+
						    dot8*M*(-5.*M2+5.*th+uh))))/(a1*a2*sqr(M2-sh)*sqr(-(a12*M2)+uh));

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

