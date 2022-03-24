// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2BC3P2Jet class.
//

#include "MEPP2BC3P2Jet.h"
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
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include <numeric>

using namespace Herwig;

double MEPP2BC3P2Jet::me2() const {
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
  TensorWaveFunction      Bcw(rescaledMomenta()[2],mePartonData()[2],outgoing);
  vector<TensorWaveFunction> ten3;
  for(unsigned int ix=0;ix<5;++ix) {
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
	      complex<Energy> dot1 = rescaledMomenta()[1]*g1[ih1].wave();
	      complex<Energy> dot2 = rescaledMomenta()[3]*g1[ih1].wave();
	      complex<Energy> dot3=u2[ih2].dimensionedWave().scalar(ubar4[ih4].dimensionedWave());
	      LorentzPolarizationVectorE vec1 = u2[ih2].dimensionedWave().vectorCurrent(ubar4[ih4].dimensionedWave());
	      LorentzPolarizationVectorE vec2 = u2[ih2].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(ubar4[ih4].dimensionedWave());
	      LorentzVector<complex<Energy2> > vec3 = u2[ih2].dimensionedWave().slash(rescaledMomenta()[0]).vectorCurrent(ubar4[ih4].dimensionedWave());
	      LorentzVector<complex<Energy2> > vec4 = u2[ih2].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[0]).vectorCurrent(ubar4[ih4].dimensionedWave());
	      LorentzPolarizationVectorE vec5 = ten3[ih3].wave().preDot(rescaledMomenta()[0]);
	      LorentzPolarizationVectorE vec6 = ten3[ih3].wave().preDot(rescaledMomenta()[1]);
	      LorentzPolarizationVector  vec7 = ten3[ih3].wave().preDot(g1[ih1].wave());
	      complex<Energy2> dot4 = vec1*vec5;
	      complex<Energy2> dot5 = vec1*vec6;
	      complex<Energy3> dot6 = vec4*vec5;
	      complex<Energy3> dot7 = vec4*vec6;
	      complex<Energy>  dot8 = vec1*vec7;
	      complex<Energy2> dot9 = vec2*vec5;
	      complex<Energy2> dot10 = vec2*vec6;
	      complex<Energy3> dot11 = vec3*vec5;
	      complex<Energy3> dot12 = vec3*vec6;
	      complex<Energy2> dot13 = vec5*rescaledMomenta()[0];
	      complex<Energy2> dot14 = vec2*rescaledMomenta()[0];
	      complex<Energy2> dot15 = vec5*rescaledMomenta()[1];
	      complex<Energy>  dot16 = vec5*g1[ih1].wave();
	      complex<Energy2> dot17 = vec6*rescaledMomenta()[1];
	      complex<Energy>  dot18 = vec6*g1[ih1].wave();
	      complex<Energy>  dot19 = vec1*g1[ih1].wave();
	      complex<Energy2> dot20 = vec1*rescaledMomenta()[0];
	      // diagrams
	      diag[0]=(16.*a1*(2.*dot1*(dot4+dot5)+dot6+dot7)*pow<3,1>(M))/(a2*pow<3,1>(a12*M2-sh));
	      diag[1]=(8.*M2*(6.*pow(a1,5)*dot16*dot3*pow<4,1>(M)-pow(a1,4)*pow<3,1>(M)*(-6.*dot1*dot4+6.*dot2*dot4+M*(6.*dot16*dot3+6.*dot9+dot8*M)-dot8*th)+
			      sh*(-2.*dot14*(dot13+2.*dot15)-2.*dot11*dot2-2.*dot1*(dot11-2.*dot13*dot3)+dot6*M-2.*dot16*dot3*sh+dot9*(M2+sh-uh))+
			      pow(a1,3)*M2*(-6.*dot13*dot14-16.*dot14*dot15-6.*dot11*dot2+6.*dot2*dot4*M+3.*dot6*M+dot8*pow<3,1>(M)-
					    6.*dot1*(dot11-2.*dot13*dot3+dot4*M)+2.*dot16*dot3*M2-2.*dot18*dot3*M2+9.*dot9*M2-10.*dot16*dot3*sh+
					    3.*dot9*sh+2.*dot18*dot3*th-dot8*M*th-2.*dot16*dot3*uh-3.*dot9*uh)+
			      a12*M*(6.*dot13*dot14*M+16.*dot14*dot15*M+6.*dot11*dot2*M+2.*dot10*pow<3,1>(M)-5.*dot9*pow<3,1>(M)-
				     2.*dot2*dot4*M2+2.*dot2*dot5*M2-3.*dot6*M2+4.*dot2*dot4*sh+8.*dot16*dot3*M*sh+dot9*M*sh+dot8*M2*sh-2.*dot2*dot5*th-
				     2.*dot10*M*th-dot8*sh*th+2.*dot2*dot4*uh+5.*dot9*M*uh+
				     2.*dot1*(3.*dot11*M-6.*dot13*dot3*M+dot4*M2-dot5*M2-2.*dot4*sh+dot5*th-dot4*uh))+
			      a1*(dot6*pow<3,1>(M)-dot7*pow<3,1>(M)-dot10*pow<4,1>(M)+dot9*pow<4,1>(M)-6.*dot14*dot15*M2+4.*dot14*dot17*M2-
				  2.*dot11*dot2*M2+2.*dot12*dot2*M2+10.*dot14*dot15*sh+4.*dot11*dot2*sh-2.*dot2*dot4*M*sh-2.*dot6*M*sh-
				  dot8*pow<3,1>(M)*sh-dot10*M2*sh-2.*dot16*dot3*M2*sh+2.*dot18*dot3*M2*sh-3.*dot9*M2*sh+4.*dot16*dot3*sqr(sh)-
				  2.*dot9*sqr(sh)-4.*dot14*dot17*th-2.*dot12*dot2*th+dot7*M*th+dot10*M2*th+dot10*sh*th-2.*dot18*dot3*sh*th+
				  dot8*M*sh*th+(6.*dot14*dot15+2.*dot11*dot2-dot6*M+dot10*M2-2.*dot9*M2+2.*dot16*dot3*sh+dot9*sh-dot10*th)*uh+
				  dot9*sqr(uh)+2.*dot13*dot14*(-M2+2.*sh+uh)+
				  2.*dot1*(dot12*M2+2.*dot13*dot3*M2-2.*dot15*dot3*M2-4.*dot13*dot3*sh+dot4*M*sh-dot12*th+2.*dot15*dot3*th-2.*dot13*dot3*uh+dot11*(-M2+2.*sh+uh)))))/
		(a1*a2*sqr(-(a12*M2)+sh)*sqr(M2-th));
	      diag[2]=(8.*M2*(dot9/a2+(2.*dot15*dot19*pow<3,1>(M)-2.*dot18*dot20*pow<3,1>(M)-2.*dot2*dot5*pow<3,1>(M)-
				       dot7*pow<3,1>(M)+4.*dot14*dot17*M2+2.*dot12*dot2*M2+2.*dot14*dot15*sh-dot10*M2*sh+
				       2.*dot18*dot3*M2*sh+a12*M2*(-2.*dot14*dot15+(dot10-2.*dot18*dot3)*(M2-th))+
				       2.*a1*dot2*dot5*M*(M2-th)+2.*dot1*(dot12-2.*dot15*dot3+a2*dot5*M)*(M2-th)+
				       (-4.*dot14*dot17-2.*dot12*dot2-2.*dot15*dot19*M+(2.*dot18*dot20+2.*dot2*dot5+dot7)*M+(dot10-2.*dot18*dot3)*sh)*th)/(a1*sqr(-(a22*M2)+uh))+
			      (2.*dot14*dot15+pow(a1,3)*(-2.*dot16*dot3+dot9)*M2+a12*M*(-2.*dot1*dot4+2.*dot2*dot4+dot8*(M2-th))+
			       dot10*(M2-th)+a1*(2.*dot14*(dot13+dot15)+2.*dot11*dot2+2.*dot13*dot19*M-2.*dot16*dot20*M+
						 2.*dot1*(dot11-2.*dot13*dot3+dot4*M)-M*(2.*dot2*dot4+dot6+M*(dot10+dot8*M))+
						 2.*dot16*dot3*sh-dot9*sh+(dot10+dot8*M)*th))/(-(a1*pow(a2,3)*M2)+a1*a2*uh)))/sqr(M2-th);
	      diag[3]=(-16.*a2*(2.*dot15*dot19-2.*dot18*dot20-2.*dot2*dot5-dot7)*pow<3,1>(M))/(a1*pow<3,1>(a22*M2-uh));
	      diag[4]=(8.*M2*(-((dot10+dot9)/a2)+((a12*M2-sh)*(2.*dot14*(dot15+2.*dot17)+2.*dot12*dot2+
							       a1*M*(2.*dot2*dot5+a1*(dot10-2.*dot18*dot3)*M)+2.*dot1*(dot12-2.*dot15*dot3+dot5*M-a1*dot5*M)-
							       (dot10-2.*dot18*dot3)*sh))/(a1*sqr(-(a22*M2)+uh))+
			      (pow(a1,4)*dot8*pow<3,1>(M)-2.*a1*(dot14*(dot13+3.*dot15+2.*dot17)+(dot11+dot12)*dot2+
								 dot1*(dot11+dot12-2.*(dot13+dot15)*dot3+(dot4+dot5)*M))-
			       pow(a1,3)*(2.*dot10-2.*(dot16+dot18)*dot3+dot9+dot8*M)*M2-dot10*sh+
			       a1*(2.*dot10-2.*(dot16+dot18)*dot3+dot9+dot8*M)*sh+a12*M*(2.*(dot1-dot2)*(dot4+dot5)+dot10*M-dot8*sh))/(-(a1*pow(a2,3)*M2)+a1*a2*uh)))/sqr(-(a12*M2)+sh);
	      Complex aSum =  (-0.5*diag[0] + diag[1] + diag[2] - 0.5*diag[3])/3. + 1.5*(diag[0] + diag[3] + 2.*diag[4]);
	      // diagram weights
	      save[0]+=norm(diag[3]);
	      save[1]+=norm(diag[1]);
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
	      complex<Energy> dot1 = rescaledMomenta()[1]*g1[ih1].wave();
	      complex<Energy> dot2 = rescaledMomenta()[3]*g1[ih1].wave();
	      complex<Energy> dot3=v4[ih4].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	      LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	      LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	      LorentzVector<complex<Energy2> > vec3 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[0]).vectorCurrent(vbar2[ih2].dimensionedWave());
	      LorentzVector<complex<Energy2> > vec4 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[0]).vectorCurrent(vbar2[ih2].dimensionedWave());
	      LorentzPolarizationVectorE vec5 = ten3[ih3].wave().preDot(rescaledMomenta()[0]);
	      LorentzPolarizationVectorE vec6 = ten3[ih3].wave().preDot(rescaledMomenta()[1]);
	      LorentzPolarizationVector vec7 = ten3[ih3].wave().preDot(g1[ih1].wave());
	      complex<Energy2> dot4 = vec1*vec5;
	      complex<Energy2> dot5 = vec1*vec6;
	      complex<Energy3> dot6 = vec4*vec5;
	      complex<Energy3> dot7 = vec4*vec6;
	      complex<Energy>  dot8 = vec1*g1[ih1].wave();
	      complex<Energy2> dot9 = vec5*rescaledMomenta()[0];
	      complex<Energy2> dot10 = vec5*rescaledMomenta()[1];
	      complex<Energy2> dot11 = vec1*rescaledMomenta()[0];
	      complex<Energy>  dot12 = vec5*g1[ih1].wave();
	      complex<Energy>  dot13 = vec6*g1[ih1].wave();
	      complex<Energy>  dot14 = vec1*vec7;
	      complex<Energy2> dot15 = vec2*vec5;
	      complex<Energy2> dot16 = vec2*vec6;
	      complex<Energy3> dot17 = vec3*vec5;
	      complex<Energy3> dot18 = vec3*vec6;
	      complex<Energy2> dot19 = vec2*rescaledMomenta()[0];
	      complex<Energy2> dot20 = vec6*rescaledMomenta()[1];
	      // diagrams
	      diag[0]=(16.*a2*(2.*dot11*(dot12+dot13)-2.*dot1*(dot4+dot5)+dot6+dot7-2.*dot8*(dot10+dot9))*pow<3,1>(M))/(a1*pow<3,1>(a22*M2-sh));
	      diag[1]=(8.*M2*(2.*pow(a1,5)*(dot15-dot12*dot3)*pow<4,1>(M)+
			      pow(a1,4)*pow<3,1>(M)*(-2.*dot1*dot4+2.*dot2*dot4-2.*dot15*M+2.*dot12*dot3*M+dot14*M2-dot14*th)-
			      pow(a1,3)*M2*(2.*dot17*dot2+2.*dot19*dot9-4.*dot2*dot3*dot9+2.*dot2*dot4*M+dot6*M+dot14*pow<3,1>(M)+
					    2.*dot1*(dot17-dot4*M)-dot15*M2-2.*dot16*M2+2.*dot13*dot3*M2+dot15*th+2.*dot16*th-
					    2.*dot13*dot3*th-dot14*M*th+4.*dot15*uh-4.*dot12*dot3*uh)+
			      uh*(-2.*dot1*dot17-2.*dot17*dot2+4.*dot2*dot3*dot9-2.*dot19*(2.*dot10+dot9)-dot6*M+dot15*(M2-th-2.*uh)+2.*dot12*dot3*uh)+
			      a12*M*(2.*dot17*dot2*M+2.*dot19*dot9*M-4.*dot2*dot3*dot9*M-dot15*pow<3,1>(M)+2.*dot2*dot5*M2+dot6*M2-2.*dot2*dot5*th+
				     dot15*M*th+(-2.*dot2*dot4-M*(-4.*dot15+4.*dot12*dot3+dot14*M)+dot14*th)*uh+2.*dot1*(M*(dot17-dot5*M)+dot5*th+dot4*uh))+
			      a1*(-((M2-th)*(2.*dot18*dot2+4.*dot19*dot20+M*(dot7-dot16*M)+dot16*th))+
				  (2.*dot17*dot2+2.*dot19*dot9+dot2*(-4.*dot3*dot9+2.*dot4*M)+M*(dot6+M*(-dot15-2.*dot16+2.*dot13*dot3+dot14*M))+
				   (dot15+2.*dot16-2.*dot13*dot3-dot14*M)*th)*uh+2.*(dot15-dot12*dot3)*sqr(uh)+2.*dot1*(dot18*(-M2+th)+(dot17-dot4*M)*uh)+
				  2.*dot10*(2.*dot2*dot3*(M2-th)+dot19*(M2-sh+uh)))))/(a1*a2*sqr(M2-th)*sqr(-(a12*M2)+uh));
	      diag[2]=(8.*M2*(-((dot15+dot16-2.*(dot12+dot13)*dot3)/a1)-((a22*M2-sh)*(4.*dot10*dot19+2.*dot17*dot2+2.*dot19*dot9-4.*dot2*dot3*dot9+
										      2.*dot11*dot12*M+2.*dot2*dot4*M-2.*a1*dot2*dot4*M+dot6*M-2.*dot8*dot9*M+
										      2.*dot1*(dot17-a2*dot4*M)+2.*dot15*M2-4.*a1*dot15*M2+2.*a12*dot15*M2-
										      2.*dot12*dot3*M2+4.*a1*dot12*dot3*M2-2.*a12*dot12*dot3*M2-2.*dot15*sh+
										      2.*dot12*dot3*sh))/(a2*sqr(M2-th))-
			      (2.*a2*dot1*(dot17+dot18)+6.*dot10*dot19+2.*dot17*dot2+2.*dot18*dot2+4.*dot19*dot20-4.*dot10*dot2*dot3+2.*dot19*dot9-
			       4.*dot2*dot3*dot9+2.*dot11*dot12*M+2.*dot11*dot13*M+2.*dot2*dot4*M+2.*dot2*dot5*M-2.*a22*dot1*(dot4+dot5)*M+dot6*M+
			       dot7*M-2.*dot10*dot8*M-2.*dot8*dot9*M-pow(a1,4)*dot14*pow<3,1>(M)+2.*dot15*M2+2.*dot16*M2-2.*dot12*dot3*M2-
			       2.*dot13*dot3*M2+pow(a1,3)*(-dot15-2.*dot16+2.*dot13*dot3+3.*dot14*M)*M2-2.*(dot15+dot16-(dot12+dot13)*dot3)*sh+
			       a12*M*(2.*dot2*(dot4+dot5)+M*(4.*dot15+6.*dot16-2.*(dot12+3.*dot13)*dot3-3.*dot14*M)+dot14*sh)-
			       a1*(6.*dot10*dot19+2.*(dot17+dot18)*dot2+4.*dot19*dot20-4.*dot10*dot2*dot3+2.*dot19*dot9-4.*dot2*dot3*dot9+
				   2.*dot11*(dot12+dot13)*M+4.*dot2*dot4*M+4.*dot2*dot5*M+dot6*M+dot7*M-2.*dot10*dot8*M-2.*dot8*dot9*M-
				   dot14*pow<3,1>(M)+5.*dot15*M2+6.*dot16*M2-4.*dot12*dot3*M2-6.*dot13*dot3*M2-dot15*sh-2.*dot16*sh+
				   2.*dot13*dot3*sh+dot14*M*sh))/(a1*a2*(M2-th))))/sqr(-(a22*M2)+sh);
	      diag[3]=(-16.*a1*(2.*dot2*dot5-dot7)*pow<3,1>(M))/(a2*pow<3,1>(a12*M2-uh));
	      diag[4]=(8.*M2*(pow(a1,6)*dot14*pow<5,1>(M)+pow(a1,5)*pow<4,1>(M)*(2.*dot12*dot3-3.*dot14*M)+
			      pow(a1,4)*pow<3,1>(M)*(2.*dot1*dot4-2.*dot2*dot4-2.*dot15*M-2.*dot12*dot3*M+
						     2.*dot13*dot3*M+3.*dot14*M2-dot14*sh-dot14*uh)+
			      uh*(2.*dot1*(dot17+dot18)+6.*dot10*dot19+2.*(dot17+dot18)*dot2+4.*dot19*dot20-4.*dot10*dot2*dot3+
				  2.*dot19*dot9-4.*dot2*dot3*dot9+2.*dot2*dot4*M+2.*dot2*dot5*M+dot15*M2+dot16*M2-dot15*sh-dot16*sh+
				  (dot15+dot16-2.*(dot12+dot13)*dot3)*uh)+
			      pow(a1,3)*M2*(4.*dot10*dot19+2.*dot17*dot2+2.*dot19*dot9-4.*dot2*dot3*dot9+4.*dot2*dot4*M-
					    2.*dot2*dot5*M-dot14*pow<3,1>(M)+2.*dot1*(dot17+(-dot4+dot5)*M)+3.*dot15*M2-
					    2.*dot16*M2-2.*dot13*dot3*M2-dot15*sh+2.*dot13*dot3*sh+dot14*M*sh+
					    (dot15-2.*(2.*dot12+dot13)*dot3+3.*dot14*M)*uh)+
			      a1*(-((M2-sh)*(2.*dot1*dot18+2.*dot18*dot2+4.*dot19*dot20+2.*dot10*(dot19-2.*dot2*dot3)+2.*dot2*dot5*M+dot16*M2-dot16*sh))-
				  (6.*dot10*dot19+2.*(dot17+dot18)*dot2+4.*dot19*dot20-4.*dot10*dot2*dot3+2.*dot19*dot9-4.*dot2*dot3*dot9+4.*dot2*dot4*M+
				   4.*dot2*dot5*M-dot14*pow<3,1>(M)+2.*dot1*(dot17+dot18-(dot4+dot5)*M)+3.*dot15*M2+4.*dot16*M2-2.*dot13*dot3*M2-dot15*sh-
				   2.*dot16*sh+2.*dot13*dot3*sh+dot14*M*sh)*uh-(dot15+dot16-2.*(dot12+dot13)*dot3)*sqr(uh))+
			      a12*M*(-2.*dot17*dot2*M+2.*dot18*dot2*M+4.*dot19*dot20*M-2.*dot10*(dot19+2.*dot2*dot3)*M-2.*dot19*dot9*M+4.*dot2*dot3*dot9*M-
				     dot15*pow<3,1>(M)+3.*dot16*pow<3,1>(M)-2.*dot2*dot4*M2+4.*dot2*dot5*M2-2.*dot2*dot5*sh+dot15*M*sh-3.*dot16*M*sh+
				     (2.*dot2*(dot4+dot5)+M*(dot15+3.*dot16+4.*dot12*dot3-3.*dot14*M)+dot14*sh)*uh-
				     2.*dot1*(dot17*M-dot18*M+dot4*uh+dot5*(M2-sh+uh)))))/(a1*a2*sqr(-(a22*M2)+sh)*sqr(-(a12*M2)+uh));
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
    meSum *= 1./2160.;
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
    meSum *=  8./243.;
  }
  // save the diagram weights
  meInfo(save);
  // final factors
  return meSum*O1_*pow<3,1>(Constants::pi*standardModel()->alphaS(scale())/M)/M2/sqr(a1*a2);
}

IBPtr MEPP2BC3P2Jet::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2BC3P2Jet::fullclone() const {
  return new_ptr(*this);
}

void MEPP2BC3P2Jet::doinit() {
  MEPP2BCJetBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<1>(bcbar,principleQuantumNumber(),1,2);

}

void MEPP2BC3P2Jet::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2*GeV2);
}

void MEPP2BC3P2Jet::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2*GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEPP2BC3P2Jet,MEPP2BCJetBase>
describeHerwigMEPP2BC3P2Jet("Herwig::MEPP2BC3P2Jet", "HwOniumParameters.so HwMEHadronOnium.so");

void MEPP2BC3P2Jet::Init() {

  static ClassDocumentation<MEPP2BC3P2Jet> documentation
    ("The MEPP2BC3P2Jet class implements the matrix element for g c -> B_c(3P2) b q and "
     "q bar -> Bc(3P2) g and charged conjugate");

}

