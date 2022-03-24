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
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      // diagrams
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
	  for(unsigned int ih3=0;ih3<3;++ih3) {
	    for(unsigned int ih1=0;ih1<2;++ih1) {
	      // diagrams
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
    meSum *= 1./81.;
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
	    
	    auto dotB1 = rescaledMomenta()[0]*g4[ih4].wave();
	    auto dotB2 = rescaledMomenta()[1]*g4[ih4].wave();
	    
	    Complex         dot5 = v3[ih3].wave()*g4[ih4].wave();
	    LorentzPolarizationVectorE vec2 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),1,-1);
	    LorentzPolarizationVectorE vec4 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).slash(v3[ih3].wave()).generalCurrent(vbar2[ih2].dimensionedWave(),-1,1);
	    complex<Energy> dot9 = vec1*g4[ih4].wave();
	    complex<Energy2> dot10 = vec2*rescaledMomenta()[3];
	    complex<Energy> dot11 = vec2*v3[ih3].wave();
	    complex<Energy2> dot13 = vec4*rescaledMomenta()[3];    
	    // diagrams 1P1
	    diag[0]= 2.*mix1*(-2.*M*(-2.*pow(a1,5)*pow<5,1>(M)*(-6.*dot5*dot7+7.*dot1*dot8+5.*dot3*dot8-dot11*M+4.*dot5*dot6*M)+
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
	    diag[1]=2.*mix1*(-4.*M2*(4.*a2*dot10*dot2*M-(-1.+2.*a1)*(dot13-2.*(dot5*dot7+dot3*dot8))*(a22*M2-th)-
	    			     4.*a2*dot2*(-2.*dot3*dot6*M+a22*dot9*M2-dot9*th)))/(a1*pow<3,1>(a22*M2-th));
	    diag[2]=2.*mix1*(-2.*M*(4.*pow(a1,5)*(-dot2+dot4)*dot9*pow<5,1>(M)+
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
	    				   2.*dot2*(M*(dot10+2.*dot9*M)*(M2-sh)+2.*dot3*(4.*dot7*M2+dot6*M*(5.*M2+sh-2.*th)-2.*dot7*(2.*sh+uh))))))/(a1*a2*sqr(M2-sh)*sqr(-(a22*M2)+th));
	    diag[3]=2.*mix1*(-4.*M2*(a1*M*(-4.*dot4*(dot10-2.*dot1*dot6)+a1*(dot13-2.*a1*dot13-2.*((-1.+2.*a1)*dot1*dot8+(dot2-2.*a1*dot2+dot4)*dot9))*M)+
	    			     ((-1.+2.*a1)*(dot13+2.*dot1*dot8)+2.*(dot2-2.*a1*dot2+dot4)*dot9)*uh))/(a2*pow<3,1>(a12*M2-uh));
	    diag[4]=2.*mix1*(-2.*M*(4.*pow(a1,5)*(dot2-dot4)*dot9*pow<5,1>(M)+
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
	    complex<Energy> dotB3=u1[ih1].dimensionedWave().scalar(vbar2[ih2].dimensionedWave());
	    auto vecB1 = u1[ih1].dimensionedWave().vectorCurrent(vbar2[ih2].dimensionedWave());
	    auto vecB2 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).vectorCurrent(vbar2[ih2].dimensionedWave());
	    auto ten1 = Complex(0.,2.)*(u1[ih1].dimensionedWave().sigma(vbar2[ih2].dimensionedWave()));
	    auto vecB3 = u1[ih1].dimensionedWave().slash(rescaledMomenta()[3]).vectorCurrent(vbar2[ih2].dimensionedWave());
	    auto ten2 = Complex(0,2.)*u1[ih1].dimensionedWave().slash(g4[ih4].wave()).sigma(vbar2[ih2].dimensionedWave());
	    auto vecB4 = u1[ih1].dimensionedWave().slash(g4[ih4].wave()).slash(rescaledMomenta()[3]).vectorCurrent(vbar2[ih2].dimensionedWave());
	    auto ten3 = Complex(0,2.)*u1[ih1].dimensionedWave().slash(rescaledMomenta()[3]).sigma(vbar2[ih2].dimensionedWave());
	    auto ten4 = Complex(0,2.)*u1[ih1].dimensionedWave().slash(g4[ih4].wave()).slash(rescaledMomenta()[3]).sigma(vbar2[ih2].dimensionedWave());
	    auto dotB4 = vecB1*rescaledMomenta()[3];
	    auto dotB5 = vecB1*g4[ih4].wave();
	    auto dotB6 = vecB2*rescaledMomenta()[3];
	    auto vecB5 = epsilon(rescaledMomenta()[0],rescaledMomenta()[1],v3[ih3].wave());
	    auto dotB7 = vecB5*vecB1;
	    auto vecB6 = epsilon(rescaledMomenta()[0],rescaledMomenta()[3],v3[ih3].wave());
	    auto dotB8 = vecB6*vecB1;
	    auto vecB7 = epsilon(rescaledMomenta()[0],v3[ih3].wave(),g4[ih4].wave());
	    auto dotB9 = vecB7*vecB1;
	    auto vecB8 = epsilon(rescaledMomenta()[1],rescaledMomenta()[3],v3[ih3].wave());
	    auto dotB10 = vecB8*vecB1;
	    auto vecB9 = epsilon(rescaledMomenta()[1],v3[ih3].wave(),g4[ih4].wave());
	    auto dotB11 = vecB9*vecB1;
	    auto vecB10 = epsilon(rescaledMomenta()[3],v3[ih3].wave(),g4[ih4].wave());
	    auto dotB12 = vecB10*vecB1;
	    auto dotB13 = vecB5*vecB2;
	    auto dotB14 = vecB6*vecB2;
	    auto dotB15 = vecB8*vecB2;
	    auto dotB16 = vecB5*vecB3;
	    auto dotB17 = vecB6*vecB3;
	    auto dotB18 = vecB8*vecB3;
	    auto dotB19 = vecB5*vecB4;
	    auto dotB20 = vecB6*vecB4;
	    auto dotB21 = vecB8*vecB4;
	    auto dotB22 = vecB8*rescaledMomenta()[0];
	    auto dotB23 = vecB9*rescaledMomenta()[0];
	    auto dotB24 = vecB10*rescaledMomenta()[0];
	    auto dotB25 = vecB10*rescaledMomenta()[1];
	    auto ten5 = epsilon(rescaledMomenta()[0],v3[ih3].wave());
	    auto dotB26=ten5*ten1;
	    auto ten6 = epsilon(rescaledMomenta()[1],v3[ih3].wave());
	    auto dotB27=ten6*ten1;
	    auto ten7 = epsilon(rescaledMomenta()[3],v3[ih3].wave());
	    auto dotB28=ten7*ten1;
	    auto dotB29=ten5*ten2;
	    auto dotB30=ten6*ten2;
	    auto dotB31=ten7*ten2;
	    auto dotB32=ten5*ten3;
	    auto dotB33=ten6*ten3;
	    auto dotB34=ten7*ten3;
	    auto dotB35 = vecB7*vecB3;
	    auto dotB36 = vecB9*vecB3;
	    auto dotB37 = vecB10*vecB3;
	    auto dotB38=ten5*ten4;
	    auto dotB39=ten6*ten4;
	    auto dotB40=ten7*ten4;
	    diag[0]+=mix2*(-2.*(-4.*pow(a1,6)*(dotB1+dotB2)*(dotB26+dotB27-dotB28)*pow<5,1>(M)+
			       2.*pow(a1,5)*pow<4,1>(M)*((dotB1+dotB2)*(4.*dotB10+dotB32+dotB33-dotB34+4.*dotB7)+
							 (7.*dotB1+5.*dotB2)*(dotB26+dotB27-dotB28)*M-4.*(dotB24+dotB25)*dotB3*M+
							 (-2.*dotB11+2.*dotB12+dotB29+dotB30-dotB31-2.*dotB9)*M2)-
			       2.*pow(a1,3)*M2*(-4.*dotB1*dotB17*M-4.*dotB1*dotB18*M+4.*dotB17*dotB2*M+4.*dotB18*dotB2*M-4.*dotB13*pow<3,1>(M)-
						6.*dotB14*pow<3,1>(M)-10.*dotB15*pow<3,1>(M)-5.*dotB1*dotB26*pow<3,1>(M)-dotB2*dotB26*pow<3,1>(M)-
						5.*dotB1*dotB27*pow<3,1>(M)-dotB2*dotB27*pow<3,1>(M)+5.*dotB1*dotB28*pow<3,1>(M)+dotB2*dotB28*pow<3,1>(M)+
						4.*dotB23*dotB3*pow<3,1>(M)+4.*dotB25*dotB3*pow<3,1>(M)+2.*dotB11*pow<4,1>(M)-2.*dotB12*pow<4,1>(M)-
						2.*dotB29*pow<4,1>(M)-2.*dotB30*pow<4,1>(M)+2.*dotB31*pow<4,1>(M)+2.*dotB9*pow<4,1>(M)-8.*dotB1*dotB10*M2-
						8.*dotB19*M2+4.*dotB10*dotB2*M2+4.*dotB20*M2-4.*dotB21*M2-5.*dotB1*dotB32*M2+dotB2*dotB32*M2-
						5.*dotB1*dotB33*M2+dotB2*dotB33*M2+5.*dotB1*dotB34*M2-dotB2*dotB34*M2+8.*dotB23*dotB4*M2+
						8.*dotB25*dotB4*M2+8.*dotB22*dotB5*M2-8.*dotB1*dotB7*M2+4.*dotB2*dotB8*M2+4.*dotB1*dotB10*th+
						4.*dotB10*dotB2*th+dotB1*dotB32*th+dotB2*dotB32*th+dotB1*dotB33*th+dotB2*dotB33*th-dotB1*dotB34*th-
						dotB2*dotB34*th+4.*dotB1*dotB7*th+4.*dotB2*dotB7*th+2.*dotB14*M*th+2.*dotB15*M*th+3.*dotB1*dotB26*M*th+
						dotB2*dotB26*M*th+3.*dotB1*dotB27*M*th+dotB2*dotB27*M*th-3.*dotB1*dotB28*M*th-dotB2*dotB28*M*th-
						4.*dotB23*dotB3*M*th-4.*dotB25*dotB3*M*th-2.*dotB11*M2*th+2.*dotB12*M2*th+2.*dotB29*M2*th+2.*dotB30*M2*th-
						2.*dotB31*M2*th-2.*dotB9*M2*th+((dotB1+dotB2)*(4.*dotB10+dotB32+dotB33-dotB34+4.*dotB7)-
										2.*(dotB14+dotB15)*M+(7.*dotB1+5.*dotB2)*(dotB26+dotB27-dotB28)*M+
										4.*(dotB23-2.*dotB24-dotB25)*dotB3*M+2.*(dotB11-dotB12+dotB9)*M2)*uh)+
			       uh*(8.*dotB16*dotB2*M-8.*dotB17*dotB2*M+8.*dotB22*dotB6*M-4.*dotB13*pow<3,1>(M)+4.*dotB14*pow<3,1>(M)+dotB29*pow<4,1>(M)+
				   dotB30*pow<4,1>(M)-dotB31*pow<4,1>(M)+8.*dotB19*M2-4.*dotB20*M2+4.*dotB21*M2-2.*dotB2*dotB32*M2-2.*dotB2*dotB33*M2+
				   2.*dotB2*dotB34*M2-8.*dotB23*dotB4*M2-8.*dotB25*dotB4*M2-8.*dotB22*dotB5*M2+
				   2.*dotB1*M*(-4.*dotB16+4.*dotB17+8.*dotB22*dotB3+2.*(4.*dotB10+dotB32+dotB33-dotB34+4.*dotB7)*M+(dotB26+dotB27-dotB28)*M2)-
				   8.*dotB19*th+4.*dotB20*th-4.*dotB21*th+2.*dotB2*dotB32*th+2.*dotB2*dotB33*th-2.*dotB2*dotB34*th+8.*dotB23*dotB4*th+8.*dotB25*dotB4*th+
				   8.*dotB22*dotB5*th+4.*dotB13*M*th-4.*dotB14*M*th-2.*dotB1*(8.*dotB10+2.*(dotB32+dotB33-dotB34+2.*(dotB7+dotB8))+(dotB26+dotB27-dotB28)*M)*th-
				   2.*dotB29*M2*th-2.*dotB30*M2*th+2.*dotB31*M2*th+dotB29*sqr(th)+dotB30*sqr(th)-dotB31*sqr(th)+
				   M*(-4.*dotB13+4.*dotB14+8.*(-dotB23+dotB24)*dotB3+(-4.*dotB11+4.*dotB12+dotB29+dotB30-dotB31-4.*dotB9)*M)*uh+
				   (4.*dotB11-4.*dotB12-dotB29-dotB30+dotB31+4.*dotB9)*th*uh)+
			       2.*a1*(2.*M*(dotB13*pow<4,1>(M)+dotB15*pow<4,1>(M)-2.*dotB16*dotB2*M2-2.*dotB18*dotB2*M2-2.*dotB22*dotB6*M2+2.*dotB22*dotB6*sh+
					    2.*dotB1*(dotB16+dotB18-2.*dotB22*dotB3-(dotB10+dotB7)*M)*(M2-th)+2.*(dotB16+dotB18)*dotB2*th-
					    2.*(dotB13+dotB15)*M2*th+(dotB13+dotB15)*sqr(th))+
				      (-4.*dotB16*dotB2*M+4.*dotB17*dotB2*M+8.*dotB13*pow<3,1>(M)-6.*dotB14*pow<3,1>(M)+2.*dotB15*pow<3,1>(M)-dotB2*dotB26*pow<3,1>(M)-
				       dotB2*dotB27*pow<3,1>(M)+dotB2*dotB28*pow<3,1>(M)+4.*dotB23*dotB3*pow<3,1>(M)+4.*dotB25*dotB3*pow<3,1>(M)+
				       2.*dotB11*pow<4,1>(M)-2.*dotB12*pow<4,1>(M)-2.*dotB29*pow<4,1>(M)-2.*dotB30*pow<4,1>(M)+2.*dotB31*pow<4,1>(M)+
				       2.*dotB9*pow<4,1>(M)+dotB1*M*(4.*dotB16-4.*dotB17-8.*dotB22*dotB3-5.*M*(4.*dotB10+dotB32+dotB33-dotB34+4.*dotB7+(dotB26+dotB27-dotB28)*M))-
				       8.*dotB19*M2-4.*dotB10*dotB2*M2+4.*dotB20*M2-4.*dotB21*M2+dotB2*dotB32*M2+dotB2*dotB33*M2-dotB2*dotB34*M2+
				       8.*dotB23*dotB4*M2+8.*dotB25*dotB4*M2+8.*dotB22*dotB5*M2-8.*dotB2*dotB7*M2+4.*dotB2*dotB8*M2+
				       dotB1*(4.*dotB10+dotB32+dotB33-dotB34+4.*dotB7+3.*(dotB26+dotB27-dotB28)*M)*th+
				       (4.*dotB10*dotB2+dotB2*(dotB32+dotB33-dotB34+4.*dotB7+(dotB26+dotB27-dotB28)*M)+
					2.*M*(-2.*dotB13+dotB14-dotB15-2.*(dotB23+dotB25)*dotB3+(-dotB11+dotB12+dotB29+dotB30-dotB31-dotB9)*M))*th)*uh+
				      M*(2.*dotB13-2.*dotB14+4.*(dotB23-dotB24)*dotB3+(4.*dotB11-4.*dotB12-dotB29-dotB30+dotB31+4.*dotB9)*M)*sqr(uh))+
			       pow(a1,4)*pow<3,1>(M)*(-8.*dotB19*M-8.*dotB10*dotB2*M+4.*dotB20*M-4.*dotB21*M-2.*dotB2*dotB32*M-2.*dotB2*dotB33*M+
						      2.*dotB2*dotB34*M+8.*dotB23*dotB4*M+8.*dotB25*dotB4*M+8.*dotB22*dotB5*M-16.*dotB2*dotB7*M+
						      8.*dotB2*dotB8*M+8.*dotB11*pow<3,1>(M)-8.*dotB12*pow<3,1>(M)-5.*dotB29*pow<3,1>(M)-5.*dotB30*pow<3,1>(M)+
						      5.*dotB31*pow<3,1>(M)+8.*dotB9*pow<3,1>(M)-8.*dotB14*M2-8.*dotB15*M2-8.*dotB2*dotB26*M2-8.*dotB2*dotB27*M2+
						      8.*dotB2*dotB28*M2+8.*dotB23*dotB3*M2+8.*dotB24*dotB3*M2+16.*dotB25*dotB3*M2+4.*dotB2*dotB26*th+
						      4.*dotB2*dotB27*th-4.*dotB2*dotB28*th+dotB29*M*th+dotB30*M*th-dotB31*M*th+4.*dotB2*(dotB26+dotB27-dotB28)*uh+
						      (4.*dotB11-4.*dotB12-dotB29-dotB30+dotB31+4.*dotB9)*M*uh+
						      2.*dotB1*(M*(-4.*(3.*dotB10+dotB32+dotB33-dotB34+3.*dotB7)-9.*(dotB26+dotB27-dotB28)*M)+
								2.*(dotB26+dotB27-dotB28)*th+2.*(dotB26+dotB27-dotB28)*uh))+
			       a12*M*(-8.*dotB19*pow<3,1>(M)+8.*dotB10*dotB2*pow<3,1>(M)+4.*dotB20*pow<3,1>(M)-4.*dotB21*pow<3,1>(M)+2.*dotB2*dotB32*pow<3,1>(M)+
				      2.*dotB2*dotB33*pow<3,1>(M)-2.*dotB2*dotB34*pow<3,1>(M)+8.*dotB23*dotB4*pow<3,1>(M)+8.*dotB25*dotB4*pow<3,1>(M)+
				      8.*dotB22*dotB5*pow<3,1>(M)+8.*dotB2*dotB7*pow<3,1>(M)-12.*dotB13*pow<4,1>(M)-4.*dotB14*pow<4,1>(M)-
				      16.*dotB15*pow<4,1>(M)-dotB29*pow<5,1>(M)-dotB30*pow<5,1>(M)+dotB31*pow<5,1>(M)+8.*dotB16*dotB2*M2+
				      8.*dotB17*dotB2*M2+16.*dotB18*dotB2*M2-8.*dotB22*dotB6*M2+8.*dotB19*M*th-8.*dotB10*dotB2*M*th-
				      4.*dotB20*M*th+4.*dotB21*M*th-2.*dotB2*dotB32*M*th-2.*dotB2*dotB33*M*th+2.*dotB2*dotB34*M*th-
				      8.*dotB23*dotB4*M*th-8.*dotB25*dotB4*M*th-8.*dotB22*dotB5*M*th-8.*dotB2*dotB7*M*th+2.*dotB29*pow<3,1>(M)*th+
				      2.*dotB30*pow<3,1>(M)*th-2.*dotB31*pow<3,1>(M)*th+12.*dotB13*M2*th+4.*dotB14*M2*th+16.*dotB15*M2*th-
				      dotB29*M*sqr(th)-dotB30*M*sqr(th)+dotB31*M*sqr(th)+
				      2.*M*(4.*dotB19+dotB2*(8.*dotB10+dotB32+dotB33-dotB34+12.*dotB7-4.*dotB8)+4.*dotB2*(dotB26+dotB27-dotB28)*M-
					    8.*(dotB24+dotB25)*dotB3*M-2.*(dotB20-dotB21+2.*(dotB23+dotB25)*dotB4+2.*dotB22*dotB5+3.*dotB13*M-dotB14*M+2.*dotB15*M)+
					    2.*(-dotB11+dotB12+dotB29+dotB30-dotB31-dotB9)*M2)*uh-4.*(dotB2*(dotB26+dotB27-dotB28)+(dotB11-dotB12+dotB9)*M)*th*uh+
				      (-4.*dotB11+4.*dotB12+dotB29+dotB30-dotB31-4.*dotB9)*M*sqr(uh)+
				      2.*dotB1*(M2*(-4.*(dotB16+dotB17+2.*dotB18-2.*dotB22*dotB3)+2.*(2.*dotB10-dotB32-dotB33+dotB34+2.*dotB7)*M-
						    (dotB26+dotB27-dotB28)*M2)+M*(4.*dotB10+2.*(dotB32+dotB33-dotB34+2.*dotB8)+
										  (dotB26+dotB27-dotB28)*M)*th+
						M*(4.*(4.*dotB10+dotB32+dotB33-dotB34+4.*dotB7)+9.*(dotB26+dotB27-dotB28)*M)*uh-
						2.*(dotB26+dotB27-dotB28)*th*uh))))/(a1*a2*sqr(-(a22*M2)+th)*sqr(-(a12*M2)+uh));
	    diag[1]+=mix2*(-4.*(4.*dotB14*pow<3,1>(M)-16.*a1*dotB14*pow<3,1>(M)+20.*a12*dotB14*pow<3,1>(M)-8.*pow(a1,3)*dotB14*pow<3,1>(M)+
			       4.*dotB15*pow<3,1>(M)-16.*a1*dotB15*pow<3,1>(M)+20.*a12*dotB15*pow<3,1>(M)-8.*pow(a1,3)*dotB15*pow<3,1>(M)+
			       2.*dotB2*dotB26*pow<3,1>(M)-8.*a1*dotB2*dotB26*pow<3,1>(M)+10.*a12*dotB2*dotB26*pow<3,1>(M)-
			       4.*pow(a1,3)*dotB2*dotB26*pow<3,1>(M)+2.*dotB2*dotB27*pow<3,1>(M)-8.*a1*dotB2*dotB27*pow<3,1>(M)+
			       10.*a12*dotB2*dotB27*pow<3,1>(M)-4.*pow(a1,3)*dotB2*dotB27*pow<3,1>(M)-2.*dotB2*dotB28*pow<3,1>(M)+
			       8.*a1*dotB2*dotB28*pow<3,1>(M)-10.*a12*dotB2*dotB28*pow<3,1>(M)+4.*pow(a1,3)*dotB2*dotB28*pow<3,1>(M)+
			       8.*dotB24*dotB3*pow<3,1>(M)-32.*a1*dotB24*dotB3*pow<3,1>(M)+40.*a12*dotB24*dotB3*pow<3,1>(M)-
			       16.*pow(a1,3)*dotB24*dotB3*pow<3,1>(M)+8.*dotB25*dotB3*pow<3,1>(M)-32.*a1*dotB25*dotB3*pow<3,1>(M)+
			       40.*a12*dotB25*dotB3*pow<3,1>(M)-16.*pow(a1,3)*dotB25*dotB3*pow<3,1>(M)+4.*dotB35*pow<3,1>(M)-
			       16.*a1*dotB35*pow<3,1>(M)+20.*a12*dotB35*pow<3,1>(M)-8.*pow(a1,3)*dotB35*pow<3,1>(M)+4.*dotB36*pow<3,1>(M)-
			       16.*a1*dotB36*pow<3,1>(M)+20.*a12*dotB36*pow<3,1>(M)-8.*pow(a1,3)*dotB36*pow<3,1>(M)-4.*dotB37*pow<3,1>(M)+
			       16.*a1*dotB37*pow<3,1>(M)-20.*a12*dotB37*pow<3,1>(M)+8.*pow(a1,3)*dotB37*pow<3,1>(M)+dotB38*pow<3,1>(M)-
			       4.*a1*dotB38*pow<3,1>(M)+5.*a12*dotB38*pow<3,1>(M)-2.*pow(a1,3)*dotB38*pow<3,1>(M)+dotB39*pow<3,1>(M)-
			       4.*a1*dotB39*pow<3,1>(M)+5.*a12*dotB39*pow<3,1>(M)-2.*pow(a1,3)*dotB39*pow<3,1>(M)-dotB40*pow<3,1>(M)+
			       4.*a1*dotB40*pow<3,1>(M)-5.*a12*dotB40*pow<3,1>(M)+2.*pow(a1,3)*dotB40*pow<3,1>(M)-4.*dotB11*pow<4,1>(M)+
			       16.*a1*dotB11*pow<4,1>(M)-24.*a12*dotB11*pow<4,1>(M)+16.*pow(a1,3)*dotB11*pow<4,1>(M)-
			       4.*pow(a1,4)*dotB11*pow<4,1>(M)+4.*dotB12*pow<4,1>(M)-16.*a1*dotB12*pow<4,1>(M)+24.*a12*dotB12*pow<4,1>(M)-
			       16.*pow(a1,3)*dotB12*pow<4,1>(M)+4.*pow(a1,4)*dotB12*pow<4,1>(M)+dotB29*pow<4,1>(M)-4.*a1*dotB29*pow<4,1>(M)+
			       6.*a12*dotB29*pow<4,1>(M)-4.*pow(a1,3)*dotB29*pow<4,1>(M)+pow(a1,4)*dotB29*pow<4,1>(M)+dotB30*pow<4,1>(M)-
			       4.*a1*dotB30*pow<4,1>(M)+6.*a12*dotB30*pow<4,1>(M)-4.*pow(a1,3)*dotB30*pow<4,1>(M)+pow(a1,4)*dotB30*pow<4,1>(M)-
			       dotB31*pow<4,1>(M)+4.*a1*dotB31*pow<4,1>(M)-6.*a12*dotB31*pow<4,1>(M)+4.*pow(a1,3)*dotB31*pow<4,1>(M)-
			       pow(a1,4)*dotB31*pow<4,1>(M)-4.*dotB9*pow<4,1>(M)+16.*a1*dotB9*pow<4,1>(M)-24.*a12*dotB9*pow<4,1>(M)+
			       16.*pow(a1,3)*dotB9*pow<4,1>(M)-4.*pow(a1,4)*dotB9*pow<4,1>(M)-4.*(-2.+a1)*a2*dotB19*M2-16.*dotB23*dotB4*M2+
			       24.*a1*dotB23*dotB4*M2-8.*a12*dotB23*dotB4*M2+16.*dotB24*dotB4*M2-24.*a1*dotB24*dotB4*M2+
			       8.*a12*dotB24*dotB4*M2-16.*dotB22*dotB5*M2+24.*a1*dotB22*dotB5*M2-8.*a12*dotB22*dotB5*M2+
			       16.*dotB2*dotB7*M2-24.*a1*dotB2*dotB7*M2+8.*a12*dotB2*dotB7*M2-16.*dotB2*dotB8*M2+24.*a1*dotB2*dotB8*M2-
			       8.*a12*dotB2*dotB8*M2-4.*dotB19*th+
			       (8.*(dotB23*dotB4-dotB24*dotB4+dotB22*dotB5-dotB2*dotB7+dotB2*dotB8)+
				(-1.+2.*a1)*(4.*dotB14+4.*dotB15+2.*dotB2*(dotB26+dotB27-dotB28)+
					     8.*(dotB24+dotB25)*dotB3+4.*dotB35+4.*dotB36-4.*dotB37+dotB38+dotB39-dotB40)*M+
				2.*a22*(4.*dotB11-4.*dotB12-dotB29-dotB30+dotB31+4.*dotB9)*M2)*th+
			       (-4.*dotB11+4.*dotB12+dotB29+dotB30-dotB31-4.*dotB9)*sqr(th)+4.*dotB20*((-2.+a1)*a2*M2+th)))/(a1*pow<3,1>(a22*M2-th));
	    diag[2]+=mix2*(2.*(8.*dotB16*dotB2*pow<4,1>(M)-16.*a12*dotB16*dotB2*pow<4,1>(M)+8.*pow(a1,3)*dotB16*dotB2*pow<4,1>(M)-8.*dotB17*dotB2*pow<4,1>(M)+
			      20.*a1*dotB17*dotB2*pow<4,1>(M)+12.*pow(a1,3)*dotB17*dotB2*pow<4,1>(M)-4.*a1*dotB18*dotB2*pow<4,1>(M)+
			      8.*a12*dotB18*dotB2*pow<4,1>(M)-4.*pow(a1,3)*dotB18*dotB2*pow<4,1>(M)-8.*a1*dotB2*dotB22*dotB3*pow<4,1>(M)+
			      16.*a12*dotB2*dotB22*dotB3*pow<4,1>(M)-8.*pow(a1,3)*dotB2*dotB22*dotB3*pow<4,1>(M)+8.*dotB22*dotB6*pow<4,1>(M)-
			      24.*a1*dotB22*dotB6*pow<4,1>(M)-16.*pow(a1,3)*dotB22*dotB6*pow<4,1>(M)+4.*a1*dotB19*pow<5,1>(M)-4.*a12*dotB19*pow<5,1>(M)+
			      8.*dotB10*dotB2*pow<5,1>(M)-28.*a1*dotB10*dotB2*pow<5,1>(M)+40.*a12*dotB10*dotB2*pow<5,1>(M)-
			      28.*pow(a1,3)*dotB10*dotB2*pow<5,1>(M)+8.*pow(a1,4)*dotB10*dotB2*pow<5,1>(M)-4.*dotB20*pow<5,1>(M)+
			      8.*a1*dotB20*pow<5,1>(M)+4.*pow(a1,3)*dotB20*pow<5,1>(M)-4.*dotB21*pow<5,1>(M)+4.*a1*dotB21*pow<5,1>(M)+
			      4.*a12*dotB21*pow<5,1>(M)-4.*pow(a1,3)*dotB21*pow<5,1>(M)+2.*dotB2*dotB32*pow<5,1>(M)-
			      4.*a1*dotB2*dotB32*pow<5,1>(M)+2.*a12*dotB2*dotB32*pow<5,1>(M)+2.*dotB2*dotB33*pow<5,1>(M)-
			      4.*a1*dotB2*dotB33*pow<5,1>(M)+2.*a12*dotB2*dotB33*pow<5,1>(M)-2.*dotB2*dotB34*pow<5,1>(M)+
			      4.*a1*dotB2*dotB34*pow<5,1>(M)-2.*a12*dotB2*dotB34*pow<5,1>(M)-8.*dotB23*dotB4*pow<5,1>(M)+
			      8.*a1*dotB23*dotB4*pow<5,1>(M)+8.*dotB24*dotB4*pow<5,1>(M)-16.*a1*dotB24*dotB4*pow<5,1>(M)-
			      8.*pow(a1,3)*dotB24*dotB4*pow<5,1>(M)+8.*a1*dotB25*dotB4*pow<5,1>(M)-16.*a12*dotB25*dotB4*pow<5,1>(M)+
			      8.*pow(a1,3)*dotB25*dotB4*pow<5,1>(M)-8.*dotB22*dotB5*pow<5,1>(M)+8.*a1*dotB22*dotB5*pow<5,1>(M)+
			      16.*dotB2*dotB7*pow<5,1>(M)-40.*a1*dotB2*dotB7*pow<5,1>(M)+32.*a12*dotB2*dotB7*pow<5,1>(M)-
			      8.*pow(a1,3)*dotB2*dotB7*pow<5,1>(M)-8.*dotB2*dotB8*pow<5,1>(M)+12.*a1*dotB2*dotB8*pow<5,1>(M)+
			      8.*a12*dotB2*dotB8*pow<5,1>(M)-20.*pow(a1,3)*dotB2*dotB8*pow<5,1>(M)+8.*pow(a1,4)*dotB2*dotB8*pow<5,1>(M)+
			      4.*dotB13*pow<6,1>(M)-8.*a1*dotB13*pow<6,1>(M)+8.*a12*dotB13*pow<6,1>(M)-12.*pow(a1,3)*dotB13*pow<6,1>(M)+
			      4.*a1*dotB14*pow<6,1>(M)-4.*a12*dotB14*pow<6,1>(M)+16.*pow(a1,3)*dotB14*pow<6,1>(M)+8.*pow(a1,5)*dotB14*pow<6,1>(M)-
			      4.*dotB15*pow<6,1>(M)+12.*a1*dotB15*pow<6,1>(M)-12.*a12*dotB15*pow<6,1>(M)-4.*pow(a1,3)*dotB15*pow<6,1>(M)+
			      16.*pow(a1,4)*dotB15*pow<6,1>(M)-8.*pow(a1,5)*dotB15*pow<6,1>(M)+2.*dotB2*dotB26*pow<6,1>(M)-
			      10.*a1*dotB2*dotB26*pow<6,1>(M)+18.*a12*dotB2*dotB26*pow<6,1>(M)-16.*pow(a1,3)*dotB2*dotB26*pow<6,1>(M)+
			      8.*pow(a1,4)*dotB2*dotB26*pow<6,1>(M)-2.*pow(a1,5)*dotB2*dotB26*pow<6,1>(M)+2.*dotB2*dotB27*pow<6,1>(M)-
			      10.*a1*dotB2*dotB27*pow<6,1>(M)+18.*a12*dotB2*dotB27*pow<6,1>(M)-16.*pow(a1,3)*dotB2*dotB27*pow<6,1>(M)+
			      8.*pow(a1,4)*dotB2*dotB27*pow<6,1>(M)-2.*pow(a1,5)*dotB2*dotB27*pow<6,1>(M)-2.*dotB2*dotB28*pow<6,1>(M)+
			      10.*a1*dotB2*dotB28*pow<6,1>(M)-18.*a12*dotB2*dotB28*pow<6,1>(M)+16.*pow(a1,3)*dotB2*dotB28*pow<6,1>(M)-
			      8.*pow(a1,4)*dotB2*dotB28*pow<6,1>(M)+2.*pow(a1,5)*dotB2*dotB28*pow<6,1>(M)-4.*a1*dotB23*dotB3*pow<6,1>(M)+
			      8.*a12*dotB23*dotB3*pow<6,1>(M)-12.*pow(a1,3)*dotB23*dotB3*pow<6,1>(M)+4.*a1*dotB24*dotB3*pow<6,1>(M)-
			      8.*a12*dotB24*dotB3*pow<6,1>(M)+20.*pow(a1,3)*dotB24*dotB3*pow<6,1>(M)+8.*pow(a1,5)*dotB24*dotB3*pow<6,1>(M)-
			      8.*pow(a1,3)*dotB25*dotB3*pow<6,1>(M)+16.*pow(a1,4)*dotB25*dotB3*pow<6,1>(M)-8.*pow(a1,5)*dotB25*dotB3*pow<6,1>(M)-
			      4.*dotB35*pow<6,1>(M)+14.*a1*dotB35*pow<6,1>(M)-16.*a12*dotB35*pow<6,1>(M)+6.*pow(a1,3)*dotB35*pow<6,1>(M)-
			      4.*dotB36*pow<6,1>(M)+14.*a1*dotB36*pow<6,1>(M)-16.*a12*dotB36*pow<6,1>(M)+6.*pow(a1,3)*dotB36*pow<6,1>(M)+
			      4.*dotB37*pow<6,1>(M)-14.*a1*dotB37*pow<6,1>(M)+16.*a12*dotB37*pow<6,1>(M)-6.*pow(a1,3)*dotB37*pow<6,1>(M)-
			      3.*dotB38*pow<6,1>(M)+11.*a1*dotB38*pow<6,1>(M)-13.*a12*dotB38*pow<6,1>(M)+5.*pow(a1,3)*dotB38*pow<6,1>(M)-
			      3.*dotB39*pow<6,1>(M)+11.*a1*dotB39*pow<6,1>(M)-13.*a12*dotB39*pow<6,1>(M)+5.*pow(a1,3)*dotB39*pow<6,1>(M)+
			      3.*dotB40*pow<6,1>(M)-11.*a1*dotB40*pow<6,1>(M)+13.*a12*dotB40*pow<6,1>(M)-5.*pow(a1,3)*dotB40*pow<6,1>(M)+
			      2.*a1*dotB11*pow<7,1>(M)-8.*a12*dotB11*pow<7,1>(M)+10.*pow(a1,3)*dotB11*pow<7,1>(M)-
			      4.*pow(a1,4)*dotB11*pow<7,1>(M)-2.*a1*dotB12*pow<7,1>(M)+8.*a12*dotB12*pow<7,1>(M)-10.*pow(a1,3)*dotB12*pow<7,1>(M)+
			      4.*pow(a1,4)*dotB12*pow<7,1>(M)+a1*dotB29*pow<7,1>(M)-a12*dotB29*pow<7,1>(M)-pow(a1,3)*dotB29*pow<7,1>(M)+
			      pow(a1,4)*dotB29*pow<7,1>(M)+a1*dotB30*pow<7,1>(M)-a12*dotB30*pow<7,1>(M)-pow(a1,3)*dotB30*pow<7,1>(M)+
			      pow(a1,4)*dotB30*pow<7,1>(M)-a1*dotB31*pow<7,1>(M)+a12*dotB31*pow<7,1>(M)+pow(a1,3)*dotB31*pow<7,1>(M)-
			      pow(a1,4)*dotB31*pow<7,1>(M)+2.*a1*dotB9*pow<7,1>(M)-8.*a12*dotB9*pow<7,1>(M)+10.*pow(a1,3)*dotB9*pow<7,1>(M)-
			      4.*pow(a1,4)*dotB9*pow<7,1>(M)-4.*a1*dotB19*pow<3,1>(M)*sh+4.*a12*dotB19*pow<3,1>(M)*sh-
			      8.*dotB10*dotB2*pow<3,1>(M)*sh+20.*a1*dotB10*dotB2*pow<3,1>(M)*sh-16.*a12*dotB10*dotB2*pow<3,1>(M)*sh+
			      4.*pow(a1,3)*dotB10*dotB2*pow<3,1>(M)*sh+4.*dotB20*pow<3,1>(M)*sh-8.*a1*dotB20*pow<3,1>(M)*sh+
			      4.*dotB21*pow<3,1>(M)*sh-8.*a1*dotB21*pow<3,1>(M)*sh+4.*a12*dotB21*pow<3,1>(M)*sh-
			      2.*dotB2*dotB32*pow<3,1>(M)*sh+4.*a1*dotB2*dotB32*pow<3,1>(M)*sh-2.*a12*dotB2*dotB32*pow<3,1>(M)*sh-
			      2.*dotB2*dotB33*pow<3,1>(M)*sh+4.*a1*dotB2*dotB33*pow<3,1>(M)*sh-2.*a12*dotB2*dotB33*pow<3,1>(M)*sh+
			      2.*dotB2*dotB34*pow<3,1>(M)*sh-4.*a1*dotB2*dotB34*pow<3,1>(M)*sh+2.*a12*dotB2*dotB34*pow<3,1>(M)*sh+
			      8.*dotB23*dotB4*pow<3,1>(M)*sh-8.*a1*dotB23*dotB4*pow<3,1>(M)*sh-8.*dotB24*dotB4*pow<3,1>(M)*sh+
			      16.*a1*dotB24*dotB4*pow<3,1>(M)*sh+8.*dotB22*dotB5*pow<3,1>(M)*sh-8.*a1*dotB22*dotB5*pow<3,1>(M)*sh-
			      16.*dotB2*dotB7*pow<3,1>(M)*sh+40.*a1*dotB2*dotB7*pow<3,1>(M)*sh-32.*a12*dotB2*dotB7*pow<3,1>(M)*sh+
			      8.*pow(a1,3)*dotB2*dotB7*pow<3,1>(M)*sh+8.*dotB2*dotB8*pow<3,1>(M)*sh-20.*a1*dotB2*dotB8*pow<3,1>(M)*sh+
			      16.*a12*dotB2*dotB8*pow<3,1>(M)*sh-4.*pow(a1,3)*dotB2*dotB8*pow<3,1>(M)*sh-8.*dotB13*pow<4,1>(M)*sh+
			      20.*a1*dotB13*pow<4,1>(M)*sh-8.*a12*dotB13*pow<4,1>(M)*sh+16.*pow(a1,3)*dotB13*pow<4,1>(M)*sh-
			      8.*a1*dotB14*pow<4,1>(M)*sh+4.*a12*dotB14*pow<4,1>(M)*sh-16.*pow(a1,3)*dotB14*pow<4,1>(M)*sh+
			      4.*dotB15*pow<4,1>(M)*sh-16.*a1*dotB15*pow<4,1>(M)*sh+20.*a12*dotB15*pow<4,1>(M)*sh-
			      8.*pow(a1,3)*dotB15*pow<4,1>(M)*sh-2.*dotB2*dotB26*pow<4,1>(M)*sh+10.*a1*dotB2*dotB26*pow<4,1>(M)*sh-
			      18.*a12*dotB2*dotB26*pow<4,1>(M)*sh+16.*pow(a1,3)*dotB2*dotB26*pow<4,1>(M)*sh-8.*pow(a1,4)*dotB2*dotB26*pow<4,1>(M)*sh+
			      2.*pow(a1,5)*dotB2*dotB26*pow<4,1>(M)*sh-2.*dotB2*dotB27*pow<4,1>(M)*sh+10.*a1*dotB2*dotB27*pow<4,1>(M)*sh-
			      18.*a12*dotB2*dotB27*pow<4,1>(M)*sh+16.*pow(a1,3)*dotB2*dotB27*pow<4,1>(M)*sh-8.*pow(a1,4)*dotB2*dotB27*pow<4,1>(M)*sh+
			      2.*pow(a1,5)*dotB2*dotB27*pow<4,1>(M)*sh+2.*dotB2*dotB28*pow<4,1>(M)*sh-10.*a1*dotB2*dotB28*pow<4,1>(M)*sh+
			      18.*a12*dotB2*dotB28*pow<4,1>(M)*sh-16.*pow(a1,3)*dotB2*dotB28*pow<4,1>(M)*sh+8.*pow(a1,4)*dotB2*dotB28*pow<4,1>(M)*sh-
			      2.*pow(a1,5)*dotB2*dotB28*pow<4,1>(M)*sh+12.*a1*dotB23*dotB3*pow<4,1>(M)*sh-8.*a12*dotB23*dotB3*pow<4,1>(M)*sh+
			      16.*pow(a1,3)*dotB23*dotB3*pow<4,1>(M)*sh-8.*a1*dotB24*dotB3*pow<4,1>(M)*sh+8.*a12*dotB24*dotB3*pow<4,1>(M)*sh-
			      20.*pow(a1,3)*dotB24*dotB3*pow<4,1>(M)*sh-4.*a1*dotB25*dotB3*pow<4,1>(M)*sh+8.*a12*dotB25*dotB3*pow<4,1>(M)*sh-
			      4.*pow(a1,3)*dotB25*dotB3*pow<4,1>(M)*sh+4.*dotB35*pow<4,1>(M)*sh-16.*a1*dotB35*pow<4,1>(M)*sh+20.*a12*dotB35*pow<4,1>(M)*sh-
			      8.*pow(a1,3)*dotB35*pow<4,1>(M)*sh+4.*dotB36*pow<4,1>(M)*sh-16.*a1*dotB36*pow<4,1>(M)*sh+20.*a12*dotB36*pow<4,1>(M)*sh-
			      8.*pow(a1,3)*dotB36*pow<4,1>(M)*sh-4.*dotB37*pow<4,1>(M)*sh+16.*a1*dotB37*pow<4,1>(M)*sh-20.*a12*dotB37*pow<4,1>(M)*sh+
			      8.*pow(a1,3)*dotB37*pow<4,1>(M)*sh+3.*dotB38*pow<4,1>(M)*sh-12.*a1*dotB38*pow<4,1>(M)*sh+15.*a12*dotB38*pow<4,1>(M)*sh-
			      6.*pow(a1,3)*dotB38*pow<4,1>(M)*sh+3.*dotB39*pow<4,1>(M)*sh-12.*a1*dotB39*pow<4,1>(M)*sh+15.*a12*dotB39*pow<4,1>(M)*sh-
			      6.*pow(a1,3)*dotB39*pow<4,1>(M)*sh-3.*dotB40*pow<4,1>(M)*sh+12.*a1*dotB40*pow<4,1>(M)*sh-15.*a12*dotB40*pow<4,1>(M)*sh+
			      6.*pow(a1,3)*dotB40*pow<4,1>(M)*sh+4.*a12*dotB11*pow<5,1>(M)*sh-8.*pow(a1,3)*dotB11*pow<5,1>(M)*sh+
			      4.*pow(a1,4)*dotB11*pow<5,1>(M)*sh-4.*a12*dotB12*pow<5,1>(M)*sh+8.*pow(a1,3)*dotB12*pow<5,1>(M)*sh-
			      4.*pow(a1,4)*dotB12*pow<5,1>(M)*sh-2.*a1*dotB29*pow<5,1>(M)*sh+3.*a12*dotB29*pow<5,1>(M)*sh-pow(a1,4)*dotB29*pow<5,1>(M)*sh-
			      2.*a1*dotB30*pow<5,1>(M)*sh+3.*a12*dotB30*pow<5,1>(M)*sh-pow(a1,4)*dotB30*pow<5,1>(M)*sh+
			      2.*a1*dotB31*pow<5,1>(M)*sh-3.*a12*dotB31*pow<5,1>(M)*sh+pow(a1,4)*dotB31*pow<5,1>(M)*sh+
			      4.*a12*dotB9*pow<5,1>(M)*sh-8.*pow(a1,3)*dotB9*pow<5,1>(M)*sh+4.*pow(a1,4)*dotB9*pow<5,1>(M)*sh-
			      8.*dotB16*dotB2*M2*sh+16.*a12*dotB16*dotB2*M2*sh-8.*pow(a1,3)*dotB16*dotB2*M2*sh+
			      8.*dotB17*dotB2*M2*sh-24.*a1*dotB17*dotB2*M2*sh-4.*pow(a1,3)*dotB17*dotB2*M2*sh-4.*a1*dotB18*dotB2*M2*sh+
			      8.*a12*dotB18*dotB2*M2*sh-4.*pow(a1,3)*dotB18*dotB2*M2*sh+8.*a1*dotB2*dotB22*dotB3*M2*sh-
			      16.*a12*dotB2*dotB22*dotB3*M2*sh+8.*pow(a1,3)*dotB2*dotB22*dotB3*M2*sh-8.*dotB22*dotB6*M2*sh+
			      24.*a1*dotB22*dotB6*M2*sh+4.*a1*dotB17*dotB2*sqr(sh)-2.*a1*dotB11*pow<3,1>(M)*sqr(sh)+
			      4.*a12*dotB11*pow<3,1>(M)*sqr(sh)-2.*pow(a1,3)*dotB11*pow<3,1>(M)*sqr(sh)+2.*a1*dotB12*pow<3,1>(M)*sqr(sh)-
			      4.*a12*dotB12*pow<3,1>(M)*sqr(sh)+2.*pow(a1,3)*dotB12*pow<3,1>(M)*sqr(sh)+a1*dotB29*pow<3,1>(M)*sqr(sh)-
			      2.*a12*dotB29*pow<3,1>(M)*sqr(sh)+pow(a1,3)*dotB29*pow<3,1>(M)*sqr(sh)+a1*dotB30*pow<3,1>(M)*sqr(sh)-
			      2.*a12*dotB30*pow<3,1>(M)*sqr(sh)+pow(a1,3)*dotB30*pow<3,1>(M)*sqr(sh)-a1*dotB31*pow<3,1>(M)*sqr(sh)+
			      2.*a12*dotB31*pow<3,1>(M)*sqr(sh)-pow(a1,3)*dotB31*pow<3,1>(M)*sqr(sh)-2.*a1*dotB9*pow<3,1>(M)*sqr(sh)+
			      4.*a12*dotB9*pow<3,1>(M)*sqr(sh)-2.*pow(a1,3)*dotB9*pow<3,1>(M)*sqr(sh)+4.*dotB13*M2*sqr(sh)-
			      16.*a1*dotB13*M2*sqr(sh)-4.*pow(a1,3)*dotB13*M2*sqr(sh)+4.*a1*dotB14*M2*sqr(sh)+4.*a1*dotB15*M2*sqr(sh)-
			      8.*a12*dotB15*M2*sqr(sh)+4.*pow(a1,3)*dotB15*M2*sqr(sh)-12.*a1*dotB23*dotB3*M2*sqr(sh)-
			      4.*pow(a1,3)*dotB23*dotB3*M2*sqr(sh)+4.*a1*dotB24*dotB3*M2*sqr(sh)+4.*a1*dotB25*dotB3*M2*sqr(sh)-
			      8.*a12*dotB25*dotB3*M2*sqr(sh)+4.*pow(a1,3)*dotB25*dotB3*M2*sqr(sh)+2.*a1*dotB35*M2*sqr(sh)-
			      4.*a12*dotB35*M2*sqr(sh)+2.*pow(a1,3)*dotB35*M2*sqr(sh)+2.*a1*dotB36*M2*sqr(sh)-
			      4.*a12*dotB36*M2*sqr(sh)+2.*pow(a1,3)*dotB36*M2*sqr(sh)-2.*a1*dotB37*M2*sqr(sh)+4.*a12*dotB37*M2*sqr(sh)-
			      2.*pow(a1,3)*dotB37*M2*sqr(sh)+a1*dotB38*M2*sqr(sh)-2.*a12*dotB38*M2*sqr(sh)+pow(a1,3)*dotB38*M2*sqr(sh)+
			      a1*dotB39*M2*sqr(sh)-2.*a12*dotB39*M2*sqr(sh)+pow(a1,3)*dotB39*M2*sqr(sh)-a1*dotB40*M2*sqr(sh)+
			      2.*a12*dotB40*M2*sqr(sh)-pow(a1,3)*dotB40*M2*sqr(sh)+4.*a1*dotB13*pow<3,1>(sh)+
			      4.*a1*dotB23*dotB3*pow<3,1>(sh)+4.*dotB19*pow<3,1>(M)*th-8.*dotB10*dotB2*pow<3,1>(M)*th+
			      12.*a1*dotB10*dotB2*pow<3,1>(M)*th-8.*a12*dotB10*dotB2*pow<3,1>(M)*th+4.*dotB21*pow<3,1>(M)*th+
			      4.*a1*dotB21*pow<3,1>(M)*th-2.*dotB2*dotB32*pow<3,1>(M)*th-2.*dotB2*dotB33*pow<3,1>(M)*th+
			      2.*dotB2*dotB34*pow<3,1>(M)*th-8.*a1*dotB25*dotB4*pow<3,1>(M)*th-8.*dotB2*dotB7*pow<3,1>(M)*th+
			      8.*a1*dotB2*dotB7*pow<3,1>(M)*th+4.*a1*dotB2*dotB8*pow<3,1>(M)*th-8.*a12*dotB2*dotB8*pow<3,1>(M)*th+
			      4.*dotB15*pow<4,1>(M)*th-4.*a1*dotB15*pow<4,1>(M)*th+8.*pow(a1,3)*dotB15*pow<4,1>(M)*th-
			      2.*dotB2*dotB26*pow<4,1>(M)*th+6.*a1*dotB2*dotB26*pow<4,1>(M)*th-4.*a12*dotB2*dotB26*pow<4,1>(M)*th+
			      2.*pow(a1,3)*dotB2*dotB26*pow<4,1>(M)*th-2.*dotB2*dotB27*pow<4,1>(M)*th+6.*a1*dotB2*dotB27*pow<4,1>(M)*th-
			      4.*a12*dotB2*dotB27*pow<4,1>(M)*th+2.*pow(a1,3)*dotB2*dotB27*pow<4,1>(M)*th+2.*dotB2*dotB28*pow<4,1>(M)*th-
			      6.*a1*dotB2*dotB28*pow<4,1>(M)*th+4.*a12*dotB2*dotB28*pow<4,1>(M)*th-2.*pow(a1,3)*dotB2*dotB28*pow<4,1>(M)*th+
			      8.*pow(a1,3)*dotB25*dotB3*pow<4,1>(M)*th+4.*dotB35*pow<4,1>(M)*th-6.*a1*dotB35*pow<4,1>(M)*th+
			      4.*dotB36*pow<4,1>(M)*th-6.*a1*dotB36*pow<4,1>(M)*th-4.*dotB37*pow<4,1>(M)*th+6.*a1*dotB37*pow<4,1>(M)*th+
			      3.*dotB38*pow<4,1>(M)*th-5.*a1*dotB38*pow<4,1>(M)*th+3.*dotB39*pow<4,1>(M)*th-5.*a1*dotB39*pow<4,1>(M)*th-
			      3.*dotB40*pow<4,1>(M)*th+5.*a1*dotB40*pow<4,1>(M)*th-2.*a1*dotB11*pow<5,1>(M)*th+4.*a12*dotB11*pow<5,1>(M)*th+
			      2.*a1*dotB12*pow<5,1>(M)*th-4.*a12*dotB12*pow<5,1>(M)*th-a1*dotB29*pow<5,1>(M)*th-a12*dotB29*pow<5,1>(M)*th-
			      a1*dotB30*pow<5,1>(M)*th-a12*dotB30*pow<5,1>(M)*th+a1*dotB31*pow<5,1>(M)*th+a12*dotB31*pow<5,1>(M)*th-
			      2.*a1*dotB9*pow<5,1>(M)*th+4.*a12*dotB9*pow<5,1>(M)*th-8.*a1*dotB16*dotB2*M2*th+4.*a1*dotB18*dotB2*M2*th+
			      8.*a1*dotB2*dotB22*dotB3*M2*th+8.*a1*dotB16*dotB2*sh*th+4.*a1*dotB18*dotB2*sh*th-
			      8.*a1*dotB2*dotB22*dotB3*sh*th-4.*dotB19*M*sh*th+8.*dotB10*dotB2*M*sh*th-4.*a1*dotB10*dotB2*M*sh*th-
			      4.*dotB21*M*sh*th+2.*dotB2*dotB32*M*sh*th+2.*dotB2*dotB33*M*sh*th-2.*dotB2*dotB34*M*sh*th+
			      8.*dotB2*dotB7*M*sh*th-8.*a1*dotB2*dotB7*M*sh*th+4.*a1*dotB2*dotB8*M*sh*th-4.*a12*dotB11*pow<3,1>(M)*sh*th+
			      4.*a12*dotB12*pow<3,1>(M)*sh*th+2.*a1*dotB29*pow<3,1>(M)*sh*th+a12*dotB29*pow<3,1>(M)*sh*th+
			      2.*a1*dotB30*pow<3,1>(M)*sh*th+a12*dotB30*pow<3,1>(M)*sh*th-2.*a1*dotB31*pow<3,1>(M)*sh*th-
			      a12*dotB31*pow<3,1>(M)*sh*th-4.*a12*dotB9*pow<3,1>(M)*sh*th-4.*dotB15*M2*sh*th+8.*a1*dotB15*M2*sh*th+
			      2.*dotB2*dotB26*M2*sh*th-6.*a1*dotB2*dotB26*M2*sh*th+4.*a12*dotB2*dotB26*M2*sh*th-
			      2.*pow(a1,3)*dotB2*dotB26*M2*sh*th+2.*dotB2*dotB27*M2*sh*th-6.*a1*dotB2*dotB27*M2*sh*th+
			      4.*a12*dotB2*dotB27*M2*sh*th-2.*pow(a1,3)*dotB2*dotB27*M2*sh*th-2.*dotB2*dotB28*M2*sh*th+
			      6.*a1*dotB2*dotB28*M2*sh*th-4.*a12*dotB2*dotB28*M2*sh*th+2.*pow(a1,3)*dotB2*dotB28*M2*sh*th+
			      4.*a1*dotB25*dotB3*M2*sh*th-4.*dotB35*M2*sh*th+8.*a1*dotB35*M2*sh*th-4.*dotB36*M2*sh*th+
			      8.*a1*dotB36*M2*sh*th+4.*dotB37*M2*sh*th-8.*a1*dotB37*M2*sh*th-3.*dotB38*M2*sh*th+
			      6.*a1*dotB38*M2*sh*th-3.*dotB39*M2*sh*th+6.*a1*dotB39*M2*sh*th+3.*dotB40*M2*sh*th-
			      6.*a1*dotB40*M2*sh*th-4.*a1*dotB15*sqr(sh)*th-4.*a1*dotB25*dotB3*sqr(sh)*th-2.*a1*dotB35*sqr(sh)*th-
			      2.*a1*dotB36*sqr(sh)*th+2.*a1*dotB37*sqr(sh)*th-a1*dotB38*sqr(sh)*th-a1*dotB39*sqr(sh)*th+
			      a1*dotB40*sqr(sh)*th+2.*a1*dotB11*M*sqr(sh)*th-2.*a1*dotB12*M*sqr(sh)*th-a1*dotB29*M*sqr(sh)*th-
			      a1*dotB30*M*sqr(sh)*th+a1*dotB31*M*sqr(sh)*th+2.*a1*dotB9*M*sqr(sh)*th+
			      (a12*M2*(-16.*dotB15*M2+M*(-16.*dotB25*dotB3*M-(dotB29+dotB30-dotB31)*(M2-sh))-
				       4.*dotB2*(dotB26+dotB27-dotB28)*(M2-sh))+
			       2.*pow(a1,3)*M2*((-8.*dotB14+4.*dotB15+dotB2*(dotB26+dotB27-dotB28)+
						 4.*(-2.*dotB24+dotB25)*dotB3)*M2-dotB2*(dotB26+dotB27-dotB28)*sh)+
			       M*(M2-sh)*(-(M*(8.*dotB13-4.*dotB14+8.*(dotB23-dotB24)*dotB3+(dotB29+dotB30-dotB31)*M))+
					  (dotB29+dotB30-dotB31)*th)+
			       2.*a1*(-2.*dotB20*pow<3,1>(M)+4.*dotB24*dotB4*pow<3,1>(M)+6.*dotB13*pow<4,1>(M)-8.*dotB14*pow<4,1>(M)+
				      4.*dotB15*pow<4,1>(M)+dotB2*dotB26*pow<4,1>(M)+dotB2*dotB27*pow<4,1>(M)-dotB2*dotB28*pow<4,1>(M)+
				      6.*dotB23*dotB3*pow<4,1>(M)-10.*dotB24*dotB3*pow<4,1>(M)+4.*dotB25*dotB3*pow<4,1>(M)+
				      dotB29*pow<5,1>(M)+dotB30*pow<5,1>(M)-dotB31*pow<5,1>(M)+8.*dotB22*dotB6*M2-dotB29*pow<3,1>(M)*sh-
				      dotB30*pow<3,1>(M)*sh+dotB31*pow<3,1>(M)*sh-8.*dotB13*M2*sh+8.*dotB14*M2*sh-dotB2*dotB26*M2*sh-
				      dotB2*dotB27*M2*sh+dotB2*dotB28*M2*sh-8.*dotB23*dotB3*M2*sh+10.*dotB24*dotB3*M2*sh+2.*dotB13*sqr(sh)+
				      2.*dotB23*dotB3*sqr(sh)+2.*dotB17*dotB2*(-3.*M2+sh)-
				      (4.*dotB15+dotB2*(dotB26+dotB27-dotB28)+4.*dotB25*dotB3)*M2*th+
				      dotB2*(dotB26+dotB27-dotB28)*sh*th))*uh+8.*a1*(dotB14+dotB24*dotB3)*M2*sqr(uh)+
			      2.*dotB1*(6.*a1*dotB18*pow<4,1>(M)-12.*a12*dotB18*pow<4,1>(M)+6.*pow(a1,3)*dotB18*pow<4,1>(M)+
					8.*dotB22*dotB3*pow<4,1>(M)-12.*a1*dotB22*dotB3*pow<4,1>(M)+8.*a12*dotB22*dotB3*pow<4,1>(M)-
					4.*pow(a1,3)*dotB22*dotB3*pow<4,1>(M)+4.*dotB10*pow<5,1>(M)-14.*a1*dotB10*pow<5,1>(M)+
					20.*a12*dotB10*pow<5,1>(M)-14.*pow(a1,3)*dotB10*pow<5,1>(M)+4.*pow(a1,4)*dotB10*pow<5,1>(M)+
					8.*dotB7*pow<5,1>(M)-20.*a1*dotB7*pow<5,1>(M)+16.*a12*dotB7*pow<5,1>(M)-
					4.*pow(a1,3)*dotB7*pow<5,1>(M)-4.*dotB8*pow<5,1>(M)+6.*a1*dotB8*pow<5,1>(M)+
					4.*a12*dotB8*pow<5,1>(M)-10.*pow(a1,3)*dotB8*pow<5,1>(M)+4.*pow(a1,4)*dotB8*pow<5,1>(M)+
					dotB26*pow<6,1>(M)-4.*a1*dotB26*pow<6,1>(M)+5.*a12*dotB26*pow<6,1>(M)-
					pow(a1,3)*dotB26*pow<6,1>(M)-2.*pow(a1,4)*dotB26*pow<6,1>(M)+pow(a1,5)*dotB26*pow<6,1>(M)+
					dotB27*pow<6,1>(M)-4.*a1*dotB27*pow<6,1>(M)+5.*a12*dotB27*pow<6,1>(M)-
					pow(a1,3)*dotB27*pow<6,1>(M)-2.*pow(a1,4)*dotB27*pow<6,1>(M)+pow(a1,5)*dotB27*pow<6,1>(M)-
					dotB28*pow<6,1>(M)+4.*a1*dotB28*pow<6,1>(M)-5.*a12*dotB28*pow<6,1>(M)+pow(a1,3)*dotB28*pow<6,1>(M)+
					2.*pow(a1,4)*dotB28*pow<6,1>(M)-pow(a1,5)*dotB28*pow<6,1>(M)-4.*dotB10*pow<3,1>(M)*sh+
					10.*a1*dotB10*pow<3,1>(M)*sh-8.*a12*dotB10*pow<3,1>(M)*sh+2.*pow(a1,3)*dotB10*pow<3,1>(M)*sh-
					8.*dotB7*pow<3,1>(M)*sh+20.*a1*dotB7*pow<3,1>(M)*sh-16.*a12*dotB7*pow<3,1>(M)*sh+
					4.*pow(a1,3)*dotB7*pow<3,1>(M)*sh+4.*dotB8*pow<3,1>(M)*sh-10.*a1*dotB8*pow<3,1>(M)*sh+
					8.*a12*dotB8*pow<3,1>(M)*sh-2.*pow(a1,3)*dotB8*pow<3,1>(M)*sh-dotB26*pow<4,1>(M)*sh+
					4.*a1*dotB26*pow<4,1>(M)*sh-5.*a12*dotB26*pow<4,1>(M)*sh+pow(a1,3)*dotB26*pow<4,1>(M)*sh+
					2.*pow(a1,4)*dotB26*pow<4,1>(M)*sh-pow(a1,5)*dotB26*pow<4,1>(M)*sh-dotB27*pow<4,1>(M)*sh+
					4.*a1*dotB27*pow<4,1>(M)*sh-5.*a12*dotB27*pow<4,1>(M)*sh+pow(a1,3)*dotB27*pow<4,1>(M)*sh+
					2.*pow(a1,4)*dotB27*pow<4,1>(M)*sh-pow(a1,5)*dotB27*pow<4,1>(M)*sh+dotB28*pow<4,1>(M)*sh-
					4.*a1*dotB28*pow<4,1>(M)*sh+5.*a12*dotB28*pow<4,1>(M)*sh-pow(a1,3)*dotB28*pow<4,1>(M)*sh-
					2.*pow(a1,4)*dotB28*pow<4,1>(M)*sh+pow(a1,5)*dotB28*pow<4,1>(M)*sh-2.*a1*dotB18*M2*sh+
					4.*a12*dotB18*M2*sh-2.*pow(a1,3)*dotB18*M2*sh-8.*dotB22*dotB3*M2*sh+12.*a1*dotB22*dotB3*M2*sh-
					8.*a12*dotB22*dotB3*M2*sh+4.*pow(a1,3)*dotB22*dotB3*M2*sh+
					M2*(-6.*a1*dotB18+4.*a1*dotB22*dotB3-4.*dotB10*M-4.*dotB7*M+
					    2.*a1*((3.-2.*a1)*dotB10+2.*dotB7+dotB8-2.*a1*dotB8)*M+
					    (-1.+2.*a12)*a2*(dotB26+dotB27-dotB28)*M2)*th+
					(2.*a1*(dotB18-2.*dotB22*dotB3)+4.*dotB10*M+4.*dotB7*M-2.*a1*(dotB10+2.*dotB7-dotB8)*M-
					 (-1.+2.*a12)*a2*(dotB26+dotB27-dotB28)*M2)*sh*th+a1*(dotB26+dotB27-dotB28)*(M2-sh)*sqr(th)-
					4.*dotB16*(M2-sh)*((1.+pow(a1,3))*M2-a1*(sh+uh))+
					2.*dotB17*(-((-2.+3.*a1+pow(a1,3))*pow<4,1>(M))+a1*sh*(sh+uh)+
						   M2*(-((2.-2.*a1+pow(a1,3))*sh)+a1*uh)))))/(a1*a2*M*sqr(M2-sh)*sqr(-(a22*M2)+th));
	    diag[3]+=mix2*(-4.*(-4.*a1*(dotB19+dotB21-2.*dotB1*(dotB10+dotB7))*M2+
			       (4.*(dotB19+dotB21-2.*dotB1*(dotB10+dotB7))+(-1.+2.*a1)*(2.*dotB1*(dotB26+dotB27-dotB28)-dotB38-dotB39+dotB40)*M)*(-(a12*M2)+uh)+
			       (dotB29+dotB30-dotB31)*sqr(-(a12*M2)+uh)))/(a2*pow<3,1>(-(a12*M2)+uh));
	    diag[4]+=mix2*(-2.*(8.*pow(a1,5)*(4.*dotB14+4.*dotB15+5.*dotB24*dotB3+7.*dotB25*dotB3)*pow<5,1>(M)+
			       pow(a1,4)*pow<3,1>(M)*(-8.*(dotB1+dotB2)*(dotB10+dotB8)*M+(dotB29+dotB30-dotB31)*pow<3,1>(M)-
						      4.*(10.*dotB14+6.*dotB15+(dotB1+dotB2)*(dotB26+dotB27-dotB28)+
							  2.*(7.*dotB24+9.*dotB25)*dotB3)*M2+4.*(dotB1+dotB2)*(dotB26+dotB27-dotB28)*sh-(dotB29+dotB30-dotB31)*M*sh)+
			       2.*pow(a1,3)*M2*(4.*dotB17*dotB2*M+12.*dotB18*dotB2*M-16.*dotB22*dotB6*M+10.*dotB14*pow<3,1>(M)-2.*dotB15*pow<3,1>(M)+
						dotB2*dotB26*pow<3,1>(M)+dotB2*dotB27*pow<3,1>(M)-dotB2*dotB28*pow<3,1>(M)-4.*dotB23*dotB3*pow<3,1>(M)+
						16.*dotB24*dotB3*pow<3,1>(M)+20.*dotB25*dotB3*pow<3,1>(M)-8.*dotB35*pow<3,1>(M)-8.*dotB36*pow<3,1>(M)+
						8.*dotB37*pow<3,1>(M)-3.*dotB38*pow<3,1>(M)-3.*dotB39*pow<3,1>(M)+3.*dotB40*pow<3,1>(M)+2.*dotB11*pow<4,1>(M)-
						2.*dotB12*pow<4,1>(M)-dotB29*pow<4,1>(M)-dotB30*pow<4,1>(M)+dotB31*pow<4,1>(M)+2.*dotB9*pow<4,1>(M)+
						dotB1*M*(-4.*dotB17-12.*dotB18+M*(8.*dotB10+dotB32+dotB33-dotB34+4.*(dotB7+dotB8)+(dotB26+dotB27-dotB28)*M))+
						8.*dotB10*dotB2*M2+2.*dotB20*M2+6.*dotB21*M2+dotB2*dotB32*M2+dotB2*dotB33*M2-dotB2*dotB34*M2+
						4.*dotB2*dotB7*M2+4.*dotB2*dotB8*M2-4.*dotB10*dotB2*sh-dotB2*dotB32*sh-dotB2*dotB33*sh+dotB2*dotB34*sh-
						4.*dotB2*dotB7*sh-dotB2*dotB26*M*sh-dotB2*dotB27*M*sh+dotB2*dotB28*M*sh+4.*dotB23*dotB3*M*sh+
						8.*dotB35*M*sh+8.*dotB36*M*sh-8.*dotB37*M*sh+3.*dotB38*M*sh+3.*dotB39*M*sh-3.*dotB40*M*sh-
						dotB1*(4.*dotB10+dotB32+dotB33-dotB34+4.*dotB7+(dotB26+dotB27-dotB28)*M)*sh-
						2.*dotB11*M2*sh+2.*dotB12*M2*sh+dotB29*M2*sh+dotB30*M2*sh-dotB31*M2*sh-
						2.*dotB9*M2*sh-6.*dotB14*M*th-2.*dotB15*M*th-8.*dotB24*dotB3*M*th-12.*dotB25*dotB3*M*th-
						2.*(13.*dotB14+15.*dotB15+16.*dotB24*dotB3+22.*dotB25*dotB3)*M*uh)+
			       uh*(-8.*dotB1*(dotB17+dotB18)*M+8.*dotB17*dotB2*M+8.*dotB18*dotB2*M-16.*dotB22*dotB6*M+4.*dotB14*pow<3,1>(M)+
				   4.*dotB15*pow<3,1>(M)+8.*dotB24*dotB3*pow<3,1>(M)+8.*dotB25*dotB3*pow<3,1>(M)-8.*dotB35*pow<3,1>(M)-
				   8.*dotB36*pow<3,1>(M)+8.*dotB37*pow<3,1>(M)-3.*dotB38*pow<3,1>(M)-3.*dotB39*pow<3,1>(M)+
				   3.*dotB40*pow<3,1>(M)+4.*dotB11*pow<4,1>(M)-4.*dotB12*pow<4,1>(M)-dotB29*pow<4,1>(M)-dotB30*pow<4,1>(M)+
				   dotB31*pow<4,1>(M)+4.*dotB9*pow<4,1>(M)-4.*dotB19*M2+8.*dotB10*dotB2*M2+8.*dotB20*M2+4.*dotB21*M2+
				   2.*dotB2*dotB32*M2+2.*dotB2*dotB33*M2-2.*dotB2*dotB34*M2+8.*dotB23*dotB4*M2-8.*dotB24*dotB4*M2+
				   8.*dotB22*dotB5*M2+8.*dotB2*dotB8*M2+4.*dotB19*sh-8.*dotB10*dotB2*sh-4.*dotB20*sh-2.*dotB2*dotB32*sh-
				   2.*dotB2*dotB33*sh+2.*dotB2*dotB34*sh-8.*dotB23*dotB4*sh+8.*dotB24*dotB4*sh-8.*dotB22*dotB5*sh-
				   8.*dotB2*dotB8*sh+8.*dotB35*M*sh+8.*dotB36*M*sh-8.*dotB37*M*sh+3.*dotB38*M*sh+3.*dotB39*M*sh-
				   3.*dotB40*M*sh-4.*dotB11*M2*sh+4.*dotB12*M2*sh+dotB29*M2*sh+dotB30*M2*sh-dotB31*M2*sh-4.*dotB9*M2*sh-
				   4.*dotB14*M*th-4.*dotB15*M*th-8.*dotB24*dotB3*M*th-8.*dotB25*dotB3*M*th-4.*dotB11*M2*th+
				   4.*dotB12*M2*th+dotB29*M2*th+dotB30*M2*th-dotB31*M2*th-4.*dotB9*M2*th+4.*dotB11*sh*th-4.*dotB12*sh*th-
				   dotB29*sh*th-dotB30*sh*th+dotB31*sh*th+4.*dotB9*sh*th-
				   4.*(3.*(dotB14+dotB15)+4.*(dotB24+dotB25)*dotB3)*M*uh)+
			       a12*M*(8.*dotB1*dotB10*pow<3,1>(M)+4.*dotB19*pow<3,1>(M)-8.*dotB20*pow<3,1>(M)-12.*dotB21*pow<3,1>(M)-
				      2.*dotB2*dotB32*pow<3,1>(M)-2.*dotB2*dotB33*pow<3,1>(M)+2.*dotB2*dotB34*pow<3,1>(M)-8.*dotB23*dotB4*pow<3,1>(M)+
				      8.*dotB24*dotB4*pow<3,1>(M)-8.*dotB22*dotB5*pow<3,1>(M)+8.*dotB1*dotB7*pow<3,1>(M)+8.*dotB2*dotB7*pow<3,1>(M)-
				      8.*dotB2*dotB8*pow<3,1>(M)-8.*dotB13*pow<4,1>(M)-4.*dotB14*pow<4,1>(M)+12.*dotB15*pow<4,1>(M)-
				      8.*dotB24*dotB3*pow<4,1>(M)-8.*dotB25*dotB3*pow<4,1>(M)+8.*dotB35*pow<4,1>(M)+8.*dotB36*pow<4,1>(M)-
				      8.*dotB37*pow<4,1>(M)+3.*dotB38*pow<4,1>(M)+3.*dotB39*pow<4,1>(M)-3.*dotB40*pow<4,1>(M)-
				      4.*dotB11*pow<5,1>(M)+4.*dotB12*pow<5,1>(M)+dotB29*pow<5,1>(M)+dotB30*pow<5,1>(M)-dotB31*pow<5,1>(M)-
				      4.*dotB9*pow<5,1>(M)+8.*dotB1*dotB17*M2+24.*dotB1*dotB18*M2-8.*dotB17*dotB2*M2-24.*dotB18*dotB2*M2+
				      32.*dotB22*dotB6*M2-8.*dotB1*dotB10*M*sh-4.*dotB19*M*sh+4.*dotB20*M*sh+2.*dotB2*dotB32*M*sh+
				      2.*dotB2*dotB33*M*sh-2.*dotB2*dotB34*M*sh+8.*dotB23*dotB4*M*sh-8.*dotB24*dotB4*M*sh+
				      8.*dotB22*dotB5*M*sh-8.*dotB1*dotB7*M*sh-8.*dotB2*dotB7*M*sh+8.*dotB2*dotB8*M*sh+4.*dotB11*pow<3,1>(M)*sh-
				      4.*dotB12*pow<3,1>(M)*sh-dotB29*pow<3,1>(M)*sh-dotB30*pow<3,1>(M)*sh+dotB31*pow<3,1>(M)*sh+
				      4.*dotB9*pow<3,1>(M)*sh+8.*dotB13*M2*sh-8.*dotB35*M2*sh-8.*dotB36*M2*sh+8.*dotB37*M2*sh-
				      3.*dotB38*M2*sh-3.*dotB39*M2*sh+3.*dotB40*M2*sh+4.*dotB11*pow<3,1>(M)*th-4.*dotB12*pow<3,1>(M)*th-
				      dotB29*pow<3,1>(M)*th-dotB30*pow<3,1>(M)*th+dotB31*pow<3,1>(M)*th+4.*dotB9*pow<3,1>(M)*th+
				      4.*dotB14*M2*th-12.*dotB15*M2*th+8.*dotB24*dotB3*M2*th+8.*dotB25*dotB3*M2*th-4.*dotB11*M*sh*th+
				      4.*dotB12*M*sh*th+dotB29*M*sh*th+dotB30*M*sh*th-dotB31*M*sh*th-4.*dotB9*M*sh*th+
				      M*(8.*(dotB1+dotB2)*(dotB10+dotB8)+52.*(dotB14+dotB15)*M+4.*(dotB1+dotB2)*(dotB26+dotB27-dotB28)*M+
					 8.*(9.*dotB24+11.*dotB25)*dotB3*M-(dotB29+dotB30-dotB31)*M2)*uh+
				      (-4.*(dotB1+dotB2)*(dotB26+dotB27-dotB28)+(dotB29+dotB30-dotB31)*M)*sh*uh)+
			       2.*a1*(-2.*M*(2.*dotB16*dotB2*(M2-sh)+M*(-2.*dotB18*dotB2*M+2.*dotB22*dotB6*M+
									(dotB19-dotB21-dotB13*M+dotB15*M)*M2-dotB19*sh+dotB13*M*sh)+
					     (2.*dotB18*dotB2-2.*dotB22*dotB6+M*(dotB21+dotB13*M-2.*dotB15*M)-dotB13*sh)*th+
					     dotB15*sqr(th))+4.*dotB1*M*(dotB16*(M2-sh)+2.*dotB22*dotB3*(-M2+sh)+dotB18*(-M2+th))+
				      dotB1*(M*(4.*dotB17+8.*dotB18-M*(8.*dotB10+dotB32+dotB33-dotB34+4.*(dotB7+dotB8)+(dotB26+dotB27-dotB28)*M))+
					     (4.*dotB10+dotB32+dotB33-dotB34+4.*dotB7+(dotB26+dotB27-dotB28)*M)*sh)*uh-
				      (4.*dotB17*dotB2*M+8.*dotB18*dotB2*M-12.*dotB22*dotB6*M-2.*dotB13*pow<3,1>(M)+10.*dotB14*pow<3,1>(M)+
				       10.*dotB15*pow<3,1>(M)+dotB2*dotB26*pow<3,1>(M)+dotB2*dotB27*pow<3,1>(M)-dotB2*dotB28*pow<3,1>(M)-
				       4.*dotB23*dotB3*pow<3,1>(M)+16.*dotB24*dotB3*pow<3,1>(M)+20.*dotB25*dotB3*pow<3,1>(M)-
				       8.*dotB35*pow<3,1>(M)-8.*dotB36*pow<3,1>(M)+8.*dotB37*pow<3,1>(M)-3.*dotB38*pow<3,1>(M)-
				       3.*dotB39*pow<3,1>(M)+3.*dotB40*pow<3,1>(M)+2.*dotB11*pow<4,1>(M)-2.*dotB12*pow<4,1>(M)-
				       dotB29*pow<4,1>(M)-dotB30*pow<4,1>(M)+dotB31*pow<4,1>(M)+2.*dotB9*pow<4,1>(M)+8.*dotB10*dotB2*M2+
				       2.*dotB20*M2+4.*dotB21*M2+dotB2*dotB32*M2+dotB2*dotB33*M2-dotB2*dotB34*M2+4.*dotB2*dotB7*M2+
				       4.*dotB2*dotB8*M2-4.*dotB10*dotB2*sh-dotB2*dotB32*sh-dotB2*dotB33*sh+dotB2*dotB34*sh-
				       4.*dotB2*dotB7*sh+2.*dotB13*M*sh-dotB2*dotB26*M*sh-dotB2*dotB27*M*sh+dotB2*dotB28*M*sh+
				       4.*dotB23*dotB3*M*sh+8.*dotB35*M*sh+8.*dotB36*M*sh-8.*dotB37*M*sh+3.*dotB38*M*sh+3.*dotB39*M*sh-
				       3.*dotB40*M*sh-2.*dotB11*M2*sh+2.*dotB12*M2*sh+dotB29*M2*sh+dotB30*M2*sh-dotB31*M2*sh-
				       2.*dotB9*M2*sh-2.*(3.*dotB14+3.*dotB15+4.*dotB24*dotB3+6.*dotB25*dotB3)*M*th)*uh+
				      2.*(5.*dotB14+6.*dotB15+6.*dotB24*dotB3+8.*dotB25*dotB3)*M*sqr(uh))))/(a1*a2*sqr(M2-sh)*sqr(-(a12*M2)+uh));
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

