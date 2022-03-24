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
	    diag[0]=(-4.*M*(a12*pow<3,1>(M)*(-((-11.-2.*a1+8.*a12)*(2.*dot1*dot3+dot5))+a1*(-4.+a1+4.*a12)*dot4*M)+
			    M*((-3.+2.*a1)*(1.+4.*a1)*(2.*dot1*dot3+dot5)+2.*a1*(2.+a1-4.*a12)*dot4*M)*sh+(-3.+4.*a1)*dot4*sqr(sh)))/(a2*pow<3,1>(a12*M2-sh));
	    diag[1]=2.*M*(2.*pow(a1,5)*pow<3,1>(M)*(-8.*dot1*dot3+8.*dot2*dot3+dot4*M)*(M2-th)+
			  pow(a1,4)*M2*(M2*(16.*(-dot1+dot2)*dot6+4.*(dot1*dot3-dot2*dot3+dot5)*M+3.*dot4*M2)+
					(8.*(dot1-dot2)*dot6-4.*(dot1*dot3-dot2*dot3+dot5)*M-3.*dot4*M2)*th)-
			  sh*(M2-th)*(2.*dot2*dot6+M*(dot5-dot4*M)+dot4*uh)-
			  2.*a1*(dot1*dot3*M*(4.*M2-sh-4.*th)*(M2-th)+dot2*sh*(dot3*pow<3,1>(M)-9.*dot6*M2+5.*dot6*th-dot3*M*th)+
				 dot1*dot6*(4.*pow<4,1>(M)+M2*(7.*sh-4.*th)-3.*sh*th)-2.*dot5*M*(M2-th)*(M2-uh)-dot4*sh*(M2-th)*(2.*M2+uh))+
			  a12*(4.*dot1*(4.*dot6*M2*(M2+sh)-2.*dot6*(2.*M2+sh)*th+5.*dot3*M*sh*(-M2+th))+
			       dot2*(20.*dot3*M*sh*(M2-th)+dot6*(-6.*pow<4,1>(M)+8.*sh*th+M2*(-16.*sh+6.*th)))+
			       M*(M2-th)*(dot5*(M2-4.*sh)+dot4*M*(-M2-3.*sh+uh)))+
			  2.*pow(a1,3)*M*(dot2*(-(dot3*(7.*M2+8.*sh)*(M2-th))+dot6*M*(-5.*M2+th))+
					  dot1*(dot3*(7.*M2+8.*sh)*(M2-th)+dot6*M*(3.*M2+th))-M*(M2-th)*(2.*dot5*M+dot4*(2.*M2+sh+uh))))/(a1*a2*sqr(-(a12*M2)+sh)*sqr(M2-th));
	    diag[2]=2.*M*(-4.*M*(M2-th)*(M*(2.*a1*(-(a2*dot1)-a1*dot2)*dot6+(2.*(2.*a1*a22*(dot1-dot2)+dot2)*dot3+dot5-2.*a1*a2*dot5)*M)-(2.*dot2*dot3+dot5)*th)-
			  ((2.*pow(a1,3)*M*(8.*dot1*dot3-8.*dot2*dot3+dot4*M)*(M2-th)+(2.*dot1*dot3*M-dot5*M-2.*dot2*(dot6+dot3*M)+dot4*sh)*(M2-th)+
			    a12*(M2*(16.*(dot1-dot2)*dot6-4.*(7.*dot1*dot3-7.*dot2*dot3+dot5)*M-9.*dot4*M2)+
				 (8.*(-dot1+dot2)*dot6+4.*(7.*dot1*dot3-7.*dot2*dot3+dot5)*M+9.*dot4*M2)*th)+
			    2.*a1*((2.*M*(dot5+2.*dot4*M)-dot4*sh)*(M2-th)+dot1*(5.*dot3*pow<3,1>(M)-7.*dot6*M2+3.*dot6*th-5.*dot3*M*th)+
				   dot2*(-5.*dot3*pow<3,1>(M)+9.*dot6*M2-5.*dot6*th+5.*dot3*M*th)))*(a22*M2-uh))/a2)/(a1*sqr(M2-th)*sqr(-(a22*M2)+uh));
	    diag[3]=4.*M*(a22*pow<3,1>(M)*((-5.+2.*a1*(-7.+4.*a1))*(2.*dot2*dot3+dot5)-(1.+a1*(-9.+4.*a1))*a2*dot4*M)+
			  M*(-((1.+2.*a1)*(-5.+4.*a1)*(2.*dot2*dot3+dot5))+2.*(1.+a1*(-7.+4.*a1))*a2*dot4*M)*uh+(-1.+4.*a1)*dot4*sqr(uh))/(a1*pow<3,1>(a22*M2-uh));
	    diag[4]=-2.*M*(16.*pow(a1,7)*(dot1-dot2)*dot3*pow<5,1>(M)+
			   4.*pow(a1,6)*pow<4,1>(M)*(2.*dot1*dot6-2.*dot2*dot6-15.*dot1*dot3*M+13.*dot2*dot3*M+2.*dot5*M)+
			   sh*(-(M*(2.*dot1*dot3*M2+4.*dot5*sh+dot4*M*(M2+sh)+dot2*(-2.*dot6*M+8.*dot3*sh)))+
			       (-2.*dot2*dot6+2.*M*(dot1*dot3+dot4*M)+dot4*sh)*uh-dot4*sqr(uh))+
			   pow(a1,4)*M2*(-14.*dot2*dot6*M2+4.*dot5*M*(M2-2.*(2.*sh+th))+dot4*M2*(M2+5.*sh-5.*uh)+8.*dot2*dot6*(sh+uh)-
					 4.*dot2*dot3*M*(4.*M2+17.*sh+uh)+2.*dot1*dot3*M*(17.*M2+38.*sh+6.*uh)+4.*dot1*dot6*(M2-2.*(sh+uh)))-
			   2.*pow(a1,5)*pow<3,1>(M)*(-13.*dot2*dot6*M+17.*dot2*dot3*M2+4.*dot5*M2+dot4*M*(sh-uh)-
						     8.*dot2*dot3*(sh+uh)+dot1*(11.*dot6*M-29.*dot3*M2+8.*dot3*(sh+uh)))-
			   2.*pow(a1,3)*M*(8.*dot5*M2*(M2-sh-th)+dot2*dot6*M*(5.*M2+17.*sh+uh)+
					   dot4*M*(pow<4,1>(M)-sqr(sh)+2.*M2*(sh-uh)+sqr(uh))+
					   dot2*dot3*(-7.*pow<4,1>(M)+8.*sh*uh+M2*(-37.*sh+7.*uh))+
					   dot1*(dot6*M*(-17.*M2-15.*sh+uh)+dot3*(45.*pow<4,1>(M)+M2*(65.*sh-19.*uh)-8.*sh*uh)))+
			   a12*(8.*dot5*M*sh*(sh+th)+4.*dot2*dot3*M*sh*(-6.*M2+5.*uh)+dot2*dot6*(6.*pow<4,1>(M)+M2*(38.*sh-6.*uh)-8.*sh*uh)+
				dot4*M2*(pow<4,1>(M)-5.*sqr(sh)-2.*M2*uh+4.*sh*uh+sqr(uh))+
				dot1*(2.*dot3*M*(25.*pow<4,1>(M)+M2*(47.*sh-25.*uh)-14.*sh*uh)+4.*dot6*(-8.*pow<4,1>(M)+2.*sh*uh+M2*(-9.*sh+4.*uh))))+
			   2.*a1*(-2.*dot5*M*(pow<4,1>(M)+2.*sqr(sh)+4.*sh*th+sqr(th)-2.*M2*(2.*sh+th))+
				  sh*(dot4*(pow<4,1>(M)+2.*M2*(sh-uh)+uh*(-sh+uh))+dot2*(dot3*M*(M2+4.*sh-uh)+dot6*(-7.*M2+5.*uh)))+
				  dot1*(dot6*(4.*pow<4,1>(M)+M2*(7.*sh-4.*uh)-3.*sh*uh)+dot3*M*(-4.*pow<4,1>(M)+(5.*sh-4.*uh)*uh+M2*(-11.*sh+8.*uh)))))/(a1*a2*sqr(-(a12*M2)+sh)*sqr(-(a22*M2)+uh));
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
	    diag[0]=(4.*M*(a22*pow<3,1>(M)*((-5.+2.*a1*(-7.+4.*a1))*(2.*dot1*dot3-dot5)+(1.+a1*(-9.+4.*a1))*a2*dot4*M)+
			   M*(-((1.+2.*a1)*(-5.+4.*a1)*(2.*dot1*dot3-dot5))-2.*(1.+a1*(-7.+4.*a1))*a2*dot4*M)*sh+(1.-4.*a1)*dot4*sqr(sh)))/(a1*pow<3,1>(a22*M2-sh));
	    diag[1]=(-2.*M*(2.*pow(a1,5)*pow<3,1>(M)*(-8.*dot1*dot3+8.*dot2*dot3+dot4*M)*(M2-th)+4.*a1*M*(M2-th)*(2.*dot2*M*(dot6+dot3*M)+dot5*(M2-sh)-2.*dot2*dot3*th)+
			    pow(a1,4)*M2*(M2*(16.*(-dot1+dot2)*dot6+4.*(dot1*dot3-dot2*dot3+dot5)*M+3.*dot4*M2)+(8.*(dot1-dot2)*dot6-4.*(dot1*dot3-dot2*dot3+dot5)*M-3.*dot4*M2)*th)+
			    (2.*dot1*dot6-dot5*M+dot4*(M2-sh))*(M2-th)*uh+2.*a1*(M2*(-9.*dot1*dot6+7.*dot2*dot6+dot1*dot3*M-dot2*dot3*M+2.*dot4*M2+dot4*sh)-
										 (-5.*dot1*dot6+3.*dot2*dot6+dot1*dot3*M-dot2*dot3*M+2.*dot4*M2+dot4*sh)*th)*uh
			    +a12*(M2*(6.*dot1*dot6-16.*dot2*dot6+M*(dot5-dot4*M)+dot4*sh)*(M2-th)+
				  (M2*(16.*(dot1-dot2)*dot6-4.*(5.*dot1*dot3-5.*dot2*dot3+dot5)*M-3.*dot4*M2)+
				   (8.*(-dot1+dot2)*dot6+4.*(5.*dot1*dot3-5.*dot2*dot3+dot5)*M+3.*dot4*M2)*th)*uh)+
			    2.*pow(a1,3)*M*(-(M*(M2-th)*(2.*dot5*M+dot4*(2.*M2+sh+uh)))-dot2*(dot6*M*(3.*M2+th)+dot3*(M2-th)*(7.*M2+8.*uh))+
					    dot1*(dot6*(5.*pow<3,1>(M)-M*th)+dot3*(M2-th)*(7.*M2+8.*uh)))))/(a1*a2*sqr(M2-th)*sqr(-(a12*M2)+uh));
	    diag[2]=(-2.*M*(2.*pow(a1,5)*pow<3,1>(M)*(8.*dot2*dot3+dot4*M)*(M2-th)-(2.*dot2*dot3*M*(-M2+sh)+dot5*M*(3.*M2+sh-4.*th)+dot4*(M2-sh)*(2.*M2-sh-th))*(M2-th)+
			    pow(a1,4)*M2*(M*(4.*dot5-dot4*M)*(M2-th)+8.*dot2*dot6*(2.*M2-th)+76.*dot2*dot3*M*(-M2+th))-
			    2.*a1*dot2*(-(dot3*M*(11.*M2-5.*sh)*(M2-th))+dot6*(11.*pow<4,1>(M)+3.*sh*th-7.*M2*(sh+th)))+
			    2.*a1*(M2-th)*(dot5*M*(3.*M2+2.*sh-2.*th)+dot4*(pow<4,1>(M)+sh*(sh+th)-M2*(sh+2.*th)))+
			    a12*(M*(M2-th)*(-(dot5*(3.*M2+4.*sh))+dot4*M*(5.*M2+2.*sh+5.*th))+
				 2.*dot2*(30.*dot6*pow<4,1>(M)-dot3*M*(47.*M2-14.*sh)*(M2-th)+4.*dot6*sh*th-2.*dot6*M2*(4.*sh+9.*th)))+
			    2.*a2*dot1*(dot3*M*((3.+a1*(-3.+2.*a1)*(4.+a1*(-9.+4.*a1)))*M2+sh+2.*(3.-4.*a1)*a1*sh-4.*th)*(M2-th)+
					dot6*((-1.+(-2.+a1)*a1*(-5.+8.*a1))*pow<4,1>(M)+(-1.+4.*a1)*sh*th+M2*(sh-8.*a1*sh+th+a1*(-6.+(13.-4.*a1)*a1)*th)))-
			    2.*pow(a1,3)*M*(M*(M2-th)*(2.*dot5*M+dot4*(3.*M2+2.*sh+th))+dot2*(-(dot3*(65.*M2-8.*sh)*(M2-th))+3.*dot6*(9.*pow<3,1>(M)-5.*M*th)))))/(a1*a2*sqr(-(a22*M2)+sh)*sqr(M2-th));
	    diag[3]=-4.*M*(a12*pow<3,1>(M)*(-((-11.-2.*a1+8.*a12)*(2.*dot2*dot3-dot5))-a1*(-4.+a1+4.*a12)*dot4*M)+
			   M*((-3.+2.*a1)*(1.+4.*a1)*(2.*dot2*dot3-dot5)+2.*a1*(-2.+a1*(-1.+4.*a1))*dot4*M)*uh+(3.-4.*a1)*dot4*sqr(uh))/(a2*pow<3,1>(a12*M2-uh));
	    diag[4]=(-2.*M*(16.*pow(a1,7)*(-dot1+dot2)*dot3*pow<5,1>(M)+4.*pow(a1,6)*pow<4,1>(M)*(-2.*dot1*dot6+2.*dot2*dot6+13.*dot1*dot3*M-15.*dot2*dot3*M+2.*dot5*M)+
			    uh*((2.*dot1*dot6-2.*dot2*dot3*M+dot4*(M2-sh))*(M2-sh)+(-8.*dot1*dot3*M+4.*dot5*M+dot4*(M2-sh))*uh)+
			    2.*pow(a1,3)*M*(pow<3,1>(M)*(-5.*dot1*dot6+17.*dot2*dot6+7.*dot1*dot3*M-45.*dot2*dot3*M+8.*dot5*M+dot4*M2)-
					    M*((dot1+dot2)*dot6+(7.*dot1-19.*dot2)*dot3*M+2.*dot4*M2)*sh+dot4*M*sqr(sh)+
					    (M*(-17.*dot1*dot6+15.*dot2*dot6+37.*dot1*dot3*M-65.*dot2*dot3*M+16.*dot5*M+2.*dot4*M2)+8.*(-dot1+dot2)*dot3*sh)*uh-dot4*M*sqr(uh))+
			    a12*(M2*(6.*dot1*dot6*(M2-sh)-(16.*dot5*M+dot4*(M2-sh))*(M2-sh)+dot2*(50.*dot3*M*(M2-sh)+16.*dot6*(-2.*M2+sh)))+
				 2.*((19.*dot1*dot6-18.*dot2*dot6-12.*dot1*dot3*M+47.*dot2*dot3*M-16.*dot5*M)*M2-2.*(2.*dot1*dot6-2.*dot2*dot6-5.*dot1*dot3*M+7.*dot2*dot3*M-2.*dot5*M+dot4*M2)*sh)*uh+
				 5.*dot4*M2*sqr(uh))+2.*a1*(2.*M*(M2-sh)*(dot5*(M2-sh)+2.*dot2*(M*(dot6-dot3*M)+dot3*sh))-
							    (M2*(7.*dot1*dot6-7.*dot2*dot6-dot1*dot3*M+11.*dot2*dot3*M-4.*dot5*M+dot4*M2)+
							     (-5.*dot1*dot6+3.*dot2*dot6+dot1*dot3*M-5.*dot2*dot3*M+4.*dot5*M-2.*dot4*M2)*sh+dot4*sqr(sh))*uh+
							    (4.*dot1*dot3*M-2.*M*(dot5+dot4*M)+dot4*sh)*sqr(uh))-
			    2.*pow(a1,5)*pow<3,1>(M)*(-13.*dot1*dot6*M+11.*dot2*dot6*M+17.*dot1*dot3*M2-29.*dot2*dot3*M2+
						      12.*dot5*M2+dot4*M*(sh-uh)-8.*dot1*dot3*(sh+uh)+8.*dot2*dot3*(sh+uh))-
			    pow(a1,4)*M2*(M*(-12.*dot5*M2+8.*dot5*(sh+uh)+dot4*M*(M2-5.*sh+5.*uh))+2.*dot1*(7.*dot6*M2-4.*dot6*(sh+uh)+2.*dot3*M*(4.*M2+sh+17.*uh))-
					  2.*dot2*(dot3*M*(17.*M2+6.*sh+38.*uh)+2.*dot6*(M2-2.*(sh+uh))))))/(a1*a2*sqr(-(a22*M2)+sh)*sqr(-(a12*M2)+uh));
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
	  diag[0]=(-2.*M*(16.*pow(a1,7)*(dot1+dot2)*dot3*pow<5,1>(M)-4.*pow(a1,6)*pow<4,1>(M)*(2.*(dot1+dot2)*dot4+(15.*dot1*dot3+13.*dot2*dot3+2.*dot6)*M)+
			  uh*((M2-th)*(2.*dot2*dot4-M*(2.*dot1*dot3+dot5*M)+dot5*th)+(M*(8.*dot2*dot3+4.*dot6-dot5*M)+dot5*th)*uh)+
			  2.*pow(a1,3)*M*(M*(M*(-(M*(17.*dot1*dot4+5.*dot2*dot4+45.*dot1*dot3*M+7.*dot2*dot3*M-8.*dot6*M+dot5*M2))-8.*dot6*sh)+
					     (-(dot2*dot4)+7.*dot2*dot3*M+dot1*(dot4+19.*dot3*M)+2.*dot5*M2)*th-dot5*sqr(th))-
					  M*(15.*dot1*dot4+17.*dot2*dot4+65.*dot1*dot3*M+37.*dot2*dot3*M+
					     8.*dot6*M+2.*dot5*M2)*uh+8.*(dot1+dot2)*dot3*th*uh+dot5*M*sqr(uh))+
			  a12*(M2*(2.*dot1*(25.*dot3*M*(M2-th)+8.*dot4*(2.*M2-th))+(6.*dot2*dot4+dot5*(M2-th))*(M2-th))+
			       2.*(M*(M*(18.*dot1*dot4+19.*dot2*dot4+47.*dot1*dot3*M+12.*dot2*dot3*M)-4.*dot6*sh)-
				   2.*(2.*(dot1+dot2)*dot4+(7.*dot1+5.*dot2)*dot3*M-dot5*M2)*th)*uh-M*(8.*dot6+5.*dot5*M)*sqr(uh))+
			  2.*a1*(2.*dot6*M*sqr(M2-sh)-4.*dot1*M*(M2-th)*(M*(dot4+dot3*M)-dot3*th)+
				 M*(M*(-7.*(dot1+dot2)*dot4-(11.*dot1*dot3+dot2*dot3+8.*dot6)*M+dot5*M2)+8.*dot6*sh)*uh+
				 (3.*dot1*dot4+5.*dot2*dot4+5.*dot1*dot3*M+dot2*dot3*M-2.*dot5*M2)*th*uh+dot5*sqr(th)*uh+
				 (-4.*dot2*dot3*M+4.*dot6*M+2.*dot5*M2-dot5*th)*sqr(uh))+
			  2.*pow(a1,5)*pow<3,1>(M)*(M*(4.*dot6*M+dot5*(th-uh))+dot2*(13.*dot4*M+17.*dot3*M2-8.*dot3*(th+uh))
						    +dot1*(11.*dot4*M+29.*dot3*M2-8.*dot3*(th+uh)))+
			  pow(a1,4)*M2*(M*(-4.*dot6*M2+8.*dot6*(sh+2.*uh)+dot5*M*(M2-5.*th+5.*uh))+
					2.*dot2*(-7.*dot4*M2+4.*dot4*(th+uh)+2.*dot3*M*(4.*M2+th+17.*uh))+
					dot1*(-4.*dot4*M2+8.*dot4*(th+uh)+2.*dot3*M*(17.*M2+6.*th+38.*uh)))))/(a1*a2*sqr(-(a22*M2)+th)*sqr(-(a12*M2)+uh));
	  diag[1]=(4.*M*(a22*pow<3,1>(M)*(-((-5.+2.*a1*(-7.+4.*a1))*(2.*dot2*dot3+dot6))-(1.+a1*(-9.+4.*a1))*a2*dot5*M)+
			 M*((1.+2.*a1)*(-5.+4.*a1)*(2.*dot2*dot3+dot6)+2.*(1.+a1*(-7.+4.*a1))*a2*dot5*M)*th+
			 (-1.+4.*a1)*dot5*sqr(th)))/(a1*pow<3,1>(a22*M2-th));
	  diag[2]=(-2.*M*(2.*pow(a1,5)*pow<3,1>(M)*(8.*(dot1+dot2)*dot3+dot5*M)*(M2-sh)+
			  pow(a1,4)*M2*(-(M2*(16.*(dot1+dot2)*dot4+76.*(dot1+dot2)*dot3*M-4.*dot6*M+13.*dot5*M2))+
					(8.*(dot1+dot2)*dot4+76.*(dot1+dot2)*dot3*M-4.*dot6*M+13.*dot5*M2)*sh)+
			  (M2-sh)*((-2.*dot2*dot4+2.*dot1*dot3*M+5.*(2.*dot2*dot3+dot6)*M)*(M2-th)-
				   (M*(8.*dot2*dot3+4.*dot6-dot5*M)+dot5*th)*uh)+
			  2.*pow(a1,3)*M*(dot2*(dot4*(29.*pow<3,1>(M)-17.*M*sh)+dot3*(M2-sh)*(57.*M2-8.*th))+
					  dot1*(3.*dot4*(9.*pow<3,1>(M)-5.*M*sh)+dot3*(M2-sh)*(65.*M2-8.*th))-
					  M*(M2-sh)*(6.*dot6*M+dot5*(-14.*M2+th+uh)))+
			  a12*(dot1*(-2.*dot3*M*(M2-sh)*(47.*M2-14.*th)+4.*dot4*(-15.*pow<4,1>(M)-2.*sh*th+M2*(9.*sh+4.*th)))+
			       dot2*(-2.*dot3*M*(M2-sh)*(31.*M2-14.*th)+dot4*(-62.*pow<4,1>(M)-8.*sh*th+2.*M2*(19.*sh+8.*th)))+
			       M*(M2-sh)*(dot6*(13.*M2-4.*th)+dot5*M*(-25.*M2+9.*th+5.*uh)))+
			  2.*a1*(dot1*(dot3*M*(M2-sh)*(11.*M2-5.*th)+dot4*(11.*pow<4,1>(M)+3.*sh*th-7.*M2*(sh+th)))+
				 dot2*(dot4*(11.*pow<4,1>(M)+5.*sh*th-M2*(7.*sh+9.*th))-dot3*M*(M2-sh)*(M2+th-4.*uh))+
				 (M2-sh)*(dot6*M*(-5.*M2+4.*th+2.*uh)+dot5*(4.*pow<4,1>(M)+th*uh-2.*M2*(2.*th+uh))))))/(a1*a2*sqr(M2-sh)*sqr(-(a22*M2)+th));
	  diag[3]=(-4.*M*(a12*pow<3,1>(M)*(-((-11.-2.*a1+8.*a12)*(2.*dot1*dot3-dot6))+a1*(-4.+a1+4.*a12)*dot5*M)+
			  M*((-3.+2.*a1)*(1.+4.*a1)*(2.*dot1*dot3-dot6)+2.*a1*(2.+a1-4.*a12)*dot5*M)*uh+(-3.+4.*a1)*dot5*sqr(uh)))/(a2*pow<3,1>(a12*M2-uh));
	  diag[4]=(-2.*M*(2.*pow(a1,5)*pow<3,1>(M)*(8.*(dot1+dot2)*dot3-dot5*M)*(M2-sh)+
			  pow(a1,4)*M2*(M2*(-16.*(dot1+dot2)*dot4+4.*(-((dot1+dot2)*dot3)+dot6)*M-3.*dot5*M2)+
					(8.*(dot1+dot2)*dot4+4.*(dot1+dot2)*dot3*M-4.*dot6*M+3.*dot5*M2)*sh)+
			  4.*a1*M*(M2-sh)*(-2.*dot1*(M*(dot4-dot3*M)+dot3*sh)+dot6*(M2-th))+(M2-sh)*(2.*dot2*dot4-M*(dot6+dot5*M)+dot5*th)*uh-
			  2.*a1*(dot2*(dot3*pow<3,1>(M)+9.*dot4*M2-5.*dot4*sh-dot3*M*sh)+dot1*(dot3*pow<3,1>(M)+7.*dot4*M2-3.*dot4*sh-dot3*M*sh)+
				 dot5*(M2-sh)*(2.*M2+th))*uh+
			  a12*(M2*(M2-sh)*(16.*dot1*dot4+6.*dot2*dot4+M*(dot6+dot5*M)-dot5*th)+
			       (M2*(16.*(dot1+dot2)*dot4+20.*(dot1+dot2)*dot3*M-4.*dot6*M+3.*dot5*M2)-
				(8.*(dot1+dot2)*dot4+20.*(dot1+dot2)*dot3*M-4.*dot6*M+3.*dot5*M2)*sh)*uh)+
			  2.*pow(a1,3)*M*(M*(M2-sh)*(-2.*dot6*M+dot5*(2.*M2+th+uh))+
					  dot1*(dot4*M*(3.*M2+sh)-dot3*(M2-sh)*(7.*M2+8.*uh))+
					  dot2*(dot4*(5.*pow<3,1>(M)-M*sh)-dot3*(M2-sh)*(7.*M2+8.*uh)))))/(a1*a2*sqr(M2-sh)*sqr(-(a12*M2)+uh));
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
