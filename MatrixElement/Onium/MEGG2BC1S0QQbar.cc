// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGG2BC1S0QQbar class.
//

#include "MEGG2BC1S0QQbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

IBPtr MEGG2BC1S0QQbar::clone() const {
  return new_ptr(*this);
}

IBPtr MEGG2BC1S0QQbar::fullclone() const {
  return new_ptr(*this);
}

void MEGG2BC1S0QQbar::doinit() {
  GGtoBCQQbarBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<0>(bcbar,principleQuantumNumber(),0,0);
}

void MEGG2BC1S0QQbar::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2);
}

void MEGG2BC1S0QQbar::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGG2BC1S0QQbar,GGtoBCQQbarBase>
describeHerwigMEGG2BC1S0QQbar("Herwig::MEGG2BC1S0QQbar", "HwOniumParameters.so HwMEHadronOnium.so");

void MEGG2BC1S0QQbar::Init() {

  static ClassDocumentation<MEGG2BC1S0QQbar> documentation
    ("There is no documentation for the MEGG2BC1S0QQbar class");

}

double MEGG2BC1S0QQbar::me2() const {
  using namespace ThePEG::Helicity;
  // for(unsigned int ix=0;ix<mePartonData().size();++ix) {
  //   cerr << mePartonData()[ix]->PDGName() << " " << meMomenta()[ix]/GeV << " " << meMomenta()[ix].mass()/GeV << "\n";
  // }
  // for(unsigned int ix=0;ix<mePartonData().size();++ix) {
  //   cerr << mePartonData()[ix]->PDGName() << " " << rescaledMomenta()[ix]/GeV << " " << rescaledMomenta()[ix].mass()/GeV << "\n";
  // }
  // gluon wavefunction2
  VectorWaveFunction      g1w(rescaledMomenta()[0],mePartonData()[0],incoming);
  VectorWaveFunction      g2w(rescaledMomenta()[1],mePartonData()[1],incoming);
  vector<VectorWaveFunction> g1,g2;
  for(unsigned int ix=0;ix<2;++ix) {
    //g1w.reset(10);
    g1w.reset(2*ix);
    g1.push_back(g1w);
    //g2w.reset(10);
    g2w.reset(2*ix);
    g2.push_back(g2w);
  }
  SpinorWaveFunction      q4w(rescaledMomenta()[3],mePartonData()[3],incoming);
  SpinorBarWaveFunction   q5w(rescaledMomenta()[4],mePartonData()[4],outgoing);
  vector<SpinorWaveFunction> v4;
  vector<SpinorBarWaveFunction> ubar5;
  for(unsigned int ix=0;ix<2;++ix) {
    q4w.reset(ix);
    v4.push_back(q4w);
    q5w.reset(ix);
    ubar5.push_back(q5w);
  }
  // colour matrix
  vector<vector<double > > cMatrix = {{24., 12., 12.},
				      {12., 48., -6.},
				      {12., -6., 48.}};
  // invariants
  Energy M  = rescaledMomenta()[2].mass();
  Energy2 M2 = sqr(M);
  double a1 = rescaledMomenta()[3].mass()/M, a2 = rescaledMomenta()[4].mass()/M;
  double a12(sqr(a1)),a22(sqr(a2));
  // matrix element
  //ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half);
  Complex me2Sum=0.;
  vector<double> flows(3,0.);
  for(unsigned int ih1=0;ih1<2;++ih1) {
    for(unsigned int ih2=0;ih2<2;++ih2) {
      for(unsigned int ih4=0;ih4<2;++ih4) {
	for(unsigned int ih5=0;ih5<2;++ih5) {
	  auto dot1 = rescaledMomenta()[0]*rescaledMomenta()[1];
	  auto dot2 = rescaledMomenta()[0]*rescaledMomenta()[3];
	  auto dot3 = rescaledMomenta()[0]*rescaledMomenta()[4];
	  auto dot4 = rescaledMomenta()[0]*g2[ih2].wave();
	  auto dot5 = rescaledMomenta()[1]*rescaledMomenta()[3];
	  auto dot6 = rescaledMomenta()[1]*rescaledMomenta()[4];
	  auto dot7 = rescaledMomenta()[1]*g1[ih1].wave();
	  auto dot8 = rescaledMomenta()[3]*g1[ih1].wave();
	  auto dot9 = rescaledMomenta()[3]*g2[ih2].wave();
	  auto dot10 = rescaledMomenta()[4]*g1[ih1].wave();
	  auto dot11 = rescaledMomenta()[4]*g2[ih2].wave();
	  auto dot12 = g1[ih1].wave()*g2[ih2].wave();
	  complex<Energy> dot13=v4[ih4].dimensionedWave().pseudoScalar(ubar5[ih5].dimensionedWave());
	  auto vec1 = v4[ih4].dimensionedWave().generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  auto vec2 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  auto vec3 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  auto vec4 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[1]).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  auto vec5 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[0]).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  auto vec6 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(rescaledMomenta()[0]).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  auto vec7 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(rescaledMomenta()[1]).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  auto vec8 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[0]).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  auto vec9 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[1]).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  auto vec10 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(g1[ih1].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  auto vec11 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[0]).slash(g2[ih2].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  auto vec12 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[1]).slash(g2[ih2].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  auto vec13 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(g1[ih1].wave()).slash(rescaledMomenta()[1]).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  auto dot14 = vec1*rescaledMomenta()[0];
	  auto dot15 = vec1*rescaledMomenta()[1];
	  auto dot16 = vec1*g1[ih1].wave();
	  auto dot17 = vec1*g2[ih2].wave();
	  auto dot18 = vec6*g1[ih1].wave();
	  auto dot19 = vec7*g1[ih1].wave();
	  auto dot20 = vec8*g2[ih2].wave();
	  auto dot21 = vec9*g2[ih2].wave();
	  auto dot22 = vec2*rescaledMomenta()[0];
	  auto dot23 = vec2*rescaledMomenta()[1];
	  auto dot24 = vec3*rescaledMomenta()[0];
	  auto dot25 = vec3*rescaledMomenta()[1];
	  auto dot26 = vec4*rescaledMomenta()[0];
	  auto dot27 = vec4*g1[ih1].wave();
	  auto dot28 = vec4*g2[ih2].wave();
	  auto dot29 = vec5*g1[ih1].wave();
	  auto dot30 = vec5*g2[ih2].wave();
	  auto dot31 = vec10*rescaledMomenta()[0];
	  auto dot32 = vec10*rescaledMomenta()[1];
	  auto dot33 = vec11*g1[ih1].wave();
	  auto dot34 = vec12*g1[ih1].wave();
	  auto dot35 = vec13*rescaledMomenta()[0];
	  auto dot36 = vec3*g1[ih1].wave();
	  auto dot37 = vec7*rescaledMomenta()[0];
	  auto dot38 = vec9*rescaledMomenta()[0];
	  Complex diag[38];
	  diag[0]=((dot18+dot19+dot20+dot21-2.*dot17*(dot10+dot8)-2.*dot16*(dot11+dot9)+4.*dot12*(dot14+dot15+dot13*M))*M2)/
	    (16.*a1*a2*(dot1-dot2-dot5)*(dot1-dot3-dot6));
	  diag[1]=-0.125*((4.*dot11*dot16+dot18+dot19-2.*dot17*(dot10-2.*dot8)-2.*(dot20+dot21+dot16*dot9)-2.*dot12*(dot14+dot15+dot13*M))*M2)/
	    (a1*a2*(dot1-dot2-dot5)*(dot1-dot3-dot6));
	  diag[2]=-0.125*((-2.*dot11*dot16+4.*dot10*dot17-2.*dot18-2.*dot19+dot20+dot21-2.*dot17*dot8+4.*dot16*dot9-2.*dot12*(dot14+dot15+dot13*M))*M2)/
	    (a1*a2*(dot1-dot2-dot5)*(dot1-dot3-dot6));
	  diag[3]=-0.25*((dot12*(dot14*(-(a2*dot2)-dot5+a1*(dot3-dot5-dot6)+2.*dot6)+dot15*((1.+a1)*dot2+(-2.+a1)*dot3+dot5-a1*(dot5+dot6))+
				 dot13*(-(a2*dot2)+dot5+a1*(dot3-dot5-dot6))*M)
			  +2.*(dot16*dot4*(dot2-dot3+dot5-dot6)-a2*dot11*(dot14+dot15)*dot7+dot17*(-dot2+dot3-dot5+dot6)*dot7-
			       a1*(dot14+dot15)*dot4*dot8+a1*(dot14+dot15)*dot7*dot9
			       +dot13*(a1*dot11*dot7+dot4*dot8-a1*dot4*dot8-a2*dot7*dot9)*M
			       +dot10*dot4*(a2*(dot14+dot15)-a1*dot13*M)))*M2)/(a1*a2*dot1*(dot1-dot2-dot5)*(dot1-dot3-dot6));
	  diag[4]=((2.*(-(a2*dot10*(dot14+dot15)*dot4)-dot16*dot4*(dot2+dot5)+a2*dot11*(dot14+dot15)*dot7+dot17*(dot2+dot5)*dot7
			+a1*(dot14+dot15)*dot4*dot8-a1*(dot14+dot15)*dot7*dot9)+
		    (-2.*a2*dot10*dot13*dot4+(2.*a2*dot11*dot13+dot24+dot25)*dot7-dot4*(dot22+dot23+2.*a2*dot13*dot8)+2.*a2*dot13*dot7*dot9)*M+
		    dot12*(-(dot15*((1.+a1)*dot2-a2*(dot3-dot5-dot6)))+dot14*(dot2-a1*dot2+dot3+dot5-dot6+a1*(-dot3+dot5+dot6))-
			   (dot26-a2*dot13*(dot2+dot3-dot5-dot6))*M)+
		    dot1*(2.*dot16*dot4-2.*dot17*dot7+dot12*(-dot14+dot15+dot13*M)))*M2)/
	    (4.*a1*dot1*(dot1-dot3-dot6)*((-1.+2.*a1)*dot1+a2*(dot2+dot3+dot5+dot6)));
	  diag[5]=(((dot1-dot3-dot6)*(dot12*(dot14-dot15)-2.*dot16*dot4+2.*dot17*dot7)+
		    (-(dot1*dot12*dot13)+2.*dot10*dot13*dot4-dot27*dot4-dot29*dot4+dot12*(dot26-dot13*dot3+dot13*dot6)-
		     2.*dot11*dot13*dot7+dot28*dot7+dot30*dot7)*M)*M2)/(4.*a1*dot1*sqr(-dot1+dot3+dot6));
	  diag[6]=-0.25*((2.*(-(dot10*(dot14+dot15)*dot4)+dot12*(dot15*dot3-dot14*dot6)+dot11*(dot14+dot15)*dot7+
			      (dot3+dot6)*(dot16*dot4-dot17*dot7))-(dot12*dot26-(dot27+dot29)*dot4+(dot28+dot30)*dot7)*M-
			  a1*(dot12*(dot2+dot3-dot5-dot6)-2.*dot4*(dot10+dot8)+2.*dot7*(dot11+dot9))*(dot14+dot15+dot13*M)+
			  dot1*(-2.*dot16*dot4+2.*dot17*dot7+dot12*(dot14-dot15+dot13*M)))*M2)/
	    (a2*dot1*(dot1-dot2-dot5)*((-1.+2.*a1)*dot1-a1*(dot2+dot3+dot5+dot6)));
	  diag[7]=-0.25*(((dot1-dot2-dot5)*(dot12*(dot14-dot15)-2.*dot16*dot4+2.*dot17*dot7)+
			  (dot1*dot12*dot13-dot12*(dot13*dot2+dot26-dot13*dot5)+(dot24+dot25)*dot7-dot4*(dot22+dot23-2.*dot13*dot8)-2.*dot13*dot7*dot9)*M)*M2)/
	    (a2*dot1*sqr(-dot1+dot2+dot5));
	  diag[8]=-0.0625*M2*(-2.*dot1*dot18-4.*dot11*dot16*dot3+4.*a1*dot11*dot16*dot3+2.*dot18*dot3+2.*dot11*dot38-2.*a1*dot11*dot38+
			      2.*dot18*dot5+2.*dot18*dot6-4.*dot11*dot14*dot7+8.*a1*dot11*dot14*dot7+4.*a1*dot11*dot15*dot7-2.*a1*dot37*dot7-
			      4.*a12*(dot14+dot15)*(dot7-dot8)*(dot11-dot4+dot9)+
			      2.*(-(dot3*dot32)+dot3*dot33+dot3*dot34-dot38*dot4-dot1*(dot19+dot33+dot34-2.*dot16*dot4)+
				  dot32*dot5+dot33*dot6+dot34*dot6-2.*dot16*dot4*dot6+dot19*(dot3+dot5+dot6)+
				  2.*dot17*dot3*dot7+2.*dot14*dot4*dot7-2.*dot17*dot5*dot7-2.*dot16*dot3*dot9)
			      -2.*a1*(dot19*dot3+dot3*dot33+dot3*dot34-dot38*dot4-dot1*(dot18+dot19+dot33+dot34-2.*dot16*dot4)+
				      dot19*dot6+dot33*dot6+dot34*dot6-2.*dot16*dot4*dot6+dot18*(dot3+dot6)+4.*dot14*dot4*dot7+
				      2.*dot15*dot4*dot7-2.*dot17*dot5*dot7+2.*dot11*dot14*dot8+2.*dot11*dot15*dot8-dot37*dot8-
				      2.*dot14*dot4*dot8-2.*dot15*dot4*dot8+2.*dot17*dot5*dot8+
				      (-2.*dot16*dot3+dot38-4.*dot14*dot7+2.*dot14*dot8)*dot9)+
			      (-2.*dot11*dot22-4.*dot11*dot23-2.*dot12*dot26-4.*dot11*dot27-4.*dot11*dot29+3.*dot35-
			       4.*a12*(dot22+dot23+dot27+dot29-dot13*(dot7+dot8))*(dot11-dot4+dot9)+2.*(-(dot1*dot36)+4.*dot11*dot13*dot7-2.*dot24*dot7+(dot24+dot25+dot28+dot30)*dot8+dot4*(dot22+dot23+2.*(dot27+dot29-dot13*(dot7+dot8)))-(dot22+2.*(dot23+dot27+dot29-2.*dot13*dot7))*dot9)-2.*a1*(-(dot25*dot7)+(dot24+2.*dot25+dot28+dot30)*dot8+dot4*(3.*dot22+4.*dot23+4.*dot27+4.*dot29-6.*dot13*dot7-4.*dot13*dot8)-dot11*(3.*dot22+4.*dot23+4.*dot27+4.*dot29-2.*dot13*(3.*dot7+dot8))-(3.*dot22+4.*dot23+4.*dot27+4.*dot29-2.*dot13*(3.*dot7+dot8))*dot9))*M+2.*dot10*(-dot37+2.*dot14*(-dot4+dot9)+(dot24+dot25+dot28+dot30-2.*dot13*dot4)*M-2.*a12*(dot4-dot9)*(dot14+dot15+dot13*M)+a1*(dot37+2.*dot15*dot4-2.*dot17*dot5+4.*dot14*(dot4-dot9)-(dot24+2.*dot25+dot28+dot30-4.*dot13*dot4+2.*dot13*dot9)*M)-2.*a2*dot11*(-(a2*dot14)+a1*(dot15+dot13*M)))+a2*(-2.*(-1.+2.*a1)*dot12*(dot14+dot15)-2.*dot18-2.*dot19+dot20+dot21-3.*dot33-3.*dot34+4.*dot16*dot4-2.*dot17*dot7+a1*(4.*dot18+4.*dot19-2.*dot20-2.*dot21+dot31+dot32+5.*dot33+5.*dot34-6.*dot16*dot4+2.*dot17*dot7))*M2)/(a1*a2*(dot1-dot2-dot3)*(dot1-dot5-dot6)*(-(a2*dot2)-a1*dot6));
	  diag[9]=-0.125*((2.*((dot18+dot19)*dot5+dot32*dot5-a1*dot37*dot7-2.*dot17*dot5*dot7+2.*a1*dot17*dot5*dot7+a1*dot37*dot8-2.*a1*dot17*dot5*dot8-2.*a1*dot14*(dot10-dot7+dot8)*dot9+dot10*(-(a2*dot37)-2.*a1*dot17*dot5+2.*dot14*dot9)-dot3*(dot32-2.*dot17*dot7+2.*dot16*dot9))+(dot35-2.*dot24*dot7-2.*(dot23+dot27+dot29-2.*dot13*dot7)*dot9+2.*a1*(dot10-dot7+dot8)*(dot25-2.*dot13*dot9))*M)*M2)/(a1*a2*(dot1-dot2-dot3)*(dot1-dot2-dot5)*dot5);
	  diag[10]=-0.125*((4.*dot12*dot15*(-dot1+dot3)+2.*(-(dot3*dot32)-2.*dot12*dot14*dot6-dot32*dot6+dot33*dot6+dot34*dot6+a2*(dot37+2.*dot15*dot4)*dot7+2.*dot17*dot6*dot7-2.*a1*dot17*dot6*dot7+a1*(dot37+2.*dot15*dot4+2.*dot17*dot6)*dot8+dot1*(-2.*dot11*dot16+dot32-2.*a2*dot17*dot7-2.*a1*dot17*dot8)+2.*dot11*(dot16*dot3+dot14*(dot7-a1*dot7+a1*dot8)))+(2.*dot11*(dot23+dot27+dot29)+dot35-2.*dot23*dot4-2.*(dot12*dot26+dot1*(-2.*dot12*dot13+dot36)+2.*(1.+a1)*dot11*dot13*dot7+a1*dot25*dot7-a1*(2.*dot11*dot13+dot25)*dot8))*M+2.*dot10*(2.*a2*dot1*dot17-dot37-2.*dot15*dot4+2.*dot11*(-(a2*dot14)+a1*dot13*M)+a1*(dot37+2.*dot15*dot4+2.*dot17*dot6+dot25*M)))*M2)/(a1*(dot1-dot2-dot3)*dot6*(-(a2*dot2)-a1*dot6));
	  diag[11]=(M2*(2.*(dot21*dot3-dot3*dot31-2.*dot3*dot32+dot11*dot38-dot38*dot4-dot10*(dot37+2.*(dot14+dot15)*dot4)+dot21*dot6-dot31*dot6-2.*dot32*dot6+dot33*dot6+dot34*dot6-2.*dot16*dot4*dot6+dot20*(dot3+dot6)+2.*dot12*(dot14*dot3+dot15*(2.*dot3+dot6))+dot37*dot7+2.*((dot14+dot15)*dot4+dot17*dot6)*dot7)+(-4.*dot12*dot26+3.*dot35+2.*((dot25+dot28+dot30)*dot7+dot11*(dot22+dot23+dot27+dot29-2.*dot13*dot7)+dot22*dot9+dot23*(-dot4+dot9)+(dot27+dot29-2.*dot13*dot7)*(dot4+dot9)))*M-2.*dot1*(2.*dot11*dot16-2.*dot10*dot17+dot20+dot21-dot31-2.*dot32-2.*dot16*dot4+2.*dot17*dot7+2.*dot36*M+2.*dot12*(dot14+2.*dot15-dot13*M))+(-2.*dot12*(dot14+dot15)+dot18+dot19-dot20-dot21+dot31+dot32+dot33+dot34)*M2+a12*(4.*(dot14+dot15)*(dot10-dot7+dot8)*(dot11-dot4+dot9)-4.*(dot10-dot7+dot8)*(dot24+dot25+dot28+dot30-dot13*(dot11+dot4+dot9))*M+(-10.*dot12*(dot14+dot15)+4.*dot18+4.*dot19-5.*dot20-5.*dot21+5.*dot31+5.*dot32+4.*(dot33+dot34)+2.*dot16*dot4)*M2)-a1*(-8.*dot1*dot12*dot14-8.*dot1*dot12*dot15-4.*dot1*dot11*dot16+4.*dot12*dot14*dot2+4.*dot12*dot15*dot2+4.*dot11*dot16*dot2-2.*dot1*dot20-2.*dot1*dot21+8.*dot12*dot14*dot3+8.*dot12*dot15*dot3+2.*dot20*dot3+2.*dot21*dot3+2.*dot1*dot31-2.*dot3*dot31+2.*dot1*dot32-2.*dot3*dot32+2.*dot1*dot33-2.*dot2*dot33-2.*dot3*dot33+2.*dot1*dot34-2.*dot2*dot34-2.*dot3*dot34+2.*dot11*dot38+4.*dot1*dot16*dot4-2.*dot38*dot4+4.*dot12*dot14*dot6+4.*dot12*dot15*dot6+2.*dot20*dot6+2.*dot21*dot6-2.*dot31*dot6-2.*dot32*dot6-4.*dot16*dot4*dot6-4.*dot11*dot14*dot7-4.*dot11*dot15*dot7-4.*dot1*dot17*dot7+2.*dot37*dot7+8.*dot14*dot4*dot7+8.*dot15*dot4*dot7+4.*dot17*dot6*dot7+4.*dot11*dot15*dot8+4.*dot1*dot17*dot8-2.*dot37*dot8-4.*dot14*dot4*dot8-8.*dot15*dot4*dot8-4.*dot17*dot6*dot8-4.*dot16*dot3*dot9+2.*dot38*dot9-4.*dot14*dot7*dot9+2.*((-dot1+dot2+dot3)*dot36+(3.*dot24+2.*dot25+3.*(dot28+dot30))*dot7+dot11*(2.*dot22+dot23+dot27+dot29-2.*dot13*dot7)-(2.*dot24+dot25+2.*(dot28+dot30))*dot8+dot4*(dot23+dot27+dot29-6.*dot13*dot7+2.*dot13*dot8)+2.*dot22*dot9+(dot23+dot27+dot29-2.*dot13*dot7)*dot9)*M+2.*dot10*(2.*dot11*(dot14+dot15)+4.*dot1*dot17-dot37-4.*(dot14+dot15)*dot4-2.*dot17*(dot2+dot3+dot6)+2.*dot14*dot9-(2.*dot24+dot25+2.*(dot28+dot30-dot13*dot4))*M)+(-10.*dot12*(dot14+dot15)+4.*dot18+4.*dot19-5.*dot20-5.*dot21+5.*dot31+5.*dot32+4.*(dot33+dot34)+2.*dot16*dot4)*M2)))/(16.*a1*a2*(dot1-dot2-dot3)*(dot1-dot2-dot5)*((-1.+2.*a1)*dot1-a1*(dot2+dot3+dot5+dot6)));
	  diag[12]=-0.03125*(M2*(8.*dot2*(-(dot12*dot15)+dot32)+4.*dot11*dot38+4.*(2.*dot10*dot14-2.*dot16*dot3+dot38)*dot4-2.*(2.*dot12*(-dot14+dot15)+dot20+dot21-2.*dot32+dot33+dot34)*dot5+2.*(2.*dot12*(dot14+dot15)+dot20+dot21+2.*dot31-dot33-dot34)*dot6-4.*(2.*dot11*dot14+dot37+2.*(dot14+dot15)*dot4+dot17*(2.*dot2+dot5+dot6))*dot7+8.*(dot37+dot15*dot4)*dot8-8.*(dot10*dot14-dot16*dot3+dot38-dot14*dot7)*dot9+2.*dot1*(2.*dot12*(-dot14+dot15)-dot20-dot21-2.*dot32+dot33+dot34+6.*dot17*dot7-4.*dot17*dot8+4.*dot16*dot9)-4.*(dot23+dot27+dot29-2.*dot13*dot7)*(dot4-dot9)*M-(dot18+dot19+dot31+dot32-2.*dot17*dot7)*M2+a12*(-8.*(dot14+dot15)*(dot10-dot7+dot8)*(dot11-dot4+dot9)-8.*dot13*(dot10-dot7+dot8)*(dot11-dot4+dot9)*M-(dot18+dot19+dot31+dot32-2.*dot17*dot7)*M2)+2.*a1*(2.*dot11*dot16*dot2-dot2*dot20-dot2*dot21-2.*dot11*dot16*dot3-dot20*dot3-dot21*dot3+2.*dot11*dot38+4.*dot16*dot3*dot4-2.*dot38*dot4-dot20*dot5-dot21*dot5+dot33*dot5+dot34*dot5-dot20*dot6-dot21*dot6+dot33*dot6+dot34*dot6-2.*dot12*(dot14+dot15)*(dot5+dot6)-4.*dot11*dot14*dot7-4.*dot11*dot15*dot7+4.*dot37*dot7+8.*dot14*dot4*dot7+8.*dot15*dot4*dot7-2.*dot17*dot5*dot7+6.*dot17*dot6*dot7+4.*dot11*dot15*dot8+2.*dot17*dot2*dot8+2.*dot17*dot3*dot8-4.*dot37*dot8-4.*dot14*dot4*dot8-8.*dot15*dot4*dot8+4.*dot17*dot5*dot8-4.*dot17*dot6*dot8+dot1*(2.*dot12*(dot14+dot15)-2.*dot11*dot16+2.*dot20+2.*dot21-dot33-dot34-6.*dot17*dot7+2.*dot17*dot8)-4.*dot16*dot3*dot9+2.*dot38*dot9-8.*dot14*dot7*dot9+4.*dot14*dot8*dot9-2.*dot11*(dot23+dot27+dot29-2.*dot13*dot7)*M+2.*(dot23+dot27+dot29-2.*dot13*dot8)*(dot4-dot9)*M+4.*dot10*(dot11*(dot14+dot15)-dot37-2.*dot14*dot4-2.*dot15*dot4+dot17*(dot1+dot5-dot6)+2.*dot14*dot9-dot13*dot4*M+dot13*dot9*M)+(dot18+dot19+dot31+dot32-2.*dot17*dot7)*M2)))/(a1*a2*(dot1-dot2-dot3)*(dot1-dot2-dot5)*(-(a2*dot2)-a1*dot6));
	  diag[13]=-0.0625*(M2*(2.*(dot11*dot38+dot32*(dot2-dot6)+(dot20+dot21)*dot6+2.*dot12*dot15*(-dot2+dot6)+dot10*(2.*dot11*dot15+2.*dot1*dot17-dot37-2.*dot15*dot4-2.*dot17*dot6)-2.*dot11*dot15*dot7+dot37*dot7+2.*dot15*dot4*dot7+2.*dot17*dot6*dot7-2.*dot1*(dot11*dot16+dot17*dot7)+2.*dot11*dot15*dot8-2.*dot17*dot6*dot8)+(2.*dot10*dot25+3.*dot35-4.*dot23*dot4-2.*(dot24+dot25)*dot7+4.*(dot1*dot12*dot13-dot12*dot26-dot1*dot36+dot13*dot4*dot7)+2.*dot25*dot8)*M+a1*(4.*dot12*dot14*dot2+4.*dot12*dot15*dot2-4.*dot11*dot16*dot2+2.*dot2*dot20+2.*dot2*dot21-2.*dot2*dot31-2.*dot2*dot32-2.*dot11*dot38+2.*dot38*dot4+4.*dot12*dot14*dot5+4.*dot12*dot15*dot5+2.*dot20*dot5+2.*dot21*dot5-2.*dot31*dot5-2.*dot32*dot5-4.*dot16*dot4*dot5+4.*dot11*dot14*dot7+8.*dot11*dot15*dot7-2.*dot37*dot7-4.*dot14*dot4*dot7-8.*dot15*dot4*dot7-4.*dot17*dot6*dot7-8.*dot11*dot15*dot8+2.*dot37*dot8+8.*dot15*dot4*dot8+4.*dot17*dot6*dot8-4.*dot16*dot2*dot9-2.*dot38*dot9+4.*dot14*dot7*dot9+4.*dot15*dot7*dot9-4.*dot15*dot8*dot9+2.*dot1*(-2.*dot12*(dot14+dot15)-dot20-dot21+dot31+dot32+2.*dot17*dot7-2.*dot17*dot8+2.*dot16*(dot11+dot9))-2.*((dot24+dot28+dot30)*dot7-(2.*dot24+dot25+2.*(dot28+dot30))*dot8+dot11*(dot22+2.*dot13*(-dot7+dot8))-2.*dot13*dot7*dot9+dot22*(-dot4+dot9)+2.*dot13*dot8*(dot4+dot9))*M+2.*dot10*(-2.*dot1*dot17+dot37-2.*(dot14+2.*dot15)*(dot11-dot4)+2.*dot17*dot6-2.*(dot14+dot15)*dot9+(2.*dot24+dot25+2.*(dot28+dot30)-2.*dot13*(dot11+dot4+dot9))*M)+(2.*dot12*(dot14+dot15)+dot20+dot21-dot31-dot32-2.*dot16*dot4)*M2)+a12*(4.*(dot14+dot15)*(dot10-dot7+dot8)*(dot11-dot4+dot9)-4.*(dot10-dot7+dot8)*(dot24+dot25+dot28+dot30-dot13*(dot11+dot4+dot9))*M+(-4.*dot12*(dot14+dot15)-2.*dot20-2.*dot21+dot31+dot32+dot33+dot34+2.*dot16*dot4+2.*dot17*dot7)*M2)))/(a1*a2*(dot1-dot2-dot3)*(a1*(dot3-dot5)+dot5)*(dot1-dot5-dot6));
	  diag[14]=-0.125*(M2*(4.*dot1*dot11*dot16-2.*(dot2*dot32+dot11*dot38+2.*dot12*dot15*(-dot2+dot6)+2.*dot11*dot15*dot8+dot6*(dot20+dot21-dot32-2.*dot17*dot8))-(4.*dot1*dot12*dot13-2.*dot12*dot26+dot35-2.*dot1*dot36-2.*dot23*dot4+2.*dot25*dot8)*M+a12*(2.*dot12*(dot14+dot15)+dot20+dot21-dot31-dot32-2.*dot16*dot4)*M2+a1*(2.*dot1*(2.*dot12*(dot14+dot15)+dot20+dot21-dot31-dot32-2.*dot16*(dot11+dot9))+2.*((-dot20-dot21+dot31+dot32)*dot5+dot11*(2.*dot16*dot2+dot38+2.*dot15*dot8)-dot4*(dot38-2.*dot16*dot5+2.*dot15*dot8)+(dot38+2.*dot15*dot8)*dot9+dot2*(-dot20-dot21+dot31+dot32+2.*dot16*dot9))-2.*(dot22-2.*dot13*dot8)*(dot11-dot4+dot9)*M+(-dot20-dot21+dot31+dot32+2.*dot16*dot4)*M2-2.*dot12*(dot14+dot15)*(2.*(dot2+dot5)+M2))))/(a1*a2*dot2*(-dot1+dot2+dot5)*(dot1-dot5-dot6));
	  diag[15]=(M2*(-2.*dot11*dot38+4.*dot11*dot14*dot7-dot35*M+2.*(dot38*dot4+(dot32-dot33-dot34)*dot6+2.*dot16*dot4*dot6-2.*dot14*dot4*dot7-2.*dot17*dot6*dot7+dot3*(dot32-2.*dot17*dot7)+dot1*(-dot32+dot33+dot34-2.*dot16*dot4+2.*dot17*dot7)+dot24*dot7*M)+dot10*(4.*dot17*(dot1-dot6)-4.*dot15*(-(a2*(dot11-dot4))+a1*dot9)+2.*(dot25-2.*a1*dot13*(dot11-dot4+dot9))*M)+2.*a12*(dot18+dot19+dot33+dot34-2.*dot16*dot4)*M2+a1*(2.*(dot19*dot3+dot3*dot33+dot3*dot34-dot38*dot4-dot1*(dot18+dot19+dot33+dot34-2.*dot16*dot4)+dot19*dot6+dot33*dot6+dot34*dot6-2.*dot16*dot4*dot6+dot18*(dot3+dot6)+2.*dot14*dot4*dot7+dot11*(-2.*dot16*dot3+dot38-2.*dot14*dot7)+(-2.*dot16*dot3+dot38-2.*dot14*dot7)*dot9)+2.*(dot23+dot27+dot29-2.*dot13*dot7)*(dot11-dot4+dot9)*M-(dot18+dot19+dot33+dot34-2.*dot16*dot4)*M2)))/(8.*a1*dot3*(a1*(dot3-dot5)+dot5)*(-dot1+dot5+dot6));
	  diag[16]=(M2*(2.*(dot3*dot32+dot10*dot37-dot11*dot38+dot38*dot4+dot32*dot6-dot33*dot6-dot34*dot6+2.*dot16*dot4*dot6-(2.*dot17*dot3+dot37-2.*(dot14+dot15)*(dot11-dot4)+4.*dot17*dot6)*dot7)+(2.*dot12*dot26-3.*dot35+2.*dot23*dot4+2.*(2.*dot24+dot25-2.*dot13*dot4)*dot7)*M+2.*dot1*(-dot32+dot33+dot34-2.*dot16*dot4+4.*dot17*dot7+dot36*M)+a12*(-4.*(dot14+dot15)*(dot10-dot7+dot8)*(dot11-dot4+dot9)+4.*(dot22+dot23+dot27+dot29-dot13*(dot10+dot7+dot8))*(dot11-dot4+dot9)*M+(-4.*dot12*(dot14+dot15)+3.*dot18+3.*dot19-2.*dot20-2.*dot21+dot31+dot32+4.*(dot33+dot34-dot16*dot4)+2.*dot17*dot7)*M2)+a1*(4.*dot1*dot11*dot16-2.*dot1*dot18-2.*dot1*dot19-4.*dot11*dot16*dot3+2.*dot18*dot3+2.*dot19*dot3-4.*dot1*dot33+2.*dot3*dot33-4.*dot1*dot34+2.*dot3*dot34+2.*dot11*dot38+4.*dot1*dot16*dot4-2.*dot38*dot4-4.*dot11*dot16*dot5+2.*dot33*dot5+2.*dot34*dot5-4.*dot11*dot16*dot6+2.*dot18*dot6+2.*dot19*dot6+4.*dot33*dot6+4.*dot34*dot6-4.*dot16*dot4*dot6-8.*dot11*dot14*dot7-8.*dot11*dot15*dot7-4.*dot1*dot17*dot7+2.*dot37*dot7+8.*dot14*dot4*dot7+8.*dot15*dot4*dot7+4.*dot17*dot6*dot7+4.*dot11*dot15*dot8+4.*dot1*dot17*dot8-2.*dot37*dot8-4.*dot15*dot4*dot8-4.*dot17*dot6*dot8-4.*dot16*dot3*dot9+2.*dot38*dot9-8.*dot14*dot7*dot9-4.*dot15*dot7*dot9-2.*(dot1*(2.*dot12*dot13-dot36)-2.*dot12*dot13*(dot5+dot6)+dot36*(dot5+dot6)-2.*dot11*dot13*dot7+dot25*dot7+2.*dot13*dot4*dot7-dot24*dot8-2.*dot25*dot8-dot28*dot8-dot30*dot8+2.*dot13*dot4*dot8-2.*dot13*dot7*dot9+(dot22+2.*(dot23+dot27+dot29))*(dot11-dot4+dot9))*M+2.*dot10*(2.*dot11*(dot14+dot15)-dot37-2.*(dot14+dot15)*dot4+2.*dot17*dot5+2.*dot14*dot9+(dot24+2.*dot25+dot28+dot30-2.*dot13*dot4)*M)+(2.*dot12*(dot14+dot15)-dot18-dot19+dot20+dot21-2.*(dot33+dot34-dot16*dot4+dot17*dot7))*M2)))/(16.*a1*a2*(dot1-dot2-dot5)*(dot1-dot5-dot6)*((-1.+2.*a1)*dot1-a1*(dot2+dot3+dot5+dot6)));
	  diag[17]=-0.03125*(M2*(2.*(-2.*dot2*dot32+2.*dot3*dot32+2.*dot10*dot37+4.*dot10*dot15*dot4-2.*dot38*dot4+3.*dot31*dot5-dot32*dot5+dot33*dot5+dot34*dot5+2.*dot16*dot4*dot5-4.*dot12*(dot15*dot3+dot14*dot5)+(dot18+dot19)*dot6+dot33*dot6+dot34*dot6-2.*dot16*dot4*dot6-4.*dot11*dot15*dot7+4.*dot17*dot2*dot7+2.*dot37*dot7+4.*dot14*dot4*dot7+4.*dot15*dot4*dot7+2.*dot17*dot5*dot7+4.*dot17*dot6*dot7+4.*dot11*dot15*dot8-4.*dot37*dot8-4.*dot15*dot4*dot8-4.*dot17*dot6*dot8+4.*dot38*dot9-4.*dot14*dot7*dot9-dot1*(4.*dot10*dot17+dot31-dot32+dot33+dot34-2.*dot16*dot4+6.*dot17*dot7-4.*dot17*dot8+4.*dot16*dot9)-2.*dot25*dot7*M+2.*dot25*dot8*M)+a12*(8.*(dot14+dot15)*(dot10-dot7+dot8)*(dot11-dot4+dot9)+8.*dot13*(dot10-dot7+dot8)*(dot11-dot4+dot9)*M+(2.*dot12*(dot14+dot15)-dot18-dot19+dot20+dot21-dot31-dot32-dot33-dot34)*M2)+a1*(-8.*dot11*dot16*dot2+2.*dot18*dot2+2.*dot19*dot2+8.*dot11*dot16*dot3-2.*dot20*dot3-2.*dot21*dot3+2.*dot3*dot31+2.*dot3*dot32+2.*dot2*dot33+2.*dot2*dot34-8.*dot11*dot38+4.*dot16*dot2*dot4-4.*dot16*dot3*dot4+8.*dot38*dot4+2.*dot18*dot5+2.*dot19*dot5-2.*dot20*dot6-2.*dot21*dot6+2.*dot31*dot6+2.*dot32*dot6-2.*dot33*dot6-2.*dot34*dot6+8.*dot16*dot4*dot6-4.*dot12*(dot14+dot15)*(dot3+dot6)+16.*dot11*dot14*dot7+16.*dot11*dot15*dot7-4.*dot37*dot7-16.*dot14*dot4*dot7-16.*dot15*dot4*dot7-8.*dot17*dot6*dot7-16.*dot11*dot15*dot8+4.*dot37*dot8+16.*dot15*dot4*dot8+8.*dot17*dot6*dot8-8.*dot16*dot2*dot9+8.*dot16*dot3*dot9-8.*dot38*dot9-4.*dot16*dot5*dot9-4.*dot16*dot6*dot9+16.*dot14*dot7*dot9+8.*dot15*dot7*dot9-8.*dot15*dot8*dot9+4.*dot10*(dot37-2.*(dot14+dot15)*(dot11-dot4)-dot17*(dot1+dot5-dot6)-2.*dot14*dot9)+2.*dot1*(2.*dot12*(dot14+dot15)+4.*dot11*dot16-dot18-dot19+dot20+dot21-dot31-dot32-6.*dot16*dot4+4.*dot17*dot7-4.*dot17*dot8+6.*dot16*dot9)+2.*dot10*(dot24-dot25+dot28+dot30-2.*dot13*dot4)*M+2.*(dot25*(2.*dot7-3.*dot8)+4.*dot11*dot13*(dot7-dot8)-(dot24+dot28+dot30)*dot8+2.*dot13*(-2.*dot4*dot7+3.*dot4*dot8+2.*dot7*dot9-2.*dot8*dot9))*M+(dot18+dot19+dot33+dot34-2.*dot16*dot4)*M2)))/(a1*a2*(dot1-dot2-dot5)*(a1*(dot3-dot5)+dot5)*(dot1-dot5-dot6));
	  diag[18]=-0.125*((2.*dot31*dot5+2.*a2*dot10*(dot37-2.*dot17*dot5)+2.*(a1*dot37+2.*a2*dot17*dot5)*(dot7-dot8)-4.*dot10*dot14*dot9+4.*a1*dot14*(dot10-dot7+dot8)*dot9+2.*dot1*(dot32-2.*dot17*dot7+2.*dot16*dot9)-2.*dot2*(dot32-2.*dot17*dot7+2.*dot16*dot9)-(dot35-2.*dot24*dot7-2.*a2*dot25*(dot10-dot7+dot8)+2.*(dot22+2.*a2*dot13*(dot10-dot7+dot8))*dot9)*M)*M2)/(a2*(dot1-dot2-dot3)*dot5*(a1*(dot3-dot5)+dot5));
	  diag[19]=-0.125*((-4.*dot12*dot15*dot2+2.*(-2.*dot11*dot16*dot2+dot2*dot32+(2.*dot12*dot15+dot20+dot21)*dot6-dot32*dot6-a2*(2.*dot1*dot17-dot37-2.*dot15*dot4)*dot7+2.*a2*dot17*dot6*dot7-2.*dot17*dot6*dot8+a1*(-2.*dot1*dot17+dot37+2.*dot15*dot4+2.*dot17*dot6)*dot8+2.*dot11*dot14*(dot7-a1*dot7+a1*dot8))+(-2.*dot11*dot22-2.*dot12*dot26+dot35+dot1*(4.*dot12*dot13-2.*dot36)-2.*dot23*dot4+2.*a2*(2.*dot11*dot13+dot25)*(dot7-dot8))*M-2.*a2*dot10*(-2.*dot1*dot17+dot37+2.*dot15*dot4+2.*dot17*dot6+dot25*M+2.*dot11*(dot14+dot13*M)))*M2)/(a1*a2*(dot1-dot2-dot3)*(dot1-dot3-dot6)*dot6);
	  diag[20]=-0.0625*(M2*(-4.*dot1*dot11*dot16+4.*a1*dot1*dot11*dot16-2.*a1*dot1*dot18-2.*a1*dot1*dot19+4.*dot11*dot16*dot2-4.*a1*dot11*dot16*dot2+2.*a1*dot18*dot2+2.*a1*dot19*dot2+2.*dot1*dot31-4.*a1*dot1*dot31-2.*dot2*dot31+2.*a1*dot2*dot31-2.*dot3*dot31+2.*a1*dot3*dot31-4.*a1*dot1*dot32+2.*a1*dot2*dot32-2.*dot3*dot32+2.*a1*dot3*dot32+2.*a1*dot2*dot33+2.*a1*dot2*dot34+2.*dot11*dot38-2.*a1*dot11*dot38-4.*a1*dot1*dot16*dot4+2.*a1*dot38*dot4+2.*a1*dot18*dot5+2.*a1*dot19*dot5-2.*dot31*dot5+2.*a1*dot31*dot5+2.*a1*dot32*dot5+2.*a1*dot31*dot6+2.*a1*dot32*dot6-2.*a1*dot33*dot6-2.*a1*dot34*dot6+4.*a1*dot16*dot4*dot6-4.*dot11*dot14*dot7+8.*a1*dot11*dot14*dot7-4.*a12*dot11*dot14*dot7-4.*dot11*dot15*dot7+8.*a1*dot11*dot15*dot7-4.*a12*dot11*dot15*dot7+4.*dot1*dot17*dot7+4.*a1*dot1*dot17*dot7-4.*dot17*dot2*dot7-2.*a1*dot37*dot7-4.*a1*dot14*dot4*dot7+4.*a12*dot14*dot4*dot7-4.*a1*dot15*dot4*dot7+4.*a12*dot15*dot4*dot7-4.*dot17*dot5*dot7-4.*a1*dot17*dot6*dot7-4.*a1*dot11*dot14*dot8+4.*a12*dot11*dot14*dot8+4.*dot11*dot15*dot8-8.*a1*dot11*dot15*dot8+4.*a12*dot11*dot15*dot8-4.*dot1*dot17*dot8+4.*a1*dot1*dot17*dot8+4.*dot17*dot2*dot8-4.*a1*dot17*dot2*dot8+4.*dot17*dot3*dot8-4.*a1*dot17*dot3*dot8+2.*a1*dot37*dot8-4.*a12*dot14*dot4*dot8+4.*a1*dot15*dot4*dot8-4.*a12*dot15*dot4*dot8+4.*dot17*dot5*dot8-4.*a1*dot17*dot5*dot8-4.*dot16*dot3*dot9+4.*a1*dot16*dot3*dot9-2.*a1*dot38*dot9+8.*a1*dot14*dot7*dot9-4.*a12*dot14*dot7*dot9+4.*a1*dot15*dot7*dot9-4.*a12*dot15*dot7*dot9-4.*a1*dot14*dot8*dot9+4.*a12*dot14*dot8*dot9-4.*a1*dot15*dot8*dot9+4.*a12*dot15*dot8*dot9+3.*dot35*M-2.*(-2.*a2*dot1*dot12*dot13+2.*dot12*dot13*dot2-2.*a1*dot12*dot13*dot2-3.*a1*dot11*dot22+2.*a12*dot11*dot22+dot11*dot23-4.*a1*dot11*dot23+2.*a12*dot11*dot23+dot12*dot26+dot11*dot27-4.*a1*dot11*dot27+2.*a12*dot11*dot27+dot11*dot29-4.*a1*dot11*dot29+2.*a12*dot11*dot29+2.*dot12*dot13*dot3-2.*a1*dot12*dot13*dot3-(-2.+a1)*dot1*dot36-dot2*dot36+a1*dot2*dot36-dot3*dot36+a1*dot3*dot36+a1*dot22*dot4-2.*a12*dot22*dot4+dot23*dot4+2.*a1*dot23*dot4-2.*a12*dot23*dot4+2.*a1*dot27*dot4-2.*a12*dot27*dot4+2.*a1*dot29*dot4-2.*a12*dot29*dot4+4.*a1*dot11*dot13*dot7-2.*a12*dot11*dot13*dot7+2.*dot24*dot7-a1*dot24*dot7+dot25*dot7-2.*a1*dot25*dot7-a1*dot28*dot7-a1*dot30*dot7-2.*dot13*dot4*dot7+2.*a12*dot13*dot4*dot7-2.*dot11*dot13*dot8+4.*a1*dot11*dot13*dot8-2.*a12*dot11*dot13*dot8-dot25*dot8+a1*dot25*dot8-2.*a1*dot13*dot4*dot8+2.*a12*dot13*dot4*dot8+(a1*(-3.+2.*a1)*dot22+dot23+dot27+dot29+2.*(-2.+a1)*a1*(dot23+dot27+dot29-dot13*dot7)-2.*a22*dot13*dot8)*dot9)*M-2.*a2*dot10*(dot37-2.*dot17*dot5-dot25*M-2.*dot9*(dot14+dot13*M)-2.*a2*dot11*(dot14+dot15+dot13*M)-2.*a1*(dot4-dot9)*(dot14+dot15+dot13*M))+a1*(-2.*(-1.+2.*a1)*dot12*(dot14+dot15)+dot20+dot21-dot31-dot32-2.*dot16*dot4+a1*(dot18+dot19-2.*dot20-2.*dot21+2.*dot31+2.*dot32+dot33+dot34+2.*dot16*dot4))*M2))/(a1*a2*(dot1-dot2-dot3)*(dot1-dot3-dot6)*((-1.+2.*a1)*dot1+a2*(dot2+dot3+dot5+dot6)));
	  diag[21]=-0.03125*(M2*(2.*(dot18*dot2+dot19*dot2+dot18*dot3+dot19*dot3+4.*dot3*dot32-2.*dot32*dot5+2.*dot33*dot5+2.*dot34*dot5+2.*dot18*dot6+2.*dot19*dot6+2.*dot32*dot6-4.*dot17*dot3*dot7+2.*dot37*dot7+4.*dot15*dot4*dot7-2.*(2.*dot12*dot15*dot3+dot11*dot38+2.*dot12*dot14*dot5+2.*dot11*(-dot14+dot15)*dot7)+4.*dot11*dot15*dot8-4.*dot37*dot8-4.*dot15*dot4*dot8+4.*dot17*dot5*dot8-4.*dot17*dot6*dot8+2.*dot16*dot2*dot9-2.*dot16*dot3*dot9+4.*dot38*dot9-4.*dot14*dot7*dot9-dot1*(dot18+dot19+2.*dot32-4.*dot17*dot8+6.*dot16*dot9)+2.*(dot22-2.*dot13*dot7+2.*dot13*dot8)*dot9*M+2.*dot10*(2.*dot11*dot15+dot1*dot17-dot17*(dot2+dot3-2.*dot5+2.*dot6)+2.*dot9*(dot14+dot13*M)))-2.*a1*(4.*dot11*dot16*dot2+dot18*dot2+dot19*dot2+dot18*dot3+dot19*dot3+2.*dot11*dot38-4.*dot16*dot2*dot4-2.*dot38*dot4-4.*dot11*dot14*dot7-8.*dot11*dot15*dot7+4.*dot37*dot7+4.*dot14*dot4*dot7+8.*dot15*dot4*dot7-4.*dot17*dot5*dot7+4.*dot17*dot6*dot7+8.*dot11*dot15*dot8-4.*dot37*dot8-8.*dot15*dot4*dot8+4.*dot17*dot5*dot8-4.*dot17*dot6*dot8+2.*dot16*dot2*dot9-2.*dot16*dot3*dot9+2.*dot38*dot9-8.*dot14*dot7*dot9-4.*dot15*dot7*dot9+4.*dot14*dot8*dot9+4.*dot15*dot8*dot9-dot1*(dot18+dot19+4.*dot17*dot7-4.*dot17*dot8+2.*dot16*(2.*dot11-2.*dot4+dot9))+2.*((dot11-dot4)*(dot22+2.*dot13*(-dot7+dot8))+(dot22+4.*dot13*(-dot7+dot8))*dot9)*M+dot10*(6.*dot1*dot17-2.*dot17*(dot2+dot3-2.*dot5+2.*dot6)+4.*dot11*(dot14+2.*dot15+dot13*M)-4.*(dot37+dot14*dot4+2.*dot15*dot4-2.*dot14*dot9-dot15*dot9+dot13*(dot4-2.*dot9)*M)))+a12*(8.*(dot14+dot15)*(dot10-dot7+dot8)*(dot11-dot4+dot9)+8.*dot13*(dot10-dot7+dot8)*(dot11-dot4+dot9)*M+(-2.*dot12*(dot14+dot15)-dot20-dot21+dot33+dot34+2.*dot17*dot7)*M2)))/(a1*a2*(dot1-dot2-dot3)*(a1*(dot3-dot5)+dot5)*(dot1-dot3-dot6));
	  diag[22]=(M2*(4.*dot11*dot16*dot2-4.*a1*dot11*dot16*dot2-2.*dot2*dot20+2.*a1*dot2*dot20-2.*dot2*dot21+2.*a1*dot2*dot21+2.*dot2*dot31-2.*a1*dot2*dot31+4.*dot2*dot32-2.*a1*dot2*dot32+2.*dot11*dot38-2.*a1*dot11*dot38+2.*a1*dot38*dot4-4.*dot20*dot5+2.*a1*dot20*dot5-4.*dot21*dot5+2.*a1*dot21*dot5+2.*dot31*dot5-2.*a1*dot31*dot5+4.*dot32*dot5-2.*a1*dot32*dot5+4.*dot16*dot4*dot5-4.*a1*dot16*dot4*dot5+4.*dot12*(-(a2*dot14)+(-2.+a1)*dot15)*(dot2+dot5)-4.*dot12*(dot14+dot15)*dot6-2.*dot20*dot6-2.*dot21*dot6+2.*dot31*dot6+2.*dot32*dot6+4.*dot16*dot4*dot6+4.*dot11*dot15*dot8-4.*a1*dot11*dot15*dot8+4.*a1*dot15*dot4*dot8+4.*dot17*dot5*dot8+4.*dot16*dot2*dot9-4.*a1*dot16*dot2*dot9-2.*a1*dot38*dot9-4.*a1*dot15*dot8*dot9+(-2.*dot12*dot26+dot35+2.*dot25*dot8-2.*a2*dot11*(dot22-2.*dot13*dot8)-2.*dot4*(-(a2*dot22)+dot23+2.*a2*dot13*dot8)-2.*a2*(dot22-2.*dot13*dot8)*dot9)*M-2.*dot1*(2.*dot11*dot16-2.*dot20-2.*dot21+dot31+2.*dot32+2.*dot16*dot4+2.*dot17*dot8-a1*(-dot20-dot21+dot31+dot32+2.*dot16*(dot11+dot9))+dot36*M+2.*dot12*(-(a2*dot14)+(-2.+a1)*dot15-dot13*M))+(-1.+2.*a1)*a2*(2.*dot12*(dot14+dot15)+dot20+dot21-dot31-dot32-2.*dot16*dot4)*M2))/(8.*a2*dot2*(dot1-dot5-dot6)*(dot2-a1*dot2+a1*dot6));
	  diag[23]=-0.125*(M2*(2.*dot1*dot18+4.*dot11*dot16*dot3-2.*(dot18*(dot3+dot5+dot6)+dot11*(dot38-2.*dot14*dot7))-(2.*dot11*(dot23+dot27+dot29)+dot35)*M+dot10*(4.*dot17*dot5-4.*dot15*(-(a2*(dot11-dot4))+a1*dot9)+2.*(dot25+2.*a2*dot13*(dot11-dot4+dot9))*M)+3.*a12*(dot18+dot19+dot33+dot34-2.*dot16*dot4)*M2+2.*(dot3*dot32-dot3*dot33-dot3*dot34+dot38*dot4+dot1*(dot19+dot33+dot34-2.*dot16*dot4)-dot32*dot5-dot33*dot6-dot34*dot6+2.*dot16*dot4*dot6-2.*dot17*dot3*dot7-2.*dot14*dot4*dot7+2.*dot17*dot5*dot7+2.*dot16*dot3*dot9+(2.*dot11*dot13+dot24)*dot7*M+(dot23+dot27+dot29-2.*dot13*dot7)*(dot4-dot9)*M-dot19*(dot3+dot5+dot6-M2)+(dot18+dot33+dot34-2.*dot16*dot4)*M2)+a1*(2.*(dot19*dot3+dot3*dot33+dot3*dot34-dot38*dot4-dot1*(dot18+dot19+dot33+dot34-2.*dot16*dot4)+dot19*dot6+dot33*dot6+dot34*dot6-2.*dot16*dot4*dot6+dot18*(dot3+dot6)+2.*dot14*dot4*dot7+dot11*(-2.*dot16*dot3+dot38-2.*dot14*dot7)+(-2.*dot16*dot3+dot38-2.*dot14*dot7)*dot9)+2.*(dot23+dot27+dot29-2.*dot13*dot7)*(dot11-dot4+dot9)*M-5.*(dot18+dot19+dot33+dot34-2.*dot16*dot4)*M2)))/(a1*a2*dot3*(dot1-dot5-dot6)*(-dot1+dot3+dot6));
	  diag[24]=-0.0625*(M2*(-4.*dot1*dot12*dot14+8.*a1*dot1*dot12*dot14-8.*dot1*dot12*dot15+8.*a1*dot1*dot12*dot15+4.*dot1*dot11*dot16-4.*a1*dot1*dot11*dot16+4.*dot12*dot14*dot2-4.*a1*dot12*dot14*dot2+8.*dot12*dot15*dot2-4.*a1*dot12*dot15*dot2-4.*dot11*dot16*dot2+4.*a1*dot11*dot16*dot2-2.*dot1*dot20+2.*a1*dot1*dot20+2.*dot2*dot20-2.*a1*dot2*dot20-2.*dot1*dot21+2.*a1*dot1*dot21+2.*dot2*dot21-2.*a1*dot2*dot21+2.*dot1*dot31-2.*a1*dot1*dot31-2.*dot2*dot31+2.*a1*dot2*dot31+4.*dot1*dot32-2.*a1*dot1*dot32-4.*dot2*dot32+2.*a1*dot2*dot32-2.*a1*dot1*dot33-2.*a1*dot1*dot34-2.*dot11*dot38+2.*a1*dot11*dot38+4.*a1*dot1*dot16*dot4-2.*a1*dot38*dot4+4.*dot12*dot14*dot5-8.*a1*dot12*dot14*dot5+8.*dot12*dot15*dot5-8.*a1*dot12*dot15*dot5+2.*dot20*dot5-2.*a1*dot20*dot5+2.*dot21*dot5-2.*a1*dot21*dot5-2.*dot31*dot5+2.*a1*dot31*dot5-4.*dot32*dot5+2.*a1*dot32*dot5+2.*a1*dot33*dot5+2.*a1*dot34*dot5+4.*dot12*dot14*dot6-4.*a1*dot12*dot14*dot6+4.*dot12*dot15*dot6-4.*a1*dot12*dot15*dot6-2.*dot31*dot6-2.*dot32*dot6+2.*a1*dot33*dot6+2.*a1*dot34*dot6-4.*a1*dot16*dot4*dot6-4.*a1*dot11*dot14*dot7+4.*a12*dot11*dot14*dot7-4.*a1*dot11*dot15*dot7+4.*a12*dot11*dot15*dot7-4.*a1*dot1*dot17*dot7+2.*a1*dot37*dot7+4.*a1*dot14*dot4*dot7-4.*a12*dot14*dot4*dot7+4.*a1*dot15*dot4*dot7-4.*a12*dot15*dot4*dot7+4.*a1*dot17*dot6*dot7+4.*a1*dot11*dot14*dot8-4.*a12*dot11*dot14*dot8-4.*dot11*dot15*dot8+8.*a1*dot11*dot15*dot8-4.*a12*dot11*dot15*dot8+4.*a1*dot1*dot17*dot8-2.*a1*dot37*dot8-4.*a1*dot14*dot4*dot8+4.*a12*dot14*dot4*dot8-8.*a1*dot15*dot4*dot8+4.*a12*dot15*dot4*dot8+4.*dot17*dot6*dot8-4.*a1*dot17*dot6*dot8+4.*dot1*dot16*dot9-8.*a1*dot1*dot16*dot9-4.*dot16*dot2*dot9+4.*a1*dot16*dot2*dot9+2.*a1*dot38*dot9-4.*dot16*dot5*dot9+4.*a1*dot16*dot5*dot9-4.*dot16*dot6*dot9+4.*a1*dot16*dot6*dot9-4.*a1*dot14*dot7*dot9+4.*a12*dot14*dot7*dot9+4.*a12*dot15*dot7*dot9+4.*a1*dot14*dot8*dot9-4.*a12*dot14*dot8*dot9+4.*a1*dot15*dot8*dot9-4.*a12*dot15*dot8*dot9+(-3.*dot35+2.*dot22*(-(a2*dot11)+dot4)+2.*(2.*dot12*dot26+dot1*(-2.*dot12*dot13+dot36+a1*dot36)+2.*dot23*dot4+dot36*(dot5+dot6)-(dot25+dot28+dot30)*dot7-2.*dot11*dot13*dot8+2.*dot24*dot8+(dot25+2.*(dot28+dot30))*dot8-a1*(dot22*dot4+dot36*(dot5+dot6)+(-3.*dot24-2.*dot25-3.*(dot28+dot30)+2.*dot13*(dot11+2.*dot4))*dot7+(4.*dot24+3.*dot25+4.*(dot28+dot30)-4.*dot13*(dot11+dot4))*dot8)-dot22*dot9+a1*(dot22-2.*dot13*dot7+4.*dot13*dot8)*dot9-2.*dot13*dot8*(dot4+dot9)-2.*a12*(dot7-dot8)*(dot24+dot25+dot28+dot30-dot13*(dot11+dot4+dot9))))*M-2.*a2*dot10*(-dot37+2.*a2*(dot14+dot15)*(dot11-dot4)+2.*dot17*dot5-2.*(-(a2*dot14)+a1*dot15)*dot9-(-2.*a2*dot11*dot13+2.*dot24+dot25+2.*(dot28+dot30)-2.*(dot13*(dot4+dot9)+a1*(dot24+dot25+dot28+dot30-dot13*(dot4+dot9))))*M)+a22*(2.*dot12*(dot14+dot15)+dot20+dot21-dot31-dot32-2.*dot16*dot4)*M2))/(a1*a2*(dot1-dot3-dot6)*(dot1-dot5-dot6)*((-1.+2.*a1)*dot1+a2*(dot2+dot3+dot5+dot6)));
	  diag[25]=(M2*(8.*dot1*dot12*dot14-4.*a1*dot1*dot12*dot14-4.*a1*dot1*dot12*dot15-4.*dot1*dot11*dot16-4.*a1*dot1*dot11*dot16-2.*dot1*dot18+2.*a1*dot1*dot18-2.*dot1*dot19+2.*a1*dot1*dot19-8.*dot11*dot16*dot2+8.*a1*dot11*dot16*dot2+2.*dot18*dot2-2.*a1*dot18*dot2+2.*dot19*dot2-2.*a1*dot19*dot2+2.*dot1*dot20-4.*a1*dot1*dot20+2.*dot1*dot21-4.*a1*dot1*dot21-4.*dot12*dot14*dot3+4.*a1*dot12*dot14*dot3+4.*dot12*dot15*dot3+4.*a1*dot12*dot15*dot3+8.*dot11*dot16*dot3-8.*a1*dot11*dot16*dot3-2.*dot20*dot3+2.*a1*dot20*dot3-2.*dot21*dot3+2.*a1*dot21*dot3-4.*dot1*dot31+2.*a1*dot1*dot31+2.*dot3*dot31-2.*a1*dot3*dot31+2.*a1*dot1*dot32+4.*dot2*dot32-2.*dot3*dot32-2.*a1*dot3*dot32-2.*dot1*dot33+2.*a1*dot1*dot33+2.*dot2*dot33-2.*a1*dot2*dot33-2.*dot1*dot34+2.*a1*dot1*dot34+2.*dot2*dot34-2.*a1*dot2*dot34+8.*a1*dot11*dot38-4.*dot1*dot16*dot4+8.*a1*dot1*dot16*dot4+4.*dot16*dot2*dot4-4.*a1*dot16*dot2*dot4-4.*dot16*dot3*dot4+4.*a1*dot16*dot3*dot4+4.*dot38*dot4-8.*a1*dot38*dot4+4.*dot11*dot16*dot5-4.*a1*dot11*dot16*dot5+2.*dot18*dot5-2.*a1*dot18*dot5+2.*dot19*dot5-2.*a1*dot19*dot5-2.*dot20*dot5+2.*a1*dot20*dot5-2.*dot21*dot5+2.*a1*dot21*dot5+2.*dot33*dot5-2.*a1*dot33*dot5+2.*dot34*dot5-2.*a1*dot34*dot5-4.*dot16*dot4*dot5+4.*a1*dot16*dot4*dot5-8.*dot12*dot14*dot6+4.*a1*dot12*dot14*dot6+4.*a1*dot12*dot15*dot6+4.*dot11*dot16*dot6-4.*a1*dot11*dot16*dot6-2.*dot20*dot6+4.*a1*dot20*dot6-2.*dot21*dot6+4.*a1*dot21*dot6+6.*dot31*dot6-2.*a1*dot31*dot6-2.*dot32*dot6-2.*a1*dot32*dot6+2.*dot33*dot6+2.*dot34*dot6+4.*dot16*dot4*dot6-4.*a1*dot16*dot4*dot6-16.*a1*dot11*dot14*dot7+8.*a12*dot11*dot14*dot7-8.*a1*dot11*dot15*dot7+8.*a12*dot11*dot15*dot7-8.*dot17*dot2*dot7+4.*a1*dot37*dot7-8.*dot14*dot4*dot7+16.*a1*dot14*dot4*dot7-8.*a12*dot14*dot4*dot7+8.*a1*dot15*dot4*dot7-8.*a12*dot15*dot4*dot7-8.*a1*dot17*dot5*dot7+4.*dot17*dot6*dot7+8.*a1*dot11*dot14*dot8-8.*a12*dot11*dot14*dot8-8.*dot11*dot15*dot8+16.*a1*dot11*dot15*dot8-8.*a12*dot11*dot15*dot8-4.*dot1*dot17*dot8+4.*a1*dot1*dot17*dot8+8.*dot37*dot8-4.*a1*dot37*dot8-8.*a1*dot14*dot4*dot8+8.*a12*dot14*dot4*dot8+8.*dot15*dot4*dot8-16.*a1*dot15*dot4*dot8+8.*a12*dot15*dot4*dot8-4.*dot17*dot5*dot8+4.*a1*dot17*dot5*dot8+4.*dot17*dot6*dot8-4.*a1*dot17*dot6*dot8+8.*dot1*dot16*dot9-8.*a1*dot1*dot16*dot9-8.*dot16*dot2*dot9+8.*a1*dot16*dot2*dot9+8.*dot16*dot3*dot9-8.*a1*dot16*dot3*dot9-8.*dot38*dot9+8.*a1*dot38*dot9+8.*dot14*dot7*dot9-16.*a1*dot14*dot7*dot9+8.*a12*dot14*dot7*dot9+8.*a12*dot15*dot7*dot9+8.*a1*dot14*dot8*dot9-8.*a12*dot14*dot8*dot9+8.*a1*dot15*dot8*dot9-8.*a12*dot15*dot8*dot9+4.*dot10*(2.*a2*dot11*(-(a2*dot14)+a1*dot15)-(1.+a1)*dot37+2.*(a22*dot14+(-1.-a1*a2)*dot15)*dot4+2.*dot17*(dot1+a1*dot5)-2.*(a22*dot14+a12*dot15)*dot9)+2.*dot10*(dot24+dot25+dot28+dot30-2.*dot13*dot4-4.*a12*dot13*(dot11-dot4+dot9)+a1*(4.*dot11*dot13-dot24+dot25-dot28-dot30-2.*dot13*dot4+4.*dot13*dot9))*M+2.*(4.*a12*dot13*(dot7-dot8)*(dot11-dot4+dot9)-dot8*(4.*dot11*dot13+dot24+3.*dot25+dot28+dot30-6.*dot13*dot4+4.*dot13*dot9)+a1*(-4.*dot11*dot13*(dot7-2.*dot8)+(dot24+dot28+dot30)*dot8+dot25*(-2.*dot7+3.*dot8)+2.*dot13*(2.*dot4*dot7-5.*dot4*dot8-2.*dot7*dot9+4.*dot8*dot9)))*M-a2*(-dot18-dot19-dot33-dot34+2.*dot16*dot4+a1*(2.*dot12*(dot14+dot15)+dot18+dot19+dot20+dot21-dot31-dot32+dot33+dot34-4.*dot16*dot4))*M2))/(32.*a1*a2*(dot1-dot3-dot6)*(dot1-dot5-dot6)*(dot2-a1*dot2+a1*dot6));
	  diag[26]=((4.*dot11*dot16*dot2-2.*dot2*dot32+4.*dot12*dot15*(dot2-dot6)-2.*(dot20+dot21)*dot6+2.*dot32*dot6+4.*dot17*dot6*dot8+(4.*dot1*dot12*dot13-2.*dot11*dot22-2.*dot12*dot26+dot35-2.*dot1*dot36-2.*dot23*dot4+4.*dot11*dot13*dot8+2.*dot25*dot8)*M)*M2)/(8.*dot2*dot6*(-(a2*dot2)-a1*dot6));
	  diag[27]=-0.125*(M2*(-2.*(dot1-dot2-dot5)*(2.*dot12*dot15+dot20+dot21-dot32-2.*dot17*dot8)+(-2.*dot12*(2.*dot13*dot2+dot26)+dot35+dot1*(4.*dot12*dot13-2.*dot36)+2.*dot2*dot36-2.*(dot22+dot23)*dot4+2.*dot25*dot8+4.*dot13*dot8*(dot4-dot9)+2.*dot22*dot9)*M+2.*a12*(2.*dot12*(dot14+dot15)+dot20+dot21-dot31-dot32-2.*dot16*dot4)*M2))/(a2*dot2*sqr(-dot1+dot2+dot5));
	  diag[28]=-0.0625*(M2*(-16.*dot12*dot15*dot2+4.*dot11*dot16*dot2-8.*a1*dot11*dot16*dot2-2.*dot2*dot20-2.*dot2*dot21+8.*dot2*dot32+4.*dot11*dot38-4.*a1*dot11*dot38-8.*dot16*dot2*dot4+8.*a1*dot16*dot2*dot4+4.*a1*dot38*dot4-8.*dot12*dot15*dot5-4.*dot20*dot5-4.*dot21*dot5+4.*dot32*dot5-8.*dot12*dot14*dot6+4.*dot31*dot6+8.*dot16*dot4*dot6+8.*dot11*dot15*dot8-8.*a1*dot11*dot15*dot8+4.*dot17*dot2*dot8+8.*a1*dot15*dot4*dot8+8.*dot17*dot5*dot8-8.*dot17*dot6*dot8+8.*dot16*dot2*dot9-8.*a1*dot16*dot2*dot9-4.*a1*dot38*dot9-8.*a1*dot15*dot8*dot9+4.*dot1*(2.*dot12*dot15+dot20+dot21-dot32-2.*dot17*dot8+2.*dot16*(-(a2*dot11)+a1*(-dot4+dot9)))+4.*(dot22-2.*dot13*dot8)*(dot4-dot9+a1*(dot11-dot4+dot9))*M+a22*(2.*dot12*(dot14+dot15)+dot20+dot21-dot31-dot32-2.*dot16*dot4)*M2))/(a2*dot2*(-dot1+dot2+dot5)*(-(a2*dot2)-a1*dot6));
	  diag[29]=((2.*dot5*(dot18+dot19+dot32-2.*dot17*(dot10+dot7))-2.*dot3*(dot32-2.*dot17*dot7+2.*dot16*dot9)+(2.*dot10*dot25-dot35+2.*dot24*dot7+2.*(dot23+dot27+dot29-2.*dot13*(dot10+dot7))*dot9)*M)*M2)/(8.*dot3*dot5*(a1*(dot3-dot5)+dot5));
	  diag[30]=-0.125*((2.*(dot1-dot2-dot5)*(dot32-2.*dot17*dot7+2.*dot16*dot9)+dot35*M-2.*(dot36*dot5+dot24*dot7-dot22*dot9+(dot7-dot8)*(dot25-2.*dot13*dot9))*M)*M2)/(a2*dot5*sqr(-dot1+dot2+dot5));
	  diag[31]=((2.*dot3*dot32+3.*dot31*dot5-dot32*dot5-4.*dot17*dot3*dot7+2.*a1*dot37*dot7+6.*dot17*dot5*dot7-4.*a1*dot17*dot5*dot7-2.*a1*dot37*dot8-4.*dot17*dot5*dot8+4.*a1*dot17*dot5*dot8+4.*a1*dot10*dot14*dot9+4.*dot16*dot3*dot9-2.*dot16*dot5*dot9-4.*a1*dot14*dot7*dot9+4.*a1*dot14*dot8*dot9-2.*dot10*(-(a2*dot37)+dot17*dot5-2.*a1*dot17*dot5+2.*dot14*dot9)+2.*dot1*(dot32-2.*dot17*dot7+2.*dot16*dot9)-2.*dot2*(dot32-2.*dot17*dot7+2.*dot16*dot9)-2.*(dot7-dot8+a1*(dot10-dot7+dot8))*(dot25-2.*dot13*dot9)*M)*M2)/(8.*a2*(dot1-dot2-dot5)*dot5*(a1*(dot3-dot5)+dot5));
	  diag[32]=-0.125*(M2*(-2.*(dot1-dot3-dot6)*(dot18+dot19+dot32-2.*dot17*(dot10+dot7))+(dot35-2.*dot3*dot36+2.*(dot23+dot27+dot29)*(dot11-dot4)-2.*dot10*(2.*dot11*dot13+dot25-2.*dot13*dot4)-2.*(2.*dot11*dot13+dot24-2.*dot13*dot4)*dot7)*M+2.*a22*(dot18+dot19+dot33+dot34-2.*dot16*dot4)*M2))/(a1*dot3*sqr(-dot1+dot3+dot6));
	  diag[33]=(M2*(4.*dot10*dot17*(-2.*dot1+dot3-2.*dot5+2.*dot6)+4.*dot11*(2.*a1*dot16*dot3+dot38-a1*dot38-2.*dot14*dot7)+4.*dot1*(dot18+dot19+dot32-2.*dot17*dot7)+8.*dot10*dot15*(-(a2*(dot11-dot4))+a1*dot9)-2.*(dot18*(dot3+2.*dot6)+dot19*(dot3+2.*dot6)+2.*((-dot32+dot33+dot34)*dot5+dot32*dot6+dot38*(dot4-a1*dot4+a1*dot9)+dot3*(2.*dot32-4.*dot17*dot7+dot16*(2.*a1*dot4+dot9-2.*a1*dot9))-2.*(dot16*dot4*dot5+dot17*(-dot5+dot6)*dot7+dot14*dot7*(dot4+a1*(dot11-dot4+dot9)))))-4.*(dot23+dot27+dot29-2.*dot13*(dot10+dot7))*(-dot9+a1*(dot11-dot4+dot9))*M-a12*(dot18+dot19+dot33+dot34-2.*dot16*dot4)*M2))/(16.*a1*dot3*(a1*(dot3-dot5)+dot5)*(-dot1+dot3+dot6));
	  diag[34]=((-2.*(2.*dot12*dot15+2.*dot11*dot16-dot32)*(dot1-dot3-dot6)-(4.*dot1*dot12*dot13-2.*dot10*(2.*dot11*dot13+dot25)-2.*dot12*dot26+2.*dot11*(dot23+dot27+dot29)+dot35-2.*dot1*dot36-2.*dot23*dot4-4.*dot12*dot13*dot6+2.*dot36*dot6+2.*dot25*dot7)*M)*M2)/(8.*a1*dot6*sqr(-dot1+dot3+dot6));
	  diag[35]=-0.125*((-4.*dot11*dot16*dot2-4.*dot12*dot15*(dot1+dot2-dot3)+4.*dot11*dot16*dot3+2.*dot2*dot32-2.*dot3*dot32+2.*dot12*(-3.*dot14+dot15)*dot6+2.*dot11*dot16*dot6+2.*dot31*dot6-2.*dot32*dot6+dot33*dot6+dot34*dot6+4.*dot16*dot4*dot6+4.*dot11*dot14*dot7-4.*a1*dot11*dot14*dot7+2.*dot37*dot7-2.*a1*dot37*dot7+4.*dot15*dot4*dot7-4.*a1*dot15*dot4*dot7+2.*dot17*dot6*dot7-4.*a1*dot17*dot6*dot7+4.*a1*dot11*dot14*dot8+2.*a1*dot37*dot8+4.*a1*dot15*dot4*dot8-2.*dot17*dot6*dot8+4.*a1*dot17*dot6*dot8+2.*dot1*(-2.*dot11*dot16+dot32-2.*a2*dot17*dot7-2.*a1*dot17*dot8)-2.*(2.*dot11*dot13+dot25)*(a1*(dot7-dot8)+dot8)*M+2.*dot10*(2.*a2*dot1*dot17-dot37-2.*dot15*dot4+2.*dot11*(-(a2*dot14)+a1*dot13*M)+a1*(dot37+2.*dot15*dot4+2.*dot17*dot6+dot25*M)))*M2)/(a1*(dot1-dot3-dot6)*dot6*(dot2-a1*dot2+a1*dot6));
	  diag[36]=-0.0625*(-4.*dot16*dot2*dot4+4.*dot16*dot3*dot4-2.*dot38*dot4-2.*dot11*dot16*dot5+dot20*dot5+dot21*dot5+4.*dot16*dot4*dot6-
			    4.*(dot12*dot15*dot2+dot10*dot14*dot4+dot12*dot14*dot6)+4.*dot11*dot14*dot7+4.*dot17*dot2*dot7+2.*dot37*dot7+
			    4.*dot14*dot4*dot7+4.*dot15*dot4*dot7+4.*dot11*dot15*dot8-4.*dot37*dot8-4.*dot15*dot4*dot8+2.*dot17*dot5*dot8-
			    4.*dot17*dot6*dot8+4.*dot10*dot14*dot9+4.*dot16*dot2*dot9-4.*dot16*dot3*dot9+4.*dot38*dot9-4.*dot14*dot7*dot9+
			    dot1*(dot20+dot21-4.*dot17*dot7+2.*dot17*dot8-2.*dot16*(dot11+2.*dot9))+4.*dot13*dot8*(-dot4+dot9)*M+
			    4.*a12*(dot10-dot7+dot8)*(dot11-dot4+dot9)*(dot14+dot15+dot13*M)
			    -a1*(4.*dot11*dot16*dot2-2.*dot11*dot16*dot3-dot20*dot3-dot21*dot3+4.*dot11*dot38-4.*dot16*dot2*dot4+4.*dot16*dot3*dot4-
				 4.*dot38*dot4-2.*dot11*dot16*dot5+dot20*dot5+dot21*dot5-4.*dot11*dot14*dot7-4.*dot11*dot15*dot7+4.*dot37*dot7+
				 8.*dot14*dot4*dot7+8.*dot15*dot4*dot7-4.*dot17*dot5*dot7+4.*dot17*dot6*dot7+8.*dot11*dot15*dot8+2.*dot17*dot3*dot8-
				 4.*dot37*dot8-4.*dot14*dot4*dot8-12.*dot15*dot4*dot8+2.*dot17*dot5*dot8-4.*dot17*dot6*dot8+4.*dot16*dot2*dot9-
				 4.*dot16*dot3*dot9+4.*dot38*dot9-8.*dot14*dot7*dot9+4.*dot14*dot8*dot9+4.*dot15*dot8*dot9-
				 4.*dot1*(dot17*(dot7-dot8)+dot16*(dot11-dot4+dot9))+4.*dot13*(dot11*dot8+(dot7-2.*dot8)*(dot4-dot9))*M+
				 4.*dot10*(dot11*(dot14+dot15)-dot37-2.*dot14*dot4-2.*dot15*dot4+dot17*(dot1+dot5-dot6)+
					   2.*dot14*dot9-dot13*dot4*M+dot13*dot9*M)))*M2/
	    (a1*a2*(dot1-dot2-dot5)*(dot1-dot3-dot6)*(-(a2*dot2)-a1*dot6));
	  diag[37]=(-2.*dot10*dot17*dot2+dot18*dot2+dot19*dot2+4.*dot10*dot15*dot4-2.*dot38*dot4+4.*dot16*dot4*dot5-
		    4.*dot12*(dot15*dot3+dot14*dot5)-4.*dot11*dot15*dot7+4.*dot17*dot3*dot7+2.*dot37*dot7+4.*dot14*dot4*dot7+
		    4.*dot15*dot4*dot7-4.*dot17*dot5*dot7+4.*dot17*dot6*dot7+4.*dot11*dot15*dot8-4.*dot37*dot8-4.*dot15*dot4*dot8+
		    4.*dot17*dot5*dot8-4.*dot17*dot6*dot8+4.*dot10*dot14*dot9+2.*dot16*dot2*dot9-4.*dot16*dot3*dot9+4.*dot38*dot9-
		    4.*dot14*dot7*dot9+dot1*(dot18+dot19-2.*dot17*(dot10+2.*dot7-2.*dot8)-6.*dot16*dot9)+4.*dot13*(-dot7+dot8)*dot9*M+
		    4.*a12*(dot10-dot7+dot8)*(dot11-dot4+dot9)*(dot14+dot15+dot13*M)-
		    a1*(4.*dot11*dot16*dot2+dot18*dot2+dot19*dot2-4.*dot11*dot16*dot3+4.*dot11*dot38-4.*dot16*dot2*dot4+4.*dot16*dot3*dot4-
			4.*dot38*dot4-dot18*dot6-dot19*dot6-8.*dot11*dot14*dot7-8.*dot11*dot15*dot7+4.*dot37*dot7+8.*dot14*dot4*dot7+
			8.*dot15*dot4*dot7-4.*dot17*dot5*dot7+4.*dot17*dot6*dot7+8.*dot11*dot15*dot8-4.*dot37*dot8-8.*dot15*dot4*dot8+
			4.*dot17*dot5*dot8-4.*dot17*dot6*dot8+2.*dot16*dot2*dot9-4.*dot16*dot3*dot9+4.*dot38*dot9+2.*dot16*dot6*dot9-
			12.*dot14*dot7*dot9-4.*dot15*dot7*dot9+4.*dot14*dot8*dot9+4.*dot15*dot8*dot9-4.*dot1*(dot17*(dot7-dot8)+dot16*(dot11-dot4+dot9))-
			4.*dot13*(dot7-dot8)*(dot11-dot4+2.*dot9)*M+
			dot10*(4.*dot11*(dot14+dot15)+4.*dot1*dot17-2.*dot17*(dot2-2.*dot5+dot6)-
			       4.*(dot37+dot14*dot4+dot15*dot4-2.*dot14*dot9-dot13*dot9*M))))*M2/
	    (16.*a1*a2*(dot1-dot2-dot5)*(a1*(dot3-dot5)+dot5)*(dot1-dot3-dot6));




















	  // Complex diagB[38];

	  // diagB[0]=-0.125*((dot15*(dot2+dot3)+dot14*(dot5+dot6)-2.*dot1*(dot14+dot15+dot13*M))*M2)/(a1*a2*(dot1-dot2-dot5)*(dot1-dot3-dot6)*MeV2);
	  // diagB[1]=((dot15*(dot1-2.*dot2+dot3)+dot14*(dot1+dot5-2.*dot6)+dot1*dot13*M)*M2)/(4.*a1*a2*(dot1-dot2-dot5)*(dot1-dot3-dot6)*MeV2);
	  // diagB[2]=((dot15*(dot1+dot2-2.*dot3)+dot14*(dot1-2.*dot5+dot6)+dot1*dot13*M)*M2)/(4.*a1*a2*(dot1-dot2-dot5)*(dot1-dot3-dot6)*MeV2);
	  // diagB[3]=-0.25*((-(dot15*((1.+a1)*dot2+(-2.+a1)*dot3+dot5-a1*(dot5+dot6)))+dot14*(dot2-a1*dot2+dot5-2.*dot6+a1*(-dot3+dot5+dot6))+dot13*(dot2-a1*dot2-dot5+a1*(-dot3+dot5+dot6))*M)*M2)/(a1*a2*(dot1-dot2-dot5)*(dot1-dot3-dot6)*MeV2);
	  // diagB[4]=((dot15*((1.+a1)*dot2-a2*(dot3-dot5-dot6))+dot14*(-(a2*dot2)-a2*dot3-dot5+dot6-a1*(dot5+dot6))+(dot26-a2*dot13*(dot2+dot3-dot5-dot6))*M+dot1*(dot14-dot15-dot13*M))*M2)/(4.*a1*(dot1-dot3-dot6)*((-1.+2.*a1)*dot1+a2*(dot2+dot3+dot5+dot6))*MeV2);
	  // diagB[5]=-0.25*(((dot14-dot15)*(dot1-dot3-dot6)+dot26*M-dot13*(dot1+dot3-dot6)*M)*M2)/(a1*sqr(-dot1+dot3+dot6)*MeV2);
	  // diagB[6]=-0.25*((-2.*dot15*dot3+2.*dot14*dot6+dot26*M-dot1*(dot14-dot15+dot13*M)+a1*(dot2+dot3-dot5-dot6)*(dot14+dot15+dot13*M))*M2)/(a2*(dot1-dot2-dot5)*((-1.+2.*a1)*dot1-a1*(dot2+dot3+dot5+dot6))*MeV2);
	  // diagB[7]=(((dot14-dot15)*(dot1-dot2-dot5)-dot26*M+dot13*(dot1-dot2+dot5)*M)*M2)/(4.*a2*sqr(-dot1+dot2+dot5)*MeV2);
	  // diagB[8]=((dot14+dot15+dot13*M)*M2)/((-4.*a2*dot2-4.*a1*dot6)*MeV2);
	  // diagB[9]=-0.5*((dot14+dot15+dot13*M)*M2)/(a2*(dot1-dot2-dot5)*MeV2);
	  // diagB[10]=((dot14+dot15+dot13*M)*M2)/((-2.*a2*dot2-2.*a1*dot6)*MeV2);
	  // diagB[11]=-0.125*((-(dot26*M)+2.*(dot14*(dot1-a1*dot1-dot6+a1*(dot5+dot6))+dot15*(dot3+a1*(-dot1+dot5+dot6))+dot13*(dot1-a1*dot1+a1*(dot5+dot6))*M))*M2)/(a2*(dot1-dot2-dot5)*((-1.+2.*a1)*dot1-a1*(dot2+dot3+dot5+dot6))*MeV2);
	  // diagB[12]=-0.125*((-2.*a2*dot14*dot5+(-1.+2.*a1)*dot14*dot6-dot15*(dot2+2.*dot5-2.*a1*dot5-2.*a1*dot6)+2.*dot13*(-(a2*dot5)+a1*dot6)*M+2.*a2*dot1*(dot14+dot15+dot13*M))*M2)/(a2*(dot1-dot2-dot5)*(-(a2*dot2)-a1*dot6)*MeV2);
	  // diagB[13]=((dot14+dot15+dot13*M)*M2)/(4.*(a1*(dot3-dot5)+dot5)*MeV2);
	  // diagB[14]=-0.5*((dot14+dot15+dot13*M)*M2)/(a2*(dot1-dot2-dot5)*MeV2);
	  // diagB[15]=-0.5*((dot14+dot15+dot13*M)*M2)/((a1*(dot3-dot5)+dot5)*MeV2);
	  // diagB[16]=((2.*dot1*dot15-2.*dot15*dot3+2.*dot14*dot6+dot26*M-2.*a1*(dot1-dot2-dot3)*(dot14+dot15+dot13*M))*M2)/(8.*a2*(dot1-dot2-dot5)*((-1.+2.*a1)*dot1-a1*(dot2+dot3+dot5+dot6))*MeV2);
	  // diagB[17]=((-2.*a2*dot15*dot2+(-1.+2.*a1)*dot15*dot3+dot14*(-2.*a2*dot2+2.*a1*dot3-dot5)+2.*dot13*(-(a2*dot2)+a1*dot3)*M+2.*a2*dot1*(dot14+dot15+dot13*M))*M2)/(8.*a2*(dot1-dot2-dot5)*(a1*(dot3-dot5)+dot5)*MeV2);
	  // diagB[18]=-0.5*((dot14+dot15+dot13*M)*M2)/((a1*(dot3-dot5)+dot5)*MeV2);
	  // diagB[19]=((dot14+dot15+dot13*M)*M2)/(2.*a1*(-dot1+dot3+dot6)*MeV2);
	  // diagB[20]=((2.*dot14*dot6+2.*dot15*(dot2+dot5+dot6)+(dot26+2.*dot13*(dot5+dot6))*M-2.*dot1*(dot15+dot13*M)+2.*a1*(dot1-dot5-dot6)*(dot14+dot15+dot13*M))*M2)/(8.*a1*(dot1-dot3-dot6)*((-1.+2.*a1)*dot1+a2*(dot2+dot3+dot5+dot6))*MeV2);
	  // diagB[21]=((dot15*(dot3-2.*dot5)-2.*a1*(dot1-dot5-dot6)*(dot14+dot15+dot13*M)-dot5*(dot14+2.*dot13*M))*M2)/(8.*a1*(a1*(dot3-dot5)+dot5)*(-dot1+dot3+dot6)*MeV2);
	  // diagB[22]=((dot14+dot15+dot13*M)*M2)/((-2.*a2*dot2-2.*a1*dot6)*MeV2);
	  // diagB[23]=((dot14+dot15+dot13*M)*M2)/(2.*a1*(-dot1+dot3+dot6)*MeV2);
	  // diagB[24]=((dot26*M-2.*(dot15*(a1*dot1+dot3-a1*(dot2+dot3))+dot14*(-(a2*(dot1-dot2-dot3))+dot5)+dot13*(a1*dot1+a2*(dot2+dot3))*M))*M2)/(8.*a1*(dot1-dot3-dot6)*((-1.+2.*a1)*dot1+a2*(dot2+dot3+dot5+dot6))*MeV2);
	  // diagB[25]=((2.*dot14*dot2+dot15*dot2-dot14*dot6+2.*dot13*dot2*M+2.*a1*(dot1-dot2-dot3)*(dot14+dot15+dot13*M))*M2)/(8.*a1*(dot1-dot3-dot6)*(dot2-a1*dot2+a1*dot6)*MeV2);
	  // diagB[26]=((dot14+dot15+dot13*M)*M2)/((-2.*a2*dot2-2.*a1*dot6)*MeV2);
	  // diagB[27]=-0.25*((2.*dot1*dot15-2.*dot15*(dot2+dot5)+(dot26-2.*dot13*dot5)*M)*M2)/(a2*sqr(-dot1+dot2+dot5)*MeV2);
	  // diagB[28]=-0.25*((-2.*a2*dot14*dot5+(-1.+2.*a1)*dot14*dot6-dot15*(dot2+2.*dot5-2.*a1*dot5-2.*a1*dot6)+2.*dot13*(-(a2*dot5)+a1*dot6)*M+2.*a2*dot1*(dot14+dot15+dot13*M))*M2)/(a2*(dot1-dot2-dot5)*(-(a2*dot2)-a1*dot6)*MeV2);
	  // diagB[29]=-0.5*((dot14+dot15+dot13*M)*M2)/((a1*(dot3-dot5)+dot5)*MeV2);
	  // diagB[30]=((2.*dot14*(dot2+dot5)+(2.*dot13*dot2+dot26)*M-2.*dot1*(dot14+dot13*M))*M2)/(4.*a2*sqr(-dot1+dot2+dot5)*MeV2);
	  // diagB[31]=((-2.*a2*dot15*dot2+(-1.+2.*a1)*dot15*dot3+dot14*(-2.*a2*dot2+2.*a1*dot3-dot5)+2.*dot13*(-(a2*dot2)+a1*dot3)*M+2.*a2*dot1*(dot14+dot15+dot13*M))*M2)/(4.*a2*(dot1-dot2-dot5)*(a1*(dot3-dot5)+dot5)*MeV2);
	  // diagB[32]=((2.*dot15*(dot3+dot6)+(dot26+2.*dot13*dot6)*M-2.*dot1*(dot15+dot13*M))*M2)/(4.*a1*sqr(-dot1+dot3+dot6)*MeV2);
	  // diagB[33]=((dot15*(dot3-2.*dot5)-2.*a1*(dot1-dot5-dot6)*(dot14+dot15+dot13*M)-dot5*(dot14+2.*dot13*M))*M2)/(4.*a1*(a1*(dot3-dot5)+dot5)*(-dot1+dot3+dot6)*MeV2);
	  // diagB[34]=((2.*dot14*(-dot1+dot3+dot6)-dot26*M+2.*dot13*dot3*M)*M2)/(4.*a1*sqr(-dot1+dot3+dot6)*MeV2);
	  // diagB[35]=((2.*dot14*dot2+dot15*dot2-dot14*dot6+2.*dot13*dot2*M+2.*a1*(dot1-dot2-dot3)*(dot14+dot15+dot13*M))*M2)/(4.*a1*(dot1-dot3-dot6)*(dot2-a1*dot2+a1*dot6)*MeV2);
	  // diagB[36]=((a1*(2.*sqr(dot1)*(dot14+dot15)+2.*dot15*dot3*dot5+2.*dot14*(2.*dot2+dot3)*dot5-2.*dot1*(dot14+dot15)*(2.*dot2+dot3+dot5)+dot14*(2.*dot2+dot3-dot5)*dot6+dot15*dot2*(dot3+3.*dot5+2.*dot6))+2.*a1*dot13*((dot1-2.*dot2-dot3)*(dot1-dot5)+dot2*dot6)*M-2.*a12*(dot1-dot2-dot3)*(dot1-dot5-dot6)*(dot14+dot15+dot13*M)+(dot1-dot5)*(2.*dot14*dot2+dot15*dot2-dot14*dot6+2.*dot13*dot2*M))*M2)/(8.*a1*a2*(dot1-dot2-dot5)*(dot1-dot3-dot6)*(-(a2*dot2)-a1*dot6)*MeV2);
	  // diagB[37]=((2.*a12*(dot1-dot2-dot3)*(dot1-dot5-dot6)*(dot14+dot15+dot13*M)-a1*(dot15*(-(dot2*dot3)+4.*dot2*dot5+2.*dot3*dot5+2.*dot2*dot6+dot3*dot6)+dot14*(3.*dot2*dot5+2.*dot3*dot5+2.*dot2*dot6+dot5*dot6)+2.*dot13*(2.*dot2*dot5+dot3*dot5+dot2*dot6)*M+2.*sqr(dot1)*(dot14+dot15+dot13*M)-2.*dot1*(dot2+2.*dot5+dot6)*(dot14+dot15+dot13*M))+(dot1-dot2)*(dot15*(dot3-2.*dot5)-dot5*(dot14+2.*dot13*M)))*M2)/(8.*a1*a2*(dot1-dot2-dot5)*(a1*(dot3-dot5)+dot5)*(dot1-dot3-dot6)*MeV2);

	  // for(unsigned int ix=0;ix<38;++ix)
	  //   cerr << "testing diags " << diag[ix] << " " << (diag[ix]-diagB[ix])/(diag[ix]+diagB[ix]) << "\n";



	  
	  Complex flow[3];
	  flow[0] = diag[0] + diag[8] - diag[11] + diag[12] - diag[13] + diag[16] + diag[17] - diag[20] + diag[21] + diag[24] + diag[25] + diag[36] + diag[37];
	  flow[1] = diag[1] + diag[3] + diag[5] + diag[7] + diag[29] + diag[30] + diag[31]  + diag[32] + diag[33] + 2.*diag[37] +
	    (diag[4] - diag[5] + diag[6] - diag[7] + diag[9] + 2.*diag[11] + 2.*diag[13] + diag[15] + diag[18] + diag[23] - 2.*diag[24] - diag[29] - diag[30] - diag[32])/9.;
	  flow[2] = diag[2] - diag[3] - diag[5] - diag[7] + diag[26] + diag[27] + diag[28]  + diag[34] + diag[35] + 2.*diag[36] +
	    (-diag[4] + diag[5] - diag[6] + diag[7] - 2.*diag[8] + diag[10] + diag[14] - 2.*diag[16] + diag[19] + 2.*diag[20] + diag[22] - diag[26] - diag[27] - diag[34])/9.;


	  // cerr << "testing leading " << (diag[2]-diag[3]-diag[5]-diag[7]+diag[26]+diag[27]+diag[28]+diag[34]+diag[35]+2.*diag[36]) << "\n";
	  // cerr << "testing leading " << (diagB[2]-diagB[3]-diagB[5]-diagB[7]+diagB[26]+diagB[27]+diagB[28]+diagB[34]+diagB[35]+2.*diagB[36]) << "\n";

	  // cerr << "testing subleading " << (-diag[4]+diag[5]-diag[6]+diag[7]-2.*diag[8]+diag[10]+diag[14]-2.*diag[16]+diag[19]+2.*diag[20]+diag[22]-diag[26]-diag[27]-diag[34])/3. << "\n";
	  // cerr << "testing subleading " << (-diagB[4]+diagB[5]-diagB[6]+diagB[7]-2.*diagB[8]+diagB[10]+diagB[14]-2.*diagB[16]+diagB[19]+2.*diagB[20]+diagB[22]-diagB[26]-diagB[27]-diagB[34])/3. << "\n";
	  
	  // cerr << "testing totals " << ih1 << " " << ih2 << " " << ih4 << " " << ih5 << " " << flow[0] << " " << flow[1] << " " << flow[2] << "\n";
	  for(unsigned int ix=0;ix<3;++ix) {
	    flows[ix]+=norm(flow[ix]);
	    for(unsigned int iy=0;iy<3;++iy) {
	      // cerr << "testing in loop" << ix << " " << iy << " " << cMatrix[ix][iy] << "\n";
	      me2Sum+=cMatrix[ix][iy]*flow[ix]*std::conj(flow[iy]);
	    }
	  }
	}
      }
    }
  }
  // return the answer
  meInfo(flows);
  double pre = O1_*sHat()/M/72.*pow<4,1>(Constants::pi*standardModel()->alphaS(scale())/M);
  return pre*me2Sum.real();
}

