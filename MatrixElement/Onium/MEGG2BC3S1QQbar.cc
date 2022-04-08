// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGG2BC3S1QQbar class.
//

#include "MEGG2BC3S1QQbar.h"
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

IBPtr MEGG2BC3S1QQbar::clone() const {
  return new_ptr(*this);
}

IBPtr MEGG2BC3S1QQbar::fullclone() const {
  return new_ptr(*this);
}

void MEGG2BC3S1QQbar::doinit() {
  GGtoBCQQbarBase::doinit();
  O1_ = oniumParameters()->singletMEProduction<0>(bcbar,principleQuantumNumber(),1,1);
}

void MEGG2BC3S1QQbar::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2);
}

void MEGG2BC3S1QQbar::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEGG2BC3S1QQbar,GGtoBCQQbarBase>
describeHerwigMEGG2BC3S1QQbar("Herwig::MEGG2BC3S1QQbar", "HwOniumParameters.so HwMEHadronOnium.so");

void MEGG2BC3S1QQbar::Init() {

  static ClassDocumentation<MEGG2BC3S1QQbar> documentation
    ("There is no documentation for the MEGG2BC3S1QQbar class");

}

double MEGG2BC3S1QQbar::me2() const {
  using namespace ThePEG::Helicity;
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
  SpinorWaveFunction      q4w(rescaledMomenta()[3],mePartonData()[3],outgoing);
  SpinorBarWaveFunction   q5w(rescaledMomenta()[4],mePartonData()[4],outgoing);
  vector<SpinorWaveFunction> v4;
  vector<SpinorBarWaveFunction> ubar5;
  for(unsigned int ix=0;ix<2;++ix) {
    q4w.reset(ix);
    v4.push_back(q4w);
    q5w.reset(ix);
    ubar5.push_back(q5w);
  }
  // B_C* wavefunction
  VectorWaveFunction Bcw(rescaledMomenta()[2],mePartonData()[2],outgoing);
  vector<VectorWaveFunction> v3;
  for(unsigned int ix=0;ix<3;++ix) {
    Bcw.reset(ix);
    v3.push_back(Bcw);
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
  ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);
  Complex me2Sum=0.;
  vector<double> flows(3,0.);
  for(unsigned int ih1=0;ih1<2;++ih1) {
    for(unsigned int ih2=0;ih2<2;++ih2) {
      for(unsigned int ih4=0;ih4<2;++ih4) {
  	for(unsigned int ih5=0;ih5<2;++ih5) {
	  for(unsigned int ih3=0;ih3<2;++ih3) {
	    Energy2 dot1 = rescaledMomenta()[0]*rescaledMomenta()[1];
	    Energy2 dot2 = rescaledMomenta()[0]*rescaledMomenta()[3];
	    Energy2 dot3 = rescaledMomenta()[0]*rescaledMomenta()[4];
	    complex<Energy> dot4 = rescaledMomenta()[0]*g2[ih2].wave();
	    complex<Energy> dot5 = rescaledMomenta()[0]*v3[ih3].wave();
	    complex<Energy> dot6 = rescaledMomenta()[1]*g1[ih1].wave();
	    complex<Energy> dot7 = rescaledMomenta()[1]*v3[ih3].wave();
	    Energy2 dot8 = rescaledMomenta()[2]*rescaledMomenta()[3];
	    Energy2 dot9 = rescaledMomenta()[2]*rescaledMomenta()[4];
	    complex<Energy> dot10 = rescaledMomenta()[3]*g1[ih1].wave();
	    complex<Energy> dot11 = rescaledMomenta()[3]*g2[ih2].wave();
	    complex<Energy> dot12 = rescaledMomenta()[4]*g1[ih1].wave();
	    complex<Energy> dot13 = rescaledMomenta()[4]*g2[ih2].wave();
	    complex<Energy> dot14 = rescaledMomenta()[4]*v3[ih3].wave();
	    Complex dot15 = g1[ih1].wave()*g2[ih2].wave();
	    Complex dot16 = g1[ih1].wave()*v3[ih3].wave();
	    Complex dot17 = g2[ih2].wave()*v3[ih3].wave();
	    complex<Energy> dot18=v4[ih4].dimensionedWave().scalar(ubar5[ih5].dimensionedWave());
	    LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzPolarizationVectorE vec3 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzPolarizationVectorE vec4 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzVector<complex<Energy2> > vec5 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[1]).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzPolarizationVectorE vec6 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g1[ih1].wave()).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzPolarizationVectorE vec7 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g2[ih2].wave()).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzVector<complex<Energy2> > vec8 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(rescaledMomenta()[1]).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzPolarizationVectorE vec9 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(g1[ih1].wave()).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzVector<complex<Energy2> > vec10 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(rescaledMomenta()[1]).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzVector<complex<Energy2> > vec11 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[1]).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzPolarizationVectorE vec12 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g2[ih2].wave()).slash(g1[ih1].wave()).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzVector<complex<Energy2> > vec13 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(g1[ih1].wave()).slash(rescaledMomenta()[1]).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzVector<complex<Energy2> > vec14 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g2[ih2].wave()).slash(rescaledMomenta()[1]).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzVector<complex<Energy2> > vec15 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g1[ih1].wave()).slash(rescaledMomenta()[1]).vectorCurrent(ubar5[ih5].dimensionedWave());
	    LorentzVector<complex<Energy2> > vec16 = v4[ih4].dimensionedWave().slash(v3[ih3].wave()).slash(g2[ih2].wave()).slash(g1[ih1].wave()).slash(rescaledMomenta()[1]).vectorCurrent(ubar5[ih5].dimensionedWave());
	    complex<Energy> dot19 = vec1*v3[ih3].wave();
	    complex<Energy2> dot20 = vec12*rescaledMomenta()[0];
	    complex<Energy2> dot21 = vec12*rescaledMomenta()[1];
	    complex<Energy2> dot22 = vec2*rescaledMomenta()[0];
	    complex<Energy2> dot23 = vec2*rescaledMomenta()[1];
	    complex<Energy> dot24 = vec2*g1[ih1].wave();
	    complex<Energy> dot25 = vec2*g2[ih2].wave();
	    complex<Energy2> dot26 = vec3*rescaledMomenta()[0];
	    complex<Energy2> dot27 = vec3*rescaledMomenta()[1];
	    complex<Energy> dot28 = vec3*g1[ih1].wave();
	    complex<Energy2> dot29 = vec4*rescaledMomenta()[0];
	    complex<Energy2> dot30 = vec4*rescaledMomenta()[1];
	    complex<Energy3> dot31 = vec5*rescaledMomenta()[0];
	    complex<Energy2> dot32 = vec6*rescaledMomenta()[0];
	    complex<Energy2> dot33 = vec6*rescaledMomenta()[1];
	    complex<Energy2> dot34 = vec7*rescaledMomenta()[0];
	    complex<Energy2> dot35 = vec7*rescaledMomenta()[1];
	    complex<Energy3> dot36 = vec8*rescaledMomenta()[0];
	    complex<Energy2> dot37 = vec1*rescaledMomenta()[0];
	    complex<Energy2> dot38 = vec1*rescaledMomenta()[1];
	    complex<Energy> dot39 = vec1*g1[ih1].wave();
	    complex<Energy> dot40 = vec1*g2[ih2].wave();
	    complex<Energy3> dot41 = vec10*rescaledMomenta()[0];
	    complex<Energy3> dot42 = vec13*rescaledMomenta()[0];
	    complex<Energy3> dot43 = vec14*rescaledMomenta()[0];
	    complex<Energy3> dot44 = vec15*rescaledMomenta()[0];
	    complex<Energy3> dot45 = vec16*rescaledMomenta()[0];
	    complex<Energy> dot46 = vec7*g1[ih1].wave();
	    complex<Energy2> dot47 = vec9*rescaledMomenta()[1];
	    complex<Energy3> dot48 = vec11*rescaledMomenta()[0];
	    complex<Energy2> dot49 = vec9*rescaledMomenta()[0];
	    // diagrams
	    Complex diag[38];
	    diag[0]=-0.125*(-2.*dot14*dot15*dot18+2.*dot12*dot17*dot18-dot15*(dot22+dot23)+
	    		    dot13*(2.*dot16*dot18-dot24)-dot11*dot24-dot10*dot25-dot12*dot25+
	    		    dot16*(dot26+dot27)+dot17*dot29+dot17*dot30-2.*dot16*dot18*dot4+dot24*dot4+
	    		    2.*dot15*dot18*dot5-2.*dot17*dot18*dot6+dot25*dot6+2.*dot15*dot18*dot7+
	    		    2.*dot15*dot19*M)*M2/(a1*a2*(dot8+a1*M2)*(dot9+a2*M2));
	    diag[1]=(-3.*dot20-3.*dot21+4.*dot15*(dot22+dot23)-2.*dot11*dot24+4.*dot13*(-2.*dot16*dot18+dot24)+
	    	     4.*dot10*dot25-4.*dot16*(dot26+dot27)+dot14*(8.*dot15*dot18-6.*dot28)+
	    	     2.*(dot12*(2.*dot17*dot18-dot25)+4.*dot16*dot18*dot4-2.*dot24*dot4-
	    		 4.*dot15*dot18*dot5+3.*dot28*dot5+dot25*dot6+dot17*(dot29+dot30-2.*dot18*dot6)-
	    		 4.*dot15*dot18*dot7+3.*dot28*dot7-dot15*dot19*M))*M2/(8.*a1*a2*(dot8+a1*M2)*(dot9+a2*M2));
	    diag[2]=(3.*dot20+3.*dot21+4.*dot11*dot24-2.*dot10*dot25-
	    	     2.*(4.*dot12*dot17*dot18+dot15*(dot22+dot23)+dot13*(-2.*dot16*dot18+dot24)-
	    		 2.*dot12*dot25-dot16*(dot26+dot27)+dot14*(2.*dot15*dot18-3.*dot28)+2.*dot17*dot29+
	    		 2.*dot17*dot30+2.*dot16*dot18*dot4-dot24*dot4-2.*dot15*dot18*dot5+3.*dot28*dot5-
	    		 4.*dot17*dot18*dot6+2.*dot25*dot6-2.*dot15*dot18*dot7+
	    		 3.*dot28*dot7+dot15*dot19*M))*M2/(8.*a1*a2*(dot8+a1*M2)*(dot9+a2*M2));
	    diag[3]=M2*(dot1*dot15*(4.*dot14*dot18-4.*dot18*(dot5+dot7)+(1.-2.*a1)*dot19*M)+
			dot15*(-2.*dot23*dot3-4.*dot14*dot31+4.*dot18*dot3*dot5+2.*dot31*dot5+
			       4.*dot18*dot3*dot7+2.*dot31*dot7+4.*dot18*dot5*dot8-2.*dot22*(dot3+dot8-dot9)-
			       dot19*dot9*M+a1*dot19*(2.*dot3+dot8+dot9)*M+(-1.+2.*a1)*dot19*pow<3,1>(M)+
			       2.*dot2*(dot22+dot23-a2*dot19*M)+2.*(dot22-2.*a1*dot22+2.*a1*dot18*dot5)*M2)+
			2.*(-2.*dot14*dot29*dot4-2.*dot14*dot30*dot4+dot29*dot4*dot5+dot30*dot4*dot5+
			    (dot11-dot13)*(dot22+dot23)*dot6+2.*dot14*dot26*dot6+2.*dot14*dot27*dot6+
			    2.*dot13*dot18*dot5*dot6-dot26*dot5*dot6-dot27*dot5*dot6+dot29*dot4*dot7+
			    dot30*dot4*dot7+2.*dot13*dot18*dot6*dot7-dot26*dot6*dot7-dot27*dot6*dot7-
			    2.*dot16*dot18*dot4*dot8+dot24*dot4*dot8+2.*dot17*dot18*dot6*dot8-
			    dot25*dot6*dot8-dot24*dot4*dot9+dot25*dot6*dot9+(-(a2*dot11)+a1*dot13)*dot19*dot6*M+
			    dot12*dot4*(dot22+dot23-2.*dot18*(dot5+dot7)-a1*dot19*M)-
			    dot10*dot4*(dot22+dot23-a2*dot19*M)+
			    (-((2.*a1*dot16*dot18+dot24-2.*a1*dot24)*dot4)+
			     (2.*a1*dot17*dot18+dot25-2.*a1*dot25)*dot6)*M2))/(4.*a1*a2*dot1*(dot8+a1*M2)*(dot9+a2*M2));
	    diag[4]=M2*(-2.*dot10*dot22*dot4+2.*dot12*dot22*dot4-2.*dot10*dot23*dot4+2.*dot12*dot23*dot4-4.*dot14*dot29*dot4-
			4.*dot14*dot30*dot4-4.*dot12*dot18*dot4*dot5+2.*dot29*dot4*dot5+2.*dot30*dot4*dot5+
			2.*dot11*dot22*dot6-2.*dot13*dot22*dot6+2.*dot11*dot23*dot6-2.*dot13*dot23*dot6+
			4.*dot14*dot26*dot6+4.*dot14*dot27*dot6+4.*dot13*dot18*dot5*dot6-2.*dot26*dot5*dot6-
			2.*dot27*dot5*dot6-4.*dot12*dot18*dot4*dot7+2.*dot29*dot4*dot7+2.*dot30*dot4*dot7+
			4.*dot13*dot18*dot6*dot7-2.*dot26*dot6*dot7-2.*dot27*dot6*dot7-4.*dot16*dot18*dot4*dot8+
			2.*dot24*dot4*dot8+4.*dot17*dot18*dot6*dot8-2.*dot25*dot6*dot8-2.*dot24*dot4*dot9+2.*dot25*dot6*dot9+
			(2.*a2*dot10*dot19*dot4+2.*a2*dot12*dot19*dot4+(dot32+dot33)*dot4+(-2.*a2*(dot11+dot13)*dot19-dot34-dot35)*dot6)*M+
			dot1*dot15*(4.*dot14*dot18-4.*dot18*(dot5+dot7)+(1.-2.*a1)*dot19*M)+
			2.*(-(dot24*dot4)+dot25*dot6+2.*a1*(-(dot16*dot18*dot4)+dot24*dot4+dot17*dot18*dot6-dot25*dot6))*M2+
			dot15*(2.*(dot2*(dot22+dot23)-dot23*dot3-2.*dot14*dot31+2.*dot18*dot3*dot5+dot31*dot5+
				   2.*dot18*dot3*dot7+dot31*dot7+2.*dot18*dot5*dot8-dot22*(dot3+dot8-dot9))+
			       (dot36-a2*dot19*(2.*dot2+2.*dot3+dot8+dot9))*M-a2*dot19*pow<3,1>(M)+
			       2.*(dot22-2.*a1*dot22+2.*a1*dot18*dot5)*M2))/(4.*a1*dot1*(dot8+a1*M2)*(dot1-a2*(dot8+dot9+M2)));
	    diag[5]=pow<3,1>(M)*(-(dot15*dot36)-(2.*dot12*dot19+dot32+dot33)*dot4+
				 (2.*dot13*dot19+dot34+dot35)*dot6+dot15*dot19*(2.*dot3+dot8+a1*M2))/(4.*a1*dot1*sqr(dot8+a1*M2));
	    diag[6]=M2*(2.*dot15*dot2*dot22+2.*dot15*dot2*dot23-2.*dot15*dot22*dot3-2.*dot15*dot23*dot3-
			4.*dot14*dot15*dot31-2.*dot10*dot22*dot4+2.*dot12*dot22*dot4-2.*dot10*dot23*dot4+
			2.*dot12*dot23*dot4-4.*dot14*dot29*dot4-4.*dot14*dot30*dot4+4.*dot15*dot18*dot3*dot5+
			2.*dot15*dot31*dot5-4.*dot12*dot18*dot4*dot5+2.*dot29*dot4*dot5+2.*dot30*dot4*dot5+
			2.*dot11*dot22*dot6-2.*dot13*dot22*dot6+2.*dot11*dot23*dot6-2.*dot13*dot23*dot6+
			4.*dot14*dot26*dot6+4.*dot14*dot27*dot6+4.*dot13*dot18*dot5*dot6-2.*dot26*dot5*dot6-
			2.*dot27*dot5*dot6+4.*dot15*dot18*dot3*dot7+2.*dot15*dot31*dot7-4.*dot12*dot18*dot4*dot7+
			2.*dot29*dot4*dot7+2.*dot30*dot4*dot7+4.*dot13*dot18*dot6*dot7-2.*dot26*dot6*dot7-
			2.*dot27*dot6*dot7-2.*dot15*dot22*dot8-4.*dot16*dot18*dot4*dot8+2.*dot24*dot4*dot8+
			4.*dot15*dot18*dot5*dot8+4.*dot17*dot18*dot6*dot8-2.*dot25*dot6*dot8+2.*dot15*dot22*dot9-
			2.*dot24*dot4*dot9+2.*dot25*dot6*dot9+
			((dot32+dot33-2.*dot16*(dot37+dot38))*dot4-(dot34+dot35-2.*dot17*dot37)*dot6+
			 2.*dot17*dot38*dot6+2.*(dot39*dot4-dot40*dot6)*(dot5+dot7)+
			 dot15*(dot36+2.*dot38*dot5-2.*dot37*dot7)+
			 a1*dot19*(-2.*(dot10+dot12)*dot4+2.*(dot11+dot13)*dot6+dot15*(2.*dot2+2.*dot3+dot8+dot9)))*M+
			a1*dot15*dot19*pow<3,1>(M)+
			dot1*dot15*(4.*dot14*dot18-4.*dot18*(dot5+dot7)-(1.+2.*a1)*dot19*M)+
			2.*(-(dot24*dot4)+dot15*(dot22-2.*a1*dot22+2.*a1*dot18*dot5)+dot25*dot6+
			    2.*a1*(-(dot16*dot18*dot4)+dot24*dot4+dot17*dot18*dot6-dot25*dot6))*M2)
	      /(4.*a2*dot1*(dot9+a2*M2)*(dot1-a1*(dot8+dot9+M2)));
	    diag[7]=-0.25*pow<3,1>(M)*(-2.*dot1*dot15*dot19-2.*dot10*dot19*dot4+dot32*dot4+dot33*dot4-
				       2.*dot16*dot37*dot4-2.*dot16*dot38*dot4+2.*dot39*dot4*dot5+
				       2.*dot11*dot19*dot6-dot34*dot6-dot35*dot6+2.*dot17*dot37*dot6+
				       2.*dot17*dot38*dot6-2.*dot40*dot5*dot6+2.*dot39*dot4*dot7-
				       2.*dot40*dot6*dot7+dot15*(dot36+2.*dot38*dot5-2.*dot37*dot7+
								 dot19*(2.*dot2+dot9+M2-a1*M2)))/(a2*dot1*sqr(dot9+a2*M2));
	    diag[8]=-0.0625*M2*(-2.*dot2*dot20-2.*dot2*dot21+4.*(dot10-dot12)*(dot11+dot13)*dot22-4.*dot11*dot2*dot24-
				4.*dot13*dot2*dot24+8.*dot11*dot14*dot29+8.*dot13*dot14*dot29+4.*dot17*dot2*dot29-
				8.*dot11*dot16*dot18*dot3-8.*dot13*dot16*dot18*dot3-2.*dot20*dot3-2.*dot21*dot3+
				4.*dot15*dot22*dot3+4.*dot15*dot23*dot3+4.*dot11*dot24*dot3+4.*dot13*dot24*dot3-
				4.*dot16*dot26*dot3-4.*dot16*dot27*dot3+4.*dot17*dot2*dot30+4.*dot14*dot15*dot31+
				4.*dot11*dot16*dot31-4.*dot12*dot17*dot31-4.*dot10*dot22*dot4+4.*dot12*dot22*dot4+
				4.*dot12*dot23*dot4-8.*dot14*dot29*dot4+8.*dot16*dot18*dot3*dot4-4.*dot24*dot3*dot4-
				4.*dot14*dot30*dot4+2.*dot10*dot43+2.*dot12*dot43-2.*dot11*dot44-2.*dot13*dot44+
				8.*dot11*dot12*dot18*dot5+8.*dot12*dot13*dot18*dot5-4.*dot12*dot27*dot5+4.*dot2*dot28*dot5-
				4.*dot11*dot29*dot5-4.*dot13*dot29*dot5+4.*dot13*dot30*dot5-4.*dot15*dot31*dot5-
				8.*dot12*dot18*dot4*dot5+4.*dot29*dot4*dot5+2.*dot42*dot5-8.*dot17*dot18*dot2*dot6+
				4.*dot13*dot22*dot6+4.*dot2*dot25*dot6-4.*dot14*dot26*dot6+8.*dot14*dot18*dot4*dot6-
				8.*dot13*dot18*dot5*dot6+4.*dot12*dot26*dot7+4.*dot2*dot28*dot7-4.*dot11*dot29*dot7-
				8.*dot13*dot29*dot7-4.*dot15*dot31*dot7-8.*dot12*dot18*dot4*dot7+4.*dot29*dot4*dot7+
				2.*dot42*dot7-2.*dot20*dot8+4.*dot15*dot22*dot8-4.*dot16*dot26*dot8-2.*dot20*dot9+
				4.*dot17*dot29*dot9-4.*dot24*dot4*dot9+4.*dot28*dot5*dot9+
				(-2.*dot15*dot36+2.*dot32*dot4+2.*dot33*dot4-4.*dot16*dot37*dot4+2.*dot16*dot41+
				 dot45+4.*dot39*dot4*dot5-2.*dot47*dot5-2.*dot11*(dot32-2.*dot16*dot37+2.*dot39*dot5)-
				 2.*dot13*(dot32-2.*dot16*dot37+2.*dot39*dot5)-4.*dot19*dot4*dot6+
				 4.*dot40*dot5*dot6+4.*dot39*dot4*dot7)*M-
				4.*a12*(dot11+dot13-dot4)*M*(dot19*(dot10+dot12-dot6)+2.*(dot16*dot18-dot24)*M)+
				2.*dot1*(4.*dot12*dot17*dot18+2.*dot13*dot24-2.*dot12*dot25-2.*dot14*dot28-2.*dot17*dot39*M+dot46*M)-
				2.*(dot20-2.*dot17*dot29+2.*dot24*dot4-2.*dot28*dot5)*M2+
				2.*a1*(4.*dot13*dot16*dot18*dot2-4.*dot12*dot17*dot18*dot2+2.*dot2*dot20+2.*dot2*dot21-
				       2.*dot10*dot11*dot22+2.*dot11*dot12*dot22-2.*dot10*dot13*dot22+2.*dot12*dot13*dot22-
				       2.*dot15*dot2*dot22-2.*dot10*dot11*dot23+2.*dot11*dot12*dot23-2.*dot10*dot13*dot23+
				       2.*dot12*dot13*dot23-2.*dot15*dot2*dot23+2.*dot11*dot2*dot24-2.*dot13*dot2*dot24-
				       2.*dot10*dot2*dot25+2.*dot12*dot2*dot25+2.*dot16*dot2*dot26+2.*dot16*dot2*dot27-
				       2.*dot17*dot2*dot29+4.*dot13*dot16*dot18*dot3-4.*dot12*dot17*dot18*dot3+2.*dot20*dot3+
				       2.*dot21*dot3-2.*dot15*dot22*dot3-2.*dot15*dot23*dot3+2.*dot11*dot24*dot3-
				       2.*dot13*dot24*dot3-2.*dot10*dot25*dot3+2.*dot12*dot25*dot3+2.*dot16*dot26*dot3+
				       2.*dot16*dot27*dot3-2.*dot17*dot29*dot3-2.*dot17*dot2*dot30-2.*dot17*dot3*dot30-
				       4.*dot16*dot18*dot2*dot4+2.*dot10*dot22*dot4-2.*dot12*dot22*dot4+2.*dot10*dot23*dot4-
				       2.*dot12*dot23*dot4+2.*dot2*dot24*dot4-4.*dot16*dot18*dot3*dot4+2.*dot24*dot3*dot4-
				       4.*dot11*dot12*dot18*dot5-4.*dot12*dot13*dot18*dot5+4.*dot15*dot18*dot2*dot5-
				       4.*dot2*dot28*dot5+2.*dot11*dot29*dot5+2.*dot13*dot29*dot5+4.*dot15*dot18*dot3*dot5-
				       4.*dot28*dot3*dot5+2.*dot11*dot30*dot5+2.*dot13*dot30*dot5+4.*dot12*dot18*dot4*dot5-
				       2.*dot29*dot4*dot5-2.*dot30*dot4*dot5+4.*dot17*dot18*dot2*dot6-2.*dot2*dot25*dot6+
				       4.*dot17*dot18*dot3*dot6-2.*dot25*dot3*dot6-4.*dot11*dot12*dot18*dot7-
				       4.*dot12*dot13*dot18*dot7+4.*dot15*dot18*dot2*dot7-4.*dot2*dot28*dot7+
				       2.*dot11*dot29*dot7+2.*dot13*dot29*dot7+4.*dot15*dot18*dot3*dot7-4.*dot28*dot3*dot7+
				       2.*dot11*dot30*dot7+2.*dot13*dot30*dot7+4.*dot12*dot18*dot4*dot7-2.*dot29*dot4*dot7-
				       2.*dot30*dot4*dot7-4.*dot11*dot16*dot18*dot8-4.*dot12*dot17*dot18*dot8+2.*dot20*dot8+
				       2.*dot21*dot8-2.*dot15*dot22*dot8-2.*dot15*dot23*dot8+4.*dot11*dot24*dot8-
				       2.*dot10*dot25*dot8+2.*dot12*dot25*dot8+2.*dot16*dot26*dot8+2.*dot16*dot27*dot8-
				       2.*dot17*dot29*dot8-2.*dot17*dot30*dot8+4.*dot15*dot18*dot5*dot8-4.*dot28*dot5*dot8+
				       4.*dot17*dot18*dot6*dot8-2.*dot25*dot6*dot8+4.*dot15*dot18*dot7*dot8-4.*dot28*dot7*dot8+
				       4.*dot13*dot16*dot18*dot9-4.*dot12*dot17*dot18*dot9+2.*dot20*dot9+2.*dot21*dot9-
				       2.*dot15*dot22*dot9-2.*dot15*dot23*dot9-4.*dot13*dot24*dot9-2.*dot10*dot25*dot9+
				       2.*dot12*dot25*dot9+2.*dot16*dot26*dot9+2.*dot16*dot27*dot9-2.*dot17*dot29*dot9-
				       2.*dot17*dot30*dot9-4.*dot16*dot18*dot4*dot9+4.*dot24*dot4*dot9+4.*dot15*dot18*dot5*dot9-
				       4.*dot28*dot5*dot9+4.*dot17*dot18*dot6*dot9-2.*dot25*dot6*dot9+
				       4.*dot15*dot18*dot7*dot9-4.*dot28*dot7*dot9+
				       (2.*dot12*dot13*dot19+2.*dot15*dot19*dot2+2.*dot15*dot19*dot3+dot13*dot32+2.*dot13*dot33-
					dot12*dot35-2.*dot13*dot16*dot37-2.*dot13*dot16*dot38+2.*dot12*dot17*dot38+
					2.*dot17*dot2*dot39+2.*dot17*dot3*dot39-2.*dot12*dot19*dot4-dot32*dot4-
					2.*dot33*dot4+2.*dot16*dot37*dot4+2.*dot16*dot38*dot4-2.*dot16*dot2*dot40-
					2.*dot16*dot3*dot40-2.*dot2*dot46-2.*dot3*dot46+2.*dot13*dot39*dot5-
					2.*dot39*dot4*dot5-4.*dot13*dot19*dot6+dot35*dot6-2.*dot17*dot38*dot6+
					4.*dot19*dot4*dot6+2.*dot13*dot39*dot7-2.*dot39*dot4*dot7-2.*dot12*dot40*dot7+
					2.*dot40*dot6*dot7+
					dot10*(-dot35+2.*dot17*dot38+2.*dot19*(dot11+dot13-dot4)-2.*dot40*dot7)+
					dot11*(dot32+2.*dot19*(dot12-2.*dot6)+2.*(dot33-dot16*(dot37+dot38)+dot39*(dot5+dot7)))+
					2.*dot15*dot19*dot8+2.*dot17*dot39*dot8-2.*dot16*dot40*dot8-2.*dot46*dot8+
					2.*(dot15*dot19+dot17*dot39-dot16*dot40-dot46)*dot9)*M+
				       2.*(dot15*dot19+dot17*dot39-dot16*dot40-dot46)*pow<3,1>(M)-
				       2.*dot1*(-2.*dot12*dot17*dot18+dot20+dot21-dot15*(dot22+dot23)+
						dot13*(2.*dot16*dot18-dot24)+dot11*dot24-dot10*dot25+dot12*dot25+
						dot16*(dot26+dot27)+2.*dot14*(-(dot15*dot18)+dot28)-dot17*dot29-
						dot17*dot30-2.*dot16*dot18*dot4+dot24*dot4+2.*dot15*dot18*dot5-
						2.*dot28*dot5+2.*dot17*dot18*dot6-dot25*dot6+2.*dot15*dot18*dot7-
						2.*dot28*dot7+(dot15*dot19+dot17*dot39-dot16*dot40-dot46)*M)+
				       2.*(2.*dot13*dot16*dot18-2.*dot12*dot17*dot18+dot20-2.*dot13*dot24-
					   dot10*dot25+dot12*dot25-2.*dot17*dot29-2.*dot16*dot18*dot4+
					   3.*dot24*dot4+2.*dot15*dot18*dot5-3.*dot28*dot5+dot17*dot18*dot6+
					   dot15*dot18*dot7-dot28*dot7)*M2-
				       4.*dot14*((dot11+dot13-dot4)*(dot29+dot30-dot18*dot6)+
						 dot15*dot18*(dot2+dot3+dot8+dot9+M2)-
						 dot28*(dot2+dot3+dot8+dot9+M2))))
	      /(a1*a2*(dot1-dot2-dot3)*(-dot1+dot2+dot3+dot8+dot9+M2)*(-dot2+a1*(-dot1+dot2+dot3+dot8+a1*M2)));
	    diag[9]=-0.125*M2*(4.*dot2*dot20+2.*dot2*dot21+4.*dot11*(-dot10+dot12)*dot22+4.*dot11*dot2*dot24+
			       4.*dot16*dot2*dot26-8.*dot11*dot14*dot29-8.*dot17*dot2*dot29+8.*dot11*dot16*dot18*dot3-
			       2.*dot21*dot3-4.*dot11*dot24*dot3-4.*dot16*dot27*dot3-4.*dot17*dot2*dot30+
			       4.*dot17*dot3*dot30-4.*dot10*dot17*dot31+4.*dot12*dot17*dot31-4.*dot14*dot42+
			       2.*dot10*dot43-2.*dot12*dot43-8.*dot11*dot12*dot18*dot5+4.*dot12*dot27*dot5-
			       4.*dot2*dot28*dot5+4.*dot11*dot29*dot5+2.*dot42*dot5+8.*dot17*dot18*dot2*dot6+
			       4.*dot11*dot22*dot6-4.*dot2*dot25*dot6+8.*dot14*dot26*dot6-
			       8.*dot17*dot18*dot3*dot6+4.*dot25*dot3*dot6+4.*dot17*dot31*dot6-
			       2.*dot43*dot6-4.*dot26*dot5*dot6+4.*dot10*dot26*dot7-4.*dot12*dot26*dot7-
			       4.*dot2*dot28*dot7+4.*dot28*dot3*dot7-4.*dot26*dot6*dot7-
			       2.*dot1*(2.*dot20+dot21+2.*dot11*dot24+2.*dot16*dot26-2.*dot25*dot6-
					2.*dot17*(2.*dot29+dot30-2.*dot18*dot6)-2.*dot28*(dot5+dot7))+
			       4.*dot20*dot9+4.*dot16*dot26*dot9-8.*dot17*dot29*dot9-4.*dot28*dot5*dot9+
			       (dot45+2.*(dot16*dot41+dot47*dot5+dot11*(dot32-2.*dot16*dot37+2.*dot39*dot5)+
					  2.*dot17*dot37*dot6-2.*dot40*dot5*dot6)-
				2.*a1*(dot10+dot12-dot6)*(2.*dot11*dot19-dot35+2.*dot17*dot38-2.*dot40*dot7)-
				2.*(dot17*dot48+dot34*dot6+dot49*dot7))*M+
			       4.*a2*(dot20+dot16*dot26-2.*dot17*dot29-dot28*dot5)*M2)
	      /(a1*a2*(dot1-dot2-dot3)*(-dot9-a2*M2)*(dot1-dot2-dot9-a2*M2));
	    diag[10]=M2*(-2.*dot2*dot21+4.*(-dot10+dot12)*dot13*dot22+4.*dot15*dot2*dot23+4.*dot13*dot2*dot24-
			 8.*dot13*dot14*dot29+8.*dot13*dot16*dot18*dot3+4.*dot20*dot3+2.*dot21*dot3-
			 8.*dot15*dot22*dot3-4.*dot15*dot23*dot3-4.*dot13*dot24*dot3+4.*dot16*dot26*dot3+
			 4.*dot16*dot27*dot3-8.*dot14*dot15*dot31-4.*dot10*dot23*dot4+4.*dot12*dot23*dot4-
			 8.*dot16*dot18*dot3*dot4+8.*dot24*dot3*dot4-8.*dot14*dot30*dot4+4.*dot14*dot42-
			 2.*dot10*dot43+2.*dot12*dot43-8.*dot12*dot13*dot18*dot5-4.*dot12*dot27*dot5+
			 4.*dot13*dot29*dot5+8.*dot15*dot18*dot3*dot5-4.*dot28*dot3*dot5+4.*dot15*dot31*dot5+
			 4.*dot30*dot4*dot5-2.*dot42*dot5-4.*dot13*dot22*dot6-4.*dot23*dot4*dot6-
			 2.*dot43*dot6+8.*dot13*dot18*dot5*dot6+4.*dot27*dot5*dot6+8.*dot13*dot29*dot7+
			 8.*dot15*dot31*dot7+8.*dot30*dot4*dot7-4.*dot42*dot7+4.*dot20*dot8-8.*dot15*dot22*dot8+
			 4.*dot16*dot26*dot8-8.*dot16*dot18*dot4*dot8+
			 8.*dot24*dot4*dot8+8.*dot15*dot18*dot5*dot8-4.*dot28*dot5*dot8-
			 (-2.*dot15*dot36-2.*dot33*dot4+4.*dot16*dot38*dot4+2.*dot16*dot41+dot45-
			  4.*dot15*dot38*dot5+2.*dot47*dot5-2.*dot13*(dot32-2.*dot16*dot37+2.*dot39*dot5)+
			  2.*a1*(2.*dot13*dot19+dot35)*(dot10+dot12-dot6))*M+
			 2.*dot1*(-2.*dot20-dot21+4.*dot15*dot22+2.*dot15*dot23+2.*dot13*(-2.*dot16*dot18+dot24)-
				  2.*dot16*(dot26+dot27)+dot14*(8.*dot15*dot18-4.*dot28)+
				  4.*(dot16*dot18-dot24)*dot4-8.*dot15*dot18*dot5+4.*dot28*dot5+
				  2.*dot25*(dot10-dot12+dot6)-8.*dot15*dot18*dot7+4.*dot28*dot7-2.*dot15*dot19*M+
				  2.*dot16*dot40*M+dot46*M)+
			 4.*a1*(dot20-2.*dot15*dot22+dot16*dot26-2.*dot16*dot18*dot4+
				2.*dot24*dot4+2.*dot15*dot18*dot5-dot28*dot5)*M2)
	      /(8.*a1*(dot1-dot2-dot3)*(-dot1+dot3+dot8+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot8+a1*M2)));
	    diag[11]=-0.0625*M2*(-4.*dot16*dot2*dot26-4.*dot16*dot2*dot27+4.*dot17*dot2*dot29-2.*(dot20+dot21)*(dot2+dot3)+
				 4.*dot17*dot2*dot30-4.*dot11*dot16*dot31+4.*dot10*dot17*dot31+
				 4.*dot15*((dot22+dot23)*dot3+dot14*dot31)-4.*(dot12*dot23+dot24*dot3)*dot4+
				 4.*dot14*dot30*dot4-2.*dot10*dot43-2.*dot12*dot43+2.*dot11*dot44+2.*dot13*dot44+
				 4.*dot10*dot26*dot5-4.*dot12*dot26*dot5+4.*dot10*dot27*dot5-4.*dot11*dot29*dot5+
				 4.*dot13*dot29*dot5-8.*dot15*dot18*dot3*dot5+4.*dot28*dot3*dot5-4.*dot11*dot30*dot5+
				 8.*dot12*dot18*dot4*dot5-2.*dot42*dot5-4.*dot13*dot22*dot6+4.*dot14*dot26*dot6+
				 4.*dot25*dot3*dot6-8.*dot14*dot18*dot4*dot6-4.*dot12*dot26*dot7+4.*dot13*dot29*dot7-
				 8.*dot15*dot18*dot3*dot7+4.*dot28*dot3*dot7+8.*dot12*dot18*dot4*dot7-2.*dot42*dot7+
				 4.*dot1*(2.*dot11*dot16*dot18-2.*dot10*dot17*dot18-dot11*dot24+dot10*dot25+
					  dot28*(-dot14+dot5+dot7))+
				 4.*dot16*dot26*dot8-4.*dot17*dot29*dot8-4.*dot28*dot5*dot8-4.*dot15*dot22*dot9+
				 4.*dot24*dot4*dot9+2.*dot20*(dot8+dot9)+
				 2.*dot1*(-2.*a2*dot15*dot19-2.*dot17*dot39+2.*dot16*dot40+2.*dot46-a1*dot46)*M+
				 (2.*dot32*dot4-2.*dot16*dot41-dot45+2.*dot17*dot48+4.*dot15*dot37*(dot5+dot7)-
				  2.*(dot47*dot5+2.*dot49*dot5+dot34*dot6-2.*dot19*dot4*dot6+dot49*dot7+
				      2.*dot39*dot4*(dot5+dot7)))*M+
				 4.*a12*(dot10+dot12-dot6)*M*(dot19*(dot11+dot13-dot4)+2.*(dot17*dot18-dot25)*M)+
				 2.*(dot20-2.*dot15*dot22+2.*dot24*dot4)*M2+
				 2.*a1*(-2.*dot12*dot13*dot22-2.*dot12*dot13*dot23+4.*dot12*dot14*dot26+
					4.*dot12*dot14*dot27-4.*dot12*dot14*dot18*dot4+4.*dot12*dot13*dot18*dot5-
					2.*dot12*dot26*dot5-2.*dot12*dot27*dot5+2.*dot11*(dot22+dot23)*(dot12-dot6)+
					2.*dot13*dot22*dot6+2.*dot13*dot23*dot6-4.*dot14*dot26*dot6-4.*dot14*dot27*dot6+
					4.*dot14*dot18*dot4*dot6-4.*dot13*dot18*dot5*dot6+2.*dot26*dot5*dot6+
					2.*dot27*dot5*dot6+4.*dot12*dot13*dot18*dot7-2.*dot12*dot26*dot7-2.*dot12*dot27*dot7-
					4.*dot13*dot18*dot6*dot7+2.*dot26*dot6*dot7+2.*dot27*dot6*dot7+
					4.*dot12*dot17*dot18*dot8-2.*dot12*dot25*dot8-4.*dot17*dot18*dot6*dot8+
					2.*dot25*dot6*dot8+2.*dot12*dot25*dot9-2.*dot25*dot6*dot9+dot11*dot32*M+
					(-2.*dot15*dot19*(dot2+dot3)+dot13*dot32+4.*dot12*dot19*dot4-dot32*dot4+
					 (dot2+dot3)*dot46-(2.*dot34+dot35-2.*dot17*(dot37+dot38))*(dot12-dot6)-
					 4.*dot19*dot4*dot6+2.*dot40*dot5*dot6+2.*dot40*dot6*dot7-2.*dot12*dot40*(dot5+dot7))*M-
					2.*(dot20-2.*dot15*dot22-dot12*dot25-dot16*dot18*dot4+2.*dot24*dot4+dot15*dot18*dot5+dot25*dot6)*M2+
					dot10*(2.*dot11*(dot22+dot23)-
					       2.*(-2.*dot14*(dot26+dot27-dot18*dot4)+(dot26+dot27)*(dot5+dot7)+
						   dot13*(dot22+dot23-2.*dot18*(dot5+dot7))-
						   2.*dot17*dot18*dot8+dot25*dot8-dot25*dot9)-
					       (2.*dot34+dot35-2.*dot17*(dot37+dot38)-4.*dot19*dot4+2.*dot40*(dot5+dot7))*M+
					       2.*dot25*M2)))/(a1*a2*(dot1-dot2-dot3)*(-dot9-a2*M2)*(-dot1+a1*(dot8+dot9+M2)));
	    diag[12]=M2*(4.*a12*(dot10+dot12-dot6)*M*(dot19*(dot11+dot13-dot4)+2.*(dot17*dot18-dot25)*M)+
			 2.*(-(dot2*dot20)+2.*dot11*(-dot10+dot12)*dot22-4.*dot11*dot14*dot29+4.*dot11*dot16*dot18*dot3+
			     dot20*dot3-2.*dot11*dot24*dot3+2.*dot16*dot26*dot3-2.*dot17*dot29*dot3-2.*dot14*dot15*dot31-
			     4.*dot11*dot16*dot31+2.*dot13*dot16*dot31+2.*dot10*dot17*dot31+2.*dot10*dot22*dot4-
			     2.*dot12*dot22*dot4-2.*dot10*dot23*dot4+4.*dot14*dot29*dot4-4.*dot16*dot18*dot3*dot4+
			     2.*dot24*dot3*dot4-2.*dot14*dot30*dot4+2.*dot16*dot31*dot4+2.*dot14*dot42-
			     2.*dot10*dot43+3.*dot11*dot44-dot13*dot44-2.*dot4*dot44-4.*dot11*dot12*dot18*dot5+
			     4.*dot10*dot27*dot5+2.*dot11*dot29*dot5-2.*dot28*dot3*dot5-4.*dot11*dot30*dot5+
			     2.*dot13*dot30*dot5+4.*dot12*dot18*dot4*dot5-2.*dot29*dot4*dot5+2.*dot30*dot4*dot5-
			     2.*dot11*dot22*dot6+2.*dot13*dot22*dot6-2.*dot14*dot26*dot6+4.*dot14*dot18*dot4*dot6+
			     2.*dot22*dot4*dot6+2.*dot23*dot4*dot6+dot43*dot6+4.*dot11*dot18*dot5*dot6-
			     4.*dot13*dot18*dot5*dot6-2.*dot27*dot5*dot6-4.*dot18*dot4*dot5*dot6-2.*dot10*dot26*dot7+
			     4.*dot11*dot29*dot7-2.*dot13*dot29*dot7+4.*dot10*dot18*dot4*dot7-4.*dot29*dot4*dot7-
			     dot42*dot7+2.*dot26*dot6*dot7-4.*dot18*dot4*dot6*dot7+
			     2.*dot2*(-dot21+dot17*dot30+dot24*(dot11-2.*dot4)-dot16*(dot26+2.*dot27-2.*dot18*dot4)+
				      dot25*dot6+dot28*(dot5+dot7)+dot15*(dot22+dot23-2.*dot18*(dot5+dot7)))+
			     dot1*(dot21+4.*dot11*(dot16*dot18-dot24)+2.*dot10*(-2.*dot17*dot18+dot25)-
				   2.*(-(dot16*dot27)+dot14*dot28-dot17*dot29+2.*dot16*dot18*dot4-2.*dot24*dot4+
				       dot25*dot6+dot15*(dot22+dot23-2.*dot18*(dot5+dot7))))+dot20*dot8+
			     2.*dot16*dot26*dot8-2.*dot17*dot29*dot8-2.*dot28*dot5*dot8-dot20*dot9+
			     2.*dot15*dot22*dot9-2.*dot16*dot26*dot9+4.*dot16*dot18*dot4*dot9-2.*dot24*dot4*dot9-
			     4.*dot15*dot18*dot5*dot9+2.*dot28*dot5*dot9+
			     (dot11-dot4)*(dot32-2.*dot16*dot37+2.*dot39*dot5)*M-
			     (dot20-2.*dot15*dot22+2.*dot16*dot26-4.*dot16*dot18*dot4+
			      2.*dot24*dot4+4.*dot15*dot18*dot5-2.*dot28*dot5)*M2)+
			 a1*(4.*dot13*dot16*dot18*dot2+dot2*dot20+dot2*dot21+4.*dot11*(dot10-dot12)*dot22+4.*dot10*dot13*dot22-
			     4.*dot12*dot13*dot22-2.*dot15*dot2*dot22+4.*dot10*dot11*dot23+4.*dot11*dot12*dot23-
			     4.*dot10*dot13*dot23-4.*dot12*dot13*dot23-2.*dot15*dot2*dot23-4.*dot11*dot2*dot24-
			     6.*dot13*dot2*dot24+2.*dot10*dot2*dot25+4.*dot12*dot2*dot25+2.*dot16*dot2*dot26+
			     2.*dot16*dot2*dot27-8.*dot11*dot16*dot18*dot3-4.*dot13*dot16*dot18*dot3+8.*dot10*dot17*dot18*dot3+
			     8.*dot12*dot17*dot18*dot3+dot20*dot3+dot21*dot3-2.*dot15*dot22*dot3-2.*dot15*dot23*dot3+
			     4.*dot11*dot24*dot3+2.*dot13*dot24*dot3-6.*dot10*dot25*dot3-4.*dot12*dot25*dot3+
			     2.*dot16*dot26*dot3+2.*dot16*dot27*dot3+4.*dot11*dot16*dot31+4.*dot13*dot16*dot31-
			     4.*dot10*dot17*dot31-4.*dot12*dot17*dot31-4.*dot16*dot18*dot2*dot4-4.*dot10*dot22*dot4+
			     4.*dot12*dot22*dot4+4.*dot10*dot23*dot4+4.*dot12*dot23*dot4+6.*dot2*dot24*dot4+
			     4.*dot16*dot18*dot3*dot4-2.*dot24*dot3*dot4-4.*dot16*dot31*dot4+4.*dot10*dot43+
			     4.*dot12*dot43-4.*dot11*dot44-4.*dot13*dot44+4.*dot4*dot44+8.*dot11*dot12*dot18*dot5+
			     8.*dot12*dot13*dot18*dot5+4.*dot15*dot18*dot2*dot5-8.*dot10*dot27*dot5-8.*dot12*dot27*dot5-
			     2.*dot2*dot28*dot5-4.*dot11*dot29*dot5-4.*dot13*dot29*dot5+4.*dot15*dot18*dot3*dot5-
			     2.*dot28*dot3*dot5+4.*dot11*dot30*dot5+4.*dot13*dot30*dot5-8.*dot12*dot18*dot4*dot5+
			     4.*dot29*dot4*dot5-4.*dot30*dot4*dot5+
			     2.*dot14*(-2.*dot15*dot18*(dot2+dot3)+dot28*(dot2+dot3)+4.*dot29*(dot11+dot13-dot4)+4.*dot27*(dot10+dot12-dot6))+
			     4.*dot11*dot22*dot6+4.*dot13*dot22*dot6-4.*dot11*dot23*dot6+4.*dot13*dot23*dot6-4.*dot2*dot25*dot6-
			     8.*dot17*dot18*dot3*dot6+4.*dot25*dot3*dot6+4.*dot17*dot31*dot6-4.*dot22*dot4*dot6-
			     4.*dot23*dot4*dot6-4.*dot43*dot6-8.*dot11*dot18*dot5*dot6-8.*dot13*dot18*dot5*dot6+
			     8.*dot27*dot5*dot6+8.*dot18*dot4*dot5*dot6+8.*dot10*dot13*dot18*dot7+8.*dot12*dot13*dot18*dot7+
			     4.*dot15*dot18*dot2*dot7+4.*dot10*dot26*dot7+4.*dot12*dot26*dot7-4.*dot10*dot27*dot7-
			     4.*dot12*dot27*dot7-2.*dot2*dot28*dot7-8.*dot11*dot29*dot7-8.*dot13*dot29*dot7+
			     4.*dot15*dot18*dot3*dot7-2.*dot28*dot3*dot7-8.*dot10*dot18*dot4*dot7-8.*dot12*dot18*dot4*dot7+
			     8.*dot29*dot4*dot7-8.*dot13*dot18*dot6*dot7-4.*dot26*dot6*dot7+4.*dot27*dot6*dot7+
			     8.*dot18*dot4*dot6*dot7-
			     dot1*(dot20+dot21-4.*dot11*dot24+2.*dot25*(dot10+2.*dot12-2.*dot6)+
				   2.*(2.*dot13*dot16*dot18-dot15*(dot22+dot23)-3.*dot13*dot24+dot16*dot26+
				       dot16*dot27+dot14*(-2.*dot15*dot18+dot28)-2.*dot16*dot18*dot4+
				       3.*dot24*dot4+2.*dot15*dot18*dot5-dot28*dot5+2.*dot15*dot18*dot7-dot28*dot7))+
			     8.*dot10*dot17*dot18*dot8+8.*dot12*dot17*dot18*dot8-4.*dot10*dot25*dot8-
			     4.*dot12*dot25*dot8-8.*dot17*dot18*dot6*dot8+4.*dot25*dot6*dot8+4.*dot10*dot25*dot9+
			     4.*dot12*dot25*dot9-4.*dot25*dot6*dot9-
			     2.*(2.*dot10*dot19*(dot11-dot4)-2.*dot12*dot19*dot4+
				 (dot13-dot4)*(dot32-2.*dot16*dot37+2.*dot39*dot5)+
				 2.*dot19*dot4*dot6+dot11*(2.*dot12*dot19+dot32-2.*dot16*dot37+2.*dot39*dot5-2.*dot19*dot6))*M+
			     4.*(dot20-dot15*dot22-dot17*dot29+dot24*dot4+2.*dot16*(dot26-dot18*dot4)+
				 2.*dot15*dot18*dot5-2.*dot28*dot5+dot25*(dot10+dot12-dot6))*M2))
	      /(16.*a1*a2*(dot1-dot2-dot3)*(-dot9-a2*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot8+a1*M2)));
	    diag[13]=-0.0625*M2*(-2.*dot2*dot20-2.*dot2*dot21+4.*dot15*dot2*dot22+4.*(dot10+dot12)*(dot11-dot13)*dot23+4.*dot15*dot2*dot23+
				 4.*dot10*dot2*dot25+4.*dot12*dot2*dot25-4.*dot16*dot2*dot26+8.*dot10*dot14*dot27+8.*dot12*dot14*dot27-
				 4.*dot16*dot2*dot27+8.*dot10*dot17*dot18*dot3+8.*dot12*dot17*dot18*dot3-2.*dot20*dot3-2.*dot21*dot3-
				 4.*dot10*dot25*dot3-4.*dot12*dot25*dot3+4.*dot17*dot29*dot3+4.*dot17*dot3*dot30-4.*dot14*dot15*dot31+
				 4.*dot13*dot16*dot31-4.*dot10*dot17*dot31+8.*dot16*dot18*dot2*dot4+4.*dot12*dot23*dot4-4.*dot2*dot24*dot4-
				 4.*dot14*dot30*dot4+2.*dot10*dot43+2.*dot12*dot43-2.*dot11*dot44-2.*dot13*dot44-8.*dot15*dot18*dot2*dot5-
				 4.*dot10*dot27*dot5-8.*dot12*dot27*dot5+4.*dot2*dot28*dot5+4.*dot13*dot30*dot5+2.*dot42*dot5+
				 4.*dot13*dot22*dot6-4.*dot11*dot23*dot6+4.*dot13*dot23*dot6-4.*dot14*dot26*dot6-8.*dot14*dot27*dot6-
				 8.*dot17*dot18*dot3*dot6+4.*dot25*dot3*dot6+8.*dot14*dot18*dot4*dot6-8.*dot13*dot18*dot5*dot6+
				 4.*dot27*dot5*dot6+8.*dot10*dot13*dot18*dot7+8.*dot12*dot13*dot18*dot7-8.*dot15*dot18*dot2*dot7+
				 4.*dot12*dot26*dot7-4.*dot10*dot27*dot7-4.*dot12*dot27*dot7+4.*dot2*dot28*dot7-4.*dot13*dot29*dot7-
				 8.*dot12*dot18*dot4*dot7+2.*dot42*dot7-8.*dot13*dot18*dot6*dot7+4.*dot27*dot6*dot7+
				 8.*dot10*dot17*dot18*dot8+8.*dot12*dot17*dot18*dot8-2.*dot20*dot8-4.*dot10*dot25*dot8-
				 4.*dot12*dot25*dot8+4.*dot17*dot29*dot8-8.*dot17*dot18*dot6*dot8+4.*dot25*dot6*dot8-
				 2.*dot20*dot9+4.*dot15*dot22*dot9+4.*dot10*dot25*dot9+4.*dot12*dot25*dot9-4.*dot16*dot26*dot9+
				 8.*dot16*dot18*dot4*dot9-4.*dot24*dot4*dot9-8.*dot15*dot18*dot5*dot9+4.*dot28*dot5*dot9-
				 4.*dot25*dot6*dot9+
				 (-2.*dot10*dot35-2.*dot12*dot35+4.*dot10*dot17*dot38+4.*dot12*dot17*dot38+dot45-2.*dot17*dot48+
				  2.*dot34*dot6+2.*dot35*dot6-4.*dot17*dot38*dot6-4.*dot19*dot4*dot6+4.*dot40*dot5*dot6+
				  2.*(-2.*dot15*dot37+2.*dot39*dot4+dot49-2.*dot40*(dot10+dot12-dot6))*dot7)*M-
				 4.*a12*(dot10+dot12-dot6)*M*(dot19*(dot11+dot13-dot4)+2.*(dot17*dot18-dot25)*M)-
				 4.*dot1*(2.*dot12*dot17*dot18-dot20-dot21+dot10*dot25-dot14*dot28+dot17*(dot29+dot30)-
					  dot24*(dot11+dot4)-dot16*(dot26+dot27-2.*dot18*dot4)-2.*dot17*dot18*dot6+
					  dot25*dot6+2.*dot28*(dot5+dot7)-dot17*dot39*M+dot16*dot40*M+dot46*M+
					  dot15*(dot22+dot23-2.*dot18*(dot5+dot7)-dot19*M))-
				 2.*(dot20-2.*(-(dot16*dot26)+2.*dot16*dot18*dot4-dot24*dot4+dot28*dot5+dot15*(dot22-2.*dot18*dot5)+dot25*(dot10+dot12-dot6)))*M2+
				 2.*a1*(4.*dot13*dot16*dot18*dot2-4.*dot12*dot17*dot18*dot2+2.*dot2*dot20+2.*dot2*dot21-2.*dot10*dot11*dot22-2.*dot11*dot12*dot22+
					2.*dot10*dot13*dot22+2.*dot12*dot13*dot22-2.*dot15*dot2*dot22-2.*dot10*dot11*dot23-2.*dot11*dot12*dot23+2.*dot10*dot13*dot23+
					2.*dot12*dot13*dot23-2.*dot15*dot2*dot23+2.*dot11*dot2*dot24-2.*dot13*dot2*dot24-2.*dot10*dot2*dot25+2.*dot12*dot2*dot25+
					2.*dot16*dot2*dot26+2.*dot16*dot2*dot27-2.*dot17*dot2*dot29+4.*dot13*dot16*dot18*dot3-4.*dot12*dot17*dot18*dot3+
					2.*dot20*dot3+2.*dot21*dot3-2.*dot15*dot22*dot3-2.*dot15*dot23*dot3+2.*dot11*dot24*dot3-2.*dot13*dot24*dot3-
					2.*dot10*dot25*dot3+2.*dot12*dot25*dot3+2.*dot16*dot26*dot3+2.*dot16*dot27*dot3-2.*dot17*dot29*dot3-2.*dot17*dot2*dot30-
					2.*dot17*dot3*dot30-4.*dot16*dot18*dot2*dot4+2.*dot2*dot24*dot4-4.*dot16*dot18*dot3*dot4+2.*dot24*dot3*dot4-
					4.*dot10*dot13*dot18*dot5-4.*dot12*dot13*dot18*dot5+4.*dot15*dot18*dot2*dot5+2.*dot10*dot26*dot5+2.*dot12*dot26*dot5+
					2.*dot10*dot27*dot5+2.*dot12*dot27*dot5-4.*dot2*dot28*dot5+4.*dot15*dot18*dot3*dot5-4.*dot28*dot3*dot5+
					4.*dot14*(-(dot15*dot18*(dot2+dot3))+dot28*(dot2+dot3)-(dot26+dot27-dot18*dot4)*(dot10+dot12-dot6))+
					4.*dot17*dot18*dot2*dot6+2.*dot11*dot22*dot6-2.*dot13*dot22*dot6+2.*dot11*dot23*dot6-2.*dot13*dot23*dot6-
					2.*dot2*dot25*dot6+4.*dot17*dot18*dot3*dot6-2.*dot25*dot3*dot6+4.*dot13*dot18*dot5*dot6-2.*dot26*dot5*dot6-
					2.*dot27*dot5*dot6-4.*dot10*dot13*dot18*dot7-4.*dot12*dot13*dot18*dot7+4.*dot15*dot18*dot2*dot7+
					2.*dot10*dot26*dot7+2.*dot12*dot26*dot7+2.*dot10*dot27*dot7+2.*dot12*dot27*dot7-4.*dot2*dot28*dot7+
					4.*dot15*dot18*dot3*dot7-4.*dot28*dot3*dot7+4.*dot13*dot18*dot6*dot7-2.*dot26*dot6*dot7-2.*dot27*dot6*dot7-
					4.*dot10*dot17*dot18*dot8-4.*dot12*dot17*dot18*dot8+2.*dot10*dot25*dot8+2.*dot12*dot25*dot8+4.*dot17*dot18*dot6*dot8-
					2.*dot25*dot6*dot8-2.*dot10*dot25*dot9-2.*dot12*dot25*dot9+2.*dot25*dot6*dot9+
					(2.*dot12*dot13*dot19+2.*dot15*dot19*dot2+2.*dot15*dot19*dot3-dot13*dot32+2.*dot12*dot34+
					 dot12*dot35+2.*dot13*dot16*dot37-2.*dot12*dot17*dot37-2.*dot12*dot17*dot38+
					 2.*dot17*dot2*dot39+2.*dot17*dot3*dot39-4.*dot12*dot19*dot4+dot32*dot4-
					 2.*dot16*dot37*dot4-2.*dot16*dot2*dot40-2.*dot16*dot3*dot40-2.*dot2*dot46-
					 2.*dot3*dot46-2.*dot13*dot39*dot5+2.*dot39*dot4*dot5+2.*dot12*dot40*dot5-
					 2.*dot13*dot19*dot6-2.*dot34*dot6-dot35*dot6+2.*dot17*dot37*dot6+2.*dot17*dot38*dot6+
					 4.*dot19*dot4*dot6-2.*dot40*dot5*dot6+
					 dot11*(2.*dot12*dot19-dot32+2.*dot16*dot37-2.*dot39*dot5-2.*dot19*dot6)+
					 2.*dot40*(dot12-dot6)*dot7+
					 dot10*(2.*dot34+dot35-2.*dot17*(dot37+dot38)+2.*dot19*(dot11+dot13-2.*dot4)+2.*dot40*(dot5+dot7)))*M-
					2.*dot1*(-2.*dot12*dot17*dot18+dot20+dot21-dot15*(dot22+dot23)+dot13*(2.*dot16*dot18-dot24)+
						 dot11*dot24-dot10*dot25+dot12*dot25+dot16*(dot26+dot27)+2.*dot14*(-(dot15*dot18)+dot28)-
						 dot17*dot29-dot17*dot30-2.*dot16*dot18*dot4+dot24*dot4+2.*dot15*dot18*dot5-2.*dot28*dot5+
						 2.*dot17*dot18*dot6-dot25*dot6+2.*dot15*dot18*dot7-2.*dot28*dot7+
						 (dot15*dot19+dot17*dot39-dot16*dot40-dot46)*M)+
					2.*(2.*dot10*dot17*dot18+2.*dot12*dot17*dot18+dot20-2.*dot15*dot22-3.*dot10*dot25-
					    3.*dot12*dot25+2.*dot16*dot26-3.*dot16*dot18*dot4+2.*dot24*dot4+
					    3.*dot15*dot18*dot5-2.*dot28*dot5-2.*dot17*dot18*dot6+3.*dot25*dot6)*M2))
	      /(a1*a2*(dot1-dot2-dot3)*(-dot1+dot2+dot3+dot8+dot9+M2)*(-(a2*dot1)+dot2+dot9-a1*(dot2+dot3+dot9)+M2+(-2.+a1)*a1*M2));
	    // diag[14]=-0.125*(M2*(2.*dot2*dot20+4.*dot2*dot21-4.*dot15*dot2*dot22+4.*dot10*(-dot11+dot13)*dot23-8.*dot15*dot2*dot23-4.*dot10*dot2*dot25+4.*dot16*dot2*dot26-8.*dot10*dot14*dot27+8.*dot16*dot2*dot27-8.*dot10*dot17*dot18*dot3-2.*dot20*dot3+4.*dot15*dot22*dot3+4.*dot10*dot25*dot3-4.*dot16*dot26*dot3+4.*dot17*dot29*dot3-4.*dot17*dot2*dot30+8.*dot14*dot15*dot31+4.*dot11*dot16*dot31-4.*dot13*dot16*dot31-8.*dot16*dot18*dot2*dot4+4.*dot10*dot23*dot4+4.*dot2*dot24*dot4+8.*dot16*dot18*dot3*dot4-4.*dot24*dot3*dot4+8.*dot14*dot30*dot4-4.*dot16*dot31*dot4-4.*dot14*dot42-2.*dot11*dot44+2.*dot13*dot44+2.*dot4*dot44+8.*dot15*dot18*dot2*dot5-4.*dot2*dot28*dot5-8.*dot15*dot18*dot3*dot5+4.*dot28*dot3*dot5+4.*dot11*dot30*dot5-4.*dot13*dot30*dot5-4.*dot30*dot4*dot5-8.*dot10*dot13*dot18*dot7+8.*dot15*dot18*dot2*dot7+4.*dot10*dot27*dot7-4.*dot2*dot28*dot7+4.*dot13*dot29*dot7-4.*dot15*dot31*dot7-4.*dot30*dot4*dot7+2.*dot42*dot7-8.*dot10*dot17*dot18*dot8-2.*dot20*dot8+4.*dot15*dot22*dot8+4.*dot10*dot25*dot8-4.*dot16*dot26*dot8+4.*dot17*dot29*dot8+8.*dot16*dot18*dot4*dot8-4.*dot24*dot4*dot8-8.*dot15*dot18*dot5*dot8+4.*dot28*dot5*dot8+2.*dot20*dot9-4.*dot15*dot22*dot9-4.*dot10*dot25*dot9+4.*dot16*dot26*dot9-8.*dot16*dot18*dot4*dot9+4.*dot24*dot4*dot9+8.*dot15*dot18*dot5*dot9-4.*dot28*dot5*dot9+(-2.*dot15*dot36+dot45-2.*(dot33*dot4-2.*dot16*dot38*dot4-dot16*dot41+dot17*dot48+2.*dot15*dot38*dot5-dot47*dot5+a1*(dot11+dot13-dot4)*(2.*dot10*dot19-dot32+2.*dot16*dot37-2.*dot39*dot5)-2.*dot15*dot37*dot7+2.*dot39*dot4*dot7+dot49*dot7-dot10*(dot35-2.*dot17*dot38+2.*dot40*dot7)))*M+2.*dot1*(dot20+4.*dot14*(-2.*dot15*dot18+dot28)-2.*(2.*dot11*dot16*dot18-2.*dot13*dot16*dot18-2.*dot10*dot17*dot18+dot15*dot22-dot11*dot24+dot13*dot24+dot10*dot25-dot16*dot26+dot17*dot29-2.*dot15*dot18*dot5+dot28*dot5-2.*dot15*dot18*dot7+dot28*dot7)+(2.*dot15*dot19+2.*dot17*dot39-2.*dot16*dot40-dot46)*M)+2.*(dot20-2.*a1*(2.*dot10*dot17*dot18+dot20-2.*dot15*dot22-2.*dot10*dot25+2.*dot16*dot26-dot17*dot29-4.*dot16*dot18*dot4+2.*dot24*dot4+4.*dot15*dot18*dot5-2.*dot28*dot5)-2.*(dot15*dot22+dot10*dot25-dot16*dot26+2.*dot16*dot18*dot4-dot24*dot4-2.*dot15*dot18*dot5+dot28*dot5))*M2))/(a1*a2*dot2*(-dot1+dot2+dot3+dot8+dot9+M2)*(-dot9-a2*M2));
	    // diag[15]=(M2*(2.*dot2*dot20+4.*dot12*dot2*dot25-2.*dot1*(dot20+2.*dot12*dot25)+8.*dot12*dot14*dot27+8.*dot12*dot17*dot18*dot3-2.*dot20*dot3-4.*dot21*dot3-4.*dot12*dot25*dot3+4.*dot17*dot29*dot3+4.*dot17*dot3*dot30+4.*dot12*dot23*dot4-4.*dot14*dot42-2.*dot4*dot44-8.*dot12*dot27*dot5+4.*dot42*dot5+8.*dot14*dot26*dot6-8.*dot17*dot18*dot3*dot6+8.*dot25*dot3*dot6+4.*dot22*dot4*dot6-8.*dot26*dot5*dot6+2.*(dot11-dot13)*(2.*dot12*dot23-dot44+2.*dot22*dot6)+8.*dot12*dot13*dot18*dot7-4.*dot12*dot27*dot7+4.*dot13*dot29*dot7+4.*dot28*dot3*dot7-8.*dot12*dot18*dot4*dot7-4.*dot29*dot4*dot7+2.*dot42*dot7-4.*dot26*dot6*dot7+8.*dot12*dot17*dot18*dot8-2.*dot20*dot8-4.*dot12*dot25*dot8+4.*dot17*dot29*dot8+2.*dot20*dot9+4.*dot12*dot25*dot9+(2.*a1*(2.*dot12*dot19+dot32)*(dot11+dot13-dot4)+dot45-2.*dot34*dot6-2.*(dot17*dot48-2.*dot17*dot37*dot6+dot49*dot7+dot12*(dot35-2.*dot17*dot38+2.*dot40*dot7)))*M+2.*(dot20+2.*dot12*dot25+2.*a1*(2.*dot12*dot17*dot18-dot20-2.*dot12*dot25+dot17*dot29))*M2))/(8.*a1*dot3*(-dot1+dot2+dot3+dot8+dot9+M2)*(-(a2*dot1)+dot2+dot9-a1*(dot2+dot3+dot9)+M2+(-2.+a1)*a1*M2));
	    // diag[16]=(M2*(-2.*dot2*dot20-2.*dot2*dot21+4.*dot15*dot2*dot22+4.*dot15*dot2*dot23-4.*dot16*dot2*dot26-4.*dot16*dot2*dot27+4.*dot17*dot2*dot29-2.*dot20*dot3-2.*dot21*dot3+4.*dot17*dot2*dot30+4.*a1*(dot11+dot13)*((dot10-dot12)*(dot22+dot23)+2.*dot14*(dot29+dot30))-4.*dot14*dot15*dot31-4.*dot11*dot16*dot31+4.*dot10*dot17*dot31-4.*a1*dot10*dot22*dot4+4.*a1*dot12*dot22*dot4-4.*a1*dot10*dot23*dot4-4.*dot12*dot23*dot4+4.*a1*dot12*dot23*dot4-8.*a1*dot14*dot29*dot4-4.*dot24*dot3*dot4+4.*dot14*dot30*dot4-8.*a1*dot14*dot30*dot4-2.*dot10*dot43-2.*dot12*dot43+2.*dot11*dot44+2.*dot13*dot44+8.*a1*dot11*dot12*dot18*dot5+8.*a1*dot12*dot13*dot18*dot5+4.*dot12*dot27*dot5-4.*a1*dot11*dot29*dot5-4.*a1*dot13*dot29*dot5+4.*dot28*dot3*dot5-4.*a1*dot11*dot30*dot5-4.*dot13*dot30*dot5-4.*a1*dot13*dot30*dot5+4.*dot15*dot31*dot5-8.*a1*dot12*dot18*dot4*dot5+4.*a1*dot29*dot4*dot5+4.*a1*dot30*dot4*dot5-2.*dot42*dot5-8.*a1*dot11*dot14*dot18*dot6-8.*a1*dot13*dot14*dot18*dot6-4.*dot13*dot22*dot6+4.*dot14*dot26*dot6+4.*dot25*dot3*dot6-8.*dot14*dot18*dot4*dot6+8.*a1*dot14*dot18*dot4*dot6+8.*dot13*dot18*dot5*dot6+8.*a1*dot11*dot12*dot18*dot7+8.*a1*dot12*dot13*dot18*dot7-4.*dot10*dot26*dot7-4.*dot10*dot27*dot7+4.*dot12*dot27*dot7+4.*dot11*dot29*dot7-4.*a1*dot11*dot29*dot7-4.*a1*dot13*dot29*dot7+4.*dot28*dot3*dot7+4.*dot11*dot30*dot7-4.*a1*dot11*dot30*dot7-4.*dot13*dot30*dot7-4.*a1*dot13*dot30*dot7+4.*dot15*dot31*dot7-8.*a1*dot12*dot18*dot4*dot7+4.*a1*dot29*dot4*dot7+4.*a1*dot30*dot4*dot7-2.*dot42*dot7+8.*dot13*dot18*dot6*dot7+4.*dot1*(dot20+dot21-dot15*(dot22+dot23)+dot16*(dot26+dot27)+dot14*dot28-dot17*(dot29+dot30)+dot24*(-dot13+dot4)+dot25*(dot12-dot6)-dot28*(dot5+dot7))+8.*a1*dot11*dot16*dot18*dot8+8.*a1*dot13*dot16*dot18*dot8-2.*dot20*dot8-4.*dot21*dot8+4.*dot15*dot23*dot8-4.*a1*dot11*dot24*dot8-4.*a1*dot13*dot24*dot8-4.*dot16*dot27*dot8+4.*dot17*dot30*dot8-8.*a1*dot16*dot18*dot4*dot8-4.*dot24*dot4*dot8+4.*a1*dot24*dot4*dot8+4.*dot28*dot5*dot8+4.*dot25*dot6*dot8-8.*dot15*dot18*dot7*dot8+8.*dot28*dot7*dot8-2.*dot20*dot9-4.*dot21*dot9+4.*dot15*dot22*dot9+4.*dot15*dot23*dot9+4.*a1*dot11*dot24*dot9+4.*a1*dot13*dot24*dot9-4.*dot16*dot26*dot9-4.*dot16*dot27*dot9+4.*dot17*dot29*dot9+4.*dot17*dot30*dot9-4.*a1*dot24*dot4*dot9+4.*dot25*dot6*dot9-2.*(1.+a1)*dot1*dot46*M+(2.*dot15*dot36-dot45+2.*dot17*dot48+2.*(-(dot33*dot4)-dot16*dot41+dot47*dot5+2.*a12*dot19*(dot11+dot13-dot4)*(dot10+dot12-dot6)+dot35*dot6+2.*dot19*dot4*dot6+(-2.*dot15*(dot37+dot38)+2.*dot47+dot49)*dot7-2.*dot40*dot6*(dot5+dot7)+a1*((dot10+dot12)*dot35-dot35*dot6-dot11*(dot32+2.*(dot33-dot16*(dot37+dot38)-2.*dot19*dot6+dot39*(dot5+dot7)))-dot13*(dot32+2.*(dot33-dot16*(dot37+dot38)-2.*dot19*dot6+dot39*(dot5+dot7)))+dot4*(dot32+2.*(dot33-dot16*(dot37+dot38)-2.*dot19*dot6+dot39*(dot5+dot7)))+dot46*(dot2+dot3+dot8+dot9))))*M+2.*a1*dot46*pow<3,1>(M)+(-2.*dot20+8.*a12*(dot16*dot18-dot24)*(dot11+dot13-dot4)+4.*(-dot21+dot15*(dot22+dot23)-dot16*(dot26+dot27)+dot17*(dot29+dot30)+dot25*dot6)+4.*a1*(dot21+dot16*(dot26+dot27)-dot17*(dot29+dot30)+dot24*(dot11+dot13-2.*dot4)+(dot17*dot18-dot25)*dot6+dot28*(dot5+dot7)-dot15*(dot22+dot23+dot18*dot7)))*M2))/(16.*a1*a2*(-dot1+dot2+dot3+dot8+dot9+M2)*(-dot9-a2*M2)*(-dot1+a1*(dot8+dot9+M2)));
	    // diag[17]=(M2*(dot1*(4.*(2.+a1)*dot12*dot17*dot18-2.*(2.+3.*a1)*dot12*dot25+2.*(dot20-2.*(dot15*dot22+dot11*dot24-dot10*dot25-dot16*dot26+dot14*dot28))-4.*dot17*dot29+4.*dot24*dot4-4.*dot25*dot6-a1*(dot20+dot21-2.*(-2.*dot10*dot25+dot24*(dot11+2.*dot13-2.*dot4)+3.*dot25*dot6+dot17*(dot29+dot30-2.*dot18*dot6)+dot28*(-dot14+dot5+dot7))))+4.*a12*(dot10+dot12-dot6)*M*(dot19*(dot11+dot13-dot4)+2.*(dot17*dot18-dot25)*M)+2.*(dot21*dot3-2.*dot15*dot23*dot3+2.*dot16*dot27*dot3-2.*dot17*dot3*dot30-2.*dot14*dot15*dot31-2.*dot11*dot16*dot31-2.*dot12*dot17*dot31+2.*dot12*dot23*dot4-2.*dot14*dot30*dot4+2.*dot14*dot42+dot12*dot43+2.*dot11*dot44-dot4*dot44-2.*dot12*dot27*dot5-2.*dot11*dot30*dot5+2.*dot15*dot31*dot5+2.*dot30*dot4*dot5-dot42*dot5-2.*dot11*dot22*dot6+2.*dot11*dot23*dot6-2.*dot13*dot23*dot6-2.*dot14*dot26*dot6+4.*dot14*dot27*dot6+4.*dot17*dot18*dot3*dot6-2.*dot25*dot3*dot6-2.*dot17*dot31*dot6+4.*dot14*dot18*dot4*dot6+2.*dot22*dot4*dot6+2.*dot23*dot4*dot6+2.*dot43*dot6+4.*dot11*dot18*dot5*dot6-4.*dot27*dot5*dot6-4.*dot18*dot4*dot5*dot6+2.*dot12*dot26*dot7+4.*dot11*dot29*dot7+4.*dot15*dot18*dot3*dot7-2.*dot28*dot3*dot7-4.*dot12*dot18*dot4*dot7-2.*dot29*dot4*dot7+4.*dot13*dot18*dot6*dot7+2.*dot26*dot6*dot7-2.*dot27*dot6*dot7-4.*dot18*dot4*dot6*dot7-dot2*(2.*dot20+dot21-2.*dot15*dot22+2.*dot16*dot26+2.*dot24*dot4+2.*dot25*(dot10-2.*dot6)-2.*dot17*(2.*dot29+dot30-2.*dot18*dot6)-2.*dot28*(dot5+dot7))+4.*dot17*dot18*dot6*dot8-2.*dot25*dot6*dot8-2.*dot20*dot9+2.*dot15*dot22*dot9-2.*dot16*dot26*dot9+4.*dot17*dot29*dot9-2.*dot24*dot4*dot9+2.*dot28*dot5*dot9+2.*dot25*dot6*dot9-dot6*(dot35-2.*dot17*dot38+2.*dot40*dot7)*M+2.*(-dot20+dot15*dot22-dot16*dot26+2.*dot17*dot29-dot24*dot4+dot28*dot5+dot25*dot6)*M2+dot10*(-4.*dot14*dot27-4.*dot17*dot18*dot3+2.*dot25*dot3+4.*dot17*dot31-3.*dot43+4.*dot27*dot5+dot35*M-2.*(dot11*dot23-dot13*dot23+dot23*dot4+2.*dot13*dot18*dot7+2.*dot26*dot7-dot27*dot7-2.*dot18*dot4*dot7+2.*dot17*dot18*dot8-dot25*dot8+dot25*dot9+dot17*dot38*M-dot40*dot7*M+dot25*M2)))+a1*(4.*dot10*(dot13*(dot22-dot23)+dot11*(dot22+dot23))+8.*dot10*dot14*dot27+8.*dot11*dot14*dot29+8.*dot13*dot14*dot29-8.*dot11*dot16*dot18*dot3-8.*dot13*dot16*dot18*dot3+8.*dot10*dot17*dot18*dot3+dot20*dot3+dot21*dot3+6.*dot11*dot24*dot3+4.*dot13*dot24*dot3-4.*dot10*dot25*dot3+2.*dot14*dot28*dot3-2.*dot17*dot29*dot3-2.*dot17*dot3*dot30+4.*dot11*dot16*dot31+4.*dot13*dot16*dot31-4.*dot10*dot17*dot31-4.*dot10*dot22*dot4+4.*dot10*dot23*dot4-8.*dot14*dot29*dot4+8.*dot16*dot18*dot3*dot4-4.*dot24*dot3*dot4-4.*dot16*dot31*dot4+4.*dot10*dot43-4.*dot11*dot44-4.*dot13*dot44+4.*dot4*dot44-8.*dot10*dot27*dot5-4.*dot11*dot29*dot5-4.*dot13*dot29*dot5-2.*dot28*dot3*dot5+4.*dot11*dot30*dot5+4.*dot13*dot30*dot5+4.*dot29*dot4*dot5-4.*dot30*dot4*dot5+4.*dot11*dot22*dot6+4.*dot13*dot22*dot6-4.*dot11*dot23*dot6+4.*dot13*dot23*dot6-8.*dot14*dot27*dot6-4.*dot17*dot18*dot3*dot6+2.*dot25*dot3*dot6+4.*dot17*dot31*dot6-4.*dot22*dot4*dot6-4.*dot23*dot4*dot6-4.*dot43*dot6-8.*dot11*dot18*dot5*dot6-8.*dot13*dot18*dot5*dot6+8.*dot27*dot5*dot6+8.*dot18*dot4*dot5*dot6+8.*dot10*dot13*dot18*dot7+4.*dot10*dot26*dot7-4.*dot10*dot27*dot7-8.*dot11*dot29*dot7-8.*dot13*dot29*dot7-2.*dot28*dot3*dot7-8.*dot10*dot18*dot4*dot7+8.*dot29*dot4*dot7-8.*dot13*dot18*dot6*dot7-4.*dot26*dot6*dot7+4.*dot27*dot6*dot7+8.*dot18*dot4*dot6*dot7+dot2*(dot20+dot21-2.*(-2.*dot10*dot25+dot24*(dot11+2.*dot13-2.*dot4)+3.*dot25*dot6+dot17*(dot29+dot30-2.*dot18*dot6)+dot28*(-dot14+dot5+dot7)))+8.*dot10*dot17*dot18*dot8+dot20*dot8+dot21*dot8+2.*dot11*dot24*dot8-4.*dot10*dot25*dot8+2.*dot14*dot28*dot8-2.*dot17*dot29*dot8-2.*dot17*dot30*dot8-2.*dot28*dot5*dot8-4.*dot17*dot18*dot6*dot8+2.*dot25*dot6*dot8-2.*dot28*dot7*dot8+dot20*dot9+dot21*dot9+2.*dot11*dot24*dot9+4.*dot10*dot25*dot9+2.*dot14*dot28*dot9-2.*dot17*dot29*dot9-2.*dot17*dot30*dot9-2.*dot28*dot5*dot9+4.*dot17*dot18*dot6*dot9-6.*dot25*dot6*dot9-2.*dot28*dot7*dot9-2.*(dot10-dot6)*(dot35-2.*dot17*dot38+2.*dot19*(dot11+dot13-dot4)+2.*dot40*dot7)*M+(5.*dot20+dot21-4.*dot15*dot22+4.*dot16*dot26+2.*dot24*(dot11+2.*dot4)-2.*(4.*dot10*dot17*dot18-6.*dot10*dot25+5.*dot17*dot29+dot17*dot30-6.*dot17*dot18*dot6+7.*dot25*dot6)+2.*dot28*(dot14-3.*dot5-dot7))*M2-2.*dot12*(2.*dot11*(dot22-dot23)-4.*dot14*dot27+2.*(dot22+dot23)*(dot13-dot4)-2.*dot43-4.*dot11*dot18*dot5-4.*dot13*dot18*dot5+4.*dot27*dot5+4.*dot18*dot4*dot5-4.*dot13*dot18*dot7-2.*dot26*dot7+2.*dot27*dot7+4.*dot18*dot4*dot7+dot35*M+2.*dot40*dot7*M+dot25*(-3.*dot2+dot3+dot8-3.*(dot9+M2))+2.*dot17*(dot31-dot38*M+dot18*(dot2-dot3-dot8+dot9+M2))))))/(16.*a1*a2*(-dot1+dot2+dot3+dot8+dot9+M2)*(-(a2*dot1)+dot2+dot9-a1*(dot2+dot3+dot9)+M2+(-2.+a1)*a1*M2)*(-dot9-a2*M2));
	    // diag[18]=-0.125*(M2*(-4.*dot2*dot20-2.*dot2*dot21+4.*dot11*(dot10-dot12)*dot22-4.*dot11*dot2*dot24-4.*dot16*dot2*dot26+8.*dot11*dot14*dot29+8.*dot17*dot2*dot29-8.*dot11*dot16*dot18*dot3+2.*dot21*dot3+4.*dot11*dot24*dot3+4.*dot16*dot27*dot3+4.*dot17*dot2*dot30-4.*dot17*dot3*dot30+4.*dot10*dot17*dot31-4.*dot12*dot17*dot31+4.*dot14*dot42-2.*dot10*dot43+2.*dot12*dot43+8.*dot11*dot12*dot18*dot5-4.*dot12*dot27*dot5+4.*dot2*dot28*dot5-4.*dot11*dot29*dot5-2.*dot42*dot5-8.*dot17*dot18*dot2*dot6-4.*dot11*dot22*dot6+4.*dot2*dot25*dot6-8.*dot14*dot26*dot6+8.*dot17*dot18*dot3*dot6-4.*dot25*dot3*dot6-4.*dot17*dot31*dot6+2.*dot43*dot6+4.*dot26*dot5*dot6-4.*dot10*dot26*dot7+4.*dot12*dot26*dot7+4.*dot2*dot28*dot7-4.*dot28*dot3*dot7+4.*dot26*dot6*dot7+2.*dot1*(2.*dot20+dot21+2.*dot11*dot24+2.*dot16*dot26-2.*dot25*dot6-2.*dot17*(2.*dot29+dot30-2.*dot18*dot6)-2.*dot28*(dot5+dot7))-4.*dot20*dot9-4.*dot16*dot26*dot9+8.*dot17*dot29*dot9+4.*dot28*dot5*dot9+(-dot45-2.*a2*dot10*(2.*dot11*dot19-dot35+2.*dot17*dot38-2.*dot40*dot7)+2.*(-(dot11*dot32)+dot17*dot48-2.*a2*dot11*dot19*(dot12-dot6)+dot34*dot6-a2*dot35*dot6-2.*dot17*(dot37-a2*dot38)*dot6+dot49*dot7-2.*a2*dot40*dot6*dot7+a2*dot12*(dot35-2.*dot17*dot38+2.*dot40*dot7)))*M-4.*a2*(dot20+dot16*dot26-2.*dot17*dot29-dot28*dot5)*M2))/(a2*(dot1-dot2-dot3)*(-(a2*dot1)+dot2+dot9-a1*(dot2+dot3+dot9)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot9-a2*M2));
	    // diag[19]=(M2*(-2.*dot2*dot21+4.*(-dot10+dot12)*dot13*dot22+4.*dot15*dot2*dot23+4.*dot13*dot2*dot24-8.*dot13*dot14*dot29+8.*dot13*dot16*dot18*dot3+4.*dot20*dot3+2.*dot21*dot3-8.*dot15*dot22*dot3-4.*dot15*dot23*dot3-4.*dot13*dot24*dot3+4.*dot16*dot26*dot3+4.*dot16*dot27*dot3-8.*dot14*dot15*dot31-4.*dot10*dot23*dot4+4.*dot12*dot23*dot4-8.*dot16*dot18*dot3*dot4+8.*dot24*dot3*dot4-8.*dot14*dot30*dot4+4.*dot14*dot42-2.*dot10*dot43+2.*dot12*dot43-8.*dot12*dot13*dot18*dot5-4.*dot12*dot27*dot5+4.*dot13*dot29*dot5+8.*dot15*dot18*dot3*dot5-4.*dot28*dot3*dot5+4.*dot15*dot31*dot5+4.*dot30*dot4*dot5-2.*dot42*dot5-4.*dot13*dot22*dot6-4.*dot23*dot4*dot6-2.*dot43*dot6+8.*dot13*dot18*dot5*dot6+4.*dot27*dot5*dot6+8.*dot13*dot29*dot7+8.*dot15*dot31*dot7+8.*dot30*dot4*dot7-4.*dot42*dot7+4.*dot20*dot8-8.*dot15*dot22*dot8+4.*dot16*dot26*dot8-8.*dot16*dot18*dot4*dot8+8.*dot24*dot4*dot8+8.*dot15*dot18*dot5*dot8-4.*dot28*dot5*dot8+(2.*dot13*dot32+2.*a2*dot10*(2.*dot13*dot19+dot35)+2.*a2*dot12*(2.*dot13*dot19+dot35)+2.*dot15*dot36+2.*dot33*dot4-dot45-2.*a2*(2.*dot13*dot19+dot35)*dot6)*M+2.*dot1*(-2.*dot20-dot21+4.*dot15*dot22+2.*dot15*dot23+2.*dot13*(-2.*dot16*dot18+dot24)-2.*dot16*(dot26+dot27)+dot14*(8.*dot15*dot18-4.*dot28)+4.*(dot16*dot18-dot24)*dot4-8.*dot15*dot18*dot5+4.*dot28*dot5+2.*dot25*(dot10-dot12+dot6)-8.*dot15*dot18*dot7+4.*dot28*dot7-2.*dot15*dot19*M+dot46*M)+4.*a1*(dot20-2.*dot15*dot22+dot16*dot26-2.*dot16*dot18*dot4+2.*dot24*dot4+2.*dot15*dot18*dot5-dot28*dot5)*M2))/(8.*a1*a2*(dot1-dot2-dot3)*(dot8+a1*M2)*(-dot1+dot3+dot8+a1*M2));
	    // diag[20]=-0.0625*(M2*(-2.*dot2*dot20-2.*dot2*dot21+4.*dot15*dot2*dot22+4.*dot15*dot2*dot23-2.*dot20*dot3-2.*dot21*dot3-4.*dot16*dot26*dot3-4.*dot16*dot27*dot3+4.*dot17*dot29*dot3+4.*dot17*dot3*dot30-4.*a2*(dot11+dot13)*((dot10-dot12)*(dot22+dot23)+2.*dot14*(dot29+dot30))-4.*dot14*dot15*dot31-4.*dot13*dot16*dot31+4.*dot12*dot17*dot31-4.*a1*dot10*dot22*dot4+4.*a1*dot12*dot22*dot4-4.*a1*dot10*dot23*dot4-4.*dot12*dot23*dot4+4.*a1*dot12*dot23*dot4-8.*a1*dot14*dot29*dot4+8.*dot16*dot18*dot3*dot4-4.*dot24*dot3*dot4+4.*dot14*dot30*dot4-8.*a1*dot14*dot30*dot4-2.*dot10*dot43-2.*dot12*dot43+2.*dot11*dot44+2.*dot13*dot44-8.*dot11*dot12*dot18*dot5+8.*a1*dot11*dot12*dot18*dot5-8.*dot12*dot13*dot18*dot5+8.*a1*dot12*dot13*dot18*dot5+4.*dot12*dot27*dot5+4.*dot11*dot29*dot5-4.*a1*dot11*dot29*dot5+4.*dot13*dot29*dot5-4.*a1*dot13*dot29*dot5+4.*dot28*dot3*dot5+4.*dot11*dot30*dot5-4.*a1*dot11*dot30*dot5-4.*a1*dot13*dot30*dot5+4.*dot15*dot31*dot5-8.*a1*dot12*dot18*dot4*dot5+4.*a1*dot29*dot4*dot5+4.*a1*dot30*dot4*dot5-2.*dot42*dot5+8.*dot11*dot14*dot18*dot6-8.*a1*dot11*dot14*dot18*dot6+8.*dot13*dot14*dot18*dot6-8.*a1*dot13*dot14*dot18*dot6-4.*dot13*dot22*dot6+4.*dot14*dot26*dot6-8.*dot17*dot18*dot3*dot6+4.*dot25*dot3*dot6-8.*dot14*dot18*dot4*dot6+8.*a1*dot14*dot18*dot4*dot6+8.*dot13*dot18*dot5*dot6+4.*dot1*(-(dot11*dot24)+dot10*dot25+(2.*dot15*dot18-dot28)*(dot14-dot5-dot7))-8.*dot11*dot12*dot18*dot7+8.*a1*dot11*dot12*dot18*dot7-8.*dot12*dot13*dot18*dot7+8.*a1*dot12*dot13*dot18*dot7-4.*dot12*dot26*dot7+4.*dot11*dot29*dot7-4.*a1*dot11*dot29*dot7+8.*dot13*dot29*dot7-4.*a1*dot13*dot29*dot7+4.*dot28*dot3*dot7+4.*dot11*dot30*dot7-4.*a1*dot11*dot30*dot7+4.*dot13*dot30*dot7-4.*a1*dot13*dot30*dot7+4.*dot15*dot31*dot7+8.*dot12*dot18*dot4*dot7-8.*a1*dot12*dot18*dot4*dot7+4.*a1*dot29*dot4*dot7+4.*a1*dot30*dot4*dot7-2.*dot42*dot7-8.*dot14*dot15*dot18*dot8+8.*a1*dot14*dot15*dot18*dot8-8.*dot11*dot16*dot18*dot8+8.*a1*dot11*dot16*dot18*dot8-8.*dot12*dot17*dot18*dot8+8.*a1*dot12*dot17*dot18*dot8+2.*dot20*dot8-4.*a1*dot20*dot8+4.*dot21*dot8-4.*a1*dot21*dot8-4.*dot15*dot22*dot8+4.*a1*dot15*dot22*dot8-4.*dot15*dot23*dot8+4.*a1*dot15*dot23*dot8+8.*dot11*dot24*dot8-8.*a1*dot11*dot24*dot8-4.*dot10*dot25*dot8+4.*a1*dot10*dot25*dot8+4.*dot12*dot25*dot8-4.*a1*dot12*dot25*dot8-4.*a1*dot16*dot26*dot8+4.*dot16*dot27*dot8-4.*a1*dot16*dot27*dot8+8.*dot14*dot28*dot8-8.*a1*dot14*dot28*dot8+4.*a1*dot17*dot29*dot8-4.*dot17*dot30*dot8+4.*a1*dot17*dot30*dot8+8.*dot15*dot18*dot5*dot8-8.*a1*dot15*dot18*dot5*dot8-4.*dot28*dot5*dot8+8.*a1*dot28*dot5*dot8+8.*dot17*dot18*dot6*dot8-8.*a1*dot17*dot18*dot6*dot8-4.*dot25*dot6*dot8+4.*a1*dot25*dot6*dot8+8.*dot15*dot18*dot7*dot8-8.*a1*dot15*dot18*dot7*dot8-8.*dot28*dot7*dot8+8.*a1*dot28*dot7*dot8-8.*dot14*dot15*dot18*dot9+8.*a1*dot14*dot15*dot18*dot9+8.*dot13*dot16*dot18*dot9-8.*a1*dot13*dot16*dot18*dot9-8.*dot12*dot17*dot18*dot9+8.*a1*dot12*dot17*dot18*dot9+2.*dot20*dot9-4.*a1*dot20*dot9+4.*dot21*dot9-4.*a1*dot21*dot9+4.*a1*dot15*dot22*dot9-4.*dot15*dot23*dot9+4.*a1*dot15*dot23*dot9-8.*dot13*dot24*dot9+8.*a1*dot13*dot24*dot9-4.*dot10*dot25*dot9+4.*a1*dot10*dot25*dot9+4.*dot12*dot25*dot9-4.*a1*dot12*dot25*dot9+4.*dot16*dot26*dot9-4.*a1*dot16*dot26*dot9+4.*dot16*dot27*dot9-4.*a1*dot16*dot27*dot9+8.*dot14*dot28*dot9-8.*a1*dot14*dot28*dot9-4.*dot17*dot29*dot9+4.*a1*dot17*dot29*dot9-4.*dot17*dot30*dot9+4.*a1*dot17*dot30*dot9-8.*dot16*dot18*dot4*dot9+8.*a1*dot16*dot18*dot4*dot9+4.*dot24*dot4*dot9-8.*a1*dot24*dot4*dot9+8.*dot15*dot18*dot5*dot9-8.*a1*dot15*dot18*dot5*dot9-8.*dot28*dot5*dot9+8.*a1*dot28*dot5*dot9+8.*dot17*dot18*dot6*dot9-8.*a1*dot17*dot18*dot6*dot9-4.*dot25*dot6*dot9+4.*a1*dot25*dot6*dot9+8.*dot15*dot18*dot7*dot9-8.*a1*dot15*dot18*dot7*dot9-8.*dot28*dot7*dot9+8.*a1*dot28*dot7*dot9-2.*dot1*(2.*dot15*dot19+2.*a2*(dot17*dot39-dot16*dot40)+(-2.+a1)*dot46)*M+(2.*dot15*dot36-dot45-2.*a2*dot10*(-2.*a2*(dot11+dot13)*dot19+dot35-2.*dot17*dot38-2.*a1*dot19*dot4+2.*dot40*dot7)+2.*(dot13*dot32-a1*dot13*dot32+2.*dot13*dot33-2.*a1*dot13*dot33-2.*dot13*dot16*dot38+2.*a1*dot13*dot16*dot38+2.*dot17*dot2*dot39-2.*a1*dot17*dot2*dot39+2.*dot17*dot3*dot39-2.*a1*dot17*dot3*dot39+a1*dot32*dot4-dot33*dot4+2.*a1*dot33*dot4+2.*dot16*dot38*dot4-2.*a1*dot16*dot38*dot4-2.*dot16*dot2*dot40+2.*a1*dot16*dot2*dot40-2.*dot16*dot3*dot40+2.*a1*dot16*dot3*dot40-dot2*dot46+a1*dot2*dot46-dot3*dot46+a1*dot3*dot46-4.*dot13*dot19*dot6+6.*a1*dot13*dot19*dot6-2.*a12*dot13*dot19*dot6+dot35*dot6-a1*dot35*dot6-2.*dot17*dot38*dot6+2.*a1*dot17*dot38*dot6+2.*dot19*dot4*dot6-4.*a1*dot19*dot4*dot6+2.*a12*dot19*dot4*dot6+2.*dot13*dot39*dot7-2.*a1*dot13*dot39*dot7-2.*dot39*dot4*dot7+2.*a1*dot39*dot4*dot7+2.*dot40*dot6*dot7-2.*a1*dot40*dot6*dot7-a2*dot12*(-2.*a2*dot13*dot19+dot35-2.*dot17*dot38-2.*a1*dot19*dot4+2.*dot40*dot7)-a2*dot11*(-2.*a2*dot12*dot19-dot32-2.*(dot33-dot16*dot38+(-2.+a1)*dot19*dot6+dot39*dot7))+2.*dot15*dot19*dot8-2.*a1*dot15*dot19*dot8+2.*dot17*dot39*dot8-2.*a1*dot17*dot39*dot8-2.*dot16*dot40*dot8+2.*a1*dot16*dot40*dot8-2.*dot46*dot8+2.*a1*dot46*dot8+2.*a2*(dot15*dot19+dot17*dot39-dot16*dot40-dot46)*dot9))*M+4.*a2*(dot15*dot19+dot17*dot39-dot16*dot40-dot46)*pow<3,1>(M)+2.*(dot20-2.*dot10*dot25-4.*a2*dot14*(dot15*dot18-dot28)+4.*a12*(dot16*dot18-dot24)*(dot11+dot13-dot4)-2.*a1*(2.*dot11*dot16*dot18+4.*dot13*dot16*dot18-2.*dot12*dot17*dot18+dot20-2.*dot11*dot24-4.*dot13*dot24-dot10*dot25+dot12*dot25+2.*dot16*dot26-2.*dot17*dot29-4.*dot16*dot18*dot4+3.*dot24*dot4+2.*dot15*dot18*dot5-3.*dot28*dot5+dot17*dot18*dot6+dot15*dot18*dot7-dot28*dot7)+2.*(2.*dot13*dot16*dot18-2.*dot12*dot17*dot18-2.*dot13*dot24+dot12*dot25+dot16*dot26-dot17*dot29-2.*dot16*dot18*dot4+dot24*dot4+2.*dot15*dot18*dot5-2.*dot28*dot5+dot17*dot18*dot6+dot15*dot18*dot7-dot28*dot7))*M2))/(a1*a2*(dot1-dot2-dot3)*(dot8+a1*M2)*(dot1-a2*(dot8+dot9+M2)));
	    // diag[21]=-0.0625*(M2*(4.*dot12*dot17*dot18*dot2-4.*a1*dot12*dot17*dot18*dot2-3.*dot2*dot20+a1*dot2*dot20-dot2*dot21+a1*dot2*dot21-4.*dot10*dot11*dot22+4.*a1*dot10*dot11*dot22+4.*dot11*dot12*dot22-4.*a1*dot11*dot12*dot22+4.*a1*dot10*dot13*dot22-4.*a1*dot12*dot13*dot22+4.*dot15*dot2*dot22-4.*dot10*dot11*dot23+4.*a1*dot10*dot11*dot23-4.*dot11*dot12*dot23+4.*a1*dot11*dot12*dot23+4.*dot10*dot13*dot23-4.*a1*dot10*dot13*dot23+4.*dot12*dot13*dot23-4.*a1*dot12*dot13*dot23+2.*dot11*dot2*dot24-2.*a1*dot11*dot2*dot24-4.*a1*dot13*dot2*dot24-4.*dot10*dot2*dot25+4.*a1*dot10*dot2*dot25-6.*dot12*dot2*dot25+6.*a1*dot12*dot2*dot25-8.*dot10*dot14*dot27+8.*a1*dot10*dot14*dot27-8.*dot12*dot14*dot27+8.*a1*dot12*dot14*dot27-2.*dot14*dot2*dot28+2.*a1*dot14*dot2*dot28-8.*dot11*dot14*dot29+8.*a1*dot11*dot14*dot29+8.*a1*dot13*dot14*dot29+2.*dot17*dot2*dot29-2.*a1*dot17*dot2*dot29+8.*dot11*dot16*dot18*dot3-8.*a1*dot11*dot16*dot18*dot3-8.*a1*dot13*dot16*dot18*dot3-8.*dot10*dot17*dot18*dot3+8.*a1*dot10*dot17*dot18*dot3-4.*dot12*dot17*dot18*dot3+4.*a1*dot12*dot17*dot18*dot3+dot20*dot3+a1*dot20*dot3+3.*dot21*dot3+a1*dot21*dot3-4.*dot15*dot23*dot3-6.*dot11*dot24*dot3+6.*a1*dot11*dot24*dot3+4.*a1*dot13*dot24*dot3+4.*dot10*dot25*dot3-4.*a1*dot10*dot25*dot3+2.*dot12*dot25*dot3-2.*a1*dot12*dot25*dot3-2.*dot14*dot28*dot3+2.*a1*dot14*dot28*dot3-2.*dot17*dot29*dot3-2.*a1*dot17*dot29*dot3+2.*dot17*dot2*dot30-2.*a1*dot17*dot2*dot30-2.*dot17*dot3*dot30-2.*a1*dot17*dot3*dot30-4.*dot14*dot15*dot31-4.*dot11*dot16*dot31+4.*a1*dot11*dot16*dot31+4.*a1*dot13*dot16*dot31+4.*dot10*dot17*dot31-4.*a1*dot10*dot17*dot31-4.*a1*dot12*dot17*dot31-4.*a1*dot10*dot22*dot4+4.*a1*dot12*dot22*dot4-4.*dot10*dot23*dot4+4.*a1*dot10*dot23*dot4+4.*a1*dot12*dot23*dot4-4.*dot2*dot24*dot4+4.*a1*dot2*dot24*dot4-8.*a1*dot14*dot29*dot4+8.*a1*dot16*dot18*dot3*dot4-4.*a1*dot24*dot3*dot4-4.*dot14*dot30*dot4-4.*a1*dot16*dot31*dot4+4.*dot14*dot42-4.*dot10*dot43+4.*a1*dot10*dot43+4.*a1*dot12*dot43+6.*dot11*dot44-4.*a1*dot11*dot44-2.*dot13*dot44-4.*a1*dot13*dot44+4.*a1*dot4*dot44-8.*dot11*dot12*dot18*dot5+8.*a1*dot11*dot12*dot18*dot5+8.*a1*dot12*dot13*dot18*dot5+8.*dot10*dot27*dot5-8.*a1*dot10*dot27*dot5+8.*dot12*dot27*dot5-8.*a1*dot12*dot27*dot5+2.*dot2*dot28*dot5-2.*a1*dot2*dot28*dot5+4.*dot11*dot29*dot5-4.*a1*dot11*dot29*dot5-4.*a1*dot13*dot29*dot5+2.*dot28*dot3*dot5-2.*a1*dot28*dot3*dot5-4.*dot11*dot30*dot5+4.*a1*dot11*dot30*dot5+4.*a1*dot13*dot30*dot5+4.*dot15*dot31*dot5-8.*a1*dot12*dot18*dot4*dot5+4.*a1*dot29*dot4*dot5+4.*dot30*dot4*dot5-4.*a1*dot30*dot4*dot5-4.*dot42*dot5-4.*dot17*dot18*dot2*dot6+4.*a1*dot17*dot18*dot2*dot6-4.*dot11*dot22*dot6+4.*a1*dot11*dot22*dot6+4.*dot13*dot22*dot6+4.*a1*dot13*dot22*dot6+4.*dot11*dot23*dot6-4.*a1*dot11*dot23*dot6-4.*dot13*dot23*dot6+4.*a1*dot13*dot23*dot6+6.*dot2*dot25*dot6-6.*a1*dot2*dot25*dot6-4.*dot14*dot26*dot6+8.*dot14*dot27*dot6-8.*a1*dot14*dot27*dot6+4.*dot17*dot18*dot3*dot6-4.*a1*dot17*dot18*dot3*dot6-6.*dot25*dot3*dot6+2.*a1*dot25*dot3*dot6+4.*a1*dot17*dot31*dot6+8.*dot14*dot18*dot4*dot6-4.*a1*dot22*dot4*dot6+4.*dot23*dot4*dot6-4.*a1*dot23*dot4*dot6+2.*dot43*dot6-4.*a1*dot43*dot6+8.*dot11*dot18*dot5*dot6-8.*a1*dot11*dot18*dot5*dot6-8.*a1*dot13*dot18*dot5*dot6+4.*dot26*dot5*dot6-8.*dot27*dot5*dot6+8.*a1*dot27*dot5*dot6-8.*dot18*dot4*dot5*dot6+8.*a1*dot18*dot4*dot5*dot6-8.*dot10*dot13*dot18*dot7+8.*a1*dot10*dot13*dot18*dot7-8.*dot12*dot13*dot18*dot7+8.*a1*dot12*dot13*dot18*dot7-4.*dot10*dot26*dot7+4.*a1*dot10*dot26*dot7+4.*a1*dot12*dot26*dot7+4.*dot10*dot27*dot7-4.*a1*dot10*dot27*dot7+4.*dot12*dot27*dot7-4.*a1*dot12*dot27*dot7+2.*dot2*dot28*dot7-2.*a1*dot2*dot28*dot7+8.*dot11*dot29*dot7-8.*a1*dot11*dot29*dot7-4.*dot13*dot29*dot7-8.*a1*dot13*dot29*dot7+8.*dot15*dot18*dot3*dot7-2.*dot28*dot3*dot7-2.*a1*dot28*dot3*dot7+8.*dot10*dot18*dot4*dot7-8.*a1*dot10*dot18*dot4*dot7-8.*a1*dot12*dot18*dot4*dot7+8.*a1*dot29*dot4*dot7-2.*dot42*dot7+8.*dot13*dot18*dot6*dot7-8.*a1*dot13*dot18*dot6*dot7+4.*dot26*dot6*dot7-4.*a1*dot26*dot6*dot7-4.*dot27*dot6*dot7+4.*a1*dot27*dot6*dot7-8.*dot18*dot4*dot6*dot7+8.*a1*dot18*dot4*dot6*dot7+dot1*(dot20-dot21-4.*dot15*dot22-6.*dot11*dot24+4.*dot10*dot25+2.*dot12*(2.*(1.+a1)*dot17*dot18+dot25-3.*a1*dot25)-2.*dot14*dot28+2.*(2.*dot24*dot4-dot25*dot6+dot17*(dot29+dot30-2.*dot18*dot6)+dot28*(dot5+dot7))-a1*(dot20+dot21-2.*(-2.*dot10*dot25+dot24*(dot11+2.*dot13-2.*dot4)+3.*dot25*dot6+dot17*(dot29+dot30-2.*dot18*dot6)+dot28*(-dot14+dot5+dot7))))-8.*dot10*dot17*dot18*dot8+8.*a1*dot10*dot17*dot18*dot8-8.*dot12*dot17*dot18*dot8+8.*a1*dot12*dot17*dot18*dot8+2.*dot20*dot8+4.*dot10*dot25*dot8-4.*a1*dot10*dot25*dot8+4.*dot12*dot25*dot8-4.*a1*dot12*dot25*dot8-4.*dot17*dot29*dot8+8.*dot17*dot18*dot6*dot8-8.*a1*dot17*dot18*dot6*dot8-4.*dot25*dot6*dot8+4.*a1*dot25*dot6*dot8-2.*dot20*dot9+4.*dot15*dot22*dot9-4.*dot10*dot25*dot9+4.*a1*dot10*dot25*dot9-4.*dot12*dot25*dot9+4.*a1*dot12*dot25*dot9-4.*dot24*dot4*dot9+4.*dot25*dot6*dot9-4.*a1*dot25*dot6*dot9+2.*(-(a2*dot11)+a1*(dot13-dot4))*(-dot32-2.*a2*dot19*(dot10+dot12-dot6))*M+2.*(-dot20+2.*dot15*dot22-2.*dot24*dot4+4.*a12*(dot17*dot18-dot25)*(dot10+dot12-dot6)-2.*dot25*(dot10+dot12-dot6)+2.*a1*(-2.*dot10*dot17*dot18-2.*dot12*dot17*dot18+dot20-dot15*dot22+3.*dot10*dot25+3.*dot12*dot25-dot17*dot29+dot24*dot4+2.*dot17*dot18*dot6-3.*dot25*dot6))*M2))/(a1*a2*(dot1-dot2-dot3)*(dot8+a1*M2)*(-(a2*dot1)+dot2+dot9-a1*(dot2+dot3+dot9)+M2+(-2.+a1)*a1*M2));
	    // diag[22]=-0.125*(M2*(16.*dot1*dot14*dot15*dot18+8.*dot1*dot11*dot16*dot18-8.*dot1*dot13*dot16*dot18-8.*dot1*dot10*dot17*dot18-2.*dot1*dot20-2.*dot2*dot20-4.*dot2*dot21+4.*dot1*dot15*dot22+4.*dot15*dot2*dot22+4.*dot10*dot11*dot23-4.*dot10*dot13*dot23+8.*dot15*dot2*dot23-4.*dot1*dot11*dot24+4.*dot1*dot13*dot24+4.*dot1*dot10*dot25+4.*dot10*dot2*dot25-4.*dot1*dot16*dot26-4.*dot16*dot2*dot26+8.*dot10*dot14*dot27-8.*dot16*dot2*dot27-8.*dot1*dot14*dot28+4.*dot1*dot17*dot29+8.*dot10*dot17*dot18*dot3+2.*dot20*dot3-4.*dot15*dot22*dot3-4.*dot10*dot25*dot3+4.*dot16*dot26*dot3-4.*dot17*dot29*dot3+4.*dot17*dot2*dot30-8.*dot14*dot15*dot31-4.*dot11*dot16*dot31+4.*dot13*dot16*dot31+8.*dot16*dot18*dot2*dot4-4.*dot10*dot23*dot4-4.*dot2*dot24*dot4-8.*dot16*dot18*dot3*dot4+4.*dot24*dot3*dot4-8.*dot14*dot30*dot4+4.*dot16*dot31*dot4+4.*dot14*dot42+2.*dot11*dot44-2.*dot13*dot44-2.*dot4*dot44-8.*dot1*dot15*dot18*dot5-8.*dot15*dot18*dot2*dot5+4.*dot1*dot28*dot5+4.*dot2*dot28*dot5+8.*dot15*dot18*dot3*dot5-4.*dot28*dot3*dot5-4.*dot11*dot30*dot5+4.*dot13*dot30*dot5+4.*dot30*dot4*dot5+8.*dot10*dot13*dot18*dot7-8.*dot1*dot15*dot18*dot7-8.*dot15*dot18*dot2*dot7-4.*dot10*dot27*dot7+4.*dot1*dot28*dot7+4.*dot2*dot28*dot7-4.*dot13*dot29*dot7+4.*dot15*dot31*dot7+4.*dot30*dot4*dot7-2.*dot42*dot7+8.*dot10*dot17*dot18*dot8+2.*dot20*dot8-4.*dot15*dot22*dot8-4.*dot10*dot25*dot8+4.*dot16*dot26*dot8-4.*dot17*dot29*dot8-8.*dot16*dot18*dot4*dot8+4.*dot24*dot4*dot8+8.*dot15*dot18*dot5*dot8-4.*dot28*dot5*dot8-2.*dot20*dot9+4.*dot15*dot22*dot9+4.*dot10*dot25*dot9-4.*dot16*dot26*dot9+8.*dot16*dot18*dot4*dot9-4.*dot24*dot4*dot9-8.*dot15*dot18*dot5*dot9+4.*dot28*dot5*dot9+(2.*a2*(dot11+dot13)*(dot32-2.*dot16*dot37)+dot10*(-2.*dot35-4.*a2*dot19*(dot11+dot13-dot4))+2.*(-(a2*dot32)+dot33-2.*dot16*(-(a2*dot37)+dot38))*dot4-dot45+2.*(-(dot16*dot41)+dot1*(-2.*dot15*dot19+2.*dot16*dot40+dot46)-(-2.*a2*dot39*(dot11+dot13-dot4)+dot47)*dot5+dot15*(dot36+2.*dot38*dot5)))*M+2.*(-dot20+2.*a1*(2.*dot10*dot17*dot18+dot20-2.*dot15*dot22-2.*dot10*dot25+2.*dot16*dot26-dot17*dot29-4.*dot16*dot18*dot4+2.*dot24*dot4+4.*dot15*dot18*dot5-2.*dot28*dot5)+2.*(dot15*dot22+dot10*dot25-dot16*dot26+2.*dot16*dot18*dot4-dot24*dot4-2.*dot15*dot18*dot5+dot28*dot5))*M2))/(a2*dot2*(-dot1+dot2+dot3+dot8+dot9+M2)*(-dot2+a1*(-dot1+dot2+dot3+dot8+a1*M2)));
	    // diag[23]=-0.125*(M2*(-2.*dot2*dot20-4.*dot12*dot2*dot25+2.*dot1*(dot20+2.*dot12*dot25)-8.*dot12*dot14*dot27-8.*dot12*dot17*dot18*dot3+2.*dot20*dot3+4.*dot21*dot3+4.*dot12*dot25*dot3-4.*dot17*dot29*dot3-4.*dot17*dot3*dot30-4.*dot12*dot23*dot4+4.*dot14*dot42+2.*dot4*dot44+8.*dot12*dot27*dot5-4.*dot42*dot5-8.*dot14*dot26*dot6+8.*dot17*dot18*dot3*dot6-8.*dot25*dot3*dot6-4.*dot22*dot4*dot6+8.*dot26*dot5*dot6-2.*(dot11-dot13)*(2.*dot12*dot23-dot44+2.*dot22*dot6)-8.*dot12*dot13*dot18*dot7+4.*dot12*dot27*dot7-4.*dot13*dot29*dot7-4.*dot28*dot3*dot7+8.*dot12*dot18*dot4*dot7+4.*dot29*dot4*dot7-2.*dot42*dot7+4.*dot26*dot6*dot7-8.*dot12*dot17*dot18*dot8+2.*dot20*dot8+4.*dot12*dot25*dot8-4.*dot17*dot29*dot8-2.*dot20*dot9-4.*dot12*dot25*dot9+(2.*a2*dot11*(2.*dot12*dot19+dot32)+2.*dot12*(dot35+2.*a2*dot19*(dot13-dot4))+2.*a2*dot32*(dot13-dot4)-dot45+2.*dot34*dot6)*M-2.*(dot20+2.*dot12*dot25+2.*a1*(2.*dot12*dot17*dot18-dot20-2.*dot12*dot25+dot17*dot29))*M2))/(a1*a2*dot3*(-dot1+dot2+dot3+dot8+dot9+M2)*(dot8+a1*M2));
	    // diag[24]=-0.0625*(M2*(-2.*dot2*dot20-2.*dot2*dot21-4.*a2*(dot10+dot12)*((dot11-dot13)*(dot22+dot23)+2.*dot14*(dot26+dot27))-2.*dot20*dot3-2.*dot21*dot3+4.*dot15*dot22*dot3+4.*dot15*dot23*dot3-4.*dot16*dot26*dot3-4.*dot16*dot27*dot3+4.*dot17*dot29*dot3+4.*dot17*dot3*dot30+4.*dot14*dot15*dot31-4.*dot13*dot16*dot31+4.*dot12*dot17*dot31+8.*dot10*dot14*dot18*dot4-8.*a1*dot10*dot14*dot18*dot4+8.*dot12*dot14*dot18*dot4-8.*a1*dot12*dot14*dot18*dot4-4.*dot12*dot23*dot4+8.*dot16*dot18*dot3*dot4-4.*dot24*dot3*dot4+4.*dot14*dot30*dot4-2.*dot10*dot43-2.*dot12*dot43+2.*dot11*dot44+2.*dot13*dot44-8.*dot10*dot13*dot18*dot5+8.*a1*dot10*dot13*dot18*dot5-8.*dot12*dot13*dot18*dot5+8.*a1*dot12*dot13*dot18*dot5+4.*dot10*dot26*dot5-4.*a1*dot10*dot26*dot5+4.*dot12*dot26*dot5-4.*a1*dot12*dot26*dot5+4.*dot10*dot27*dot5-4.*a1*dot10*dot27*dot5+8.*dot12*dot27*dot5-4.*a1*dot12*dot27*dot5-8.*dot15*dot18*dot3*dot5+4.*dot28*dot3*dot5-4.*dot13*dot30*dot5-2.*dot42*dot5-4.*a1*dot11*dot22*dot6-4.*dot13*dot22*dot6+4.*a1*dot13*dot22*dot6-4.*a1*dot11*dot23*dot6+4.*a1*dot13*dot23*dot6+4.*dot14*dot26*dot6-8.*a1*dot14*dot26*dot6-8.*a1*dot14*dot27*dot6-8.*dot17*dot18*dot3*dot6+4.*dot25*dot3*dot6-8.*dot14*dot18*dot4*dot6+8.*a1*dot14*dot18*dot4*dot6+8.*dot13*dot18*dot5*dot6-8.*a1*dot13*dot18*dot5*dot6+4.*a1*dot26*dot5*dot6+4.*a1*dot27*dot5*dot6-8.*dot10*dot13*dot18*dot7+8.*a1*dot10*dot13*dot18*dot7-8.*dot12*dot13*dot18*dot7+8.*a1*dot12*dot13*dot18*dot7+4.*dot10*dot26*dot7-4.*a1*dot10*dot26*dot7-4.*a1*dot12*dot26*dot7+4.*dot10*dot27*dot7-4.*a1*dot10*dot27*dot7+4.*dot12*dot27*dot7-4.*a1*dot12*dot27*dot7+4.*dot13*dot29*dot7-8.*dot15*dot18*dot3*dot7+4.*dot28*dot3*dot7+8.*dot12*dot18*dot4*dot7-2.*dot42*dot7-8.*a1*dot13*dot18*dot6*dot7+4.*a1*dot26*dot6*dot7+4.*a1*dot27*dot6*dot7+4.*dot1*(-2.*dot12*dot17*dot18+dot20+dot21-dot15*(dot22+dot23)+dot13*(2.*dot16*dot18-dot24)+dot12*dot25+dot16*(dot26+dot27)+dot14*(-2.*dot15*dot18+dot28)-dot17*dot29-dot17*dot30-2.*dot16*dot18*dot4+dot24*dot4+2.*dot15*dot18*dot5-dot28*dot5+2.*dot17*dot18*dot6-dot25*dot6+2.*dot15*dot18*dot7-dot28*dot7)+8.*dot14*dot15*dot18*dot8-8.*a1*dot14*dot15*dot18*dot8-8.*dot13*dot16*dot18*dot8+8.*a1*dot13*dot16*dot18*dot8-8.*dot10*dot17*dot18*dot8+8.*a1*dot10*dot17*dot18*dot8-6.*dot20*dot8+4.*a1*dot20*dot8-4.*dot21*dot8+4.*a1*dot21*dot8+8.*dot15*dot22*dot8-4.*a1*dot15*dot22*dot8+4.*dot15*dot23*dot8-4.*a1*dot15*dot23*dot8-4.*dot11*dot24*dot8+4.*a1*dot11*dot24*dot8+4.*dot13*dot24*dot8-4.*a1*dot13*dot24*dot8+8.*dot10*dot25*dot8-8.*a1*dot10*dot25*dot8-8.*dot16*dot26*dot8+4.*a1*dot16*dot26*dot8-4.*dot16*dot27*dot8+4.*a1*dot16*dot27*dot8-8.*dot14*dot28*dot8+8.*a1*dot14*dot28*dot8+8.*dot17*dot29*dot8-4.*a1*dot17*dot29*dot8+4.*dot17*dot30*dot8-4.*a1*dot17*dot30*dot8+16.*dot16*dot18*dot4*dot8-8.*a1*dot16*dot18*dot4*dot8-8.*dot24*dot4*dot8+4.*a1*dot24*dot4*dot8-16.*dot15*dot18*dot5*dot8+8.*a1*dot15*dot18*dot5*dot8+12.*dot28*dot5*dot8-8.*a1*dot28*dot5*dot8-8.*dot17*dot18*dot6*dot8+4.*dot25*dot6*dot8-8.*dot15*dot18*dot7*dot8+8.*a1*dot15*dot18*dot7*dot8+8.*dot28*dot7*dot8-8.*a1*dot28*dot7*dot8+8.*dot14*dot15*dot18*dot9-8.*a1*dot14*dot15*dot18*dot9-8.*dot13*dot16*dot18*dot9+8.*a1*dot13*dot16*dot18*dot9+8.*dot12*dot17*dot18*dot9-8.*a1*dot12*dot17*dot18*dot9-6.*dot20*dot9+4.*a1*dot20*dot9-4.*dot21*dot9+4.*a1*dot21*dot9+4.*dot15*dot22*dot9-4.*a1*dot15*dot22*dot9+4.*dot15*dot23*dot9-4.*a1*dot15*dot23*dot9-4.*dot11*dot24*dot9+4.*a1*dot11*dot24*dot9+4.*dot13*dot24*dot9-4.*a1*dot13*dot24*dot9-8.*dot12*dot25*dot9+8.*a1*dot12*dot25*dot9-4.*dot16*dot26*dot9+4.*a1*dot16*dot26*dot9-4.*dot16*dot27*dot9+4.*a1*dot16*dot27*dot9-8.*dot14*dot28*dot9+8.*a1*dot14*dot28*dot9+4.*dot17*dot29*dot9-4.*a1*dot17*dot29*dot9+4.*dot17*dot30*dot9-4.*a1*dot17*dot30*dot9+8.*dot16*dot18*dot4*dot9-8.*a1*dot16*dot18*dot4*dot9-4.*dot24*dot4*dot9+4.*a1*dot24*dot4*dot9-8.*dot15*dot18*dot5*dot9+8.*a1*dot15*dot18*dot5*dot9+8.*dot28*dot5*dot9-8.*a1*dot28*dot5*dot9-8.*dot17*dot18*dot6*dot9+8.*a1*dot17*dot18*dot6*dot9+4.*dot25*dot6*dot9-8.*a1*dot25*dot6*dot9-8.*dot15*dot18*dot7*dot9+8.*a1*dot15*dot18*dot7*dot9+8.*dot28*dot7*dot9-8.*a1*dot28*dot7*dot9-2.*a2*dot1*(2.*dot15*dot19+2.*dot17*dot39-2.*dot16*dot40-dot46)*M+(-dot45-2.*a2*dot10*(-2.*a2*(dot11+dot13)*dot19-2.*dot34-dot35+2.*dot17*dot37-2.*(-2.+a1)*dot19*dot4-2.*dot40*dot5)+2.*(2.*dot15*dot19*dot2-2.*a1*dot15*dot19*dot2+2.*dot15*dot19*dot3-2.*a1*dot15*dot19*dot3-dot13*dot32+a1*dot13*dot32+2.*dot13*dot16*dot37-2.*a1*dot13*dot16*dot37+2.*dot17*dot2*dot39-2.*a1*dot17*dot2*dot39+2.*dot17*dot3*dot39-2.*a1*dot17*dot3*dot39+dot32*dot4-a1*dot32*dot4-2.*dot16*dot37*dot4+2.*a1*dot16*dot37*dot4-2.*dot16*dot2*dot40+2.*a1*dot16*dot2*dot40-2.*dot16*dot3*dot40+2.*a1*dot16*dot3*dot40-dot2*dot46+a1*dot2*dot46-dot3*dot46+a1*dot3*dot46-2.*dot13*dot39*dot5+2.*a1*dot13*dot39*dot5+2.*dot39*dot4*dot5-2.*a1*dot39*dot4*dot5-a2*dot12*(-2.*a2*dot13*dot19-2.*dot34-dot35+2.*dot17*dot37-2.*(-2.+a1)*dot19*dot4-2.*dot40*dot5)+2.*a1*dot13*dot19*dot6-2.*a12*dot13*dot19*dot6-dot34*dot6+2.*a1*dot34*dot6+a1*dot35*dot6+2.*dot17*dot37*dot6-2.*a1*dot17*dot37*dot6+2.*dot19*dot4*dot6-4.*a1*dot19*dot4*dot6+2.*a12*dot19*dot4*dot6-2.*dot40*dot5*dot6+2.*a1*dot40*dot5*dot6-a2*dot11*(-2.*a2*dot12*dot19+dot32-2.*dot16*dot37+2.*dot39*dot5-2.*a1*dot19*dot6)+dot46*dot8-a1*dot46*dot8+dot46*dot9-a1*dot46*dot9))*M+2.*a2*dot46*pow<3,1>(M)-2.*(dot20+2.*dot21+2.*dot11*dot24-4.*a2*dot14*(dot15*dot18-dot28)-2.*(2.*dot12*dot17*dot18+dot15*dot23-a2*dot13*(2.*dot16*dot18-dot24)-2.*dot12*dot25-dot16*dot27+dot17*dot30+dot16*dot18*dot4-dot15*dot18*dot5+dot28*dot5+2.*a12*(dot17*dot18-dot25)*(dot10+dot12-dot6)-2.*dot17*dot18*dot6+dot25*dot6-2.*dot15*dot18*dot7+2.*dot28*dot7+a1*(dot21+dot15*dot22-dot15*dot23+2.*dot10*(-(dot17*dot18)+dot25)+4.*dot12*(-(dot17*dot18)+dot25)+dot17*(dot29-dot30)+dot24*(dot11-dot4)+dot16*(-dot26+dot27+dot18*dot4)-dot15*dot18*dot5+2.*dot17*dot18*dot6-2.*dot25*dot6+2.*dot15*dot18*dot7-2.*dot28*dot7)))*M2))/(a1*a2*(-dot1+dot2+dot3+dot8+dot9+M2)*(dot8+a1*M2)*(dot1-a2*(dot8+dot9+M2)));
	    // diag[25]=-0.0625*(M2*(-4.*dot13*dot16*dot18*dot2+4.*a1*dot13*dot16*dot18*dot2-dot2*dot20+a1*dot2*dot20-3.*dot2*dot21+a1*dot2*dot21-4.*dot10*dot11*dot22+4.*a1*dot10*dot11*dot22+4.*dot11*dot12*dot22-4.*a1*dot11*dot12*dot22-4.*dot10*dot13*dot22+4.*a1*dot10*dot13*dot22+4.*dot12*dot13*dot22-4.*a1*dot12*dot13*dot22+2.*dot15*dot2*dot22-2.*a1*dot15*dot2*dot22-4.*dot10*dot11*dot23+4.*a1*dot10*dot11*dot23+4.*a1*dot11*dot12*dot23+4.*dot10*dot13*dot23-4.*a1*dot10*dot13*dot23-4.*a1*dot12*dot13*dot23+2.*dot15*dot2*dot23-2.*a1*dot15*dot2*dot23+4.*dot11*dot2*dot24-4.*a1*dot11*dot2*dot24+6.*dot13*dot2*dot24-6.*a1*dot13*dot2*dot24-2.*dot10*dot2*dot25+2.*a1*dot10*dot2*dot25+4.*a1*dot12*dot2*dot25-2.*dot16*dot2*dot26+2.*a1*dot16*dot2*dot26-2.*dot16*dot2*dot27+2.*a1*dot16*dot2*dot27+8.*dot11*dot16*dot18*dot3-8.*a1*dot11*dot16*dot18*dot3+4.*dot13*dot16*dot18*dot3-4.*a1*dot13*dot16*dot18*dot3-8.*dot10*dot17*dot18*dot3+8.*a1*dot10*dot17*dot18*dot3+8.*a1*dot12*dot17*dot18*dot3+3.*dot20*dot3+a1*dot20*dot3+dot21*dot3+a1*dot21*dot3-2.*dot15*dot22*dot3-2.*a1*dot15*dot22*dot3-2.*dot15*dot23*dot3-2.*a1*dot15*dot23*dot3-4.*dot11*dot24*dot3+4.*a1*dot11*dot24*dot3-2.*dot13*dot24*dot3+2.*a1*dot13*dot24*dot3+6.*dot10*dot25*dot3-6.*a1*dot10*dot25*dot3-4.*a1*dot12*dot25*dot3+2.*dot16*dot26*dot3+2.*a1*dot16*dot26*dot3+2.*dot16*dot27*dot3+2.*a1*dot16*dot27*dot3-4.*dot11*dot16*dot31+4.*a1*dot11*dot16*dot31+4.*a1*dot13*dot16*dot31+4.*dot10*dot17*dot31-4.*a1*dot10*dot17*dot31-4.*a1*dot12*dot17*dot31+4.*dot16*dot18*dot2*dot4-4.*a1*dot16*dot18*dot2*dot4+4.*dot10*dot22*dot4-4.*a1*dot10*dot22*dot4-4.*dot12*dot22*dot4+4.*a1*dot12*dot22*dot4-4.*dot10*dot23*dot4+4.*a1*dot10*dot23*dot4+4.*dot12*dot23*dot4+4.*a1*dot12*dot23*dot4-6.*dot2*dot24*dot4+6.*a1*dot2*dot24*dot4-4.*dot16*dot18*dot3*dot4+4.*a1*dot16*dot18*dot3*dot4+6.*dot24*dot3*dot4-2.*a1*dot24*dot3*dot4-4.*a1*dot16*dot31*dot4-6.*dot10*dot43+4.*a1*dot10*dot43+2.*dot12*dot43+4.*a1*dot12*dot43+4.*dot11*dot44-4.*a1*dot11*dot44-4.*a1*dot13*dot44-2.*dot4*dot44+4.*a1*dot4*dot44-8.*dot11*dot12*dot18*dot5+8.*a1*dot11*dot12*dot18*dot5-8.*dot12*dot13*dot18*dot5+8.*a1*dot12*dot13*dot18*dot5-4.*dot15*dot18*dot2*dot5+4.*a1*dot15*dot18*dot2*dot5+8.*dot10*dot27*dot5-8.*a1*dot10*dot27*dot5-4.*dot12*dot27*dot5-8.*a1*dot12*dot27*dot5+2.*dot2*dot28*dot5-2.*a1*dot2*dot28*dot5+4.*dot11*dot29*dot5-4.*a1*dot11*dot29*dot5+4.*dot13*dot29*dot5-4.*a1*dot13*dot29*dot5-4.*dot15*dot18*dot3*dot5+4.*a1*dot15*dot18*dot3*dot5-2.*dot28*dot3*dot5-2.*a1*dot28*dot3*dot5-4.*dot11*dot30*dot5+4.*a1*dot11*dot30*dot5+4.*a1*dot13*dot30*dot5+4.*dot15*dot31*dot5+8.*dot12*dot18*dot4*dot5-8.*a1*dot12*dot18*dot4*dot5-4.*dot29*dot4*dot5+4.*a1*dot29*dot4*dot5+4.*dot30*dot4*dot5-4.*a1*dot30*dot4*dot5-2.*dot42*dot5-4.*dot11*dot22*dot6+4.*a1*dot11*dot22*dot6+4.*a1*dot13*dot22*dot6-4.*a1*dot11*dot23*dot6+4.*a1*dot13*dot23*dot6+4.*dot2*dot25*dot6-4.*a1*dot2*dot25*dot6-8.*a1*dot17*dot18*dot3*dot6+4.*a1*dot25*dot3*dot6+4.*a1*dot17*dot31*dot6+4.*dot22*dot4*dot6-4.*a1*dot22*dot4*dot6-4.*a1*dot23*dot4*dot6-4.*a1*dot43*dot6+8.*dot11*dot18*dot5*dot6-8.*a1*dot11*dot18*dot5*dot6-8.*a1*dot13*dot18*dot5*dot6+8.*a1*dot27*dot5*dot6-8.*dot18*dot4*dot5*dot6+8.*a1*dot18*dot4*dot5*dot6-8.*dot10*dot13*dot18*dot7+8.*a1*dot10*dot13*dot18*dot7+8.*a1*dot12*dot13*dot18*dot7-4.*dot15*dot18*dot2*dot7+4.*a1*dot15*dot18*dot2*dot7-4.*dot10*dot26*dot7+4.*a1*dot10*dot26*dot7+4.*a1*dot12*dot26*dot7+4.*dot10*dot27*dot7-4.*a1*dot10*dot27*dot7-4.*a1*dot12*dot27*dot7+2.*dot2*dot28*dot7-2.*a1*dot2*dot28*dot7+8.*dot11*dot29*dot7-8.*a1*dot11*dot29*dot7+8.*dot13*dot29*dot7-8.*a1*dot13*dot29*dot7-4.*dot15*dot18*dot3*dot7+4.*a1*dot15*dot18*dot3*dot7+2.*dot28*dot3*dot7-2.*a1*dot28*dot3*dot7+4.*dot15*dot31*dot7+8.*dot10*dot18*dot4*dot7-8.*a1*dot10*dot18*dot4*dot7-8.*a1*dot12*dot18*dot4*dot7-8.*dot29*dot4*dot7+8.*a1*dot29*dot4*dot7+4.*dot30*dot4*dot7-4.*dot42*dot7-8.*a1*dot13*dot18*dot6*dot7+4.*dot26*dot6*dot7-4.*a1*dot26*dot6*dot7+4.*a1*dot27*dot6*dot7-8.*dot18*dot4*dot6*dot7+8.*a1*dot18*dot4*dot6*dot7+dot1*(-4.*a2*dot14*dot15*dot18-dot20-a1*dot20+dot21-a1*dot21-2.*dot15*dot22+2.*a1*dot15*dot22-2.*dot15*dot23+2.*a1*dot15*dot23-4.*dot11*dot24+4.*a1*dot11*dot24-2.*dot13*(-2.*a2*dot16*dot18+dot24-3.*a1*dot24)+2.*dot10*dot25-2.*a1*dot10*dot25-4.*dot12*dot25-4.*a1*dot12*dot25+2.*dot16*dot26-2.*a1*dot16*dot26+2.*dot16*dot27-2.*a1*dot16*dot27-2.*(1.+a1)*dot14*dot28-4.*dot16*dot18*dot4+4.*a1*dot16*dot18*dot4+2.*dot24*dot4-6.*a1*dot24*dot4+4.*dot15*dot18*dot5-4.*a1*dot15*dot18*dot5+2.*dot28*dot5+2.*a1*dot28*dot5+4.*a1*dot25*dot6+4.*dot15*dot18*dot7-4.*a1*dot15*dot18*dot7+2.*dot28*dot7+2.*a1*dot28*dot7)-4.*dot13*dot16*dot18*dot8+4.*a1*dot13*dot16*dot18*dot8-8.*dot10*dot17*dot18*dot8+8.*a1*dot10*dot17*dot18*dot8+8.*a1*dot12*dot17*dot18*dot8+3.*dot20*dot8+a1*dot20*dot8-dot21*dot8+a1*dot21*dot8-2.*dot15*dot22*dot8-2.*a1*dot15*dot22*dot8+2.*dot15*dot23*dot8-2.*a1*dot15*dot23*dot8+2.*dot13*dot24*dot8-2.*a1*dot13*dot24*dot8+6.*dot10*dot25*dot8-6.*a1*dot10*dot25*dot8-4.*a1*dot12*dot25*dot8+2.*dot16*dot26*dot8+2.*a1*dot16*dot26*dot8-2.*dot16*dot27*dot8+2.*a1*dot16*dot27*dot8+4.*dot16*dot18*dot4*dot8-4.*a1*dot16*dot18*dot4*dot8+2.*dot24*dot4*dot8+2.*a1*dot24*dot4*dot8-4.*dot15*dot18*dot5*dot8+4.*a1*dot15*dot18*dot5*dot8-2.*dot28*dot5*dot8-2.*a1*dot28*dot5*dot8-8.*a1*dot17*dot18*dot6*dot8+4.*a1*dot25*dot6*dot8-4.*dot15*dot18*dot7*dot8+4.*a1*dot15*dot18*dot7*dot8+2.*dot28*dot7*dot8-2.*a1*dot28*dot7*dot8-4.*dot13*dot16*dot18*dot9+4.*a1*dot13*dot16*dot18*dot9-dot20*dot9+a1*dot20*dot9-dot21*dot9+a1*dot21*dot9+2.*dot15*dot22*dot9-2.*a1*dot15*dot22*dot9+2.*dot15*dot23*dot9-2.*a1*dot15*dot23*dot9+2.*dot13*dot24*dot9-2.*a1*dot13*dot24*dot9-2.*dot10*dot25*dot9+2.*a1*dot10*dot25*dot9+4.*a1*dot12*dot25*dot9-2.*dot16*dot26*dot9+2.*a1*dot16*dot26*dot9-2.*dot16*dot27*dot9+2.*a1*dot16*dot27*dot9+4.*dot16*dot18*dot4*dot9-4.*a1*dot16*dot18*dot4*dot9-2.*dot24*dot4*dot9+2.*a1*dot24*dot4*dot9-4.*dot15*dot18*dot5*dot9+4.*a1*dot15*dot18*dot5*dot9+2.*dot28*dot5*dot9-2.*a1*dot28*dot5*dot9-4.*a1*dot25*dot6*dot9-4.*dot15*dot18*dot7*dot9+4.*a1*dot15*dot18*dot7*dot9+2.*dot28*dot7*dot9-2.*a1*dot28*dot7*dot9+2.*(-dot35-2.*a2*dot19*(dot11+dot13-dot4))*(-(a2*dot10)+a1*(dot12-dot6))*M+(-dot20-dot21-2.*a2*dot13*(2.*dot16*dot18-dot24)+8.*a12*(dot17*dot18-dot25)*(dot10+dot12-dot6)-2.*(dot10*dot25+dot24*dot4+dot16*(dot26+dot27-2.*dot18*dot4)-dot28*(dot5+dot7)-dot15*(dot22+dot23-2.*dot18*(dot5+dot7)))+a1*(5.*dot20+dot21+dot10*(-8.*dot17*dot18+10.*dot25)+6.*dot24*dot4+2.*dot16*(3.*dot26+dot27-2.*dot18*dot4)+4.*dot25*(dot12-dot6)-2.*(dot28*(3.*dot5+dot7)+dot15*(3.*dot22+dot23-2.*dot18*(dot5+dot7)))))*M2-2.*dot14*(4.*dot10*dot27+4.*dot11*dot29+4.*dot13*dot29-4.*dot29*dot4+2.*dot30*dot4-2.*dot42+2.*dot26*dot6-4.*dot18*dot4*dot6+dot28*(dot2+dot3+dot8+dot9+M2)+2.*dot15*(dot31-a2*dot18*(dot2+dot3+dot8+dot9+M2))-a1*(4.*dot29*(dot11+dot13-dot4)+4.*dot27*(dot10+dot12-dot6)+dot28*(dot2+dot3+dot8+dot9+M2)))))/(a1*a2*(-dot1+dot2+dot3+dot8+dot9+M2)*(dot8+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot8+a1*M2)));
	    // diag[26]=((-2.*dot13*dot32+2.*dot10*(2.*dot13*dot19+dot35)+4.*dot13*dot16*dot37-2.*dot33*dot4+4.*dot16*dot38*dot4+2.*dot16*dot41+dot45-2.*dot1*(-2.*dot15*dot19+2.*dot16*dot40+dot46)-4.*dot13*dot39*dot5+2.*dot47*dot5-2.*dot15*(dot36+2.*dot38*dot5))*pow<3,1>(M))/(8.*dot2*(-dot1+dot3+dot8+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot8+a1*M2)));
	    // diag[27]=-0.125*((-2.*dot11*(dot32-2.*dot16*dot37)+2.*(dot32+dot33-2.*dot16*(dot37+dot38))*dot4-dot45+2.*dot1*(-2.*dot15*dot19-2.*dot17*dot39+2.*dot16*dot40+dot46)+2.*(-(dot16*(2.*dot2*dot40+dot41))-dot2*dot46+dot17*(2.*dot2*dot39+dot48)-(2.*dot11*dot39-2.*dot39*dot4+dot47)*dot5+(2.*dot39*dot4+dot49)*dot7+dot15*(2.*dot19*dot2+dot36+2.*dot38*dot5-2.*dot37*dot7))+dot10*(4.*dot11*dot19-2.*(dot35-2.*dot17*dot38+2.*dot19*dot4+2.*dot40*dot7)))*pow<3,1>(M))/(a2*dot2*sqr(dot9+a2*M2));
	    // diag[28]=-0.125*(M2*(-4.*dot14*dot15*dot18*dot2+4.*dot13*dot16*dot18*dot2-2.*dot1*dot20-dot2*dot20-3.*dot2*dot21+2.*dot15*dot2*dot22+4.*dot10*dot11*dot23-4.*dot10*dot13*dot23+6.*dot15*dot2*dot23-2.*dot13*dot2*dot24+2.*dot10*dot2*dot25-2.*dot16*dot2*dot26+8.*dot10*dot14*dot27-6.*dot16*dot2*dot27+2.*dot14*dot2*dot28+8.*dot10*dot17*dot18*dot3+2.*dot20*dot3-4.*dot15*dot22*dot3-4.*dot10*dot25*dot3+4.*dot16*dot26*dot3-4.*dot17*dot29*dot3+4.*dot17*dot2*dot30-8.*dot14*dot15*dot31-4.*dot11*dot16*dot31+4.*dot13*dot16*dot31+4.*dot16*dot18*dot2*dot4-4.*dot10*dot23*dot4-2.*dot2*dot24*dot4-8.*dot16*dot18*dot3*dot4+4.*dot24*dot3*dot4-8.*dot14*dot30*dot4+4.*dot16*dot31*dot4+4.*dot14*dot42+2.*dot11*dot44-2.*dot13*dot44-2.*dot4*dot44-4.*dot15*dot18*dot2*dot5+2.*dot2*dot28*dot5+8.*dot15*dot18*dot3*dot5-4.*dot28*dot3*dot5-4.*dot11*dot30*dot5+4.*dot13*dot30*dot5+4.*dot30*dot4*dot5+8.*dot10*dot13*dot18*dot7-4.*dot15*dot18*dot2*dot7-4.*dot10*dot27*dot7+2.*dot2*dot28*dot7-4.*dot13*dot29*dot7+4.*dot15*dot31*dot7+4.*dot30*dot4*dot7-2.*dot42*dot7+4.*dot1*(4.*dot14*dot15*dot18+2.*dot11*dot16*dot18-2.*dot13*dot16*dot18-2.*dot10*dot17*dot18+dot15*dot22-dot11*dot24+dot13*dot24+dot10*dot25-dot16*dot26-2.*dot14*dot28+dot17*dot29-2.*dot15*dot18*dot5+dot28*dot5-2.*dot15*dot18*dot7+dot28*dot7)+8.*dot10*dot17*dot18*dot8+2.*dot20*dot8-4.*dot15*dot22*dot8-4.*dot10*dot25*dot8+4.*dot16*dot26*dot8-4.*dot17*dot29*dot8-8.*dot16*dot18*dot4*dot8+4.*dot24*dot4*dot8+8.*dot15*dot18*dot5*dot8-4.*dot28*dot5*dot8-2.*dot20*dot9+4.*dot15*dot22*dot9+4.*dot10*dot25*dot9-4.*dot16*dot26*dot9+8.*dot16*dot18*dot4*dot9-4.*dot24*dot4*dot9-8.*dot15*dot18*dot5*dot9+4.*dot28*dot5*dot9+2.*(-(a2*dot11)+a1*dot13+dot4-a1*dot4)*(2.*dot10*dot19-dot32+2.*dot16*dot37-2.*dot39*dot5)*M+2.*(-dot20+2.*dot15*dot22+2.*dot10*dot25-2.*dot16*dot26+4.*dot16*dot18*dot4-2.*dot24*dot4-4.*dot15*dot18*dot5+2.*dot28*dot5+2.*a1*(2.*dot10*dot17*dot18+dot20-2.*dot15*dot22-2.*dot10*dot25+2.*dot16*dot26-dot17*dot29-4.*dot16*dot18*dot4+2.*dot24*dot4+4.*dot15*dot18*dot5-2.*dot28*dot5))*M2))/(a2*dot2*(-dot9-a2*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot8+a1*M2)));
	   diag[29]=-0.125*(2.*dot11*(2.*dot12*dot19+dot32)+dot45-2.*dot12*(dot35-2.*dot17*dot38+2.*dot40*dot7)-
			    2.*(dot17*dot48+dot34*dot6-2.*dot17*dot37*dot6+dot49*dot7))*pow<3,1>(M)
	     /(dot3*(-(a2*dot1)+dot2+dot9-a1*(dot2+dot3+dot9)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot9-a2*M2));
	   diag[30]=-0.125*pow<3,1>(M)*(-dot45+2.*dot17*(-2.*dot1*dot39+2.*dot2*dot39+dot48)+2.*dot49*dot7+
					dot10*(4.*dot11*dot19-2.*dot35+4.*dot17*dot38-4.*dot40*dot7)-
					2.*(2.*dot16*dot2*dot40+dot16*dot41+dot2*dot46-dot1*(2.*dot16*dot40+dot46)+
					    dot47*dot5-(dot34+dot35-2.*dot17*(dot37+dot38))*dot6-2.*dot40*dot5*dot6+
					    dot11*(dot32-2.*dot16*dot37+2.*dot39*dot5+2.*dot19*dot6)-
					    2.*dot40*dot6*dot7-2.*dot17*dot39*dot9+2.*dot16*dot40*dot9+dot46*dot9)-
					2.*a2*(-2.*dot17*dot39+2.*dot16*dot40+dot46)*M2)/(a2*(dot1-dot2-dot9-a2*M2)*sqr(dot9+a2*M2));
	    diag[31]=M2*(-4.*dot12*dot17*dot18*dot2-3.*dot2*dot20-dot2*dot21+4.*dot10*dot11*dot22-4.*dot11*dot12*dot22-
			 2.*dot11*dot2*dot24+2.*dot12*dot2*dot25-4.*dot16*dot2*dot26+2.*dot14*dot2*dot28+
			 8.*dot11*dot14*dot29+6.*dot17*dot2*dot29-8.*dot11*dot16*dot18*dot3+2.*dot21*dot3+
			 4.*dot11*dot24*dot3+4.*dot16*dot27*dot3+2.*dot17*dot2*dot30-4.*dot17*dot3*dot30+
			 4.*dot10*dot17*dot31-4.*dot12*dot17*dot31+4.*dot14*dot42-2.*dot10*dot43+2.*dot12*dot43+
			 8.*dot11*dot12*dot18*dot5-4.*dot12*dot27*dot5+2.*dot2*dot28*dot5-4.*dot11*dot29*dot5-
			 2.*dot42*dot5-4.*dot17*dot18*dot2*dot6-4.*dot11*dot22*dot6+2.*dot2*dot25*dot6-
			 8.*dot14*dot26*dot6+8.*dot17*dot18*dot3*dot6-4.*dot25*dot3*dot6-4.*dot17*dot31*dot6+
			 2.*dot43*dot6+4.*dot26*dot5*dot6-4.*dot10*dot26*dot7+4.*dot12*dot26*dot7+
			 2.*dot2*dot28*dot7-4.*dot28*dot3*dot7+4.*dot26*dot6*dot7+
			 dot1*(3.*dot20+dot21+2.*dot11*dot24+dot12*(4.*dot17*dot18-2.*dot25)-
			       2.*(-2.*dot16*dot26+dot25*dot6+dot17*(3.*dot29+dot30-2.*dot18*dot6)+dot28*(dot14+dot5+dot7)))-
			 4.*dot12*dot17*dot18*dot9-3.*dot20*dot9+dot21*dot9+2.*dot11*dot24*dot9+
			 2.*dot12*dot25*dot9-4.*dot16*dot26*dot9+2.*dot14*dot28*dot9+6.*dot17*dot29*dot9-
			 2.*dot17*dot30*dot9+2.*dot28*dot5*dot9+4.*dot17*dot18*dot6*dot9-2.*dot25*dot6*dot9-2.*dot28*dot7*dot9+
			 2.*(-(a2*dot10)+a1*dot12+dot6-a1*dot6)*(2.*dot11*dot19-dot35+2.*dot17*dot38-2.*dot40*dot7)*M-
			 a2*(3.*dot20-dot21-2.*dot11*dot24+dot12*(4.*dot17*dot18-2.*dot25)+4.*dot16*dot26-
			     6.*dot17*dot29+2.*dot17*dot30-4.*dot17*dot18*dot6+2.*dot25*dot6-2.*dot28*(dot14+dot5-dot7))*M2)
	      /(8.*a2*(-(a2*dot1)+dot2+dot9-a1*(dot2+dot3+dot9)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot9-a2*M2)*(dot9+a2*M2));
	   diag[32]=-0.125*(4.*dot12*dot19*(dot13-dot4)-dot45+
			    2.*(dot13*dot32+dot12*dot35-dot32*dot4+dot3*dot46+dot34*dot6))*pow<3,1>(M)/(a1*dot3*sqr(dot8+a1*M2));
	   diag[33]=M2*(-2.*dot2*dot20+4.*dot12*(-dot11+dot13)*dot23-4.*dot12*dot2*dot25+2.*dot1*(dot20+2.*dot12*dot25)-
			8.*dot12*dot14*dot27-4.*dot12*dot17*dot18*dot3+dot20*dot3+3.*dot21*dot3-2.*dot11*dot24*dot3+
			2.*dot12*dot25*dot3-2.*dot14*dot28*dot3-2.*dot17*dot29*dot3-2.*dot17*dot3*dot30-4.*dot12*dot23*dot4+
			4.*dot14*dot42+2.*dot11*dot44-2.*dot13*dot44+2.*dot4*dot44+8.*dot12*dot27*dot5+2.*dot28*dot3*dot5-
			4.*dot42*dot5-4.*dot11*dot22*dot6+4.*dot13*dot22*dot6-8.*dot14*dot26*dot6+4.*dot17*dot18*dot3*dot6-
			6.*dot25*dot3*dot6-4.*dot22*dot4*dot6+8.*dot26*dot5*dot6-8.*dot12*dot13*dot18*dot7+4.*dot12*dot27*dot7-
			4.*dot13*dot29*dot7-2.*dot28*dot3*dot7+8.*dot12*dot18*dot4*dot7+4.*dot29*dot4*dot7-2.*dot42*dot7+
			4.*dot26*dot6*dot7-8.*dot12*dot17*dot18*dot8+2.*dot20*dot8+4.*dot12*dot25*dot8-4.*dot17*dot29*dot8-
			2.*dot20*dot9-4.*dot12*dot25*dot9-2.*(2.*dot12*dot19+dot32)*(-(a2*dot11)+a1*(dot13-dot4))*M-
			2.*(dot20+2.*dot12*dot25+2.*a1*(2.*dot12*dot17*dot18-dot20-2.*dot12*dot25+dot17*dot29))*M2)
	     /(8.*a1*dot3*(dot8+a1*M2)*(-(a2*dot1)+dot2+dot9-a1*(dot2+dot3+dot9)+M2+(-2.+a1)*a1*M2));
	   diag[34]=-0.125*pow<3,1>(M)*(dot45+4.*dot19*(-(dot12*dot13)+dot13*dot6+dot15*(dot3+dot8+a1*M2))-
					2.*(dot13*dot32+dot12*dot35+dot15*dot36+dot33*dot4-dot35*dot6+dot46*(dot3+dot8+a1*M2)))/(a1*sqr(dot8+a1*M2)*(-dot1+dot3+dot8+a1*M2));
	   // diag[35]=(M2*(4.*(dot10-dot12)*dot13*dot22+2.*dot2*(dot21-2.*(dot15*dot23+dot13*dot24))+8.*dot13*dot14*dot29-4.*dot14*dot15*dot18*dot3-4.*dot13*dot16*dot18*dot3-3.*dot20*dot3-dot21*dot3+6.*dot15*dot22*dot3+2.*dot15*dot23*dot3+2.*dot13*dot24*dot3-2.*dot10*dot25*dot3-2.*dot16*dot26*dot3-2.*dot16*dot27*dot3+2.*dot14*dot28*dot3+8.*dot14*dot15*dot31+4.*dot10*dot23*dot4-4.*dot12*dot23*dot4+4.*dot16*dot18*dot3*dot4-6.*dot24*dot3*dot4+8.*dot14*dot30*dot4-4.*dot14*dot42+2.*dot10*dot43-2.*dot12*dot43+8.*dot12*dot13*dot18*dot5+4.*dot12*dot27*dot5-4.*dot13*dot29*dot5-4.*dot15*dot18*dot3*dot5+2.*dot28*dot3*dot5-4.*dot15*dot31*dot5-4.*dot30*dot4*dot5+2.*dot42*dot5+4.*dot13*dot22*dot6+4.*dot23*dot4*dot6+2.*dot43*dot6-8.*dot13*dot18*dot5*dot6-4.*dot27*dot5*dot6-8.*dot13*dot29*dot7+4.*dot15*dot18*dot3*dot7-2.*dot28*dot3*dot7-8.*dot15*dot31*dot7-8.*dot30*dot4*dot7+4.*dot42*dot7+dot1*(3.*dot20+dot21+dot13*(4.*dot16*dot18-2.*dot24)+6.*dot14*(-2.*dot15*dot18+dot28)+6.*dot24*dot4+2.*dot16*(dot26+dot27-2.*dot18*dot4)-2.*dot25*(dot10-2.*dot12+2.*dot6)-6.*dot28*(dot5+dot7)-2.*dot15*(3.*dot22+dot23-6.*dot18*(dot5+dot7)))-4.*dot14*dot15*dot18*dot8+4.*dot13*dot16*dot18*dot8-3.*dot20*dot8+dot21*dot8+6.*dot15*dot22*dot8-2.*dot15*dot23*dot8-2.*dot13*dot24*dot8-2.*dot10*dot25*dot8-2.*dot16*dot26*dot8+2.*dot16*dot27*dot8+2.*dot14*dot28*dot8+4.*dot16*dot18*dot4*dot8-6.*dot24*dot4*dot8-4.*dot15*dot18*dot5*dot8+2.*dot28*dot5*dot8+4.*dot15*dot18*dot7*dot8-2.*dot28*dot7*dot8+2.*(2.*dot13*dot19+dot35)*(-(a2*dot10)+a1*(dot12-dot6))*M+a1*(-3.*dot20+dot21+dot13*(4.*dot16*dot18-2.*dot24)+2.*dot14*(-2.*dot15*dot18+dot28)-2.*(-3.*dot15*dot22+dot15*dot23+dot10*dot25+dot16*dot26-dot16*dot27-2.*dot16*dot18*dot4+3.*dot24*dot4+2.*dot15*dot18*dot5-dot28*dot5-2.*dot15*dot18*dot7+dot28*dot7))*M2))/(8.*a1*(dot8+a1*M2)*(-dot1+dot3+dot8+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot8+a1*M2)));
	    // diag[36]=(M2*(4.*dot13*dot16*dot18*dot2-4.*a1*dot13*dot16*dot18*dot2+dot2*dot20-a1*dot2*dot20+dot2*dot21-a1*dot2*dot21+4.*dot10*dot11*dot22-4.*a1*dot10*dot11*dot22-4.*dot11*dot12*dot22+4.*a1*dot11*dot12*dot22-4.*a1*dot10*dot13*dot22+4.*a1*dot12*dot13*dot22-2.*dot15*dot2*dot22+2.*a1*dot15*dot2*dot22+4.*dot10*dot11*dot23-4.*a1*dot10*dot11*dot23-4.*a1*dot11*dot12*dot23-4.*dot10*dot13*dot23+4.*a1*dot10*dot13*dot23+4.*a1*dot12*dot13*dot23+2.*dot15*dot2*dot23+2.*a1*dot15*dot2*dot23-4.*dot11*dot2*dot24+4.*a1*dot11*dot2*dot24-2.*dot13*dot2*dot24+6.*a1*dot13*dot2*dot24+2.*dot10*dot2*dot25-2.*a1*dot10*dot2*dot25-4.*a1*dot12*dot2*dot25+2.*dot16*dot2*dot26-2.*a1*dot16*dot2*dot26+2.*dot16*dot2*dot27-2.*a1*dot16*dot2*dot27-8.*dot11*dot16*dot18*dot3+8.*a1*dot11*dot16*dot18*dot3+4.*a1*dot13*dot16*dot18*dot3+8.*dot10*dot17*dot18*dot3-8.*a1*dot10*dot17*dot18*dot3-8.*a1*dot12*dot17*dot18*dot3-a1*dot20*dot3-a1*dot21*dot3-4.*dot15*dot22*dot3+2.*a1*dot15*dot22*dot3+2.*a1*dot15*dot23*dot3+4.*dot11*dot24*dot3-4.*a1*dot11*dot24*dot3-2.*a1*dot13*dot24*dot3-4.*dot10*dot25*dot3+6.*a1*dot10*dot25*dot3+4.*a1*dot12*dot25*dot3-2.*a1*dot16*dot26*dot3-2.*a1*dot16*dot27*dot3+4.*dot11*dot16*dot31-4.*a1*dot11*dot16*dot31-4.*a1*dot13*dot16*dot31-4.*dot10*dot17*dot31+4.*a1*dot10*dot17*dot31+4.*a1*dot12*dot17*dot31-4.*dot16*dot18*dot2*dot4+4.*a1*dot16*dot18*dot2*dot4-4.*dot10*dot22*dot4+4.*a1*dot10*dot22*dot4+4.*dot12*dot22*dot4-4.*a1*dot12*dot22*dot4-4.*a1*dot10*dot23*dot4-4.*a1*dot12*dot23*dot4+6.*dot2*dot24*dot4-6.*a1*dot2*dot24*dot4-4.*a1*dot16*dot18*dot3*dot4+2.*a1*dot24*dot3*dot4+4.*a1*dot16*dot31*dot4+4.*dot10*dot43-4.*a1*dot10*dot43-4.*a1*dot12*dot43-4.*dot11*dot44+4.*a1*dot11*dot44+4.*a1*dot13*dot44+2.*dot4*dot44-4.*a1*dot4*dot44+8.*dot11*dot12*dot18*dot5-8.*a1*dot11*dot12*dot18*dot5-8.*a1*dot12*dot13*dot18*dot5+4.*dot15*dot18*dot2*dot5-4.*a1*dot15*dot18*dot2*dot5-8.*dot10*dot27*dot5+8.*a1*dot10*dot27*dot5+8.*a1*dot12*dot27*dot5-2.*dot2*dot28*dot5+2.*a1*dot2*dot28*dot5-4.*dot11*dot29*dot5+4.*a1*dot11*dot29*dot5+4.*a1*dot13*dot29*dot5+8.*dot15*dot18*dot3*dot5-4.*a1*dot15*dot18*dot3*dot5+2.*a1*dot28*dot3*dot5+4.*dot11*dot30*dot5-4.*a1*dot11*dot30*dot5-4.*a1*dot13*dot30*dot5-8.*dot12*dot18*dot4*dot5+8.*a1*dot12*dot18*dot4*dot5+4.*dot29*dot4*dot5-4.*a1*dot29*dot4*dot5+4.*a1*dot30*dot4*dot5+4.*dot11*dot22*dot6-4.*a1*dot11*dot22*dot6-4.*dot13*dot22*dot6-4.*a1*dot13*dot22*dot6+4.*a1*dot11*dot23*dot6-4.*a1*dot13*dot23*dot6-4.*dot2*dot25*dot6+4.*a1*dot2*dot25*dot6+8.*a1*dot17*dot18*dot3*dot6-4.*a1*dot25*dot3*dot6-4.*a1*dot17*dot31*dot6-4.*dot22*dot4*dot6+4.*a1*dot22*dot4*dot6-4.*dot23*dot4*dot6+4.*a1*dot23*dot4*dot6-2.*dot43*dot6+4.*a1*dot43*dot6-8.*dot11*dot18*dot5*dot6+8.*a1*dot11*dot18*dot5*dot6+8.*dot13*dot18*dot5*dot6+8.*a1*dot13*dot18*dot5*dot6+4.*dot27*dot5*dot6-8.*a1*dot27*dot5*dot6+8.*dot18*dot4*dot5*dot6-8.*a1*dot18*dot4*dot5*dot6+8.*dot10*dot13*dot18*dot7-8.*a1*dot10*dot13*dot18*dot7-8.*a1*dot12*dot13*dot18*dot7+4.*dot15*dot18*dot2*dot7-4.*a1*dot15*dot18*dot2*dot7+4.*dot10*dot26*dot7-4.*a1*dot10*dot26*dot7-4.*a1*dot12*dot26*dot7-4.*dot10*dot27*dot7+4.*a1*dot10*dot27*dot7+4.*a1*dot12*dot27*dot7-2.*dot2*dot28*dot7+2.*a1*dot2*dot28*dot7-8.*dot11*dot29*dot7+8.*a1*dot11*dot29*dot7+8.*a1*dot13*dot29*dot7-4.*a1*dot15*dot18*dot3*dot7+2.*a1*dot28*dot3*dot7+4.*dot15*dot31*dot7-8.*dot10*dot18*dot4*dot7+8.*a1*dot10*dot18*dot4*dot7+8.*a1*dot12*dot18*dot4*dot7+8.*dot29*dot4*dot7-8.*a1*dot29*dot4*dot7+4.*dot30*dot4*dot7+8.*a1*dot13*dot18*dot6*dot7-4.*dot26*dot6*dot7+4.*a1*dot26*dot6*dot7-4.*a1*dot27*dot6*dot7+8.*dot18*dot4*dot6*dot7-8.*a1*dot18*dot4*dot6*dot7+dot1*(-4.*(-4.+a1)*dot14*dot15*dot18+4.*(-2.+a1)*dot13*dot16*dot18+a1*dot20+a1*dot21+8.*dot15*dot22-2.*a1*dot15*dot22+4.*dot15*dot23-2.*a1*dot15*dot23-4.*a1*dot11*dot24+(4.-6.*a1)*dot13*dot24-2.*(dot20+dot21-2.*dot11*dot24)+2.*a1*dot10*dot25+4.*a1*dot12*dot25-4.*dot16*dot26+2.*a1*dot16*dot26-4.*dot16*dot27+2.*a1*dot16*dot27+2.*(-2.+a1)*dot14*dot28+8.*dot16*dot18*dot4-4.*a1*dot16*dot18*dot4-8.*dot24*dot4+6.*a1*dot24*dot4-16.*dot15*dot18*dot5+4.*a1*dot15*dot18*dot5+4.*dot28*dot5-2.*a1*dot28*dot5+4.*dot25*dot6-4.*a1*dot25*dot6+4.*(-4.+a1)*dot15*dot18*dot7-2.*(-2.+a1)*dot28*dot7)+8.*dot10*dot17*dot18*dot8-8.*a1*dot10*dot17*dot18*dot8-8.*a1*dot12*dot17*dot18*dot8-4.*dot15*dot22*dot8-4.*dot10*dot25*dot8+4.*a1*dot10*dot25*dot8+4.*a1*dot12*dot25*dot8-8.*dot16*dot18*dot4*dot8+4.*dot24*dot4*dot8+8.*dot15*dot18*dot5*dot8+8.*a1*dot17*dot18*dot6*dot8-4.*a1*dot25*dot6*dot8+4.*dot13*dot16*dot18*dot9-4.*a1*dot13*dot16*dot18*dot9+dot20*dot9-a1*dot20*dot9+dot21*dot9-a1*dot21*dot9-2.*dot15*dot22*dot9+2.*a1*dot15*dot22*dot9-2.*dot15*dot23*dot9+2.*a1*dot15*dot23*dot9-2.*dot13*dot24*dot9+2.*a1*dot13*dot24*dot9+2.*dot10*dot25*dot9-2.*a1*dot10*dot25*dot9-4.*a1*dot12*dot25*dot9+2.*dot16*dot26*dot9-2.*a1*dot16*dot26*dot9+2.*dot16*dot27*dot9-2.*a1*dot16*dot27*dot9-4.*dot16*dot18*dot4*dot9+4.*a1*dot16*dot18*dot4*dot9+2.*dot24*dot4*dot9-2.*a1*dot24*dot4*dot9+4.*dot15*dot18*dot5*dot9-4.*a1*dot15*dot18*dot5*dot9-2.*dot28*dot5*dot9+2.*a1*dot28*dot5*dot9+4.*a1*dot25*dot6*dot9+4.*dot15*dot18*dot7*dot9-4.*a1*dot15*dot18*dot7*dot9-2.*dot28*dot7*dot9+2.*a1*dot28*dot7*dot9-4.*dot19*(-(a2*dot11)+a1*(dot13-dot4)+dot4)*(-(a2*dot10)+a1*(dot12-dot6))*M+(dot20+dot21+2.*a22*dot13*(2.*dot16*dot18-dot24)+a12*(-8.*dot10*dot17*dot18-8.*dot12*dot17*dot18+dot20+dot21-2.*dot15*(dot22+dot23)+6.*dot10*dot25+8.*dot12*dot25+2.*dot16*(dot26+dot27)-4.*dot16*dot18*dot4+2.*dot24*dot4+4.*dot15*dot18*dot5-2.*dot28*dot5+8.*dot17*dot18*dot6-8.*dot25*dot6+4.*dot15*dot18*dot7-2.*dot28*dot7)-2.*a1*(dot20+dot21+4.*dot10*(-(dot17*dot18)+dot25)+2.*dot16*(dot26+dot27)+2.*dot25*(dot12-dot6)-2.*dot28*(dot5+dot7)-2.*dot15*(dot23-2.*dot18*dot7))+2.*(dot10*dot25+dot24*dot4+dot16*(dot26+dot27-2.*dot18*dot4)-dot28*(dot5+dot7)-dot15*(dot22+dot23-2.*dot18*(dot5+dot7))))*M2+2.*dot14*(4.*dot10*dot27+4.*dot11*dot29+2.*dot26*dot6-2.*dot4*(2.*dot29+dot30+2.*dot18*dot6)+a12*dot28*M2+dot28*(dot2+dot9+M2)-a1*(4.*dot29*(dot11+dot13-dot4)+4.*dot27*(dot10+dot12-dot6)+dot28*(dot2+dot3+dot9+2.*M2))+2.*dot15*(-dot31+dot18*(-(a2*dot2)+a1*dot3-dot9+a1*dot9-a22*M2)))))/(16.*a1*a2*(dot8+a1*M2)*(-dot9-a2*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot8+a1*M2)));
	    // diag[37]=(M2*(4.*dot12*dot17*dot18*dot2-4.*a1*dot12*dot17*dot18*dot2-dot2*dot20+a1*dot2*dot20-dot2*dot21+a1*dot2*dot21-4.*dot10*dot11*dot22+4.*a1*dot10*dot11*dot22+4.*dot11*dot12*dot22-4.*a1*dot11*dot12*dot22+4.*a1*dot10*dot13*dot22-4.*a1*dot12*dot13*dot22+4.*dot15*dot2*dot22-4.*dot10*dot11*dot23+4.*a1*dot10*dot11*dot23+4.*a1*dot11*dot12*dot23+4.*dot10*dot13*dot23-4.*a1*dot10*dot13*dot23-4.*a1*dot12*dot13*dot23+2.*dot11*dot2*dot24-2.*a1*dot11*dot2*dot24-4.*a1*dot13*dot2*dot24-4.*dot10*dot2*dot25+4.*a1*dot10*dot2*dot25-2.*dot12*dot2*dot25+6.*a1*dot12*dot2*dot25-8.*dot10*dot14*dot27+8.*a1*dot10*dot14*dot27+8.*a1*dot12*dot14*dot27-2.*dot14*dot2*dot28+2.*a1*dot14*dot2*dot28-8.*dot11*dot14*dot29+8.*a1*dot11*dot14*dot29+8.*a1*dot13*dot14*dot29+2.*dot17*dot2*dot29-2.*a1*dot17*dot2*dot29+8.*dot11*dot16*dot18*dot3-8.*a1*dot11*dot16*dot18*dot3-8.*a1*dot13*dot16*dot18*dot3-8.*dot10*dot17*dot18*dot3+8.*a1*dot10*dot17*dot18*dot3+4.*a1*dot12*dot17*dot18*dot3+a1*dot20*dot3+a1*dot21*dot3-4.*dot15*dot23*dot3-4.*dot11*dot24*dot3+6.*a1*dot11*dot24*dot3+4.*a1*dot13*dot24*dot3+4.*dot10*dot25*dot3-4.*a1*dot10*dot25*dot3-2.*a1*dot12*dot25*dot3+2.*a1*dot14*dot28*dot3-2.*a1*dot17*dot29*dot3+2.*dot17*dot2*dot30-2.*a1*dot17*dot2*dot30-2.*a1*dot17*dot3*dot30-4.*dot14*dot15*dot31-4.*dot11*dot16*dot31+4.*a1*dot11*dot16*dot31+4.*a1*dot13*dot16*dot31+4.*dot10*dot17*dot31-4.*a1*dot10*dot17*dot31-4.*a1*dot12*dot17*dot31-4.*a1*dot10*dot22*dot4+4.*a1*dot12*dot22*dot4-4.*dot10*dot23*dot4+4.*a1*dot10*dot23*dot4+4.*dot12*dot23*dot4+4.*a1*dot12*dot23*dot4-4.*dot2*dot24*dot4+4.*a1*dot2*dot24*dot4-8.*a1*dot14*dot29*dot4+8.*a1*dot16*dot18*dot3*dot4-4.*a1*dot24*dot3*dot4-4.*dot14*dot30*dot4-4.*a1*dot16*dot31*dot4-4.*dot10*dot43+4.*a1*dot10*dot43+4.*a1*dot12*dot43+4.*dot11*dot44-4.*a1*dot11*dot44-4.*a1*dot13*dot44-2.*dot4*dot44+4.*a1*dot4*dot44-8.*dot11*dot12*dot18*dot5+8.*a1*dot11*dot12*dot18*dot5+8.*a1*dot12*dot13*dot18*dot5+8.*dot10*dot27*dot5-8.*a1*dot10*dot27*dot5-8.*a1*dot12*dot27*dot5+2.*dot2*dot28*dot5-2.*a1*dot2*dot28*dot5+4.*dot11*dot29*dot5-4.*a1*dot11*dot29*dot5-4.*a1*dot13*dot29*dot5-2.*a1*dot28*dot3*dot5-4.*dot11*dot30*dot5+4.*a1*dot11*dot30*dot5+4.*a1*dot13*dot30*dot5+4.*dot15*dot31*dot5-8.*a1*dot12*dot18*dot4*dot5+4.*a1*dot29*dot4*dot5+4.*dot30*dot4*dot5-4.*a1*dot30*dot4*dot5-4.*dot17*dot18*dot2*dot6+4.*a1*dot17*dot18*dot2*dot6+4.*a1*dot11*dot22*dot6+4.*a1*dot13*dot22*dot6+4.*dot11*dot23*dot6-4.*a1*dot11*dot23*dot6-4.*dot13*dot23*dot6+4.*a1*dot13*dot23*dot6+6.*dot2*dot25*dot6-6.*a1*dot2*dot25*dot6+4.*dot14*dot26*dot6+8.*dot14*dot27*dot6-8.*a1*dot14*dot27*dot6-4.*a1*dot17*dot18*dot3*dot6+2.*a1*dot25*dot3*dot6+4.*a1*dot17*dot31*dot6+8.*dot14*dot18*dot4*dot6+4.*dot22*dot4*dot6-4.*a1*dot22*dot4*dot6+4.*dot23*dot4*dot6-4.*a1*dot23*dot4*dot6+2.*dot43*dot6-4.*a1*dot43*dot6+8.*dot11*dot18*dot5*dot6-8.*a1*dot11*dot18*dot5*dot6-8.*a1*dot13*dot18*dot5*dot6-4.*dot26*dot5*dot6-8.*dot27*dot5*dot6+8.*a1*dot27*dot5*dot6-8.*dot18*dot4*dot5*dot6+8.*a1*dot18*dot4*dot5*dot6-8.*dot10*dot13*dot18*dot7+8.*a1*dot10*dot13*dot18*dot7+8.*a1*dot12*dot13*dot18*dot7-4.*dot10*dot26*dot7+4.*a1*dot10*dot26*dot7+4.*a1*dot12*dot26*dot7+4.*dot10*dot27*dot7-4.*a1*dot10*dot27*dot7-4.*a1*dot12*dot27*dot7+2.*dot2*dot28*dot7-2.*a1*dot2*dot28*dot7+8.*dot11*dot29*dot7-8.*a1*dot11*dot29*dot7-8.*a1*dot13*dot29*dot7+8.*dot15*dot18*dot3*dot7-2.*a1*dot28*dot3*dot7+8.*dot10*dot18*dot4*dot7-8.*a1*dot10*dot18*dot4*dot7-8.*dot12*dot18*dot4*dot7-8.*a1*dot12*dot18*dot4*dot7-4.*dot29*dot4*dot7+8.*a1*dot29*dot4*dot7+8.*dot13*dot18*dot6*dot7-8.*a1*dot13*dot18*dot6*dot7-4.*a1*dot26*dot6*dot7-4.*dot27*dot6*dot7+4.*a1*dot27*dot6*dot7-8.*dot18*dot4*dot6*dot7+8.*a1*dot18*dot4*dot6*dot7+dot1*(4.*(1.+a1)*dot12*dot17*dot18-dot20-dot21-4.*dot15*dot22-6.*dot11*dot24+4.*dot10*dot25-2.*(1.+3.*a1)*dot12*dot25-2.*dot14*dot28+2.*(2.*dot24*dot4-dot25*dot6+dot17*(dot29+dot30-2.*dot18*dot6)+dot28*(dot5+dot7))-a1*(dot20+dot21-2.*(-2.*dot10*dot25+dot24*(dot11+2.*dot13-2.*dot4)+3.*dot25*dot6+dot17*(dot29+dot30-2.*dot18*dot6)+dot28*(-dot14+dot5+dot7))))-8.*dot10*dot17*dot18*dot8+8.*a1*dot10*dot17*dot18*dot8+4.*a1*dot12*dot17*dot18*dot8+a1*dot20*dot8+a1*dot21*dot8+2.*a1*dot11*dot24*dot8+4.*dot10*dot25*dot8-4.*a1*dot10*dot25*dot8-2.*a1*dot12*dot25*dot8+2.*a1*dot14*dot28*dot8-2.*a1*dot17*dot29*dot8-2.*a1*dot17*dot30*dot8-2.*a1*dot28*dot5*dot8+8.*dot17*dot18*dot6*dot8-4.*a1*dot17*dot18*dot6*dot8-4.*dot25*dot6*dot8+2.*a1*dot25*dot6*dot8-2.*a1*dot28*dot7*dot8+4.*dot15*dot22*dot9-4.*dot10*dot25*dot9+4.*a1*dot10*dot25*dot9+4.*a1*dot12*dot25*dot9-4.*dot24*dot4*dot9+4.*dot25*dot6*dot9-4.*a1*dot25*dot6*dot9+4.*dot19*(-(a2*dot11)+a1*(dot13-dot4))*(-(a2*dot10)+a1*(dot12-dot6)+dot6)*M+(4.*a1*(-2.*dot10*dot17*dot18-dot15*dot22+3.*dot10*dot25+dot12*dot25+dot24*dot4+2.*dot17*dot18*dot6-3.*dot25*dot6)+4.*(dot15*dot22-dot10*dot25-dot24*dot4+dot25*dot6)+a12*(8.*dot10*dot17*dot18+4.*dot12*dot17*dot18+dot20+dot21+2.*dot11*dot24-8.*dot10*dot25-6.*dot12*dot25-2.*dot17*dot29-2.*dot17*dot30-4.*dot17*dot18*dot6+6.*dot25*dot6+2.*dot28*(dot14-dot5-dot7)))*M2))/(16.*a1*a2*(dot8+a1*M2)*(-(a2*dot1)+dot2+dot9-a1*(dot2+dot3+dot9)+M2+(-2.+a1)*a1*M2)*(dot9+a2*M2));
	    // colour flows
	    Complex flow[3] = {diag[0] + diag[8] + diag[11] + diag[12] + diag[13] - diag[16] + diag[17] + diag[20] + diag[21] + diag[24] + diag[25] + diag[36] + diag[37],
			       diag[1] + diag[3] + diag[5] + diag[7] + diag[29] + diag[30] + diag[31] + diag[32] + diag[33] + 2.*diag[37] +
			       (diag[4] - diag[5] + diag[6] - diag[7] + diag[9] - 2.*diag[11] - 2.*diag[13] + diag[15] + diag[18] + diag[23] - 2.*diag[24] - diag[29] - diag[30] - diag[32])/9.,
			       diag[2] - diag[3] - diag[5] - diag[7] + diag[26] + diag[27] + diag[28] + diag[34] + diag[35] + 2.*diag[36] +
			       (-diag[4] + diag[5] - diag[6] + diag[7] - 2.*diag[8] + diag[10] + diag[14] + 2.*diag[16] + diag[19] - 2.*diag[20] + diag[22] - diag[26] - diag[27] - diag[34])/9.};
	    for(unsigned int ix=0;ix<3;++ix) {
	      flows[ix]+=norm(flow[ix]);
	      for(unsigned int iy=0;iy<3;++iy) {
		me2Sum+=cMatrix[ix][iy]*flow[ix]*std::conj(flow[iy]);
	      }
	    }
	  }
  	}
      }
    }
  }
  // return the answer
  meInfo(flows);
  return O1_*sHat()/M/216.*pow<4,1>(Constants::pi*standardModel()->alphaS(scale())/M)*me2Sum.real();
}
