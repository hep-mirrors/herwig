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
  // colour matrix
  vector<vector<double > > cMatrix = {{24., 12., 12.},
  				      {12., 48., -6.},
  				      {12., -6., 48.}};
  // invariants
  Energy M  = rescaledMomenta()[2].mass();
  Energy2 M2 = sqr(M);
  double a1 = rescaledMomenta()[3].mass()/M, a2 = rescaledMomenta()[4].mass()/M;
  double a12(sqr(a1)),a22(sqr(a2));
  Energy2 dot1 = rescaledMomenta()[0]*rescaledMomenta()[1];
  Energy2 dot2 = rescaledMomenta()[0]*rescaledMomenta()[3];
  Energy2 dot3 = rescaledMomenta()[0]*rescaledMomenta()[4];
  Energy2 dot6 = rescaledMomenta()[2]*rescaledMomenta()[3];
  Energy2 dot7 = rescaledMomenta()[2]*rescaledMomenta()[4];
  // matrix element
  ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half);
  Complex me2Sum=0.;
  vector<double> flows(3,0.);
  for(unsigned int ih1=0;ih1<2;++ih1) {
    complex<Energy> dot5 = rescaledMomenta()[1]*g1[ih1].wave();
    complex<Energy> dot8 = rescaledMomenta()[3]*g1[ih1].wave();
    complex<Energy> dot10 = rescaledMomenta()[4]*g1[ih1].wave();
    for(unsigned int ih2=0;ih2<2;++ih2) {
      complex<Energy> dot4 = rescaledMomenta()[0]*g2[ih2].wave();
      complex<Energy> dot9 = rescaledMomenta()[3]*g2[ih2].wave();
      complex<Energy> dot11 = rescaledMomenta()[4]*g2[ih2].wave();
      for(unsigned int ih4=0;ih4<2;++ih4) {
  	for(unsigned int ih5=0;ih5<2;++ih5) {
	  Complex dot12 = g1[ih1].wave()*g2[ih2].wave();
	  complex<Energy> dot13=v4[ih4].dimensionedWave().pseudoScalar(ubar5[ih5].dimensionedWave());
	  LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  LorentzPolarizationVectorE vec3 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  LorentzVector<complex<Energy2> > vec4 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[1]).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  LorentzPolarizationVectorE vec5 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(g1[ih1].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  LorentzVector<complex<Energy2> > vec6 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[1]).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  LorentzVector<complex<Energy2> > vec7 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(rescaledMomenta()[1]).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  LorentzVector<complex<Energy2> > vec8 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(g1[ih1].wave()).slash(rescaledMomenta()[1]).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  complex<Energy2> dot14 = vec1*rescaledMomenta()[0];
	  complex<Energy2> dot15 = vec1*rescaledMomenta()[1];
	  complex<Energy> dot16 = vec1*g1[ih1].wave();
	  complex<Energy> dot17 = vec1*g2[ih2].wave();
	  complex<Energy2> dot18 = vec5*rescaledMomenta()[0];
	  complex<Energy2> dot19 = vec5*rescaledMomenta()[1];
	  complex<Energy2> dot20 = vec2*rescaledMomenta()[0];
	  complex<Energy2> dot21 = vec2*rescaledMomenta()[1];
	  complex<Energy2> dot22 = vec3*rescaledMomenta()[0];
	  complex<Energy2> dot23 = vec3*rescaledMomenta()[1];
	  complex<Energy3> dot24 = vec4*rescaledMomenta()[0];
	  complex<Energy> dot25 = vec3*g1[ih1].wave();
	  complex<Energy3> dot26 = vec6*rescaledMomenta()[0];
	  complex<Energy3> dot27 = vec7*rescaledMomenta()[0];
	  complex<Energy3> dot28 = vec8*rescaledMomenta()[0];
	  // diagrams
	  Complex diag[38];
	  diag[0]=(-(dot17*(dot10-dot5+dot8))-dot16*(dot11-dot4+dot9)+dot12*(dot14+dot15+2.*dot13*M))*M2/(8.*a1*a2*(dot6+a1*M2)*(dot7+a2*M2));
	  diag[1]=(3.*(dot18+dot19)+2.*dot17*(dot10-dot5-2.*dot8)+2.*dot16*(-2.*dot11+2.*dot4+dot9)-2.*dot12*(dot14+dot15-dot13*M))*M2/(8.*a1*a2*(dot6+a1*M2)*(dot7+a2*M2));
	  diag[2]=(4.*dot12*(dot14+dot15)-3.*(dot18+dot19)+2.*dot17*(-2.*dot10+2.*dot5+dot8)+2.*dot16*(dot11-dot4-2.*dot9)+2.*dot12*dot13*M)*M2/(8.*a1*a2*(dot6+a1*M2)*(dot7+a2*M2));
	  diag[3]=M2*(2.*dot12*dot14*dot3+2.*dot12*dot15*dot3-2.*dot10*dot14*dot4-2.*dot10*dot15*dot4+2.*dot11*dot14*dot5+
		      2.*dot11*dot15*dot5+2.*dot12*dot14*dot6-2.*dot16*dot4*dot6+2.*dot17*dot5*dot6-dot12*dot14*dot7+dot12*dot15*dot7+
		      2.*dot16*dot4*dot7-2.*dot17*dot5*dot7+dot13*(2.*dot12*dot2+dot12*dot7-2.*dot4*dot8+2.*dot5*dot9)*M+dot12*dot13*pow<3,1>(M)+
		      (-1.+2.*a1)*dot1*dot12*(dot14+dot15+dot13*M)+(-(dot12*dot14)+dot12*dot15+2.*dot16*dot4-2.*dot17*dot5)*M2+
		      a1*(-((dot14+dot15)*(dot12*(2.*dot2+2.*dot3+dot6+dot7)-2.*dot4*(dot10+dot8)+2.*dot5*(dot11+dot9)))-
			  dot13*(dot12*(2.*dot2+2.*dot3+dot6+dot7)-2.*dot4*(dot10+dot8)+2.*dot5*(dot11+dot9))*M-
			  2.*dot12*dot13*pow<3,1>(M)+2.*(dot12*dot14-dot12*dot15-2.*dot16*dot4+2.*dot17*dot5)*M2))/(4.*a1*a2*dot1*(dot6+a1*M2)*(dot7+a2*M2));
	  diag[4]=M2*(2.*dot12*dot14*dot3+2.*dot12*dot15*dot3-2.*dot10*dot14*dot4-2.*dot10*dot15*dot4+2.*dot11*dot14*dot5+
		      2.*dot11*dot15*dot5+dot12*dot14*dot6+dot12*dot15*dot6-dot12*dot14*dot7+dot12*dot15*dot7+
		      2.*dot16*dot4*dot7-2.*dot17*dot5*dot7+
		      (-(dot12*dot24)-2.*dot10*dot13*dot4+(2.*dot11*dot13+dot22+dot23)*dot5+
		       dot12*dot13*(2.*dot2+2.*dot3+dot6+dot7)-dot4*(dot20+dot21+2.*dot13*dot8)+2.*dot13*dot5*dot9)*M+
		      dot12*dot13*pow<3,1>(M)+(-1.+2.*a1)*dot1*dot12*(dot14+dot15+dot13*M)+
		      (dot12*(-dot14+dot15)+2.*dot16*dot4-2.*dot17*dot5)*M2+
		      a1*(2.*(dot14+dot15)*(dot4*(dot10+dot8)-dot5*(dot11+dot9))+
			  2.*dot13*(dot4*(dot10+dot8)-dot5*(dot11+dot9))*M-dot12*dot14*(2.*dot2+2.*dot3+dot6+dot7-M2)+
			  2.*(-(dot16*dot4)+dot17*dot5)*M2-dot12*(dot15+dot13*M)*(2.*dot2+2.*dot3+dot6+dot7+M2)))/(4.*a1*dot1*(dot6+a1*M2)*(dot1-a2*(dot6+dot7+M2)));
	  diag[5]=M2*(dot12*(dot24-2.*dot13*dot3)*M+(2.*dot10*dot13+dot20+dot21)*dot4*M-(2.*dot11*dot13+dot22+dot23)*dot5*M-
		      2.*dot16*dot4*(dot6+a1*M2)+2.*dot17*dot5*(dot6+a1*M2)+dot12*(dot14-dot15-dot13*M)*(dot6+a1*M2))/(4.*a1*dot1*sqr(dot6+a1*M2));
	  diag[6]=-0.25*M2*(-2.*(dot14+dot15)*(dot12*dot3-dot10*dot4+dot11*dot5)-2.*(dot12*dot14-dot16*dot4+dot17*dot5)*dot6+
			    (dot12*dot24+(dot20+dot21)*dot4-(dot22+dot23)*dot5)*M+
			    dot1*dot12*(dot14-2.*a1*dot14+dot15-dot13*M-2.*a1*(dot15+dot13*M))+
			    a1*((dot14+dot15)*(dot12*(2.*dot2+2.*dot3+dot6+dot7)-2.*dot4*(dot10+dot8)+2.*dot5*(dot11+dot9))+
				dot13*(dot12*(2.*dot2+2.*dot3+dot6+dot7)-2.*dot4*(dot10+dot8)+2.*dot5*(dot11+dot9))*M+
				dot12*dot13*pow<3,1>(M)+(-(dot12*dot14)+dot12*dot15+2.*dot16*dot4-2.*dot17*dot5)*M2))/(a2*dot1*(dot7+a2*M2)*(dot1-a1*(dot6+dot7+M2)));
	  diag[7]=-0.25*M2*(-2.*dot16*dot4*dot7+2.*dot17*dot5*dot7+dot12*(2.*dot1*dot13-2.*dot13*dot2-dot24)*M+(dot22+dot23)*dot5*M-
			    dot4*(dot20+dot21-2.*dot13*dot8)*M-2.*dot13*dot5*dot9*M-
			    2.*a2*(dot16*dot4-dot17*dot5)*M2+dot12*(dot14-dot15-dot13*M)*(dot7+a2*M2))/(a2*dot1*sqr(dot7+a2*M2));
	  diag[8]=M2*(2.*dot18*dot2-2.*dot19*dot3-2.*a2*(dot26-2.*dot16*dot3)*dot4-4.*a1*dot17*dot2*dot5-2.*a1*dot27*dot5+
		      4.*dot17*dot3*dot5+4.*dot14*dot4*dot5-8.*a1*dot14*dot4*dot5+4.*a12*dot14*dot4*dot5-4.*a1*dot15*dot4*dot5+
		      4.*a12*dot15*dot4*dot5+2.*a2*dot11*(dot26-2.*dot16*dot3+2.*(-(a2*dot14)+a1*dot15)*dot5)+2.*dot18*dot7-
		      4.*a1*dot17*dot5*dot7-4.*a1*dot11*dot14*dot8+4.*a12*dot11*dot14*dot8-4.*a1*dot11*dot15*dot8+4.*a12*dot11*dot15*dot8+
		      4.*a1*dot17*dot2*dot8+2.*a1*dot27*dot8+4.*a1*dot14*dot4*dot8-4.*a12*dot14*dot4*dot8+4.*a1*dot15*dot4*dot8-
		      4.*a12*dot15*dot4*dot8+4.*a1*dot17*dot7*dot8-2.*a1*dot26*dot9-4.*dot16*dot3*dot9+4.*a1*dot16*dot3*dot9+
		      8.*a1*dot14*dot5*dot9-4.*a12*dot14*dot5*dot9-4.*a12*dot15*dot5*dot9-4.*a1*dot14*dot8*dot9+4.*a12*dot14*dot8*dot9+4.*a12*dot15*dot8*dot9+
		      (-2.*dot12*dot24+3.*dot28-2.*a2*dot20*dot4-2.*dot21*dot4-4.*dot22*dot5+2.*a1*dot23*dot5+
		       4.*(1.-a1*a2)*dot13*dot4*dot5+2.*a2*dot11*(dot20+2.*a1*dot13*(dot5-dot8))-2.*a1*(dot23-2.*a2*dot13*dot4)*dot8+
		       2.*a2*(dot20+2.*a1*dot13*(dot5-dot8))*dot9)*M-
		      2.*dot1*(dot18+2.*a1*dot17*(-dot5+dot8)+dot25*M)+2.*a2*(dot18+2.*a1*dot17*(-dot5+dot8))*M2+
		      2.*dot10*(-dot27-2.*dot14*dot4+2.*dot14*dot9-2.*a2*dot11*(-(a2*dot14)+a1*(dot15+dot13*M))-
				2.*a12*((dot14+dot15)*(dot4-dot9)+dot13*(dot4-dot9)*M+dot17*M2)+
				a1*(dot27+4.*dot14*dot4+2.*dot15*dot4-4.*dot14*dot9-(dot23-2.*dot13*dot4+2.*dot13*dot9)*M+
				    2.*dot17*(-dot1+dot2+dot7+M2))))
	    /(16.*a1*a2*(dot1-dot2-dot3)*(-dot1+dot2+dot3+dot6+dot7+M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	  diag[9]=-0.125*M2*(2.*(dot1*dot18+dot10*dot27+dot19*dot3-2.*dot17*dot3*dot5-dot18*(dot2+dot7)-
				 2.*dot10*dot14*dot9+2.*dot16*dot3*dot9+
				 a1*(dot10-dot5+dot8)*(2.*dot1*dot17-dot27-2.*dot17*(dot2+dot7)+2.*dot14*dot9))-
			     (dot28-2.*dot22*dot5+2.*dot20*dot9+2.*a1*(dot10-dot5+dot8)*(dot23-2.*dot13*dot9))*M-
			     2.*a2*(dot18+2.*a1*dot17*(dot10-dot5+dot8))*M2)
	    /(a1*a2*(dot1-dot2-dot3)*(-dot7-a2*M2)*(dot1-dot2-dot7-a2*M2));
	  diag[10]=M2*(2.*(-(dot18*dot3)-dot19*dot3-2.*dot16*dot3*dot4+dot27*dot5-a1*dot27*dot5+2.*a1*dot17*dot3*dot5+
			   2.*dot15*dot4*dot5-2.*a1*dot15*dot4*dot5-dot18*dot6-2.*dot16*dot4*dot6+2.*a1*dot17*dot5*dot6+
			   2.*dot12*(dot15*dot3+dot14*(dot3+dot6))+a1*(dot27+2.*dot15*dot4-2.*dot17*(dot3+dot6))*dot8+
			   2.*dot11*(dot16*dot3+dot14*(dot5-a1*dot5+a1*dot8)))+
		       (-2.*dot12*dot24+dot28-2.*dot21*dot4-2.*(dot11*dot20+a1*(2.*dot11*dot13+dot23)*dot5)+2.*a1*(2.*dot11*dot13+dot23)*dot8)*M+
		       2.*dot1*(-2.*dot11*dot16+dot18+dot19+2.*dot16*dot4-2.*dot17*dot5-dot25*M-2.*dot12*(dot14+dot15-dot13*M))+
		       2.*a1*(2.*dot12*dot14-dot18-2.*(dot16*dot4+a1*dot17*(-dot5+dot8)))*M2+
		       2.*dot10*(2.*dot1*dot17-dot27-2.*dot15*dot4+2.*dot11*(-(a2*dot14)+a1*dot13*M)+
				 a1*(dot27+2.*dot15*dot4+dot23*M-2.*dot17*(dot3+dot6+a1*M2))))
	    /(8.*a1*(dot1-dot2-dot3)*(-dot1+dot3+dot6+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	  diag[11]=-0.0625*M2*(-2.*dot11*dot26+2.*(-2.*dot12*(dot14+dot15)+dot18+dot19)*dot3+2.*dot26*dot4+
			       2.*dot10*(dot27+2.*(dot14+dot15)*dot4)-2.*dot27*dot5-4.*dot14*dot4*dot5-
			       4.*dot15*dot4*dot5-4.*dot12*dot14*dot6+2.*dot18*dot6+4.*dot16*dot4*dot6+
			       (4.*dot12*dot24-3.*dot28+2.*dot22*dot5+2.*dot4*(dot20+2.*dot21-2.*dot13*dot5))*M-
			       2.*dot1*(-2.*dot11*dot16+2.*dot10*dot17+dot18+dot19+2.*dot16*dot4-2.*dot17*dot5-2.*dot25*M-2.*dot12*(dot14+dot15-dot13*M))-
			       4.*a12*(dot10-dot5+dot8)*((dot14+dot15)*(dot11-dot4+dot9)+dot13*(dot11-dot4+dot9)*M-dot17*M2)+
			       2.*a1*(2.*dot12*dot14*dot2+2.*dot12*dot15*dot2+2.*dot11*dot16*dot2-dot18*dot2-dot19*dot2+dot11*dot26+
				      2.*dot12*dot14*dot3+2.*dot12*dot15*dot3-dot18*dot3-dot19*dot3-2.*dot16*dot2*dot4-dot26*dot4-
				      2.*dot11*dot14*dot5-2.*dot11*dot15*dot5+2.*dot17*dot2*dot5+dot27*dot5+4.*dot14*dot4*dot5+
				      4.*dot15*dot4*dot5-2.*dot17*dot5*dot6+2.*dot11*dot15*dot8-dot27*dot8+2.*dot17*dot3*dot8-
				      2.*dot14*dot4*dot8-4.*dot15*dot4*dot8+2.*dot17*dot6*dot8+dot26*dot9-2.*dot16*dot3*dot9-2.*dot14*dot5*dot9+
				      (dot11*dot20+dot25*(dot2+dot3)-(dot23-2.*dot13*dot4)*(dot5-dot8)+dot20*(-dot4+dot9))*M+
				      dot1*(-2.*dot12*(dot14+dot15)-2.*dot11*dot16+dot18+dot19+2.*dot16*dot4-2.*dot17*dot5-dot25*M)+
				      dot10*(2.*dot11*(dot14+dot15)-dot27-4.*dot14*dot4-4.*dot15*dot4+
					     2.*dot17*(dot1-dot2+dot6)+2.*dot14*dot9+dot23*M-2.*dot13*dot4*M)+
				      (-2.*dot12*dot14+dot18+2.*dot16*dot4)*M2))/(a1*a2*(dot1-dot2-dot3)*(-dot7-a2*M2)*(-dot1+a1*(dot6+dot7+M2)));
	  diag[12]=-0.0625*M2*(4.*a12*(dot10-dot5+dot8)*((dot14+dot15)*(dot11-dot4+dot9)+dot13*(dot11-dot4+dot9)*M-2.*dot17*M2)+
			       a1*(2.*dot1*dot12*dot14+2.*dot1*dot12*dot15+2.*dot1*dot11*dot16-dot1*dot18-dot1*dot19-
				   2.*dot12*dot14*dot2-2.*dot12*dot15*dot2-2.*dot11*dot16*dot2+dot18*dot2+dot19*dot2-
				   2.*dot11*dot26-2.*dot12*dot14*dot3-2.*dot12*dot15*dot3+2.*dot11*dot16*dot3+dot18*dot3+
				   dot19*dot3-2.*dot1*dot16*dot4+2.*dot16*dot2*dot4+2.*dot26*dot4-2.*dot16*dot3*dot4+
				   4.*dot11*dot14*dot5+4.*dot11*dot15*dot5+4.*dot1*dot17*dot5-4.*dot17*dot2*dot5-
				   4.*dot27*dot5+4.*dot17*dot3*dot5-8.*dot14*dot4*dot5-8.*dot15*dot4*dot5+4.*dot17*dot5*dot6-
				   4.*dot17*dot5*dot7-4.*dot11*dot15*dot8-2.*dot1*dot17*dot8+2.*dot17*dot2*dot8+4.*dot27*dot8-
				   6.*dot17*dot3*dot8+4.*dot14*dot4*dot8+8.*dot15*dot4*dot8-4.*dot17*dot6*dot8+4.*dot17*dot7*dot8-
				   2.*dot26*dot9+4.*dot16*dot3*dot9+8.*dot14*dot5*dot9-4.*dot14*dot8*dot9-2.*dot11*dot20*M+
				   2.*(dot20+2.*dot13*(-dot5+dot8))*(dot4-dot9)*M-
				   4.*dot10*(dot11*(dot14+dot15)-dot27-2.*dot14*dot4-2.*dot15*dot4+2.*dot14*dot9-
					     dot13*dot4*M+dot13*dot9*M+dot17*(dot1-dot2+dot3+dot6-dot7-M2))+
				   4.*(-(dot12*dot14)+dot18+dot16*dot4-dot17*dot5+dot17*dot8)*M2)-
			       2.*(dot18*dot2+2.*dot19*dot2+dot11*dot26-dot18*dot3+2.*dot10*dot14*dot4+2.*dot16*dot2*dot4+
				   dot26*dot4-2.*dot16*dot3*dot4-2.*dot11*dot14*dot5-2.*dot17*dot2*dot5-dot27*dot5-
				   2.*dot14*dot4*dot5-2.*dot15*dot4*dot5-dot18*dot6+dot18*dot7+2.*dot16*dot4*dot7+
				   2.*dot27*dot8+2.*dot15*dot4*dot8-2.*dot10*dot14*dot9-2.*dot26*dot9+2.*dot16*dot3*dot9+2.*dot14*dot5*dot9+
				   dot1*(2.*dot12*(dot14+dot15)-dot19-2.*dot16*dot4+2.*dot17*dot5-2.*dot17*dot8+2.*dot16*dot9)+
				   dot20*(dot4-dot9)*M+(dot18+2.*dot16*dot4)*M2-2.*dot12*(dot15*dot2+dot14*(dot2+dot7+M2))))
	    /(a1*a2*(dot1-dot2-dot3)*(-dot7-a2*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	  diag[13]=M2*(2.*(2.*dot12*dot15*dot2-dot19*dot2-dot11*dot26-2.*dot12*dot14*dot3+dot18*dot3+2.*dot16*dot3*dot4+
			   2.*dot11*dot15*dot5-dot27*dot5+2.*dot17*dot3*dot5-2.*dot15*dot4*dot5-2.*dot12*dot14*dot6+
			   dot18*dot6+2.*dot16*dot4*dot6+2.*dot17*dot5*dot6+
			   dot10*(dot27+2.*dot15*(-dot11+dot4)-2.*dot17*(dot3+dot6))-2.*(dot11*dot15+dot17*(dot3+dot6))*dot8+
			   dot1*(2.*dot12*dot14+2.*dot11*dot16-dot18-2.*dot16*dot4+2.*dot17*dot8))+
		       (4.*dot12*dot24+4.*dot1*(-(dot12*dot13)+dot25)-3.*dot28+4.*dot21*dot4+
			2.*dot22*dot5-4.*dot13*dot4*dot5-2.*dot23*(dot10-dot5+dot8))*M-
		       4.*a12*(dot10-dot5+dot8)*((dot14+dot15)*(dot11-dot4+dot9)+dot13*(dot11-dot4+dot9)*M-dot17*M2)+
		       2.*a1*(-2.*dot1*dot11*dot16+2.*dot11*dot16*dot2+dot11*dot26+2.*dot1*dot16*dot4-2.*dot16*dot2*dot4-
			      dot26*dot4-2.*dot11*dot14*dot5-4.*dot11*dot15*dot5+dot27*dot5-2.*dot17*dot3*dot5+
			      2.*dot14*dot4*dot5+4.*dot15*dot4*dot5-2.*dot17*dot5*dot6+4.*dot11*dot15*dot8-
			      dot27*dot8+2.*dot17*dot3*dot8-4.*dot15*dot4*dot8+2.*dot17*dot6*dot8-2.*dot1*dot16*dot9+
			      2.*dot16*dot2*dot9+dot26*dot9-2.*dot14*dot5*dot9-2.*dot15*dot5*dot9+2.*dot15*dot8*dot9+
			      (dot23*(-dot5+dot8)+dot11*(dot20-2.*dot13*dot5+2.*dot13*dot8)-
			       (dot20-2.*dot13*dot5+2.*dot13*dot8)*(dot4-dot9))*M+
			      (-2.*dot12*dot14+dot18+2.*dot16*dot4+2.*dot17*dot5-2.*dot17*dot8)*M2+
			      dot10*(-dot27+2.*(dot14+2.*dot15)*(dot11-dot4)+2.*dot17*(dot3+dot6)+
				     2.*(dot14+dot15)*dot9+(dot23+2.*dot13*(dot11-dot4+dot9))*M-2.*dot17*M2)))
	    /(16.*a1*a2*(dot1-dot2-dot3)*(-dot1+dot2+dot3+dot6+dot7+M2)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2));
	  diag[14]=-0.125*M2*(2.*(-(dot19*dot2)-dot11*dot26+(dot18+2.*dot16*dot4)*(dot3+dot6)+2.*dot12*(dot15*dot2-dot14*(dot3+dot6))-
				  2.*(dot11*dot15+dot17*(dot3+dot6))*dot8+a1*(2.*dot16*dot2+dot26+2.*dot15*dot8)*(dot11-dot4+dot9))-
			      (-2.*dot12*dot24+dot28-2.*dot21*dot4+2.*dot23*dot8+2.*a1*(dot20-2.*dot13*dot8)*(dot11-dot4+dot9))*M+
			      2.*dot1*(-dot18+2.*dot17*dot8-2.*dot16*(-(a2*(dot11-dot4))+a1*dot9)+dot25*M+2.*dot12*(dot14-dot13*M))+
			      2.*a1*(-2.*dot12*dot14+dot18+2.*dot16*dot4-2.*dot17*dot8)*M2)
	    /(a1*a2*dot2*(-dot1+dot2+dot3+dot6+dot7+M2)*(-dot7-a2*M2));
	  diag[15]=M2*(2.*(dot18*dot3+dot19*dot3+dot26*dot4-a1*dot26*dot4+2.*a1*dot16*dot3*dot4-2.*dot17*dot3*dot5-
			   2.*dot14*dot4*dot5+2.*a1*dot14*dot4*dot5+dot11*(-(a2*dot26)-2.*a1*dot16*dot3+2.*a2*dot14*dot5)+
			   dot18*dot6+a1*(dot26-2.*(dot16*dot3+dot14*dot5))*dot9+
			   2.*dot10*(dot17*(dot3+dot6)-dot15*(-(a2*(dot11-dot4))+a1*dot9)))-
		       (-2.*dot10*dot23+dot28-2.*dot22*dot5+2.*a1*(2.*dot10*dot13+dot20)*(dot11-dot4+dot9))*M+
		       2.*a1*(2.*dot10*dot17+dot18)*M2)
	    /(8.*a1*dot3*(-dot1+dot2+dot3+dot6+dot7+M2)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2));
	  diag[16]=-0.0625*M2*(-2.*(-(dot11*dot26)+dot10*dot27+(dot18+dot19)*dot3+dot26*dot4+2.*dot11*(dot14+dot15)*dot5-
				    (dot27+2.*(dot14+dot15)*dot4)*dot5+(dot18+2.*dot17*dot5)*dot6)-
			       (2.*dot12*dot24+2.*dot1*dot25-3.*dot28+2.*dot21*dot4+2.*(2.*dot22+dot23-2.*dot13*dot4)*dot5)*M+
			       4.*a12*(dot10-dot5+dot8)*((dot14+dot15)*(dot11-dot4+dot9)+dot13*(dot11-dot4+dot9)*M-dot17*M2)-
			       2.*a1*(2.*dot11*dot16*dot2-dot18*dot2-dot19*dot2+dot11*dot26-dot18*dot3-dot19*dot3-2.*dot16*dot2*dot4-
				      dot26*dot4-4.*dot11*dot14*dot5-4.*dot11*dot15*dot5+2.*dot17*dot2*dot5+dot27*dot5+
				      4.*dot14*dot4*dot5+4.*dot15*dot4*dot5+2.*dot11*dot16*dot6-dot18*dot6-dot19*dot6-
				      2.*dot16*dot4*dot6+2.*dot11*dot16*dot7-dot18*dot7-dot19*dot7-2.*dot16*dot4*dot7+
				      2.*dot17*dot5*dot7+2.*dot11*dot15*dot8-dot27*dot8+2.*dot17*dot3*dot8-2.*dot15*dot4*dot8+
				      2.*dot17*dot6*dot8+dot26*dot9-2.*dot16*dot3*dot9-4.*dot14*dot5*dot9-2.*dot15*dot5*dot9+
				      (dot11*(dot20-2.*dot13*dot5)-2.*dot12*dot13*(dot2+dot3+dot6+dot7)+
				       dot25*(dot2+dot3+dot6+dot7)+dot23*(-dot5+dot8)-(dot20-2.*dot13*dot5)*(dot4-dot9))*M+
				      (-2.*dot12*dot13+dot25)*pow<3,1>(M)+
				      dot1*(-2.*dot11*dot16+dot18+dot19+2.*dot16*dot4-2.*dot17*dot5+2.*dot12*dot13*M-dot25*M)+
				      (2.*dot11*dot16-dot19-2.*dot16*dot4+4.*dot17*dot5)*M2+
				      dot10*(2.*dot1*dot17-dot27+2.*(dot14+dot15)*(dot11-dot4)+
					     2.*dot14*dot9+dot23*M-2.*dot17*(dot2+dot7+M2))))
	    /(a1*a2*(-dot1+dot2+dot3+dot6+dot7+M2)*(-dot7-a2*M2)*(-dot1+a1*(dot6+dot7+M2)));
	  diag[17]=-0.0625*M2*(4.*a12*(dot10-dot5+dot8)*((dot14+dot15)*(dot11-dot4+dot9)+dot13*(dot11-dot4+dot9)*M-dot17*M2)+
			       a1*(-4.*dot11*dot16*dot2+dot18*dot2+dot19*dot2-4.*dot11*dot26+4.*dot11*dot16*dot3+dot18*dot3+dot19*dot3+
				   4.*dot16*dot2*dot4+4.*dot26*dot4-4.*dot16*dot3*dot4+8.*dot11*dot14*dot5+8.*dot11*dot15*dot5-
				   2.*dot17*dot2*dot5-2.*dot27*dot5+2.*dot17*dot3*dot5-8.*dot14*dot4*dot5-8.*dot15*dot4*dot5+
				   dot18*dot6+dot19*dot6+2.*dot17*dot5*dot6+dot18*dot7+dot19*dot7-2.*dot17*dot5*dot7-
				   8.*dot11*dot15*dot8+2.*dot27*dot8-4.*dot17*dot3*dot8+8.*dot15*dot4*dot8-4.*dot17*dot6*dot8-
				   2.*dot16*dot2*dot9-4.*dot26*dot9+6.*dot16*dot3*dot9+8.*dot14*dot5*dot9+4.*dot15*dot5*dot9+
				   2.*dot16*dot6*dot9+2.*dot16*dot7*dot9-4.*dot15*dot8*dot9-
				   dot1*(dot18+dot19-2.*dot17*dot5-2.*dot16*(2.*dot11-2.*dot4+dot9))+
				   2.*(dot5-dot8)*(dot23+2.*dot13*(dot11-dot4+dot9))*M-
				   2.*dot10*(-dot27+2.*(dot14+dot15)*(dot11-dot4)+2.*dot14*dot9+dot23*M+dot17*(dot1-dot2+dot3+dot6-dot7-M2))+
				   (-4.*dot12*dot14+5.*dot18+dot19+4.*dot16*dot4-6.*dot17*dot5+4.*dot17*dot8+2.*dot16*dot9)*M2)+
			       2.*(-2.*dot18*dot2-dot19*dot2+dot10*dot27-2.*dot12*dot15*dot3+dot19*dot3+2.*dot10*dot15*dot4-
				   2.*dot16*dot2*dot4-dot26*dot4-2.*dot11*dot15*dot5+2.*dot17*dot2*dot5+dot27*dot5-
				   2.*dot17*dot3*dot5+2.*dot14*dot4*dot5+2.*dot15*dot4*dot5-2.*dot17*dot5*dot6-
				   2.*dot18*dot7-2.*dot16*dot4*dot7+2.*dot11*dot15*dot8-2.*dot27*dot8+2.*dot17*dot3*dot8-
				   2.*dot15*dot4*dot8+2.*dot17*dot6*dot8+2.*dot26*dot9-2.*dot14*dot5*dot9+
				   dot1*(-2.*dot12*dot14-2.*dot10*dot17+dot18+2.*dot16*dot4-2.*dot16*dot9)+
				   dot23*(-dot5+dot8)*M-2.*(dot18+dot16*dot4)*M2+2.*dot12*dot14*(dot2+dot7+M2)))
	    /(a1*a2*(-dot1+dot2+dot3+dot6+dot7+M2)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2)*(-dot7-a2*M2));
	  diag[18]=-0.125*M2*(-2.*dot1*(dot18+dot19+2.*dot17*(-(a2*dot10)-a1*dot5-a2*dot8)+2.*dot16*dot9)+
			      2.*(dot19*dot2-2.*a1*dot17*dot2*dot5-a1*dot27*dot5+2.*dot17*dot5*dot7-2.*a1*dot17*dot5*dot7+
				  dot18*(dot2+dot7)-2.*dot17*dot2*dot8+2.*a1*dot17*dot2*dot8+a1*dot27*dot8-2.*dot17*dot7*dot8+
				  2.*a1*dot17*dot7*dot8+2.*(dot16*dot2+a1*dot14*(dot5-dot8))*dot9-
				  a2*dot10*(dot27+2.*dot17*(dot2+dot7)-2.*dot14*dot9))+
			      (dot28-2.*dot22*dot5-2.*a2*dot23*(dot10-dot5+dot8)+2.*(dot20+2.*a2*dot13*(dot10-dot5+dot8))*dot9)*M+
			      2.*a2*(dot18-2.*a2*dot17*(dot10-dot5+dot8))*M2)
	    /(a2*(dot1-dot2-dot3)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot7-a2*M2));
	  diag[19]=-0.125*M2*(4.*dot1*dot12*dot14-2.*dot1*dot18+4.*dot12*dot15*dot2+4.*dot11*dot16*dot2-2.*dot19*dot2-4.*dot12*dot14*dot3+
			      2.*dot18*dot3-4.*dot1*dot16*dot4+4.*dot16*dot3*dot4-4.*dot11*dot14*dot5+4.*a1*dot11*dot14*dot5-2.*dot27*dot5+
			      2.*a1*dot27*dot5+4.*dot17*dot3*dot5-4.*a1*dot17*dot3*dot5-4.*dot15*dot4*dot5+4.*a1*dot15*dot4*dot5-
			      4.*dot12*dot14*dot6+2.*dot18*dot6+4.*dot16*dot4*dot6+4.*dot17*dot5*dot6-4.*a1*dot17*dot5*dot6-
			      4.*a1*dot11*dot14*dot8+4.*dot1*dot17*dot8-2.*a1*dot27*dot8-4.*dot17*dot3*dot8+4.*a1*dot17*dot3*dot8-
			      4.*a1*dot15*dot4*dot8-4.*dot17*dot6*dot8+4.*a1*dot17*dot6*dot8+
			      (2.*dot11*dot20-dot28+2.*dot21*dot4+2.*(dot12*dot24+dot1*(-2.*dot12*dot13+dot25)-a2*(2.*dot11*dot13+dot23)*(dot5-dot8)))*M+
			      2.*a1*(-2.*dot12*dot14+dot18+2.*dot16*dot4+2.*a2*dot17*(dot5-dot8))*M2+
			      2.*a2*dot10*(dot27+2.*dot15*dot4+dot23*M+2.*dot11*(dot14+dot13*M)-2.*dot17*(dot3+dot6+a1*M2)))
	    /(a1*a2*(dot1-dot2-dot3)*(dot6+a1*M2)*(-dot1+dot3+dot6+a1*M2));
	  diag[20]=-0.0625*M2*(-4.*dot1*dot11*dot16+4.*a1*dot1*dot11*dot16-2.*a1*dot1*dot18-2.*a1*dot1*dot19+4.*dot11*dot16*dot2-
			       4.*a1*dot11*dot16*dot2+2.*a1*dot18*dot2+2.*a1*dot19*dot2+2.*dot11*dot26-2.*a1*dot11*dot26-
			       2.*dot18*dot3+2.*a1*dot18*dot3-2.*dot19*dot3+2.*a1*dot19*dot3-4.*a1*dot1*dot16*dot4+
			       4.*a1*dot16*dot2*dot4+2.*a1*dot26*dot4-4.*dot11*dot14*dot5+8.*a1*dot11*dot14*dot5-
			       4.*a12*dot11*dot14*dot5-4.*dot11*dot15*dot5+8.*a1*dot11*dot15*dot5-4.*a12*dot11*dot15*dot5+
			       4.*a1*dot1*dot17*dot5-4.*a1*dot17*dot2*dot5-2.*a1*dot27*dot5-4.*a1*dot14*dot4*dot5+
			       4.*a12*dot14*dot4*dot5-4.*a1*dot15*dot4*dot5+4.*a12*dot15*dot4*dot5+2.*dot18*dot7+
			       4.*dot17*dot5*dot7-4.*a1*dot17*dot5*dot7-4.*a1*dot11*dot14*dot8+4.*a12*dot11*dot14*dot8+
			       4.*dot11*dot15*dot8-8.*a1*dot11*dot15*dot8+4.*a12*dot11*dot15*dot8+2.*a1*dot27*dot8+
			       4.*dot17*dot3*dot8-4.*a1*dot17*dot3*dot8-4.*a12*dot14*dot4*dot8+4.*a1*dot15*dot4*dot8-
			       4.*a12*dot15*dot4*dot8-4.*dot17*dot7*dot8+4.*a1*dot17*dot7*dot8-2.*a1*dot26*dot9-
			       4.*dot16*dot3*dot9+4.*a1*dot16*dot3*dot9+8.*a1*dot14*dot5*dot9-4.*a12*dot14*dot5*dot9+
			       4.*a1*dot15*dot5*dot9-4.*a12*dot15*dot5*dot9-4.*a1*dot14*dot8*dot9+4.*a12*dot14*dot8*dot9-
			       4.*a1*dot15*dot8*dot9+4.*a12*dot15*dot8*dot9+
			       (4.*a2*dot1*dot12*dot13+2.*(-2.+a1)*dot1*dot25+2.*dot2*dot25+3.*dot28+dot12*(-2.*dot24-4.*a2*dot13*(dot2+dot3))+
				2.*(dot11*dot20-a1*dot11*dot20-a1*dot2*dot25+dot25*dot3-a1*dot25*dot3+a1*dot20*dot4-dot21*dot4-
				    (2.*a22*dot11*dot13+2.*dot22+dot23-a1*dot23-2.*(1.-a1*a2)*dot13*dot4)*dot5-
				    a2*(-2.*a2*dot11*dot13-dot23-2.*a1*dot13*dot4)*dot8+a2*(dot20-2.*a2*dot13*(dot5-dot8))*dot9))*M+
			       2.*a2*(dot18+2.*a2*dot17*(dot5-dot8))*M2-
			       2.*a2*dot10*(dot27-dot23*M-2.*dot9*(dot14+dot13*M)-2.*a2*dot11*(dot14+dot15+dot13*M)-
					    2.*a1*(dot4-dot9)*(dot14+dot15+dot13*M)+2.*dot17*(-dot1+dot2+dot7+M2-a1*M2)))
	    /(a1*a2*(dot1-dot2-dot3)*(dot6+a1*M2)*(dot1-a2*(dot6+dot7+M2)));
	  diag[21]=M2*(4.*dot10*dot11*dot15+4.*dot12*dot14*dot2-6.*dot10*dot17*dot2-3.*dot18*dot2-dot19*dot2-2.*dot11*dot26-
		       4.*dot12*dot15*dot3+2.*dot10*dot17*dot3+dot18*dot3+3.*dot19*dot3-4.*dot16*dot2*dot4+
		       4.*dot11*dot14*dot5-4.*dot11*dot15*dot5+6.*dot17*dot2*dot5+2.*dot27*dot5-6.*dot17*dot3*dot5+
		       4.*dot15*dot4*dot5+4.*dot10*dot17*dot6+2.*dot18*dot6-4.*dot17*dot5*dot6+4.*dot12*dot14*dot7-
		       4.*dot10*dot17*dot7-2.*dot18*dot7-4.*dot16*dot4*dot7+4.*dot17*dot5*dot7+4.*dot11*dot15*dot8-
		       4.*dot17*dot2*dot8-4.*dot27*dot8+4.*dot17*dot3*dot8-4.*dot15*dot4*dot8+4.*dot17*dot6*dot8-
		       4.*dot17*dot7*dot8+4.*dot10*dot14*dot9+2.*dot16*dot2*dot9+4.*dot26*dot9-2.*dot16*dot3*dot9-4.*dot14*dot5*dot9+
		       dot1*(-4.*dot12*dot14+dot18-dot19+4.*dot16*dot4+2.*dot17*(dot10-dot5+2.*dot8)-6.*dot16*dot9)+
		       2.*(dot20+2.*dot13*(dot10-dot5+dot8))*dot9*M-2.*(-2.*dot12*dot14+dot18+2.*dot16*dot4+2.*dot17*(dot10-dot5+dot8))*M2+
		       4.*a12*(dot10-dot5+dot8)*((dot14+dot15)*(dot11-dot4+dot9)+dot13*(dot11-dot4+dot9)*M-2.*dot17*M2)+
		       a1*(-4.*dot11*dot16*dot2+dot18*dot2+dot19*dot2-2.*dot11*dot26+dot18*dot3+dot19*dot3+
			   4.*dot16*dot2*dot4+2.*dot26*dot4+4.*dot11*dot14*dot5+8.*dot11*dot15*dot5-6.*dot17*dot2*dot5-
			   4.*dot27*dot5+2.*dot17*dot3*dot5-4.*dot14*dot4*dot5-8.*dot15*dot4*dot5+4.*dot17*dot5*dot6-
			   4.*dot17*dot5*dot7-8.*dot11*dot15*dot8+4.*dot17*dot2*dot8+4.*dot27*dot8-4.*dot17*dot3*dot8+
			   8.*dot15*dot4*dot8-4.*dot17*dot6*dot8+4.*dot17*dot7*dot8-2.*dot16*dot2*dot9-2.*dot26*dot9+
			   2.*dot16*dot3*dot9+8.*dot14*dot5*dot9+4.*dot15*dot5*dot9-4.*dot14*dot8*dot9-
			   4.*dot15*dot8*dot9-dot1*(dot18+dot19-6.*dot17*dot5+4.*dot17*dot8-2.*dot16*(2.*dot11-2.*dot4+dot9))-
			   2.*((dot11-dot4)*(dot20-2.*dot13*dot5+2.*dot13*dot8)+(dot20-4.*dot13*dot5+4.*dot13*dot8)*dot9)*M-
			   2.*dot10*(2.*dot11*(dot14+2.*dot15+dot13*M)-
				     2.*(dot27+dot14*dot4+2.*dot15*dot4-2.*dot14*dot9-dot15*dot9+dot13*(dot4-2.*dot9)*M)+
				     dot17*(3.*dot1-3.*dot2+dot3+2.*dot6-2.*dot7-6.*M2))+
			   4.*(-(dot12*dot14)+dot18+dot16*dot4-3.*dot17*dot5+3.*dot17*dot8)*M2))
	    /(16.*a1*a2*(dot1-dot2-dot3)*(dot6+a1*M2)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2));
	  diag[22]=-0.125*M2*(-4.*dot12*(dot15*dot2+dot14*(dot2+dot7))+
			      2.*(dot19*dot2+2.*a1*dot16*dot2*dot4+a1*dot26*dot4+2.*dot16*dot4*dot7
				  +dot18*(dot2+dot7)-2.*dot17*dot2*dot8+2.*a1*dot15*dot4*dot8-2.*dot17*dot7*dot8-
				  (-2.*a2*dot16*dot2+a1*(dot26+2.*dot15*dot8))*dot9)+
			      (-2.*dot12*dot24+dot28-2.*dot21*dot4+2.*dot8*(dot23-2.*a2*dot13*(dot4-dot9))+2.*a2*dot20*(dot4-dot9))*M+
			      dot1*(-4.*a2*dot11*dot16+4.*a1*dot16*(-dot4+dot9)+4.*dot12*dot13*M-2.*dot25*M)+
			      2.*a2*dot11*(2.*dot16*dot2+dot26-dot20*M+2.*dot8*(dot15+dot13*M))+
			      2.*a2*(-2.*dot12*dot14+dot18+2.*dot16*dot4-2.*dot17*dot8)*M2)
	    /(a2*dot2*(-dot1+dot2+dot3+dot6+dot7+M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	  diag[23]=-0.125*M2*(-2.*(dot19*dot3+dot26*dot4-a1*dot26*dot4-2.*dot16*dot3*dot4+2.*a1*dot16*dot3*dot4-
				   2.*dot17*dot3*dot5-2.*dot14*dot4*dot5+2.*a1*dot14*dot4*dot5-
				   a2*dot11*(dot26-2.*dot16*dot3-2.*dot14*dot5)+2.*dot16*dot3*dot9+a1*(dot26-2.*dot16*dot3-2.*dot14*dot5)*dot9)+
			      (dot28-2.*dot22*dot5-2.*a2*dot20*(dot11-dot4+dot9))*M+
			      2.*dot18*(-dot1+dot2+dot7+M2-a1*M2)+
			      2.*dot10*(-2.*a2*dot15*(dot11-dot4)+2.*dot17*(-dot1+dot2+dot7)+2.*a1*dot15*dot9-
					dot23*M-2.*a2*dot13*(dot11-dot4+dot9)*M+2.*a2*dot17*M2))
	    /(a1*a2*dot3*(-dot1+dot2+dot3+dot6+dot7+M2)*(dot6+a1*M2));
	  diag[24]=M2*(-4.*dot1*dot12*dot14+4.*a1*dot1*dot12*dot14-4.*dot1*dot12*dot15+4.*a1*dot1*dot12*dot15-4.*dot1*dot11*dot16+
		       4.*a1*dot1*dot11*dot16+2.*dot1*dot18-2.*a1*dot1*dot18+2.*dot1*dot19-2.*a1*dot1*dot19-4.*a1*dot12*dot14*dot2-
		       4.*a1*dot12*dot15*dot2+4.*dot11*dot16*dot2-4.*a1*dot11*dot16*dot2+2.*a1*dot18*dot2+2.*a1*dot19*dot2+
		       2.*dot11*dot26-2.*a1*dot11*dot26+4.*dot12*dot14*dot3-4.*a1*dot12*dot14*dot3+4.*dot12*dot15*dot3-
		       4.*a1*dot12*dot15*dot3-2.*dot18*dot3+2.*a1*dot18*dot3-2.*dot19*dot3+2.*a1*dot19*dot3-4.*a1*dot1*dot16*dot4+
		       4.*a1*dot16*dot2*dot4+2.*a1*dot26*dot4+4.*a1*dot11*dot14*dot5-4.*a12*dot11*dot14*dot5+
		       4.*a1*dot11*dot15*dot5-4.*a12*dot11*dot15*dot5+4.*a1*dot1*dot17*dot5-4.*a1*dot17*dot2*dot5-
		       2.*a1*dot27*dot5-4.*a1*dot14*dot4*dot5+4.*a12*dot14*dot4*dot5-4.*a1*dot15*dot4*dot5+
		       4.*a12*dot15*dot4*dot5+4.*dot12*dot14*dot6-4.*a1*dot12*dot14*dot6+4.*dot12*dot15*dot6-
		       4.*a1*dot12*dot15*dot6-2.*dot18*dot6+2.*a1*dot18*dot6-2.*dot19*dot6+2.*a1*dot19*dot6-
		       4.*a1*dot12*dot14*dot7+4.*dot12*dot15*dot7-4.*a1*dot12*dot15*dot7+2.*a1*dot18*dot7-
		       2.*dot19*dot7+2.*a1*dot19*dot7+4.*dot16*dot4*dot7-4.*a1*dot17*dot5*dot7-4.*a1*dot11*dot14*dot8+
		       4.*a12*dot11*dot14*dot8+4.*dot11*dot15*dot8-8.*a1*dot11*dot15*dot8+4.*a12*dot11*dot15*dot8-
		       4.*dot1*dot17*dot8+2.*a1*dot27*dot8+4.*dot17*dot3*dot8-4.*a1*dot17*dot3*dot8+4.*a1*dot14*dot4*dot8-
		       4.*a12*dot14*dot4*dot8+8.*a1*dot15*dot4*dot8-4.*a12*dot15*dot4*dot8+4.*dot17*dot6*dot8-
		       4.*a1*dot17*dot6*dot8+4.*dot1*dot16*dot9-2.*a1*dot26*dot9-4.*dot16*dot3*dot9+4.*a1*dot16*dot3*dot9+
		       4.*a1*dot14*dot5*dot9-4.*a12*dot14*dot5*dot9-4.*a12*dot15*dot5*dot9-4.*dot16*dot6*dot9+
		       4.*a1*dot16*dot6*dot9-4.*dot16*dot7*dot9+4.*a1*dot16*dot7*dot9-4.*a1*dot14*dot8*dot9+
		       4.*a12*dot14*dot8*dot9-4.*a1*dot15*dot8*dot9+4.*a12*dot15*dot8*dot9+
		       (4.*dot1*dot12*dot13-4.*dot12*dot24+2.*(-3.+a1)*dot1*dot25+2.*dot2*dot25+3.*dot28
			+2.*(-(dot22*dot5)+dot25*(dot3+dot6+dot7)+dot23*dot8-
			     dot4*(dot20+2.*dot21-2.*dot13*dot5+2.*dot13*dot8)+
			     a2*dot11*(dot20+2.*dot13*(a1*dot5+dot8-a1*dot8))-
			     a1*(dot25*(dot2+dot3+dot6+dot7)+dot23*(-dot5+dot8)-(dot20-2.*dot13*dot5+4.*dot13*dot8)*(dot4-dot9))+
			     2.*a12*dot13*(dot5-dot8)*(dot4-dot9)+(dot20+2.*dot13*dot8)*dot9))*M+
		       2.*a2*dot25*pow<3,1>(M)-2.*a2*(-2.*dot12*dot15+dot19+2.*a1*dot17*(dot5-dot8)+2.*dot16*(-dot4+dot9))*M2-
		       2.*a2*dot10*(dot27-2.*a2*(dot14+dot15)*(dot11-dot4)+2.*dot17*(-dot1+dot2+dot7)-
				    2.*a2*dot14*dot9+2.*a1*dot15*dot9-dot23*M-2.*a2*dot13*(dot11-dot4+dot9)*M+2.*a2*dot17*M2))
		    /(16.*a1*a2*(-dot1+dot2+dot3+dot6+dot7+M2)*(dot6+a1*M2)*(dot1-a2*(dot6+dot7+M2)));
	  diag[25]=M2*(2.*dot12*dot14*dot2-2.*a1*dot12*dot14*dot2+2.*dot12*dot15*dot2-2.*a1*dot12*dot15*dot2+6.*dot11*dot16*dot2-
		       6.*a1*dot11*dot16*dot2-dot18*dot2+a1*dot18*dot2-3.*dot19*dot2+a1*dot19*dot2-4.*a1*dot11*dot26-
		       2.*dot12*dot14*dot3-2.*a1*dot12*dot14*dot3-2.*dot12*dot15*dot3-2.*a1*dot12*dot15*dot3-
		       2.*dot11*dot16*dot3+2.*a1*dot11*dot16*dot3+3.*dot18*dot3+a1*dot18*dot3+dot19*dot3+a1*dot19*dot3-
		       6.*dot16*dot2*dot4+6.*a1*dot16*dot2*dot4-2.*dot26*dot4+4.*a1*dot26*dot4+6.*dot16*dot3*dot4-
		       2.*a1*dot16*dot3*dot4+8.*a1*dot11*dot14*dot5-4.*a12*dot11*dot14*dot5+4.*a1*dot11*dot15*dot5-
		       4.*a12*dot11*dot15*dot5+4.*dot17*dot2*dot5-4.*a1*dot17*dot2*dot5-2.*a1*dot27*dot5+4.*dot14*dot4*dot5-
		       8.*a1*dot14*dot4*dot5+4.*a12*dot14*dot4*dot5-4.*a1*dot15*dot4*dot5+4.*a12*dot15*dot4*dot5-
		       2.*dot12*dot14*dot6-2.*a1*dot12*dot14*dot6+2.*dot12*dot15*dot6-2.*a1*dot12*dot15*dot6+
		       2.*dot11*dot16*dot6-2.*a1*dot11*dot16*dot6+3.*dot18*dot6+a1*dot18*dot6-dot19*dot6+a1*dot19*dot6+
		       2.*dot16*dot4*dot6+2.*a1*dot16*dot4*dot6+2.*dot12*dot14*dot7-2.*a1*dot12*dot14*dot7+
		       2.*dot12*dot15*dot7-2.*a1*dot12*dot15*dot7+2.*dot11*dot16*dot7-2.*a1*dot11*dot16*dot7-
		       dot18*dot7+a1*dot18*dot7-dot19*dot7+a1*dot19*dot7-2.*dot16*dot4*dot7+2.*a1*dot16*dot4*dot7-
		       4.*a1*dot17*dot5*dot7-4.*a1*dot11*dot14*dot8+4.*a12*dot11*dot14*dot8+4.*dot11*dot15*dot8-
		       8.*a1*dot11*dot15*dot8+4.*a12*dot11*dot15*dot8-2.*dot17*dot2*dot8+2.*a1*dot17*dot2*dot8-
		       4.*dot27*dot8+2.*a1*dot27*dot8+2.*dot17*dot3*dot8-2.*a1*dot17*dot3*dot8+4.*a1*dot14*dot4*dot8-
		       4.*a12*dot14*dot4*dot8-4.*dot15*dot4*dot8+8.*a1*dot15*dot4*dot8-4.*a12*dot15*dot4*dot8+
		       2.*dot17*dot6*dot8-2.*a1*dot17*dot6*dot8-2.*dot17*dot7*dot8+2.*a1*dot17*dot7*dot8+
		       4.*dot16*dot2*dot9-4.*a1*dot16*dot2*dot9+4.*dot26*dot9-4.*a1*dot26*dot9-4.*dot16*dot3*dot9+
		       4.*a1*dot16*dot3*dot9-4.*dot14*dot5*dot9+8.*a1*dot14*dot5*dot9-4.*a12*dot14*dot5*dot9-
		       4.*a12*dot15*dot5*dot9-4.*a1*dot14*dot8*dot9+4.*a12*dot14*dot8*dot9-4.*a1*dot15*dot8*dot9+
		       4.*a12*dot15*dot8*dot9+
		       dot1*(-2.*a2*dot12*(dot14+dot15)-2.*dot11*dot16-dot18+dot19+2.*dot17*dot8+2.*dot16*(dot4-2.*dot9)+
			     a1*(6.*dot11*dot16-dot18-dot19-6.*dot16*dot4+4.*dot17*dot5-2.*dot17*dot8+4.*dot16*dot9))-
		       2.*(a1*(dot5-dot8)+dot8)*(-dot23-2.*a2*dot13*(dot11-dot4+dot9))*M-
		       (2.*dot12*((-1.+3.*a1)*dot14-a2*dot15)-2.*a2*dot11*dot16+dot18+dot19+
			2.*dot16*dot4-a1*(5.*dot18+dot19+6.*dot16*dot4-4.*a2*dot17*dot5)-
			2.*(-1.+2.*a1)*a2*dot17*dot8)*M2+
		       2.*dot10*(-2.*dot1*dot17+dot27-2.*dot14*dot4+2.*dot15*dot4+2.*dot14*dot9-
				 2.*a2*dot11*(-(a2*dot14)+a1*(dot15+dot13*M))-
				 2.*a12*((dot14+dot15)*(dot4-dot9)+dot13*(dot4-dot9)*M+dot17*M2)+
				 a1*(dot27+4.*dot14*dot4+2.*dot15*dot4-4.*dot14*dot9-
				     (dot23+2.*dot13*(-dot4+dot9))*M+2.*dot17*(-dot1+dot2+dot7+M2))))
	    /(16.*a1*a2*(-dot1+dot2+dot3+dot6+dot7+M2)*(dot6+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	  diag[26]=-0.125*M2*(4.*dot12*dot15*dot2+4.*dot11*dot16*dot2-2.*dot19*dot2+4.*dot16*dot3*dot4+4.*dot16*dot4*dot6-
			      4.*dot12*dot14*(dot3+dot6)+2.*dot18*(dot3+dot6)-4.*dot17*(dot3+dot6)*dot8+
			      (-2.*dot11*dot20-2.*dot12*dot24+dot28-2.*dot21*dot4+4.*dot11*dot13*dot8+2.*dot23*dot8)*M+
			      4.*dot1*dot12*(dot14+dot13*M)-2.*dot1*(dot18+2.*dot16*dot4-2.*dot17*dot8+dot25*M)+
			      2.*a1*(-2.*dot12*dot14+dot18+2.*dot16*dot4-2.*dot17*dot8)*M2)
	    /(dot2*(-dot1+dot3+dot6+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	  diag[27]=-0.125*M2*(2.*dot7*(2.*dot12*dot14-dot18-2.*dot16*dot4+2.*dot17*dot8)+
			      (-2.*dot12*(2.*dot13*dot2+dot24)+dot1*(4.*dot12*dot13-2.*dot25)+2.*dot2*dot25+dot28-
			       2.*(dot20+dot21)*dot4+2.*dot23*dot8+4.*dot13*dot8*(dot4-dot9)+2.*dot20*dot9)*M+
			      2.*a2*(2.*dot12*dot14-dot18-2.*dot16*dot4+2.*dot17*dot8)*M2)/(a2*dot2*sqr(dot7+a2*M2));
	  diag[28]=-0.125*M2*(-2.*dot12*dot14*dot2-6.*dot12*dot15*dot2+2.*dot11*dot16*dot2-4.*a1*dot11*dot16*dot2+
			      dot18*dot2+3.*dot19*dot2+2.*dot11*dot26-2.*a1*dot11*dot26+4.*dot12*dot14*dot3-
			      2.*dot18*dot3-2.*dot16*dot2*dot4+4.*a1*dot16*dot2*dot4+2.*a1*dot26*dot4-
			      4.*dot16*dot3*dot4+4.*dot12*dot14*dot6-2.*dot18*dot6-4.*dot16*dot4*dot6-
			      4.*dot12*dot14*dot7+2.*dot18*dot7+4.*dot16*dot4*dot7+4.*dot11*dot15*dot8-
			      4.*a1*dot11*dot15*dot8-2.*dot17*dot2*dot8+4.*dot17*dot3*dot8+4.*a1*dot15*dot4*dot8+
			      4.*dot17*dot6*dot8-4.*dot17*dot7*dot8+4.*dot16*dot2*dot9-4.*a1*dot16*dot2*dot9-
			      2.*a1*dot26*dot9-4.*a1*dot15*dot8*dot9+
			      2.*dot1*(-2.*dot12*dot14+dot18-2.*dot17*dot8+2.*dot16*(-(a2*(dot11-dot4))+a1*dot9))+
			      2.*(dot20-2.*dot13*dot8)*(dot4-dot9+a1*(dot11-dot4+dot9))*M+
			      2.*(-1.+2.*a1)*(2.*dot12*dot14-dot18-2.*dot16*dot4+2.*dot17*dot8)*M2)
	    /(a2*dot2*(-dot7-a2*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	  diag[29]=-0.125*M2*(-2.*dot1*(2.*dot10*dot17+dot18)+4.*dot10*dot17*(dot2+dot7)+2.*dot18*(dot2+dot7)-
			      2.*dot3*(dot19-2.*dot17*dot5+2.*dot16*dot9)+
			      (2.*dot10*dot23-dot28+2.*dot22*dot5-2.*(2.*dot10*dot13+dot20)*dot9)*M+
			      2.*a2*(2.*dot10*dot17+dot18)*M2)
	    /(dot3*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot7-a2*M2));
	  diag[30]=-0.125*M2*(2.*dot7*(dot19-2.*dot17*dot5+2.*dot16*dot9)+
			      (dot28-2.*dot22*dot5+2.*dot25*(-dot1+dot2+dot7)+
			       2.*dot20*dot9-2.*(dot5-dot8)*(dot23-2.*dot13*dot9))*M+
			      2.*a2*dot25*pow<3,1>(M)+2.*a2*(dot19-2.*dot17*dot5+2.*dot16*dot9)*M2)
	    /(a2*(dot1-dot2-dot7-a2*M2)*sqr(dot7+a2*M2));
	  diag[31]=-0.125*M2*(2.*dot10*dot17*dot2-4.*a1*dot10*dot17*dot2-3.*dot18*dot2-dot19*dot2+2.*dot10*dot27-
			      2.*a1*dot10*dot27+2.*dot19*dot3-2.*dot17*dot2*dot5+4.*a1*dot17*dot2*dot5+
			      2.*a1*dot27*dot5-4.*dot17*dot3*dot5+2.*dot10*dot17*dot7-4.*a1*dot10*dot17*dot7-
			      3.*dot18*dot7+dot19*dot7-6.*dot17*dot5*dot7+4.*a1*dot17*dot5*dot7+
			      4.*dot17*dot2*dot8-4.*a1*dot17*dot2*dot8-2.*a1*dot27*dot8+4.*dot17*dot7*dot8-
			      4.*a1*dot17*dot7*dot8-4.*dot10*dot14*dot9+4.*a1*dot10*dot14*dot9-2.*dot16*dot2*dot9+
			      4.*dot16*dot3*dot9-4.*a1*dot14*dot5*dot9+2.*dot16*dot7*dot9+4.*a1*dot14*dot8*dot9+
			      dot1*(3.*dot18+dot19+2.*(-1.+2.*a1)*dot17*(dot10-dot5)-4.*a2*dot17*dot8+2.*dot16*dot9)-
			      2.*(dot5-dot8+a1*(dot10-dot5+dot8))*(dot23-2.*dot13*dot9)*M-
			      a2*(3.*dot18-dot19+2.*dot17*(-dot10+2.*a1*dot10+3.*dot5-2.*a1*dot5-2.*a2*dot8)-2.*dot16*dot9)*M2)
	    /(a2*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot7-a2*M2)*(dot7+a2*M2));
	  diag[32]=-0.125*M2*(2.*(2.*dot10*dot17+dot18)*dot6+
			      (-2.*dot11*dot20+dot28-2.*dot25*dot3+2.*dot20*dot4-2.*dot10*(2.*dot11*dot13+dot23-2.*dot13*dot4)-2.*dot22*dot5)*M+
			      2.*a1*(2.*dot10*dot17+dot18)*M2)/(a1*dot3*sqr(dot6+a1*M2));
	  diag[33]=M2*(-2.*dot1*dot18+2.*dot18*dot2-dot18*dot3-3.*dot19*dot3-2.*a2*dot26*dot4-4.*a1*dot16*dot3*dot4+4.*a1*dot11*dot14*dot5+
		       6.*dot17*dot3*dot5+4.*dot14*dot4*dot5-4.*a1*dot14*dot4*dot5+
		       2.*dot11*(dot26-a1*dot26+2.*a1*dot16*dot3-2.*dot14*dot5)-2.*dot18*dot6-
		       2.*dot10*dot17*(2.*dot1-2.*dot2+dot3+2.*dot6-2.*dot7)+2.*dot18*dot7-2.*a1*dot26*dot9-
		       2.*dot16*dot3*dot9+4.*a1*dot16*dot3*dot9+4.*a1*dot14*dot5*dot9-2.*dot20*dot9*M+
		       2.*a1*dot20*(dot11-dot4+dot9)*M+2.*(1.-2.*a1)*dot18*M2+
		       4.*dot10*(-(a2*dot15*(dot11-dot4))+a1*dot15*dot9-dot13*dot9*M+a1*dot13*(dot11-dot4+dot9)*M+(1.-2.*a1)*dot17*M2))
	    /(8.*a1*dot3*(dot6+a1*M2)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2));
	  diag[34]=M2*(-2.*dot12*(dot24-2.*dot13*dot3)*M-2.*dot11*(2.*dot10*dot13+dot20-2.*dot13*dot5)*M+
		       (-2.*dot10*dot23+dot28-2.*dot25*dot3-2.*dot21*dot4+2.*dot23*dot5)*M+
		       4.*dot11*dot16*(dot6+a1*M2)+4.*dot12*(dot15+dot13*M)*(dot6+a1*M2)-
		       2.*(dot19+dot25*M)*(dot6+a1*M2))/(8.*a1*sqr(dot6+a1*M2)*(-dot1+dot3+dot6+a1*M2));
	  diag[35]=M2*(4.*dot12*dot15*dot2+4.*dot11*dot16*dot2-2.*dot19*dot2-6.*dot12*dot14*dot3-2.*dot12*dot15*dot3-
		       2.*dot11*dot16*dot3+3.*dot18*dot3+dot19*dot3+6.*dot16*dot3*dot4-4.*dot11*dot14*dot5+4.*a1*dot11*dot14*dot5-
		       2.*dot27*dot5+2.*a1*dot27*dot5-4.*a1*dot17*dot3*dot5-4.*dot15*dot4*dot5+4.*a1*dot15*dot4*dot5-
		       6.*dot12*dot14*dot6+2.*dot12*dot15*dot6+2.*dot11*dot16*dot6+3.*dot18*dot6-dot19*dot6+
		       6.*dot16*dot4*dot6-4.*a1*dot17*dot5*dot6-4.*a1*dot11*dot14*dot8-2.*a1*dot27*dot8-2.*dot17*dot3*dot8+
		       4.*a1*dot17*dot3*dot8-4.*a1*dot15*dot4*dot8-2.*dot17*dot6*dot8+4.*a1*dot17*dot6*dot8+
		       dot1*(6.*dot12*dot14+2.*dot12*dot15+2.*dot11*dot16-3.*dot18-dot19-6.*dot16*dot4+4.*dot17*dot5+2.*dot17*dot8)+
		       2.*(2.*dot11*dot13+dot23)*(a1*dot5+dot8-a1*dot8)*M+
		       a1*(2.*dot12*(-3.*dot14+dot15)+3.*dot18-dot19+2.*dot16*(dot11+3.*dot4)-2.*dot17*(2.*a1*dot5+dot8-2.*a1*dot8))*M2+
		       2.*dot10*(-2.*dot1*dot17+dot27+2.*dot15*dot4-2.*dot11*(-(a2*dot14)+a1*dot13*M)+
				 a1*(-dot27-2.*dot15*dot4-dot23*M+2.*dot17*(dot3+dot6+a1*M2))))
	    /(8.*a1*(dot6+a1*M2)*(-dot1+dot3+dot6+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	  diag[36]=-0.0625*M2*(-2.*dot12*dot14*dot2+2.*dot12*dot15*dot2-2.*dot11*dot16*dot2+dot18*dot2+dot19*dot2-4.*dot12*dot14*dot3+4.*dot10*dot14*dot4+
			       6.*dot16*dot2*dot4+2.*dot26*dot4-4.*dot11*dot14*dot5-4.*dot17*dot2*dot5-2.*dot27*dot5-4.*dot14*dot4*dot5-4.*dot15*dot4*dot5-
			       4.*dot12*dot14*dot6+4.*dot16*dot4*dot6-2.*dot12*dot14*dot7-2.*dot12*dot15*dot7-2.*dot11*dot16*dot7+dot18*dot7+dot19*dot7+
			       2.*dot16*dot4*dot7-4.*dot11*dot15*dot8+2.*dot17*dot2*dot8+4.*dot27*dot8-4.*dot17*dot3*dot8+4.*dot15*dot4*dot8-
			       4.*dot17*dot6*dot8+2.*dot17*dot7*dot8-4.*dot10*dot14*dot9-4.*dot16*dot2*dot9-4.*dot26*dot9+4.*dot16*dot3*dot9+
			       4.*dot14*dot5*dot9+dot1*(4.*dot12*(2.*dot14+dot15)-2.*(dot18+dot19-2.*dot17*dot5)+4.*dot16*(dot11-2.*dot4+dot9))+
			       4.*dot13*dot8*(dot4-dot9)*M+(-2.*dot12*(dot14+dot15)-2.*dot11*dot16+dot18+dot19+2.*dot16*dot4+2.*dot17*dot8)*M2+
			       a1*(2.*dot12*dot14*dot2+2.*dot12*dot15*dot2+6.*dot11*dot16*dot2-dot18*dot2-dot19*dot2+4.*dot11*dot26+2.*dot12*dot14*dot3+
				   2.*dot12*dot15*dot3-2.*dot11*dot16*dot3-dot18*dot3-dot19*dot3-6.*dot16*dot2*dot4-4.*dot26*dot4+2.*dot16*dot3*dot4-
				   4.*dot11*dot14*dot5-4.*dot11*dot15*dot5+4.*dot17*dot2*dot5+4.*dot27*dot5-4.*dot17*dot3*dot5+8.*dot14*dot4*dot5+
				   8.*dot15*dot4*dot5-4.*dot17*dot5*dot6+2.*dot12*dot14*dot7+2.*dot12*dot15*dot7+2.*dot11*dot16*dot7-dot18*dot7-
				   dot19*dot7-2.*dot16*dot4*dot7+4.*dot17*dot5*dot7+8.*dot11*dot15*dot8-2.*dot17*dot2*dot8-4.*dot27*dot8+
				   6.*dot17*dot3*dot8-4.*dot14*dot4*dot8-12.*dot15*dot4*dot8+4.*dot17*dot6*dot8-2.*dot17*dot7*dot8+
				   4.*dot16*dot2*dot9+4.*dot26*dot9-4.*dot16*dot3*dot9-8.*dot14*dot5*dot9+4.*dot14*dot8*dot9+4.*dot15*dot8*dot9+
				   dot1*(-2.*dot12*(dot14+dot15)-6.*dot11*dot16+dot18+dot19+6.*dot16*dot4-4.*dot17*dot5+2.*dot17*dot8-4.*dot16*dot9)+
				   4.*dot13*(dot11*dot8+(dot5-2.*dot8)*(dot4-dot9))*M+
				   4.*dot10*(dot11*(dot14+dot15)-dot27-2.*dot14*dot4-2.*dot15*dot4+2.*dot14*dot9-dot13*dot4*M+dot13*dot9*M+dot17*(dot1-dot2+dot3+dot6-dot7-M2))+
				   2.*(2.*dot12*dot15+2.*dot11*dot16-dot18-dot19+2.*dot17*dot5-4.*dot17*dot8)*M2)+
			       a12*(-4.*(dot14+dot15)*(dot10-dot5+dot8)*(dot11-dot4+dot9)-4.*dot13*(dot10-dot5+dot8)*(dot11-dot4+dot9)*M+
				    (-2.*dot12*(dot14+dot15)-2.*dot11*dot16+8.*dot10*dot17+dot18+dot19+2.*dot16*dot4-8.*dot17*dot5+6.*dot17*dot8)*M2))
	    /(a1*a2*(dot6+a1*M2)*(-dot7-a2*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	  diag[37]=-0.0625*M2*(4.*dot12*dot14*dot2-2.*dot10*dot17*dot2-dot18*dot2-dot19*dot2-4.*dot12*dot15*dot3+4.*dot10*dot15*dot4-4.*dot16*dot2*dot4-
			       2.*dot26*dot4-4.*dot11*dot15*dot5+6.*dot17*dot2*dot5+2.*dot27*dot5+4.*dot14*dot4*dot5+4.*dot15*dot4*dot5-4.*dot17*dot5*dot6+
			       4.*dot12*dot14*dot7-4.*dot16*dot4*dot7+4.*dot17*dot5*dot7+4.*dot11*dot15*dot8-4.*dot17*dot2*dot8-4.*dot27*dot8+4.*dot17*dot3*dot8-
			       4.*dot15*dot4*dot8+4.*dot17*dot6*dot8-4.*dot17*dot7*dot8+4.*dot10*dot14*dot9+2.*dot16*dot2*dot9+4.*dot26*dot9-4.*dot16*dot3*dot9-
			       4.*dot14*dot5*dot9-dot1*(4.*dot12*dot14+dot18+dot19-4.*dot16*dot4+2.*dot17*(dot10+dot5-2.*dot8)+6.*dot16*dot9)+
			       4.*dot13*(-dot5+dot8)*dot9*M+4.*(dot12*dot14-dot16*dot4+dot17*(dot5-dot8))*M2+
			       a12*(4.*(dot14+dot15)*(dot10-dot5+dot8)*(dot11-dot4+dot9)+4.*dot13*(dot10-dot5+dot8)*(dot11-dot4+dot9)*M+
				    (-6.*dot10*dot17+dot18+dot19+6.*dot17*dot5-8.*dot17*dot8+2.*dot16*dot9)*M2)+
			       a1*(-4.*dot11*dot16*dot2+dot18*dot2+dot19*dot2-4.*dot11*dot26+4.*dot11*dot16*dot3+dot18*dot3+dot19*dot3+
				   4.*dot16*dot2*dot4+4.*dot26*dot4-4.*dot16*dot3*dot4+8.*dot11*dot14*dot5+8.*dot11*dot15*dot5-
				   6.*dot17*dot2*dot5-4.*dot27*dot5+2.*dot17*dot3*dot5-8.*dot14*dot4*dot5-8.*dot15*dot4*dot5+
				   dot18*dot6+dot19*dot6+2.*dot17*dot5*dot6-4.*dot17*dot5*dot7-8.*dot11*dot15*dot8+4.*dot17*dot2*dot8+
				   4.*dot27*dot8-4.*dot17*dot3*dot8+8.*dot15*dot4*dot8-4.*dot17*dot6*dot8+4.*dot17*dot7*dot8-
				   2.*dot16*dot2*dot9-4.*dot26*dot9+6.*dot16*dot3*dot9+12.*dot14*dot5*dot9+4.*dot15*dot5*dot9+
				   2.*dot16*dot6*dot9-4.*dot14*dot8*dot9-4.*dot15*dot8*dot9-
				   dot1*(dot18+dot19-6.*dot17*dot5+4.*dot17*dot8-2.*dot16*(2.*dot11-2.*dot4+dot9))+
				   4.*dot13*(dot5-dot8)*(dot11-dot4+2.*dot9)*M+4.*(-(dot12*dot14)+dot16*dot4-3.*dot17*dot5+3.*dot17*dot8)*M2-
				   2.*dot10*(2.*dot11*(dot14+dot15)-2.*(dot27+dot14*dot4+dot15*dot4-2.*dot14*dot9-dot13*dot9*M)+
					     dot17*(3.*dot1-3.*dot2+dot3+dot6-2.*(dot7+M2)))))
	    /(a1*a2*(dot6+a1*M2)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2)*(dot7+a2*M2));
	  // colour flows
	  Complex flow[3] = {diag[0] - diag[8] + diag[11] + diag[12] + diag[13] - diag[16] + diag[17] - diag[20] + diag[21] + diag[24] + diag[25] + diag[36] + diag[37],
			     diag[1] + diag[3] + diag[5] + diag[7] + diag[29] + diag[30] + diag[31] + diag[32] + diag[33] + 2.*diag[37] +
			     (diag[4] - diag[5] + diag[6] - diag[7] + diag[9] - 2.*diag[11] - 2.*diag[13] + diag[15] + diag[18] + diag[23] - 2.*diag[24] - diag[29] - diag[30] - diag[32])/9.,
			     diag[2] - diag[3] - diag[5] - diag[7] + diag[26] + diag[27] + diag[28] + diag[34] + diag[35] + 2.*diag[36] +
			     (-diag[4] + diag[5] - diag[6] + diag[7] + 2.*diag[8] + diag[10] + diag[14] + 2.*diag[16] + diag[19] + 2.*diag[20] + diag[22] - diag[26] - diag[27] - diag[34])/9.};
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
  // return the answer
  meInfo(flows);
  return O1_*sHat()/M/72.*pow<4,1>(Constants::pi*standardModel()->alphaS(scale())/M)*me2Sum.real();
}
