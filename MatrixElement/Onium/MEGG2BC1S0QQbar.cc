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
    ("The MEGG2BC1S0QQbar class implements the matrix element for gg -> B_c(1S0) cbar b +cc");

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
  Energy2 dot5 = rescaledMomenta()[2]*rescaledMomenta()[3];
  Energy2 dot6 = rescaledMomenta()[2]*rescaledMomenta()[4];
  // matrix element
  ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half);
  Complex me2Sum=0.;
  vector<double> flows(3,0.);
  complex<Energy> dot4 [2] = {rescaledMomenta()[1]*g1[0].wave(),rescaledMomenta()[1]*g1[1].wave()};
  complex<Energy> dot8 [2] = {rescaledMomenta()[3]*g1[0].wave(),rescaledMomenta()[3]*g1[1].wave()};
  complex<Energy> dot10[2] = {rescaledMomenta()[4]*g1[0].wave(),rescaledMomenta()[4]*g1[1].wave()};
  complex<Energy> dot7 [2] = {rescaledMomenta()[2]*g2[0].wave(),rescaledMomenta()[2]*g2[1].wave()};
  complex<Energy> dot9 [2] = {rescaledMomenta()[3]*g2[0].wave(),rescaledMomenta()[3]*g2[1].wave()};
  complex<Energy> dot11[2] = {rescaledMomenta()[4]*g2[0].wave(),rescaledMomenta()[4]*g2[1].wave()};
  for(unsigned int ih4=0;ih4<2;++ih4) {
    for(unsigned int ih5=0;ih5<2;++ih5) {
      complex<Energy> dot13=v4[ih4].dimensionedWave().pseudoScalar(ubar5[ih5].dimensionedWave());
      LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
      complex<Energy2> dot14 = vec1*rescaledMomenta()[2];
      complex<Energy2> dot19 = vec1*rescaledMomenta()[1];
      LorentzVector<complex<Energy2> > vec4 = v4[ih4].dimensionedWave().slash(rescaledMomenta()[2]).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
      complex<Energy3> dot22 = vec4*rescaledMomenta()[1];
      for(unsigned int ih1=0;ih1<2;++ih1) {
	complex<Energy> dot15 = vec1*g1[ih1].wave();
	LorentzPolarizationVectorE vec3 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	complex<Energy2> dot21 = vec3*rescaledMomenta()[2];
	complex<Energy2> dot24 = vec3*rescaledMomenta()[1];
	LorentzVector<complex<Energy2> > vec7 = v4[ih4].dimensionedWave().slash(g1[ih1].wave()).slash(rescaledMomenta()[2]).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	complex<Energy3> dot27 = vec7*rescaledMomenta()[1];
	for(unsigned int ih2=0;ih2<2;++ih2) {
	  Complex dot12 = g1[ih1].wave()*g2[ih2].wave();
	  LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  LorentzPolarizationVectorE vec5 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(g1[ih1].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  LorentzVector<complex<Energy2> > vec6 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(rescaledMomenta()[2]).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	  LorentzVector<complex<Energy2> > vec8 = v4[ih4].dimensionedWave().slash(g2[ih2].wave()).slash(g1[ih1].wave()).slash(rescaledMomenta()[2]).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	  complex<Energy> dot16 = vec2*g1[ih1].wave();
	  complex<Energy2> dot17 = vec5*rescaledMomenta()[2];
	  complex<Energy> dot18 = vec1*g2[ih2].wave();
	  complex<Energy2> dot20 = vec2*rescaledMomenta()[2];
	  complex<Energy2> dot23 = vec2*rescaledMomenta()[1];
	  complex<Energy2> dot25 = vec5*rescaledMomenta()[1];
	  complex<Energy3> dot26 = vec6*rescaledMomenta()[1];
	  complex<Energy3> dot28 = vec8*rescaledMomenta()[1];
	  // diagrams
	  Complex diag[38];
	  diag[0]=(dot15*dot7[ih2]-dot18*(dot10[ih1]-dot4[ih1]+dot8[ih1])+dot12*(dot14+3.*dot13*M))*M2
	    /(8.*a1*a2*(dot5+a1*M2)*(dot6+a2*M2));
	  diag[1]=(-2.*dot12*dot14+3.*dot17+4.*dot15*dot7[ih2]+2.*dot18*(dot10[ih1]-dot4[ih1]+dot8[ih1])+3.*dot16*M)*M2
	    /(8.*a1*a2*(dot5+a1*M2)*(dot6+a2*M2));
	  diag[2]=(4.*dot12*dot14-3.*dot17-2.*dot15*dot7[ih2]-4.*dot18*(dot10[ih1]-dot4[ih1]+dot8[ih1])+6.*dot12*dot13*M-3.*dot16*M)*M2
	    /(8.*a1*a2*(dot5+a1*M2)*(dot6+a2*M2));
	  diag[3]=M2*(-2.*a1*dot12*dot14*dot2+2.*dot12*dot14*dot3-2.*a1*dot12*dot14*dot3+2.*dot11[ih2]*dot14*dot4[ih1]-2.*a1*dot11[ih2]*dot14*dot4[ih1]+
		      2.*dot12*dot14*dot5-a1*dot12*dot14*dot5-2.*dot11[ih2]*dot15*dot5-2.*dot12*dot19*dot5+2.*dot18*dot4[ih1]*dot5-
		      dot12*dot14*dot6-a1*dot12*dot14*dot6+2.*dot11[ih2]*dot15*dot6+2.*dot12*dot19*dot6-2.*dot18*dot4[ih1]*dot6-
		      2.*dot15*dot5*dot7[ih2]+2.*dot15*dot6*dot7[ih2]+2.*a1*dot11[ih2]*dot14*dot8[ih1]+2.*a1*dot14*dot7[ih2]*dot8[ih1]-2.*a1*dot14*dot4[ih1]*dot9[ih2]-
		      2.*dot15*dot5*dot9[ih2]+2.*dot15*dot6*dot9[ih2]+2.*a1*dot14*dot8[ih1]*dot9[ih2]-
		      2.*dot13*(dot12*((-1.+2.*a1)*dot2-dot3-dot5+a1*(2.*dot3+dot5+dot6))+(-1.+2.*a1)*(dot11[ih2]*(dot4[ih1]-dot8[ih1])+dot4[ih1]*dot9[ih2]-dot8[ih1]*(dot7[ih2]+dot9[ih2])))*M+
		      (-1.+2.*a1)*dot1*dot12*(dot14+2.*dot13*M)+2.*dot10[ih1]*(dot11[ih2]+dot7[ih2]+dot9[ih2])*(-(a2*dot14)+(-1.+2.*a1)*dot13*M)+
		      (-1.+2.*a1)*(dot12*(dot14-2.*dot19)+2.*dot18*dot4[ih1]-2.*dot15*(dot11[ih2]+dot7[ih2]+dot9[ih2]))*M2)
	    /(4.*a1*a2*dot1*(dot5+a1*M2)*(dot6+a2*M2));
	  diag[4]=M2*(-2.*a1*dot12*dot14*dot2+2.*dot12*dot14*dot3-2.*a1*dot12*dot14*dot3+2.*dot11[ih2]*dot14*dot4[ih1]-2.*a1*dot11[ih2]*dot14*dot4[ih1]+
		      dot12*dot14*dot5-a1*dot12*dot14*dot5-dot12*dot14*dot6-a1*dot12*dot14*dot6+2.*dot11[ih2]*dot15*dot6+2.*dot12*dot19*dot6-
		      2.*dot18*dot4[ih1]*dot6+2.*dot15*dot6*dot7[ih2]+2.*a1*dot11[ih2]*dot14*dot8[ih1]+2.*a1*dot14*dot7[ih2]*dot8[ih1]-2.*a1*dot14*dot4[ih1]*dot9[ih2]+2.*dot15*dot6*dot9[ih2]+
		      2.*a1*dot14*dot8[ih1]*dot9[ih2]-
		      (-(dot12*dot22)-dot20*dot4[ih1]+2.*dot12*dot13*((-1.+2.*a1)*dot2-dot3+a1*(2.*dot3+dot5+dot6))-
		       4.*a2*dot11[ih2]*dot13*(dot4[ih1]-dot8[ih1])+4.*a2*dot13*dot7[ih2]*dot8[ih1]-4.*a2*dot13*(dot4[ih1]-dot8[ih1])*dot9[ih2]+dot21*(dot11[ih2]+dot7[ih2]+dot9[ih2]))*M-
		      2.*a1*dot12*dot13*pow<3,1>(M)+(-1.+2.*a1)*dot1*dot12*(dot14+2.*dot13*M)-
		      2.*a2*dot10[ih1]*(dot11[ih2]+dot7[ih2]+dot9[ih2])*(dot14+2.*dot13*M)+(dot12*(-(a2*dot14)+dot19)-dot18*dot4[ih1]+dot15*(dot11[ih2]+dot7[ih2]+dot9[ih2]))*M2)
	    /(4.*a1*dot1*(dot5+a1*M2)*(dot1-a2*(dot5+dot6+M2)));
	  diag[5]=M2*(dot5*(dot12*(dot14-2.*dot19)+2.*dot18*dot4[ih1]-2.*dot15*(dot11[ih2]+dot7[ih2]+dot9[ih2]))+
		      (-(dot12*dot22)-dot20*dot4[ih1]+2.*dot12*dot13*dot5+2.*dot13*dot7[ih2]*dot8[ih1]+
		       dot11[ih2]*(dot21-2.*dot13*dot4[ih1]+2.*dot13*dot8[ih1])+2.*dot13*(-dot4[ih1]+dot8[ih1])*dot9[ih2]+
		       dot21*(dot7[ih2]+dot9[ih2])+2.*dot10[ih1]*dot13*(dot11[ih2]+dot7[ih2]+dot9[ih2]))*M+2.*a1*dot12*dot13*pow<3,1>(M)+
		      (dot12*dot19-dot18*dot4[ih1]+dot15*(dot11[ih2]+dot7[ih2]+dot9[ih2])+
		       a1*(dot12*(dot14-4.*dot19)+4.*dot18*dot4[ih1]-4.*dot15*(dot11[ih2]+dot7[ih2]+dot9[ih2])))*M2)
	    /(4.*a1*dot1*sqr(dot5+a1*M2));
	  diag[6]=-0.25*M2*(dot1*dot12*dot14-2.*a1*dot1*dot12*dot14+2.*a1*dot12*dot14*dot2-2.*dot12*dot14*dot3+2.*a1*dot12*dot14*dot3-
			    2.*dot11[ih2]*dot14*dot4[ih1]+2.*a1*dot11[ih2]*dot14*dot4[ih1]-2.*dot12*dot14*dot5+a1*dot12*dot14*dot5+2.*dot11[ih2]*dot15*dot5+
			    2.*dot12*dot19*dot5-2.*dot18*dot4[ih1]*dot5+a1*dot12*dot14*dot6+2.*dot15*dot5*dot7[ih2]-2.*a1*dot11[ih2]*dot14*dot8[ih1]-
			    2.*a1*dot14*dot7[ih2]*dot8[ih1]+2.*a1*dot14*dot4[ih1]*dot9[ih2]+2.*dot15*dot5*dot9[ih2]-2.*a1*dot14*dot8[ih1]*dot9[ih2]+
			    (-(dot12*dot22)-2.*dot11[ih2]*dot13*dot4[ih1]-dot20*dot4[ih1]+2.*dot11[ih2]*dot13*dot8[ih1]+2.*dot13*dot7[ih2]*dot8[ih1]-2.*dot13*dot4[ih1]*dot9[ih2]+
			     2.*dot13*dot8[ih1]*dot9[ih2]+dot21*(dot11[ih2]+dot7[ih2]+dot9[ih2])+
			     2.*a1*dot13*(dot12*(-2.*dot1+2.*dot2+2.*dot3+dot5+dot6)+2.*dot11[ih2]*(dot4[ih1]-dot8[ih1])+2.*dot4[ih1]*dot9[ih2]-2.*dot8[ih1]*(dot7[ih2]+dot9[ih2])))*M+
			    2.*a1*dot12*dot13*pow<3,1>(M)-
			    2.*dot10[ih1]*(dot11[ih2]+dot7[ih2]+dot9[ih2])*(-(a2*dot14)+(-1.+2.*a1)*dot13*M)+
			    (-(a1*dot12*dot14)+dot12*dot19-dot18*dot4[ih1]+dot15*(dot11[ih2]+dot7[ih2]+dot9[ih2]))*M2)
	    /(a2*dot1*(dot6+a2*M2)*(dot1-a1*(dot5+dot6+M2)));
	  diag[7]=-0.25*M2*(2.*dot18*dot4[ih1]*dot6-2.*dot15*dot6*(dot11[ih2]+dot7[ih2]+dot9[ih2])+dot20*dot4[ih1]*M-dot21*(dot11[ih2]+dot7[ih2]+dot9[ih2])*M+
			    (-3.+4.*a1)*(-(dot18*dot4[ih1])+dot15*(dot11[ih2]+dot7[ih2]+dot9[ih2]))*M2+
			    dot12*((dot14-2.*dot19)*dot6+(2.*dot1*dot13+dot22-2.*dot13*(dot2+dot3+dot5))*M-
				   2.*a1*dot13*pow<3,1>(M)+(dot14-a1*dot14+(-3.+4.*a1)*dot19)*M2))
	    /(a2*dot1*sqr(dot6+a2*M2));
	  diag[8]=M2*(-2.*dot17*dot2-4.*a1*dot10[ih1]*dot18*dot2+2.*dot2*dot25-2.*dot10[ih1]*dot26+2.*a1*dot10[ih1]*dot26+4.*dot10[ih1]*dot18*dot3-
		      4.*a1*dot10[ih1]*dot18*dot3+2.*dot25*dot3+4.*a1*dot18*dot2*dot4[ih1]-2.*a1*dot26*dot4[ih1]-4.*dot18*dot3*dot4[ih1]+4.*a1*dot18*dot3*dot4[ih1]+
		      4.*dot10[ih1]*dot18*dot5-4.*a1*dot10[ih1]*dot18*dot5+4.*a1*dot18*dot4[ih1]*dot5-2.*dot17*dot6-4.*a1*dot10[ih1]*dot18*dot6+2.*dot25*dot6+
		      4.*a1*dot18*dot4[ih1]*dot6+4.*dot10[ih1]*dot14*dot7[ih2]-8.*a1*dot10[ih1]*dot14*dot7[ih2]+4.*a12*dot10[ih1]*dot14*dot7[ih2]-4.*dot10[ih1]*dot19*dot7[ih2]+
		      4.*a1*dot10[ih1]*dot19*dot7[ih2]-2.*dot27*dot7[ih2]+2.*a1*dot27*dot7[ih2]-4.*dot14*dot4[ih1]*dot7[ih2]+8.*a1*dot14*dot4[ih1]*dot7[ih2]-4.*a12*dot14*dot4[ih1]*dot7[ih2]+
		      4.*dot19*dot4[ih1]*dot7[ih2]-4.*a1*dot19*dot4[ih1]*dot7[ih2]+4.*dot15*dot5*dot7[ih2]-4.*a1*dot15*dot5*dot7[ih2]-4.*dot18*dot2*dot8[ih1]-
		      4.*a1*dot18*dot2*dot8[ih1]+2.*a1*dot26*dot8[ih1]-4.*a1*dot18*dot3*dot8[ih1]-4.*a1*dot18*dot5*dot8[ih1]-4.*dot18*dot6*dot8[ih1]-
		      4.*a1*dot18*dot6*dot8[ih1]-4.*a1*dot14*dot7[ih2]*dot8[ih1]+4.*a12*dot14*dot7[ih2]*dot8[ih1]-4.*dot19*dot7[ih2]*dot8[ih1]+4.*a1*dot19*dot7[ih2]*dot8[ih1]-
		      4.*dot10[ih1]*dot19*dot9[ih2]+4.*dot15*dot2*dot9[ih2]-2.*dot27*dot9[ih2]+4.*dot15*dot3*dot9[ih2]-4.*dot14*dot4[ih1]*dot9[ih2]+
		      4.*dot19*dot4[ih1]*dot9[ih2]+4.*dot15*dot5*dot9[ih2]+4.*dot15*dot6*dot9[ih2]-4.*dot19*dot8[ih1]*dot9[ih2]+
		      (-2.*dot12*dot22+2.*dot11[ih2]*dot24+3.*dot28-4.*dot11[ih2]*dot13*dot4[ih1]+4.*dot20*dot4[ih1]+4.*dot12*dot13*(dot3+dot5)-
		       2.*dot16*(dot2+3.*(dot3+dot5)+dot6)+
		       2.*(dot10[ih1]*(dot23-2.*(-1.+2.*a1)*a2*dot13*dot7[ih2])+dot23*(-2.*dot4[ih1]+3.*dot8[ih1])+
			   dot7[ih2]*(dot21-a1*dot21+dot24-4.*dot13*dot4[ih1]+2.*dot13*(-(a1*(-3.+2.*a1)*(dot4[ih1]-dot8[ih1]))+dot8[ih1]))-dot24*dot9[ih2]))*M-
		      2.*(dot16+2.*a1*(-(dot12*dot13)+dot16))*pow<3,1>(M)+
		      2.*dot1*(dot17-dot25+2.*dot18*(dot8[ih1]+a1*(dot10[ih1]-dot4[ih1]+dot8[ih1]))-2.*dot15*dot9[ih2]+2.*dot16*M)+
		      (2.*dot12*dot19-dot25+4.*dot18*dot4[ih1]+2.*dot15*dot7[ih2]-
		       2.*(dot17-a1*dot17+a1*(2.*dot12*dot19-2.*dot25+2.*dot18*dot4[ih1]+dot15*dot7[ih2])+2.*dot18*dot8[ih1]-2.*dot15*dot9[ih2]))*M2)
	  /(16.*a1*a2*(dot1-dot2-dot3)*(-dot1+dot2+dot3+dot5+dot6+M2)*(-dot2+a1*(-dot1+dot2+dot3+dot5+a1*M2)));
	  diag[9]=-0.125*M2*(-4.*a1*dot10[ih1]*dot18*dot2+2.*dot2*dot25-2.*dot10[ih1]*dot26+2.*a1*dot10[ih1]*dot26+4.*dot10[ih1]*dot18*dot3-
			     4.*a1*dot10[ih1]*dot18*dot3+2.*dot25*dot3+4.*a1*dot18*dot2*dot4[ih1]-2.*a1*dot26*dot4[ih1]-4.*dot18*dot3*dot4[ih1]+
			     4.*a1*dot18*dot3*dot4[ih1]+4.*dot10[ih1]*dot18*dot5-4.*a1*dot10[ih1]*dot18*dot5+4.*a1*dot18*dot4[ih1]*dot5-
			     4.*a1*dot10[ih1]*dot18*dot6+2.*dot25*dot6+4.*a1*dot18*dot4[ih1]*dot6-4.*dot18*dot2*dot8[ih1]-4.*a1*dot18*dot2*dot8[ih1]+
			     2.*a1*dot26*dot8[ih1]-4.*a1*dot18*dot3*dot8[ih1]-4.*a1*dot18*dot5*dot8[ih1]-4.*dot18*dot6*dot8[ih1]-4.*a1*dot18*dot6*dot8[ih1]-
			     4.*dot10[ih1]*dot14*dot9[ih2]+4.*a1*dot10[ih1]*dot14*dot9[ih2]+4.*dot15*dot2*dot9[ih2]+4.*dot15*dot3*dot9[ih2]-4.*a1*dot14*dot4[ih1]*dot9[ih2]+
			     4.*dot15*dot6*dot9[ih2]+4.*a1*dot14*dot8[ih1]*dot9[ih2]+
			     (dot28+2.*dot20*dot4[ih1]-2.*dot16*(dot2+dot3+dot5+dot6)-2.*(-1.+2.*a1)*dot23*(dot10[ih1]-dot4[ih1]+dot8[ih1])-
			      2.*dot21*dot9[ih2]+4.*(-1.+2.*a1)*dot13*(dot10[ih1]-dot4[ih1]+dot8[ih1])*dot9[ih2])*M-2.*dot16*pow<3,1>(M)+
			     2.*dot1*(dot17-dot25+2.*dot18*(dot8[ih1]+a1*(dot10[ih1]-dot4[ih1]+dot8[ih1]))-2.*dot15*dot9[ih2]+dot16*M)+
			     (dot25+2.*dot18*dot4[ih1]-4.*dot18*dot8[ih1]+2.*dot15*dot9[ih2])*M2-2.*dot17*(dot2+dot6+M2-a1*M2))
	    /(a1*a2*(dot1-dot2-dot3)*(-dot6-a2*M2)*(dot1-dot2-dot6-a2*M2));
	  diag[10]=-0.125*M2*(2.*(dot17*dot3-a2*(2.*dot11[ih2]*dot14-dot26)*dot4[ih1]-2.*dot18*dot3*dot4[ih1]+2.*dot11[ih2]*dot15*dot5+dot17*dot5+
				  2.*dot12*dot19*dot5-dot25*dot5-2.*dot18*dot4[ih1]*dot5-2.*dot12*dot14*(dot3+dot5)+
				  2.*dot15*dot3*dot7[ih2]-2.*dot19*dot4[ih1]*dot7[ih2]+2.*a1*dot19*dot4[ih1]*dot7[ih2]+2.*dot15*dot5*dot7[ih2]+
				  2.*dot18*(dot3+dot5)*dot8[ih1]+a1*(-2.*dot11[ih2]*dot14+dot26-2.*dot19*dot7[ih2])*dot8[ih1]+
				  dot1*(2.*dot12*dot14-dot17+2.*dot18*dot4[ih1]-2.*dot15*dot7[ih2]-2.*dot18*dot8[ih1]))+
			      (-2.*dot12*dot22+dot28+2.*dot24*dot7[ih2]+2.*dot11[ih2]*(dot21+2.*(-1.+2.*a1)*dot13*(dot4[ih1]-dot8[ih1]))+
			       2.*(-1.+2.*a1)*dot23*(dot4[ih1]-dot8[ih1]))*M+
			      (2.*dot11[ih2]*dot15+2.*dot12*dot19-dot25+2.*a1*(-2.*dot12*dot14+dot17-2.*dot18*dot4[ih1]+2.*dot15*dot7[ih2]+2.*dot18*dot8[ih1]))*M2+
			      2.*dot10[ih1]*(-(a2*(dot26-2.*dot19*dot7[ih2]))+(1.-2.*a1)*dot23*M+
					2.*dot11[ih2]*(dot14-a1*dot14+(1.-2.*a1)*dot13*M)+2.*dot18*(-dot1+dot3+dot5+a1*M2)))
	    /(a1*(dot1-dot2-dot3)*(-dot1+dot3+dot5+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot5+a1*M2)));
	  diag[11]=-0.0625*M2*(4.*dot10[ih1]*dot11[ih2]*dot14-4.*a1*dot10[ih1]*dot11[ih2]*dot14+4.*dot1*dot12*dot14-4.*a1*dot1*dot12*dot14-2.*dot1*dot17+2.*a1*dot1*dot17-
			       4.*dot1*dot10[ih1]*dot18+4.*a1*dot1*dot10[ih1]*dot18+4.*a1*dot12*dot14*dot2-2.*a1*dot17*dot2-4.*a1*dot10[ih1]*dot18*dot2-
			       2.*dot10[ih1]*dot26+2.*a1*dot10[ih1]*dot26-4.*dot12*dot14*dot3+4.*a1*dot12*dot14*dot3+2.*dot17*dot3-2.*a1*dot17*dot3+
			       4.*dot10[ih1]*dot18*dot3-4.*a1*dot10[ih1]*dot18*dot3-4.*dot11[ih2]*dot14*dot4[ih1]+4.*a1*dot11[ih2]*dot14*dot4[ih1]+4.*dot1*dot18*dot4[ih1]-
			       4.*a1*dot1*dot18*dot4[ih1]+4.*a1*dot18*dot2*dot4[ih1]+2.*dot26*dot4[ih1]-2.*a1*dot26*dot4[ih1]-4.*dot18*dot3*dot4[ih1]+4.*a1*dot18*dot3*dot4[ih1]-
			       4.*dot12*dot14*dot5+4.*dot11[ih2]*dot15*dot5+2.*dot17*dot5+4.*dot10[ih1]*dot18*dot5+4.*dot12*dot19*dot5-2.*dot25*dot5-
			       4.*dot18*dot4[ih1]*dot5+4.*dot10[ih1]*dot14*dot7[ih2]-8.*a1*dot10[ih1]*dot14*dot7[ih2]+4.*a12*dot10[ih1]*dot14*dot7[ih2]-4.*dot1*dot15*dot7[ih2]+
			       4.*a1*dot1*dot15*dot7[ih2]-4.*a1*dot15*dot2*dot7[ih2]-2.*dot27*dot7[ih2]+2.*a1*dot27*dot7[ih2]+4.*dot15*dot3*dot7[ih2]-4.*a1*dot15*dot3*dot7[ih2]-
			       4.*dot14*dot4[ih1]*dot7[ih2]+8.*a1*dot14*dot4[ih1]*dot7[ih2]-4.*a12*dot14*dot4[ih1]*dot7[ih2]+8.*dot15*dot5*dot7[ih2]-4.*a1*dot15*dot5*dot7[ih2]-
			       4.*a1*dot11[ih2]*dot14*dot8[ih1]-4.*dot1*dot18*dot8[ih1]+4.*a1*dot1*dot18*dot8[ih1]-4.*a1*dot18*dot2*dot8[ih1]+2.*a1*dot26*dot8[ih1]+
			       4.*dot18*dot3*dot8[ih1]-4.*a1*dot18*dot3*dot8[ih1]+4.*dot18*dot5*dot8[ih1]-4.*a1*dot14*dot7[ih2]*dot8[ih1]+4.*a12*dot14*dot7[ih2]*dot8[ih1]-
			       4.*dot19*dot7[ih2]*dot8[ih1]+4.*dot10[ih1]*dot14*dot9[ih2]-4.*a1*dot10[ih1]*dot14*dot9[ih2]-4.*dot10[ih1]*dot19*dot9[ih2]-2.*dot27*dot9[ih2]-
			       4.*dot14*dot4[ih1]*dot9[ih2]+4.*a1*dot14*dot4[ih1]*dot9[ih2]+4.*dot19*dot4[ih1]*dot9[ih2]+
			       4.*dot15*dot5*dot9[ih2]-4.*a1*dot14*dot8[ih1]*dot9[ih2]-4.*dot19*dot8[ih1]*dot9[ih2]+
			       (3.*dot28+2.*(dot1*dot16+dot11[ih2]*dot21-2.*dot12*dot22+dot11[ih2]*dot24+2.*dot12*dot13*dot3-4.*dot11[ih2]*dot13*dot4[ih1]+
					     dot20*dot4[ih1]-2.*dot23*dot4[ih1]+2.*dot12*dot13*dot5-2.*dot16*(dot3+dot5)+dot21*dot7[ih2]-a1*dot21*dot7[ih2]+
					     2.*dot24*dot7[ih2]-4.*dot13*dot4[ih1]*dot7[ih2]+2.*dot11[ih2]*dot13*dot8[ih1]+3.*dot23*dot8[ih1]+2.*dot13*dot7[ih2]*dot8[ih1]+
					     4.*a12*dot13*dot7[ih2]*(-dot4[ih1]+dot8[ih1])+(dot21-dot24-2.*dot13*dot4[ih1]+2.*dot13*dot8[ih1])*dot9[ih2]+
					     dot10[ih1]*(dot23-2.*(-1.+2.*a1)*dot13*(dot11[ih2]+dot7[ih2]-a1*dot7[ih2]+dot9[ih2]))+
					     2.*a1*dot13*(dot12*(-dot1+dot2+dot3)+(dot4[ih1]-dot8[ih1])*(2.*dot11[ih2]+3.*dot7[ih2]+2.*dot9[ih2]))))*M+
			       4.*a1*(dot12*dot13-dot16)*pow<3,1>(M)+
			       (4.*dot12*dot19-3.*dot25+2.*dot18*dot4[ih1]+
				2.*a1*(dot17-2.*dot12*(dot14+dot19)+2.*dot25+dot15*dot7[ih2]+2.*dot18*(dot10[ih1]-2.*dot4[ih1]+dot8[ih1]))+2.*dot15*(dot11[ih2]+dot7[ih2]+dot9[ih2]))*M2)
	    /(a1*a2*(dot1-dot2-dot3)*(-dot6-a2*M2)*(-dot1+a1*(dot5+dot6+M2)));
	  diag[12]=-0.0625*M2*(-4.*dot1*dot12*dot14+2.*a1*dot1*dot12*dot14+4.*dot1*dot11[ih2]*dot15-a1*dot1*dot17+4.*dot12*dot14*dot2-2.*a1*dot12*dot14*dot2-
			       4.*dot11[ih2]*dot15*dot2-2.*dot17*dot2+a1*dot17*dot2+2.*dot1*dot25-2.*dot2*dot25+4.*dot11[ih2]*dot27-2.*a1*dot12*dot14*dot3-
			       4.*dot11[ih2]*dot15*dot3+2.*dot17*dot3+a1*dot17*dot3-2.*dot25*dot3+8.*dot11[ih2]*dot14*dot4[ih1]-4.*a1*dot11[ih2]*dot14*dot4[ih1]-
			       4.*dot1*dot18*dot4[ih1]+4.*a1*dot1*dot18*dot4[ih1]-4.*dot11[ih2]*dot19*dot4[ih1]+4.*dot18*dot2*dot4[ih1]-4.*a1*dot18*dot2*dot4[ih1]-
			       2.*dot26*dot4[ih1]+4.*a1*dot26*dot4[ih1]+4.*dot18*dot3*dot4[ih1]-4.*a1*dot18*dot3*dot4[ih1]-8.*dot11[ih2]*dot15*dot5+2.*dot17*dot5-
			       2.*dot25*dot5+4.*dot18*dot4[ih1]*dot5-4.*a1*dot18*dot4[ih1]*dot5+4.*dot12*dot14*dot6-4.*dot11[ih2]*dot15*dot6-2.*dot17*dot6-
			       4.*dot12*dot19*dot6+2.*dot25*dot6-4.*a1*dot18*dot4[ih1]*dot6+4.*dot1*dot15*dot7[ih2]-2.*a1*dot1*dot15*dot7[ih2]-4.*dot15*dot2*dot7[ih2]+
			       2.*a1*dot15*dot2*dot7[ih2]+2.*dot27*dot7[ih2]-2.*a1*dot27*dot7[ih2]+2.*a1*dot15*dot3*dot7[ih2]+4.*dot14*dot4[ih1]*dot7[ih2]-
			       8.*a1*dot14*dot4[ih1]*dot7[ih2]+4.*a12*dot14*dot4[ih1]*dot7[ih2]-4.*dot15*dot5*dot7[ih2]+4.*a1*dot15*dot5*dot7[ih2]-
			       4.*dot15*dot6*dot7[ih2]+4.*a1*dot11[ih2]*dot14*dot8[ih1]+4.*dot1*dot18*dot8[ih1]-4.*a1*dot1*dot18*dot8[ih1]+4.*dot11[ih2]*dot19*dot8[ih1]-
			       4.*dot18*dot2*dot8[ih1]+4.*a1*dot18*dot2*dot8[ih1]+4.*dot26*dot8[ih1]-4.*a1*dot26*dot8[ih1]-4.*dot18*dot3*dot8[ih1]+
			       4.*a1*dot18*dot3*dot8[ih1]-4.*dot18*dot5*dot8[ih1]+4.*a1*dot18*dot5*dot8[ih1]-4.*dot18*dot6*dot8[ih1]+4.*a1*dot18*dot6*dot8[ih1]+
			       4.*a1*dot14*dot7[ih2]*dot8[ih1]-4.*a12*dot14*dot7[ih2]*dot8[ih1]-2.*dot27*dot9[ih2]+
			       (2.*(-2.+a1)*dot1*dot12*dot13-a1*dot1*dot16-2.*dot16*dot2+a1*dot16*dot2+2.*dot16*dot3+a1*dot16*dot3-
				8.*a1*dot11[ih2]*dot13*dot4[ih1]+2.*dot23*dot4[ih1]-4.*a1*dot23*dot4[ih1]+2.*dot16*dot5-
				2.*dot12*dot13*((-2.+a1)*dot2+a1*dot3-2.*dot6)-2.*dot16*dot6+2.*a1*dot21*dot7[ih2]+4.*dot13*dot4[ih1]*dot7[ih2]-
				12.*a1*dot13*dot4[ih1]*dot7[ih2]+8.*a12*dot13*dot4[ih1]*dot7[ih2]+
				4.*((-1.+2.*a1)*dot11[ih2]*dot13-a2*(dot23+dot13*dot7[ih2]-2.*a1*dot13*dot7[ih2]))*dot8[ih1]-
				2.*(dot11[ih2]*(dot21+dot24-4.*dot13*dot4[ih1])+dot21*dot7[ih2]-dot24*dot9[ih2]))*M-
			       2.*(-2.*a2*dot12*dot13+dot16-2.*a1*dot16)*pow<3,1>(M)-
			       2.*(dot17-2.*a1*dot17-2.*a2*dot12*(dot14-dot19)-dot25+2.*a1*dot25+3.*dot15*(dot11[ih2]+dot7[ih2]-a1*dot7[ih2])+
				   2.*dot18*dot8[ih1]-2.*a1*dot18*dot8[ih1])*M2+
			       4.*dot10[ih1]*(-(dot7[ih2]*(dot14-dot19+dot13*M))-a12*dot7[ih2]*(dot14+2.*dot13*M)+
					      dot11[ih2]*(-(a2*dot14)+dot19+(-1.+2.*a1)*dot13*M)+
					      a1*(-dot26+2.*dot14*dot7[ih2]+dot23*M+3.*dot13*dot7[ih2]*M+dot18*(-dot1+dot2+dot3+dot5+dot6+M2))))
	    /(a1*a2*(dot1-dot2-dot3)*(-dot6-a2*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot5+a1*M2)));
	  diag[13]=M2*(2.*dot11[ih2]*dot27+2.*dot17*(dot3+dot5)-2.*dot25*(dot2+dot3+dot5)+
		       2.*(-2.*a12*dot14*dot4[ih1]-2.*dot19*dot4[ih1]+2.*dot15*(dot3+dot5)+
			   a1*(dot27+2.*(dot14+dot19)*dot4[ih1]-2.*dot15*(dot2+dot3+dot5)))*dot7[ih2]-
		       2.*a2*(-(dot26*dot4[ih1])+dot10[ih1]*(dot26+2.*a1*dot14*dot7[ih2]-2.*dot19*dot7[ih2]))+
		       2.*a1*(dot26+2.*a1*dot14*dot7[ih2]-2.*dot19*dot7[ih2])*dot8[ih1]+2.*dot11[ih2]*dot24*M+3.*dot28*M-
		       2.*(2.*dot11[ih2]*dot13*dot4[ih1]+2.*dot16*(dot3+dot5)-4.*a12*dot10[ih1]*dot13*dot7[ih2]+
			   a1*(4.*dot10[ih1]*dot13+dot21)*dot7[ih2]+dot4[ih1]*(-dot20+dot23+2.*(1.-2.*a1*a2)*dot13*dot7[ih2])-
			   2.*(dot23-2.*a1*a2*dot13*dot7[ih2])*dot8[ih1]+dot24*(-2.*dot7[ih2]+dot9[ih2]))*M-
		       4.*a1*dot16*pow<3,1>(M)+2.*dot1*(-dot17+2.*dot12*(dot14-dot19)+dot25-2.*a2*dot15*dot7[ih2]+dot16*M)+
		       (-3.*dot25+2.*dot18*dot4[ih1]+2.*a1*(dot17+2.*dot25-2.*dot18*dot4[ih1]+dot15*dot7[ih2]))*M2+
		       4.*dot12*(-(dot22*M)+dot19*(dot2+dot3+dot5+M2-a1*M2)-(dot14-dot13*M)*(dot3+dot5+a1*M2)))
	    /(16.*a1*a2*(dot1-dot2-dot3)*(-dot1+dot2+dot3+dot5+dot6+M2)*(-(a2*dot1)+dot2+dot6-a1*(dot2+dot3+dot6)+M2+(-2.+a1)*a1*M2));
	  diag[14]=-0.125*M2*(2.*(dot11[ih2]*dot27-2.*dot12*dot14*(dot3+dot5)+dot17*(dot3+dot5)+2.*dot12*dot19*(dot2+dot3+dot5)-
				  dot25*(dot2+dot3+dot5)+2.*dot15*(dot3+dot5)*dot7[ih2]+a1*(dot27-2.*dot15*(dot2+dot3+dot5))*dot7[ih2]+
				  dot1*(-dot17+2.*dot12*(dot14-dot19)+dot25-2.*a2*dot15*dot7[ih2]))+
			      (-2.*dot12*dot22+dot28+2.*(a1*dot21+dot24-2.*a1*dot24)*dot7[ih2])*M-
			      (4.*a1*dot12*dot14-2.*dot12*dot19+dot25+8.*a12*dot15*dot7[ih2]-2.*a1*(dot17+3.*dot15*dot7[ih2]))*M2)
	    /(a1*a2*dot2*(-dot1+dot2+dot3+dot5+dot6+M2)*(-dot6-a2*M2));
	  diag[15]=M2*(2.*(-(dot25*dot5)+dot17*(dot3+dot5)-a2*dot27*dot7[ih2]+
			   2.*(-(dot18*dot3*dot4[ih1])+dot15*(dot3+dot5-a1*dot5)*dot7[ih2]+dot18*(dot3+dot5)*dot8[ih1]+
			       (dot14*dot4[ih1]-dot19*dot4[ih1]+dot19*dot8[ih1])*(-(a2*dot7[ih2])-dot9[ih2]))-
			   dot27*dot9[ih2]+2.*dot10[ih1]*(dot18*(dot3+dot5)-a2*dot19*dot7[ih2]-dot19*dot9[ih2]))+
		       (dot28+2.*dot10[ih1]*(dot23+2.*a1*dot13*dot7[ih2])+
			2.*(dot20*dot4[ih1]+dot23*(-dot4[ih1]+dot8[ih1])+
			    dot7[ih2]*(dot24-2.*dot13*dot4[ih1]+a1*(dot21-2.*dot24+2.*dot13*(dot4[ih1]+dot8[ih1])))))*M+
		       (-dot25+2.*dot18*dot4[ih1]-8.*a12*dot15*dot7[ih2]+
			2.*a1*(dot17+3.*dot15*dot7[ih2]+2.*dot18*(dot10[ih1]-dot4[ih1]+dot8[ih1])))*M2)
	    /(8.*a1*dot3*(-dot1+dot2+dot3+dot5+dot6+M2)*(-(a2*dot1)+dot2+dot6-a1*(dot2+dot3+dot6)+M2+(-2.+a1)*a1*M2));
	  diag[16]=M2*(2.*(-(dot26*dot4[ih1])+dot25*dot5-dot17*(dot3+dot5)+dot27*(dot7[ih2]+dot9[ih2])+
			   dot10[ih1]*(dot26-2.*dot18*(dot3+dot5)+2.*dot19*dot9[ih2])+
			   2.*(dot18*dot3*dot4[ih1]+dot14*dot4[ih1]*dot7[ih2]-dot15*(dot3+dot5)*dot7[ih2]-
			       dot18*(dot3+dot5)*dot8[ih1]+dot19*dot7[ih2]*dot8[ih1]+
			       dot14*dot4[ih1]*dot9[ih2]-dot19*dot4[ih1]*dot9[ih2]+dot19*dot8[ih1]*dot9[ih2]))-
		       (2.*dot1*dot16-2.*dot12*dot22+2.*dot10[ih1]*dot23+3.*dot28+4.*dot12*dot13*(dot3+dot5)-4.*dot16*(dot3+dot5)-
			4.*dot4[ih1]*(dot11[ih2]*dot13-dot20+dot23+2.*dot13*dot7[ih2])+6.*dot23*dot8[ih1]+
			2.*dot24*(dot11[ih2]+2.*dot7[ih2]-dot9[ih2]))*M-
		       4.*a12*dot7[ih2]*(dot10[ih1]-dot4[ih1]+dot8[ih1])*(dot14+2.*dot13*M)+
		       (-2.*dot12*dot19+3.*dot25-4.*dot18*dot4[ih1])*M2+
		       2.*a1*(dot17*dot2+2.*dot10[ih1]*dot18*dot2-dot10[ih1]*dot26+dot17*dot3+2.*dot10[ih1]*dot18*dot3-
			      2.*dot18*dot2*dot4[ih1]+dot26*dot4[ih1]-2.*dot18*dot3*dot4[ih1]+dot17*dot5+2.*dot10[ih1]*dot18*dot5-
			      2.*dot18*dot4[ih1]*dot5+dot17*dot6+2.*dot10[ih1]*dot18*dot6-2.*dot18*dot4[ih1]*dot6+
			      2.*dot10[ih1]*dot14*dot7[ih2]+2.*dot15*dot2*dot7[ih2]-dot27*dot7[ih2]+2.*dot15*dot3*dot7[ih2]-
			      4.*dot14*dot4[ih1]*dot7[ih2]+4.*dot15*dot5*dot7[ih2]+2.*dot15*dot6*dot7[ih2]+
			      2.*dot18*dot2*dot8[ih1]-dot26*dot8[ih1]+2.*dot18*dot3*dot8[ih1]+2.*dot18*dot5*dot8[ih1]+
			      2.*dot18*dot6*dot8[ih1]+2.*dot12*dot13*(dot2+dot3+dot5+dot6)*M+
			      dot7[ih2]*(dot21+2.*dot13*(dot10[ih1]-3.*dot4[ih1]+dot8[ih1]))*M+
			      2.*dot16*pow<3,1>(M)-dot1*(dot17+2.*(dot15*dot7[ih2]+dot18*(dot10[ih1]-dot4[ih1]+dot8[ih1])+dot12*dot13*M))+
			      (2.*dot12*dot19-2.*dot25+2.*dot18*dot4[ih1]+dot15*dot7[ih2])*M2))
	    /(16.*a1*a2*(-dot1+dot2+dot3+dot5+dot6+M2)*(-dot6-a2*M2)*(-dot1+a1*(dot5+dot6+M2)));
	  diag[17]=-0.0625*M2*(4.*dot12*dot14*dot2-4.*dot11[ih2]*dot15*dot2-4.*dot17*dot2+a1*dot17*dot2-4.*dot12*dot19*dot2+
			       2.*dot2*dot25+2.*dot11[ih2]*dot27-4.*dot11[ih2]*dot15*dot3+a1*dot17*dot3-4.*dot12*dot19*dot3+
			       2.*dot25*dot3+4.*dot11[ih2]*dot14*dot4[ih1]-4.*dot11[ih2]*dot19*dot4[ih1]+4.*dot18*dot2*dot4[ih1]-
			       2.*a1*dot18*dot2*dot4[ih1]-2.*dot26*dot4[ih1]+2.*a1*dot26*dot4[ih1]-2.*a1*dot18*dot3*dot4[ih1]-
			       4.*dot11[ih2]*dot15*dot5+a1*dot17*dot5-2.*a1*dot18*dot4[ih1]*dot5+4.*dot12*dot14*dot6-
			       4.*dot11[ih2]*dot15*dot6-4.*dot17*dot6+a1*dot17*dot6-4.*dot12*dot19*dot6+4.*dot25*dot6-
			       2.*a1*dot18*dot4[ih1]*dot6-4.*dot15*dot2*dot7[ih2]+4.*a1*dot15*dot2*dot7[ih2]+2.*dot27*dot7[ih2]-
			       4.*a1*dot27*dot7[ih2]-4.*dot15*dot3*dot7[ih2]+4.*a1*dot15*dot3*dot7[ih2]+4.*dot14*dot4[ih1]*dot7[ih2]-
			       8.*a1*dot14*dot4[ih1]*dot7[ih2]+4.*a12*dot14*dot4[ih1]*dot7[ih2]-4.*dot15*dot5*dot7[ih2]+
			       8.*a1*dot15*dot5*dot7[ih2]-4.*dot15*dot6*dot7[ih2]+4.*dot11[ih2]*dot19*dot8[ih1]-8.*dot18*dot2*dot8[ih1]+
			       2.*a1*dot18*dot2*dot8[ih1]+4.*dot26*dot8[ih1]-2.*a1*dot26*dot8[ih1]-4.*dot18*dot3*dot8[ih1]+
			       2.*a1*dot18*dot3*dot8[ih1]-4.*dot18*dot5*dot8[ih1]+2.*a1*dot18*dot5*dot8[ih1]-8.*dot18*dot6*dot8[ih1]+
			       2.*a1*dot18*dot6*dot8[ih1]-4.*a12*dot14*dot7[ih2]*dot8[ih1]+4.*dot15*dot2*dot9[ih2]-2.*dot27*dot9[ih2]+
			       4.*dot15*dot3*dot9[ih2]+4.*dot15*dot5*dot9[ih2]+4.*dot15*dot6*dot9[ih2]+
			       ((-4.+a1)*dot16*dot2-2.*dot11[ih2]*(dot24-2.*dot13*dot4[ih1])-4.*dot16*dot6+
				4.*dot12*dot13*(dot2+dot6)+a1*dot16*(dot3+dot5+dot6)-2.*dot23*dot8[ih1]+
				2.*(-1.+2.*a1)*dot7[ih2]*(dot24-2.*a2*dot13*dot4[ih1]-2.*a1*dot13*dot8[ih1])+2.*dot24*dot9[ih2])*M+
			       (4.*a2*dot12*dot13+(-4.+5.*a1)*dot16)*pow<3,1>(M)-
			       dot1*(-4.*dot11[ih2]*dot15+(-2.+a1)*dot17+2.*a1*dot18*dot8[ih1]+
				     2.*((2.+a1)*dot10[ih1]*dot18+dot25-dot18*(a1*dot4[ih1]+2.*dot8[ih1]))+
				     4.*dot15*(-(a2*dot7[ih2])+dot9[ih2])+(-2.+a1)*dot16*M+4.*dot12*(dot14-dot19+dot13*M))+
			       (4.*a2*dot12*(dot14-dot19)+8.*a12*dot15*dot7[ih2]-4.*(dot17-dot25+dot15*(dot11[ih2]+dot7[ih2]))-
				8.*dot18*dot8[ih1]+a1*(5.*dot17-4.*dot25-2.*dot18*dot4[ih1]+6.*dot18*dot8[ih1])+4.*dot15*dot9[ih2])*M2+
			       2.*dot10[ih1]*(2.*dot11[ih2]*dot19-dot26+2.*dot18*(dot3+dot5)+2.*dot19*dot7[ih2]+dot23*M-
					      2.*a12*dot7[ih2]*(dot14+2.*dot13*M)+
					      a1*(-dot26+2.*dot7[ih2]*(dot14+dot13*M)+dot18*(dot2+dot3+dot5+dot6+3.*M2))))
	    /(a1*a2*(-dot1+dot2+dot3+dot5+dot6+M2)*(-(a2*dot1)+dot2+dot6-a1*(dot2+dot3+dot6)+M2+(-2.+a1)*a1*M2)*(-dot6-a2*M2));
	  diag[18]=-0.125*M2*(2.*dot17*dot2-4.*dot10[ih1]*dot18*dot2+4.*a1*dot10[ih1]*dot18*dot2+2.*dot10[ih1]*dot26-
			      2.*a1*dot10[ih1]*dot26-4.*dot10[ih1]*dot18*dot3+4.*a1*dot10[ih1]*dot18*dot3-4.*a1*dot18*dot2*dot4[ih1]+
			      2.*a1*dot26*dot4[ih1]-4.*a1*dot18*dot3*dot4[ih1]-4.*dot10[ih1]*dot18*dot5+4.*a1*dot10[ih1]*dot18*dot5-
			      4.*a1*dot18*dot4[ih1]*dot5+2.*dot17*dot6-4.*dot10[ih1]*dot18*dot6+4.*a1*dot10[ih1]*dot18*dot6-
			      2.*dot25*dot6+4.*dot18*dot4[ih1]*dot6-4.*a1*dot18*dot4[ih1]*dot6+4.*a1*dot18*dot2*dot8[ih1]-
			      2.*a1*dot26*dot8[ih1]+4.*a1*dot18*dot3*dot8[ih1]+4.*a1*dot18*dot5*dot8[ih1]+4.*a1*dot18*dot6*dot8[ih1]+
			      4.*dot10[ih1]*dot14*dot9[ih2]-4.*a1*dot10[ih1]*dot14*dot9[ih2]+4.*a1*dot14*dot4[ih1]*dot9[ih2]-
			      4.*dot15*dot6*dot9[ih2]-4.*a1*dot14*dot8[ih1]*dot9[ih2]+
			      (-dot28-2.*dot20*dot4[ih1]+2.*dot16*(dot2+dot3+dot5+dot6)-4.*a2*dot23*(dot10[ih1]-dot4[ih1]+dot8[ih1])+
			       2.*(dot21+4.*a2*dot13*(dot10[ih1]-dot4[ih1]+dot8[ih1]))*dot9[ih2])*M+2.*dot16*pow<3,1>(M)-
			      2.*dot1*(dot17+2.*dot18*(-(a2*dot10[ih1])+a1*(-dot4[ih1]+dot8[ih1]))+dot16*M)-
			      (4.*dot10[ih1]*dot18+dot25-2.*(dot17-a1*dot17+dot18*(dot4[ih1]+2.*a1*(dot10[ih1]-dot4[ih1]+dot8[ih1])))+
			       2.*dot15*dot9[ih2])*M2)
	    /(a2*(dot1-dot2-dot3)*(-(a2*dot1)+dot2+dot6-a1*(dot2+dot3+dot6)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot6-a2*M2));
	  diag[19]=-0.125*M2*(4.*dot1*dot12*dot14-4.*dot1*dot11[ih2]*dot15-2.*dot1*dot17-4.*dot1*dot12*dot19+
			      4.*dot11[ih2]*dot15*dot2+4.*dot12*dot19*dot2+2.*dot1*dot25-2.*dot2*dot25-
			      4.*dot12*dot14*dot3+4.*dot11[ih2]*dot15*dot3+2.*dot17*dot3+4.*dot12*dot19*dot3-
			      2.*dot25*dot3-4.*dot11[ih2]*dot14*dot4[ih1]+4.*a1*dot11[ih2]*dot14*dot4[ih1]+
			      2.*dot26*dot4[ih1]-2.*a1*dot26*dot4[ih1]-4.*dot12*dot14*dot5+4.*dot11[ih2]*dot15*dot5+
			      2.*dot17*dot5+4.*dot12*dot19*dot5-2.*dot25*dot5-4.*dot1*dot15*dot7[ih2]+
			      4.*dot15*dot3*dot7[ih2]-4.*dot19*dot4[ih1]*dot7[ih2]+4.*a1*dot19*dot4[ih1]*dot7[ih2]+
			      4.*dot15*dot5*dot7[ih2]-4.*a1*dot11[ih2]*dot14*dot8[ih1]+2.*a1*dot26*dot8[ih1]-
			      4.*a1*dot19*dot7[ih2]*dot8[ih1]+
			      (-2.*dot12*dot22+dot28+2.*dot24*dot7[ih2]+
			       2.*dot11[ih2]*(dot21-4.*a2*dot13*(dot4[ih1]-dot8[ih1]))-4.*a2*dot23*(dot4[ih1]-dot8[ih1]))*M-
			      2.*a2*dot10[ih1]*(dot26-2.*dot19*dot7[ih2]-2.*dot23*M-2.*dot11[ih2]*(dot14+2.*dot13*M))+
			      (2.*dot11[ih2]*dot15+2.*dot12*dot19-dot25+2.*a1*(-2.*dot12*dot14+dot17+2.*dot15*dot7[ih2]))*M2)
	    /(a1*a2*(dot1-dot2-dot3)*(dot5+a1*M2)*(-dot1+dot3+dot5+a1*M2));
	  diag[20]=-0.0625*M2*(4.*dot1*dot11[ih2]*dot15+2.*a1*dot1*dot17-4.*dot11[ih2]*dot15*dot2-2.*a1*dot17*dot2+
			       2.*dot11[ih2]*dot27-4.*dot11[ih2]*dot15*dot3+2.*dot17*dot3-2.*a1*dot17*dot3+
			       4.*dot11[ih2]*dot14*dot4[ih1]-4.*a1*dot11[ih2]*dot14*dot4[ih1]-4.*a1*dot1*dot18*dot4[ih1]+
			       4.*a1*dot18*dot2*dot4[ih1]-2.*a1*dot26*dot4[ih1]+4.*a1*dot18*dot3*dot4[ih1]-
			       4.*dot11[ih2]*dot15*dot5+4.*a1*dot18*dot4[ih1]*dot5-2.*dot17*dot6+2.*dot25*dot6-
			       4.*dot18*dot4[ih1]*dot6+4.*a1*dot18*dot4[ih1]*dot6+4.*a1*dot1*dot15*dot7[ih2]-
			       4.*a1*dot15*dot2*dot7[ih2]+2.*a1*dot27*dot7[ih2]-4.*a1*dot15*dot3*dot7[ih2]+
			       4.*a1*dot14*dot4[ih1]*dot7[ih2]-4.*a12*dot14*dot4[ih1]*dot7[ih2]-4.*a1*dot15*dot5*dot7[ih2]+
			       4.*a1*dot11[ih2]*dot14*dot8[ih1]+4.*a1*dot1*dot18*dot8[ih1]-4.*a1*dot18*dot2*dot8[ih1]+
			       2.*a1*dot26*dot8[ih1]-4.*a1*dot18*dot3*dot8[ih1]-4.*a1*dot18*dot5*dot8[ih1]-
			       4.*a1*dot18*dot6*dot8[ih1]+4.*a12*dot14*dot7[ih2]*dot8[ih1]-4.*a1*dot14*dot4[ih1]*dot9[ih2]+
			       4.*dot15*dot6*dot9[ih2]+4.*a1*dot14*dot8[ih1]*dot9[ih2]+
			       (4.*dot1*(-(a2*dot12*dot13)+dot16)+2.*dot11[ih2]*(-dot21+dot24)+3.*dot28-
				2.*((2.*(-1.+2.*a1)*dot11[ih2]*dot13-2.*dot20+dot23)*dot4[ih1]+
				    dot12*(dot22+2.*dot13*(-(a2*dot2)+(-2.+a1)*dot3-dot5))+
				    dot16*(dot2+3.*(dot3+dot5)+dot6)+(a1*dot21-dot24+2.*(1.-2.*a1*a2)*dot13*dot4[ih1])*dot7[ih2])+
				4.*(dot23-2.*a2*dot13*(dot11[ih2]+a1*dot7[ih2]))*dot8[ih1]-
				2.*(dot21+dot24-4.*a2*dot13*(dot4[ih1]-dot8[ih1]))*dot9[ih2])*M+
			       2.*(2.*a1*dot12*dot13-dot16-2.*a1*dot16)*pow<3,1>(M)-
			       (2.*dot11[ih2]*dot15-2.*dot12*dot19+dot25+
				2.*(dot17-a1*dot17+a1*(2.*dot12*dot19-2.*dot25+dot15*dot7[ih2]+2.*dot18*dot8[ih1])-dot15*dot9[ih2]))*M2-
			       2.*a2*dot10[ih1]*(2.*dot1*dot18+dot26+2.*dot11[ih2]*(dot14+2.*dot13*M)+
						 2.*(a1*dot7[ih2]+dot9[ih2])*(dot14+2.*dot13*M)-
						 2.*dot18*(dot2+dot3+dot5+dot6+M2)))
	    /(a1*a2*(dot1-dot2-dot3)*(dot5+a1*M2)*(dot1-a2*(dot5+dot6+M2)));
	  diag[21]=M2*(4.*dot10[ih1]*dot11[ih2]*dot19+4.*dot12*dot14*dot2-4.*dot11[ih2]*dot15*dot2-3.*dot17*dot2+a1*dot17*dot2-
		       6.*dot10[ih1]*dot18*dot2+6.*a1*dot10[ih1]*dot18*dot2-4.*dot12*dot19*dot2+2.*dot2*dot25-
		       4.*a1*dot10[ih1]*dot26+2.*dot11[ih2]*dot27-4.*dot11[ih2]*dot15*dot3+dot17*dot3+a1*dot17*dot3+
		       2.*dot10[ih1]*dot18*dot3+6.*a1*dot10[ih1]*dot18*dot3-4.*dot12*dot19*dot3+2.*dot25*dot3+
		       4.*dot11[ih2]*dot14*dot4[ih1]-4.*dot11[ih2]*dot19*dot4[ih1]+6.*dot18*dot2*dot4[ih1]-
		       6.*a1*dot18*dot2*dot4[ih1]-2.*dot26*dot4[ih1]+4.*a1*dot26*dot4[ih1]-2.*dot18*dot3*dot4[ih1]-
		       6.*a1*dot18*dot3*dot4[ih1]-4.*dot11[ih2]*dot15*dot5+2.*dot17*dot5+4.*dot10[ih1]*dot18*dot5+
		       4.*a1*dot10[ih1]*dot18*dot5-2.*dot25*dot5-4.*a1*dot18*dot4[ih1]*dot5+4.*dot12*dot14*dot6-
		       4.*dot11[ih2]*dot15*dot6-2.*dot17*dot6-4.*dot10[ih1]*dot18*dot6+4.*a1*dot10[ih1]*dot18*dot6-
		       4.*dot12*dot19*dot6+2.*dot25*dot6+4.*dot18*dot4[ih1]*dot6-4.*a1*dot18*dot4[ih1]*dot6+
		       4.*a1*dot10[ih1]*dot14*dot7[ih2]-4.*a12*dot10[ih1]*dot14*dot7[ih2]+4.*a1*dot10[ih1]*dot19*dot7[ih2]-
		       4.*dot15*dot2*dot7[ih2]+4.*a1*dot15*dot2*dot7[ih2]-2.*a1*dot27*dot7[ih2]+4.*a1*dot15*dot3*dot7[ih2]-
		       4.*a1*dot14*dot4[ih1]*dot7[ih2]+4.*a12*dot14*dot4[ih1]*dot7[ih2]+4.*dot19*dot4[ih1]*dot7[ih2]-
		       4.*a1*dot19*dot4[ih1]*dot7[ih2]+4.*a1*dot15*dot5*dot7[ih2]-4.*dot15*dot6*dot7[ih2]+
		       4.*dot11[ih2]*dot19*dot8[ih1]-10.*dot18*dot2*dot8[ih1]+6.*a1*dot18*dot2*dot8[ih1]+
		       4.*dot26*dot8[ih1]-4.*a1*dot26*dot8[ih1]-2.*dot18*dot3*dot8[ih1]+6.*a1*dot18*dot3*dot8[ih1]+
		       4.*a1*dot18*dot5*dot8[ih1]-8.*dot18*dot6*dot8[ih1]+4.*a1*dot18*dot6*dot8[ih1]-
		       4.*a12*dot14*dot7[ih2]*dot8[ih1]-4.*dot19*dot7[ih2]*dot8[ih1]+4.*a1*dot19*dot7[ih2]*dot8[ih1]+
		       4.*dot10[ih1]*dot14*dot9[ih2]-4.*a1*dot10[ih1]*dot14*dot9[ih2]-4.*dot10[ih1]*dot19*dot9[ih2]+
		       4.*dot15*dot2*dot9[ih2]-4.*dot27*dot9[ih2]+4.*dot15*dot3*dot9[ih2]-4.*dot14*dot4[ih1]*dot9[ih2]+
		       4.*a1*dot14*dot4[ih1]*dot9[ih2]+4.*dot19*dot4[ih1]*dot9[ih2]+4.*dot15*dot5*dot9[ih2]-
		       4.*a1*dot14*dot8[ih1]*dot9[ih2]-4.*dot19*dot8[ih1]*dot9[ih2]+
		       (dot16*((-3.+a1)*dot2+dot3+a1*dot3+2.*dot5-2.*dot6)+4.*dot12*dot13*(dot2+dot6)+
			2.*(-(dot11[ih2]*dot24)+2.*dot11[ih2]*dot13*dot4[ih1]+dot23*dot4[ih1]-2.*dot23*dot8[ih1]-
			    4.*a12*dot13*dot7[ih2]*(dot10[ih1]-dot4[ih1]+dot8[ih1])+
			    (dot21+dot24+4.*dot13*(dot10[ih1]-dot4[ih1]+dot8[ih1]))*dot9[ih2]+
			    a1*(dot21*dot7[ih2]+2.*dot10[ih1]*(dot23+2.*dot13*dot7[ih2]-2.*dot13*dot9[ih2])-
				2.*(dot4[ih1]-dot8[ih1])*(dot23+2.*dot13*dot7[ih2]-2.*dot13*dot9[ih2]))))*M-
		       2.*(-2.*a2*dot12*dot13+dot16-2.*a1*dot16)*pow<3,1>(M)+
		       dot1*(4.*dot11[ih2]*dot15+dot17-a1*dot17-
			     2.*((-1.+3.*a1)*dot10[ih1]*dot18+dot25+
				 dot18*(dot4[ih1]-3.*a1*dot4[ih1]-3.*a2*dot8[ih1])+
				 2.*dot15*(-(a2*dot7[ih2])+dot9[ih2]))+
			     a2*dot16*M-4.*dot12*(dot14-dot19+dot13*M))+
		       2.*(-dot17-2.*dot10[ih1]*dot18+2.*a2*dot12*(dot14-dot19)+dot25+2.*dot18*dot4[ih1]-
			   2.*dot15*(dot11[ih2]+dot7[ih2])-4.*dot18*dot8[ih1]+
			   a1*(2.*dot17+6.*dot10[ih1]*dot18-2.*dot25-4.*dot18*dot4[ih1]+
			       3.*dot15*dot7[ih2]+6.*dot18*dot8[ih1])+dot15*dot9[ih2])*M2)
	    /(16.*a1*a2*(dot1-dot2-dot3)*(dot5+a1*M2)*(-(a2*dot1)+dot2+dot6-a1*(dot2+dot3+dot6)+M2+(-2.+a1)*a1*M2));
	  diag[22]=-0.125*M2*(-4.*dot12*dot14*(dot2+dot6)+
			      2.*(-(dot11[ih2]*dot27)+dot17*(dot2+dot6)+2.*dot11[ih2]*dot15*(dot2+dot3+dot5+dot6)-
				  a1*dot27*dot7[ih2]+2.*a1*dot15*(dot2+dot3+dot5)*dot7[ih2]+
				  dot6*(2.*dot12*dot19-dot25+2.*dot15*dot7[ih2]))-
			      (-2.*dot12*dot22+dot28+4.*dot12*dot13*(dot2+dot3+dot5+dot6)-
			       2.*dot16*(dot2+dot3+dot5+dot6)-2.*a2*(dot21-2.*dot24)*dot7[ih2])*M+
			      2.*(-2.*dot12*dot13+dot16)*pow<3,1>(M)-
			      2.*dot1*(2.*dot15*(dot11[ih2]+a1*dot7[ih2])+(-2.*dot12*dot13+dot16)*M)+
			      (4.*dot11[ih2]*dot15+2.*dot17-2.*a1*dot17+2.*dot12*(-2.*a2*dot14+dot19)-
			       dot25+2.*(3.+a1*(-5.+4.*a1))*dot15*dot7[ih2])*M2)
	    /(a2*dot2*(-dot1+dot2+dot3+dot5+dot6+M2)*(-dot2+a1*(-dot1+dot2+dot3+dot5+a1*M2)));
	  diag[23]=-0.125*(M2*(-2.*dot1*dot17-4.*dot1*dot10[ih1]*dot18+2.*dot17*dot2+4.*dot10[ih1]*dot18*dot2+
			       2.*dot1*dot25-2.*dot2*dot25-2.*dot25*dot3+4.*dot18*dot3*dot4[ih1]+2.*dot17*dot6+
			       4.*dot10[ih1]*dot18*dot6-2.*dot25*dot6+4.*dot10[ih1]*dot19*dot7[ih2]-
			       4.*a1*dot10[ih1]*dot19*dot7[ih2]+2.*dot27*dot7[ih2]-2.*a1*dot27*dot7[ih2]+
			       4.*dot14*dot4[ih1]*dot7[ih2]-4.*a1*dot14*dot4[ih1]*dot7[ih2]-4.*dot19*dot4[ih1]*dot7[ih2]+
			       4.*a1*dot19*dot4[ih1]*dot7[ih2]-4.*dot15*dot5*dot7[ih2]+4.*a1*dot15*dot5*dot7[ih2]-
			       4.*dot1*dot18*dot8[ih1]+4.*dot18*dot2*dot8[ih1]+4.*dot18*dot6*dot8[ih1]+
			       4.*dot19*dot7[ih2]*dot8[ih1]-4.*a1*dot19*dot7[ih2]*dot8[ih1]+4.*dot1*dot15*dot9[ih2]+
			       4.*dot10[ih1]*dot19*dot9[ih2]-4.*dot15*dot2*dot9[ih2]+2.*dot27*dot9[ih2]-
			       4.*dot15*dot3*dot9[ih2]+4.*dot14*dot4[ih1]*dot9[ih2]-4.*dot19*dot4[ih1]*dot9[ih2]-
			       4.*dot15*dot5*dot9[ih2]-4.*dot15*dot6*dot9[ih2]+4.*dot19*dot8[ih1]*dot9[ih2]-
			       (2.*dot1*dot16+dot28-2.*dot16*(dot2+dot3+dot5+dot6)+
				2.*dot10[ih1]*(dot23-2.*a2*dot13*dot7[ih2])+
				2.*(dot20*dot4[ih1]+dot23*(-dot4[ih1]+dot8[ih1])-
				    a2*dot7[ih2]*(dot21-2.*dot24+2.*dot13*(dot4[ih1]+dot8[ih1]))))*M+
			       2.*dot16*pow<3,1>(M)-
			       (-4.*dot10[ih1]*dot18+dot25-2.*dot15*dot7[ih2]+
				2.*(-(a2*dot17)+dot18*dot4[ih1]+5.*a1*dot15*dot7[ih2]-4.*a12*dot15*dot7[ih2]-
				    2.*dot18*dot8[ih1]+2.*a1*dot18*(dot10[ih1]-dot4[ih1]+dot8[ih1])+
				    2.*dot15*dot9[ih2]))*M2))
	    /(a1*a2*dot3*(-dot1+dot2+dot3+dot5+dot6+M2)*(dot5+a1*M2));
	  diag[24]=-0.0625*M2*(4.*dot12*dot14*(dot3+dot5)-
			       2.*(dot17*(dot3+dot5)+dot11[ih2]*(dot27-2.*dot15*(dot2+dot3+dot5+dot6))+
				   dot10[ih1]*(-dot26+2.*dot18*(dot2+dot3+dot5+dot6)+2.*dot14*dot7[ih2])+
				   dot6*(-2.*dot12*dot19+dot25-2.*dot15*dot7[ih2]))+
			       (-3.*dot28+4.*dot12*(dot22-dot13*(dot3+dot5))+2.*dot16*(dot2+3.*(dot3+dot5)+dot6)-
				2.*(dot20*dot4[ih1]+dot11[ih2]*(dot24-2.*dot13*dot4[ih1])+
				    (4.*dot10[ih1]*dot13+dot21+dot24)*dot7[ih2]-
				    (dot23+2.*dot13*dot7[ih2])*(dot4[ih1]-2.*dot8[ih1])-dot24*dot9[ih2]))*M+
			       2.*dot16*pow<3,1>(M)-4.*a12*dot7[ih2]*(dot10[ih1]-dot4[ih1]+dot8[ih1])*(dot14+2.*dot13*M)+
			       2.*dot1*(-2.*dot11[ih2]*dot15+dot17+2.*dot10[ih1]*dot18-
					a1*(dot17+2.*dot15*dot7[ih2]+2.*dot18*(dot10[ih1]-dot4[ih1]+dot8[ih1]))-
					2.*dot16*M+2.*dot12*(-(a2*dot14)+a1*dot13*M))+
			       (4.*dot11[ih2]*dot15+dot25-2.*dot18*(2.*dot10[ih1]+dot4[ih1])+2.*dot15*dot7[ih2])*M2+
			       2.*a1*(2.*dot10[ih1]*dot18*dot2-dot10[ih1]*dot26+2.*dot10[ih1]*dot18*dot3-2.*dot18*dot2*dot4[ih1]+
				      dot26*dot4[ih1]-2.*dot18*dot3*dot4[ih1]+2.*dot10[ih1]*dot18*dot5-2.*dot18*dot4[ih1]*dot5+
				      2.*dot10[ih1]*dot18*dot6-2.*dot18*dot4[ih1]*dot6+dot17*(dot2+dot3+dot5+dot6)+
				      4.*dot10[ih1]*dot14*dot7[ih2]+2.*dot15*dot2*dot7[ih2]-dot27*dot7[ih2]+
				      2.*dot15*dot3*dot7[ih2]-2.*dot14*dot4[ih1]*dot7[ih2]+2.*dot15*dot5*dot7[ih2]+
				      2.*dot18*dot2*dot8[ih1]-dot26*dot8[ih1]+2.*dot18*dot3*dot8[ih1]+
				      2.*dot18*dot5*dot8[ih1]+2.*dot18*dot6*dot8[ih1]+2.*dot14*dot7[ih2]*dot8[ih1]+
				      dot7[ih2]*(8.*dot10[ih1]*dot13+dot21-4.*dot13*dot4[ih1]+8.*dot13*dot8[ih1])*M+
				      2.*dot16*pow<3,1>(M)+(-2.*dot25+dot15*dot7[ih2]+2.*dot18*(dot10[ih1]+dot8[ih1]))*M2+
				      2.*dot12*(-(dot14*(dot2+dot3+dot5+dot6))+M*(dot19*M-dot13*(dot2+dot3+dot5+dot6+2.*M2)))))
	    /(a1*a2*(-dot1+dot2+dot3+dot5+dot6+M2)*(dot5+a1*M2)*(dot1-a2*(dot5+dot6+M2)));
	  diag[25]=M2*(4.*dot10[ih1]*dot11[ih2]*dot19-a2*(-2.*dot12*dot14+dot17)*dot2+4.*a1*dot10[ih1]*dot18*dot2-
		       2.*dot2*dot25-2.*dot10[ih1]*dot26-2.*a1*dot10[ih1]*dot26+2.*dot11[ih2]*dot27-2.*dot12*dot14*dot3-
		       2.*a1*dot12*dot14*dot3+3.*dot17*dot3+a1*dot17*dot3+4.*dot10[ih1]*dot18*dot3+
		       4.*a1*dot10[ih1]*dot18*dot3-2.*dot25*dot3+4.*dot11[ih2]*dot14*dot4[ih1]-
		       4.*dot11[ih2]*dot19*dot4[ih1]+4.*dot18*dot2*dot4[ih1]-4.*a1*dot18*dot2*dot4[ih1]+
		       2.*a1*dot26*dot4[ih1]-4.*a1*dot18*dot3*dot4[ih1]-2.*dot12*dot14*dot5-2.*a1*dot12*dot14*dot5+
		       3.*dot17*dot5+a1*dot17*dot5+4.*dot10[ih1]*dot18*dot5+4.*a1*dot10[ih1]*dot18*dot5+
		       4.*dot12*dot19*dot5-4.*dot25*dot5-4.*a1*dot18*dot4[ih1]*dot5+2.*dot12*dot14*dot6-
		       2.*a1*dot12*dot14*dot6-dot17*dot6+a1*dot17*dot6+4.*a1*dot10[ih1]*dot18*dot6-
		       4.*a1*dot18*dot4[ih1]*dot6-4.*dot10[ih1]*dot14*dot7[ih2]+8.*a1*dot10[ih1]*dot14*dot7[ih2]-
		       4.*a12*dot10[ih1]*dot14*dot7[ih2]+8.*dot10[ih1]*dot19*dot7[ih2]-4.*a1*dot10[ih1]*dot19*dot7[ih2]-
		       6.*dot15*dot2*dot7[ih2]+6.*a1*dot15*dot2*dot7[ih2]+2.*dot27*dot7[ih2]-4.*a1*dot27*dot7[ih2]+
		       2.*dot15*dot3*dot7[ih2]+6.*a1*dot15*dot3*dot7[ih2]+4.*dot14*dot4[ih1]*dot7[ih2]-
		       8.*a1*dot14*dot4[ih1]*dot7[ih2]+4.*a12*dot14*dot4[ih1]*dot7[ih2]-4.*dot19*dot4[ih1]*dot7[ih2]+
		       4.*a1*dot19*dot4[ih1]*dot7[ih2]-2.*dot15*dot5*dot7[ih2]+10.*a1*dot15*dot5*dot7[ih2]-
		       2.*dot15*dot6*dot7[ih2]+2.*a1*dot15*dot6*dot7[ih2]+4.*dot11[ih2]*dot19*dot8[ih1]-
		       4.*dot18*dot2*dot8[ih1]+4.*a1*dot18*dot2*dot8[ih1]+4.*dot26*dot8[ih1]-2.*a1*dot26*dot8[ih1]+
		       4.*a1*dot18*dot3*dot8[ih1]+4.*a1*dot18*dot5*dot8[ih1]-4.*dot18*dot6*dot8[ih1]+
		       4.*a1*dot18*dot6*dot8[ih1]+4.*a1*dot14*dot7[ih2]*dot8[ih1]-4.*a12*dot14*dot7[ih2]*dot8[ih1]-
		       4.*a1*dot19*dot7[ih2]*dot8[ih1]-2.*dot27*dot9[ih2]+
		       (2.*dot10[ih1]*dot23+4.*dot11[ih2]*dot13*dot4[ih1]+
			dot16*(-(a2*dot2)+(3.+a1)*(dot3+dot5)-a2*dot6)-
			2.*dot12*dot13*(-(a2*dot2)+dot3+dot5-dot6+a1*(dot3+dot5+dot6))-
			2.*(dot11[ih2]*dot24+dot23*dot8[ih1])-
			2.*(-1.+2.*a1)*dot7[ih2]*(-dot24-2.*a2*dot13*(dot10[ih1]-dot4[ih1]+dot8[ih1]))+
			2.*dot24*dot9[ih2])*M-(2.*(-1.+3.*a1)*dot12*dot13+dot16-5.*a1*dot16)*pow<3,1>(M)+
		       dot1*(-dot17+2.*(-2.*dot10[ih1]*dot18+dot25+dot15*dot7[ih2])-
			     a1*(dot17+6.*dot15*dot7[ih2]+4.*dot18*(dot10[ih1]-dot4[ih1]+dot8[ih1]))-
			     (1.+a1)*dot16*M-2.*a2*dot12*(dot14+dot13*M))+
		       (-dot17+2.*dot12*(dot14-3.*a1*dot14+2.*a1*dot19)-2.*dot15*dot7[ih2]+
			a1*(5.*dot17+8.*dot10[ih1]*dot18-4.*(dot25+dot18*dot4[ih1])+
			    2.*(1.+4.*a1)*dot15*dot7[ih2])+4.*(-1.+2.*a1)*dot18*dot8[ih1])*M2)
	    /(16.*a1*a2*(-dot1+dot2+dot3+dot5+dot6+M2)*(dot5+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot5+a1*M2)));
	  diag[26]=-0.125*M2*(-4.*dot12*dot14*(dot3+dot5)+2.*dot17*(dot3+dot5)+4.*dot11[ih2]*dot15*(dot2+dot3+dot5)+
			      4.*dot12*dot19*(dot2+dot3+dot5)-2.*dot25*(dot2+dot3+dot5)+4.*dot15*(dot3+dot5)*dot7[ih2]-
			      (2.*dot11[ih2]*dot21-2.*dot12*dot22+dot28+
			       8.*dot12*dot13*(dot3+dot5)-4.*dot16*(dot3+dot5)+2.*dot24*dot7[ih2])*M+
			      4.*a1*(-2.*dot12*dot13+dot16)*pow<3,1>(M)+
			      2.*dot1*(-dot17+dot25-2.*dot15*(dot11[ih2]+dot7[ih2])-2.*dot16*M+2.*dot12*(dot14-dot19+2.*dot13*M))+
			      (-2.*dot11[ih2]*dot15-2.*dot12*dot19+dot25+
			       2.*a1*(-2.*dot12*dot14+4.*dot11[ih2]*dot15+dot17+4.*dot12*dot19-2.*dot25+2.*dot15*dot7[ih2]))*M2)
	    /(dot2*(-dot1+dot3+dot5+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot5+a1*M2)));
	  diag[27]=M2*(2.*dot6*(dot17+2.*dot12*(-dot14+dot19)-dot25+2.*dot15*(dot11[ih2]+dot7[ih2]))+
		       (2.*dot1*(-2.*dot12*dot13+dot16)-2.*dot12*dot22+dot28+4.*dot12*dot13*(dot2+dot3+dot5-dot6)-
			2.*dot16*(dot2+dot3+dot5-dot6)+2.*dot21*(dot11[ih2]+dot7[ih2]))*M+
		       2.*(-1.+2.*a1)*(2.*dot12*dot13-dot16)*pow<3,1>(M)+
		       (6.*dot11[ih2]*dot15+2.*dot17+dot12*(-4.*a2*dot14+(6.-8.*a1)*dot19)-
			3.*dot25+6.*dot15*dot7[ih2]-2.*a1*(dot17-2.*dot25+4.*dot15*(dot11[ih2]+dot7[ih2])))*M2)
	    /(8.*a2*dot2*sqr(dot6+a2*M2));
	  diag[28]=-0.125*M2*(-2.*dot12*dot14*dot2+dot17*dot2-4.*dot12*dot19*dot2+2.*dot2*dot25-2.*dot11[ih2]*dot27+
			      4.*dot12*dot14*dot3-2.*dot17*dot3-4.*dot12*dot19*dot3+2.*dot25*dot3+4.*dot12*dot14*dot5-
			      2.*dot17*dot5-4.*dot12*dot19*dot5+2.*dot25*dot5-4.*dot12*dot14*dot6+4.*dot11[ih2]*dot15*dot6+
			      2.*dot17*dot6+4.*dot12*dot19*dot6-2.*dot25*dot6-2.*dot15*dot2*dot7[ih2]+
			      4.*a1*dot15*dot2*dot7[ih2]-2.*a1*dot27*dot7[ih2]-4.*dot15*dot3*dot7[ih2]+
			      4.*a1*dot15*dot3*dot7[ih2]-4.*dot15*dot5*dot7[ih2]+4.*a1*dot15*dot5*dot7[ih2]+4.*dot15*dot6*dot7[ih2]+
			      (2.*dot11[ih2]*dot21-(2.*dot12*dot13-dot16)*(dot2-2.*(dot3+dot5-dot6))-
			       2.*(-(a2*dot21)+dot24-2.*a1*dot24)*dot7[ih2])*M+2.*(-1.+2.*a1)*(2.*dot12*dot13-dot16)*pow<3,1>(M)+
			      2.*dot1*(dot17-dot25+2.*a2*dot15*dot7[ih2]+dot16*M-2.*dot12*(dot14-dot19+dot13*M))+
			      2.*(3.*dot11[ih2]*dot15+dot17+2.*(-1.+2.*a1)*dot12*(dot14-dot19)-
				  2.*a1*(2.*dot11[ih2]*dot15+dot17-dot25)-dot25-(-3.+4.*a1)*a2*dot15*dot7[ih2])*M2)
	    /(a2*dot2*(-dot6-a2*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot5+a1*M2)));
	  diag[29]=-0.125*M2*(2.*(2.*dot10[ih1]*dot18*dot2-dot2*dot25-dot25*dot3+2.*dot18*dot3*dot4[ih1]+
				  2.*dot10[ih1]*dot18*dot6-dot25*dot6+dot17*(dot2+dot6)+2.*dot18*dot2*dot8[ih1]+
				  2.*dot18*dot6*dot8[ih1]-2.*dot15*(dot2+dot3+dot6)*dot9[ih2]+
				  dot1*(-dot17+dot25-2.*dot18*(dot10[ih1]+dot8[ih1])+2.*dot15*dot9[ih2]))+
			      (dot28+2.*dot20*dot4[ih1]-2.*dot16*(dot1-dot2+dot3+dot5-dot6)+
			       2.*dot23*(dot10[ih1]-dot4[ih1]+dot8[ih1])-
			       2.*(dot21+2.*dot13*(dot10[ih1]-dot4[ih1]+dot8[ih1]))*dot9[ih2])*M+
			      2.*(1.-2.*a1)*dot16*pow<3,1>(M)+
			      (4.*dot10[ih1]*dot18-3.*(dot25+2.*dot15*dot9[ih2])+
			       2.*(dot17-a1*dot17+dot18*(dot4[ih1]+2.*dot8[ih1])+
				   2.*a1*(dot25-dot18*(dot10[ih1]+dot4[ih1]+dot8[ih1])+2.*dot15*dot9[ih2])))*M2)
	    /(dot3*(-(a2*dot1)+dot2+dot6-a1*(dot2+dot3+dot6)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot6-a2*M2));
	  diag[30]=-0.125*M2*(2.*dot6*(dot25-2.*dot18*dot4[ih1]+2.*dot15*dot9[ih2])+
			      (-dot28-2.*dot20*dot4[ih1]+2.*dot16*(-dot1+dot2+dot3+dot5+dot6)+2.*dot21*dot9[ih2])*M+
			      2.*dot16*pow<3,1>(M)-(-3.+4.*a1)*(dot25-2.*dot18*dot4[ih1]+2.*dot15*dot9[ih2])*M2)
	    /(a2*(dot1-dot2-dot6-a2*M2)*sqr(dot6+a2*M2));
	  diag[31]=M2*(3.*dot17*dot2-2.*dot10[ih1]*dot18*dot2+4.*a1*dot10[ih1]*dot18*dot2-2.*dot2*dot25+2.*dot10[ih1]*dot26-
		       2.*a1*dot10[ih1]*dot26-4.*dot10[ih1]*dot18*dot3+4.*a1*dot10[ih1]*dot18*dot3-2.*dot25*dot3+2.*dot18*dot2*dot4[ih1]-
		       4.*a1*dot18*dot2*dot4[ih1]+2.*a1*dot26*dot4[ih1]+4.*dot18*dot3*dot4[ih1]-4.*a1*dot18*dot3*dot4[ih1]-
		       4.*dot10[ih1]*dot18*dot5+4.*a1*dot10[ih1]*dot18*dot5-4.*a1*dot18*dot4[ih1]*dot5+3.*dot17*dot6-
		       2.*dot10[ih1]*dot18*dot6+4.*a1*dot10[ih1]*dot18*dot6-4.*dot25*dot6+6.*dot18*dot4[ih1]*dot6-
		       4.*a1*dot18*dot4[ih1]*dot6+2.*dot18*dot2*dot8[ih1]+4.*a1*dot18*dot2*dot8[ih1]-2.*a1*dot26*dot8[ih1]+
		       4.*a1*dot18*dot3*dot8[ih1]+4.*a1*dot18*dot5*dot8[ih1]+2.*dot18*dot6*dot8[ih1]+4.*a1*dot18*dot6*dot8[ih1]+
		       4.*dot10[ih1]*dot14*dot9[ih2]-4.*a1*dot10[ih1]*dot14*dot9[ih2]-4.*dot15*dot2*dot9[ih2]-
		       4.*dot15*dot3*dot9[ih2]+4.*a1*dot14*dot4[ih1]*dot9[ih2]-8.*dot15*dot6*dot9[ih2]-
		       4.*a1*dot14*dot8[ih1]*dot9[ih2]+
		       (3.*dot16*(dot2+dot6)+2.*(-1.+2.*a1)*(dot10[ih1]-dot4[ih1]+dot8[ih1])*(dot23-2.*dot13*dot9[ih2]))*M+
		       3.*a2*dot16*pow<3,1>(M)-
		       dot1*(3.*dot17-2.*dot25+
			     2.*dot18*(-dot10[ih1]+2.*a1*dot10[ih1]+dot4[ih1]-2.*a1*dot4[ih1]+dot8[ih1]+2.*a1*dot8[ih1])-
			     4.*dot15*dot9[ih2]+3.*dot16*M)+
		       (3.*a2*dot17+2.*a1*dot18*(dot10[ih1]-5.*dot4[ih1]+dot8[ih1])+
			2.*dot18*(3.*dot4[ih1]+dot8[ih1])+4.*a1*(dot25+2.*dot15*dot9[ih2])-
			2.*(dot10[ih1]*dot18+2.*dot25+4.*dot15*dot9[ih2]))*M2)
	    /(8.*a2*(-(a2*dot1)+dot2+dot6-a1*(dot2+dot3+dot6)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot6-a2*M2)*(dot6+a2*M2));
	  diag[32]=-0.125*M2*(2.*dot5*(dot17-dot25+2.*dot18*(dot10[ih1]+dot8[ih1])-2.*dot15*dot9[ih2])+
			      (-2.*dot10[ih1]*dot23-dot28-2.*dot20*dot4[ih1]+4.*dot16*dot5-2.*dot24*dot7[ih2]+
			       2.*dot23*(dot4[ih1]-dot8[ih1])-4.*dot13*dot4[ih1]*dot9[ih2]+4.*dot10[ih1]*dot13*(dot7[ih2]+dot9[ih2])+
			       2.*dot21*(dot7[ih2]+dot9[ih2])+4.*dot13*dot8[ih1]*(dot7[ih2]+dot9[ih2]))*M+
			      4.*a1*dot16*pow<3,1>(M)+
			      (dot25-2.*dot18*dot4[ih1]+2.*dot15*(dot7[ih2]+dot9[ih2])+
			       2.*a1*(dot17+2.*dot18*(dot10[ih1]+dot4[ih1]+dot8[ih1])-
				      2.*(dot25+dot15*dot7[ih2]+2.*dot15*dot9[ih2])))*M2)
	    /(a1*dot3*sqr(dot5+a1*M2));
	  diag[33]=-0.125*M2*(-2.*dot17*dot2-4.*dot10[ih1]*dot18*dot2+2.*dot2*dot25+dot17*dot3+2.*dot10[ih1]*dot18*dot3+
			      2.*dot25*dot3-6.*dot18*dot3*dot4[ih1]+2.*dot17*dot5+4.*dot10[ih1]*dot18*dot5-2.*dot25*dot5-
			      2.*dot17*dot6-4.*dot10[ih1]*dot18*dot6+2.*dot25*dot6-4.*dot10[ih1]*dot19*dot7[ih2]+
			      4.*a1*dot10[ih1]*dot19*dot7[ih2]-2.*dot27*dot7[ih2]+2.*a1*dot27*dot7[ih2]+4.*dot15*dot3*dot7[ih2]-
			      4.*dot14*dot4[ih1]*dot7[ih2]+4.*a1*dot14*dot4[ih1]*dot7[ih2]+4.*dot19*dot4[ih1]*dot7[ih2]-
			      4.*a1*dot19*dot4[ih1]*dot7[ih2]+4.*dot15*dot5*dot7[ih2]-4.*a1*dot15*dot5*dot7[ih2]-
			      4.*dot18*dot2*dot8[ih1]+2.*dot18*dot3*dot8[ih1]+4.*dot18*dot5*dot8[ih1]-4.*dot18*dot6*dot8[ih1]-
			      4.*dot19*dot7[ih2]*dot8[ih1]+4.*a1*dot19*dot7[ih2]*dot8[ih1]-4.*dot10[ih1]*dot19*dot9[ih2]+
			      4.*dot15*dot2*dot9[ih2]-2.*dot27*dot9[ih2]+4.*dot15*dot3*dot9[ih2]-4.*dot14*dot4[ih1]*dot9[ih2]+
			      4.*dot19*dot4[ih1]*dot9[ih2]+4.*dot15*dot6*dot9[ih2]-4.*dot19*dot8[ih1]*dot9[ih2]+
			      (dot16*(-2.*dot2+dot3+2.*dot5-2.*dot6)+
			       2.*dot7[ih2]*(dot24-2.*dot13*dot4[ih1]+a1*(dot21-2.*dot24+2.*dot13*(dot10[ih1]+dot4[ih1]+dot8[ih1])))+
			       2.*(dot21+2.*dot13*(dot10[ih1]-dot4[ih1]+dot8[ih1]))*dot9[ih2])*M+
			      2.*(-1.+2.*a1)*dot16*pow<3,1>(M)+2.*dot1*(dot17-dot25+2.*dot18*(dot10[ih1]+dot8[ih1])-2.*dot15*dot9[ih2]+dot16*M)+
			      2.*((-1.+2.*a1)*dot17-2.*dot10[ih1]*dot18+dot25-4.*a12*dot15*dot7[ih2]-2.*dot18*dot8[ih1]+3.*dot15*dot9[ih2]+
				  a1*(-2.*dot25+3.*dot15*dot7[ih2]+4.*dot18*(dot10[ih1]+dot8[ih1])-4.*dot15*dot9[ih2]))*M2)
			      /(a1*dot3*(dot5+a1*M2)*(-(a2*dot1)+dot2+dot6-a1*(dot2+dot3+dot6)+M2+(-2.+a1)*a1*M2));
	  diag[34]=-0.125*M2*(2.*(-2.*dot11[ih2]*dot15-2.*dot12*dot19+dot25)*dot5+
			      (-2.*dot12*dot22+2.*dot10[ih1]*(2.*dot11[ih2]*dot13+dot23)+dot28-2.*dot23*dot4[ih1]+2.*dot24*dot7[ih2]+
			       2.*dot23*dot8[ih1]+2.*dot11[ih2]*(dot21-2.*dot13*dot4[ih1]+2.*dot13*dot8[ih1]))*M-
			      (-1.+4.*a1)*(2.*dot11[ih2]*dot15+2.*dot12*dot19-dot25)*M2)
	    /(a1*sqr(dot5+a1*M2)*(-dot1+dot3+dot5+a1*M2));
	  diag[35]=M2*(4.*dot11[ih2]*dot15*dot2+4.*dot12*dot19*dot2-2.*dot2*dot25-6.*dot12*dot14*dot3+4.*dot11[ih2]*dot15*dot3+3.*dot17*dot3+
		       4.*dot12*dot19*dot3-2.*dot25*dot3-4.*dot11[ih2]*dot14*dot4[ih1]+4.*a1*dot11[ih2]*dot14*dot4[ih1]+2.*dot26*dot4[ih1]-
		       2.*a1*dot26*dot4[ih1]-4.*dot18*dot3*dot4[ih1]-6.*dot12*dot14*dot5+8.*dot11[ih2]*dot15*dot5+3.*dot17*dot5+
		       8.*dot12*dot19*dot5-4.*dot25*dot5-4.*dot18*dot4[ih1]*dot5+6.*dot15*dot3*dot7[ih2]-4.*dot19*dot4[ih1]*dot7[ih2]+
		       4.*a1*dot19*dot4[ih1]*dot7[ih2]+6.*dot15*dot5*dot7[ih2]-4.*a1*dot11[ih2]*dot14*dot8[ih1]+2.*a1*dot26*dot8[ih1]+
		       4.*dot18*dot3*dot8[ih1]+4.*dot18*dot5*dot8[ih1]-4.*a1*dot19*dot7[ih2]*dot8[ih1]+
		       (-6.*dot12*dot13*(dot3+dot5)+3.*dot16*(dot3+dot5)+2.*(-1.+2.*a1)*(2.*dot11[ih2]*dot13+dot23)*(dot4[ih1]-dot8[ih1]))*M+
		       3.*a1*(-2.*dot12*dot13+dot16)*pow<3,1>(M)-
		       dot1*(4.*dot11[ih2]*dot15+3.*dot17-2.*dot25-4.*dot18*dot4[ih1]+6.*dot15*dot7[ih2]+
			     4.*dot18*dot8[ih1]+3.*dot16*M+dot12*(-6.*dot14+4.*dot19-6.*dot13*M))+
		       a1*(-6.*dot12*dot14+8.*dot11[ih2]*dot15+3.*dot17+8.*dot12*dot19-4.*dot25-
			   4.*dot18*dot4[ih1]+6.*dot15*dot7[ih2]+4.*dot18*dot8[ih1])*M2+
		       2.*dot10[ih1]*(-(a2*(dot26-2.*dot19*dot7[ih2]))+(1.-2.*a1)*dot23*M+
				      2.*dot11[ih2]*(dot14-a1*dot14+dot13*M-2.*a1*dot13*M)+2.*dot18*(-dot1+dot3+dot5+a1*M2)))
	    /(8.*a1*(dot5+a1*M2)*(-dot1+dot3+dot5+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot5+a1*M2)));
	  diag[36]=-0.0625*M2*(-2.*dot12*dot14*dot2+2.*a1*dot12*dot14*dot2+4.*dot11[ih2]*dot15*dot2+dot17*dot2-a1*dot17*dot2+
			       4.*dot12*dot19*dot2-2.*dot11[ih2]*dot27-4.*dot12*dot14*dot3+2.*a1*dot12*dot14*dot3+
			       4.*dot11[ih2]*dot15*dot3-a1*dot17*dot3+4.*dot12*dot19*dot3-8.*dot11[ih2]*dot14*dot4[ih1]+
			       4.*a1*dot11[ih2]*dot14*dot4[ih1]+4.*dot11[ih2]*dot19*dot4[ih1]-4.*dot18*dot2*dot4[ih1]+
			       4.*a1*dot18*dot2*dot4[ih1]+2.*dot26*dot4[ih1]-4.*a1*dot26*dot4[ih1]-4.*dot18*dot3*dot4[ih1]+
			       4.*a1*dot18*dot3*dot4[ih1]-4.*dot12*dot14*dot5+8.*dot11[ih2]*dot15*dot5+4.*dot12*dot19*dot5-
			       4.*dot18*dot4[ih1]*dot5+4.*a1*dot18*dot4[ih1]*dot5-2.*dot12*dot14*dot6+2.*a1*dot12*dot14*dot6+
			       dot17*dot6-a1*dot17*dot6+4.*a1*dot18*dot4[ih1]*dot6+6.*dot15*dot2*dot7[ih2]-
			       6.*a1*dot15*dot2*dot7[ih2]-2.*dot27*dot7[ih2]+4.*a1*dot27*dot7[ih2]+4.*dot15*dot3*dot7[ih2]-
			       6.*a1*dot15*dot3*dot7[ih2]-4.*dot14*dot4[ih1]*dot7[ih2]+8.*a1*dot14*dot4[ih1]*dot7[ih2]-
			       4.*a12*dot14*dot4[ih1]*dot7[ih2]+8.*dot15*dot5*dot7[ih2]-8.*a1*dot15*dot5*dot7[ih2]+
			       2.*dot15*dot6*dot7[ih2]-2.*a1*dot15*dot6*dot7[ih2]-4.*a1*dot11[ih2]*dot14*dot8[ih1]-
			       4.*dot11[ih2]*dot19*dot8[ih1]+4.*dot18*dot2*dot8[ih1]-4.*a1*dot18*dot2*dot8[ih1]-
			       4.*dot26*dot8[ih1]+4.*a1*dot26*dot8[ih1]+4.*dot18*dot3*dot8[ih1]-4.*a1*dot18*dot3*dot8[ih1]+
			       4.*dot18*dot5*dot8[ih1]-4.*a1*dot18*dot5*dot8[ih1]+4.*dot18*dot6*dot8[ih1]-
			       4.*a1*dot18*dot6*dot8[ih1]-4.*a1*dot14*dot7[ih2]*dot8[ih1]+4.*a12*dot14*dot7[ih2]*dot8[ih1]+2.*dot27*dot9[ih2]+
			       (-8.*dot11[ih2]*dot13*dot4[ih1]+dot16*(dot2-a1*dot2+dot6-a1*(dot3+dot6))+
				2.*dot12*dot13*(-(a2*dot2)-2.*(dot3+dot5)-dot6+a1*(dot3+dot6))+
				2.*dot24*(dot11[ih2]+dot7[ih2])-
				2.*(dot4[ih1]*(dot23+2.*dot13*dot7[ih2])+
				    2.*a1*(dot24*dot7[ih2]-(2.*dot11[ih2]*dot13+dot23+3.*dot13*dot7[ih2])*(dot4[ih1]-dot8[ih1]))+
				    4.*a12*dot13*dot7[ih2]*(dot4[ih1]-dot8[ih1])-2.*(dot23+dot13*(dot11[ih2]+dot7[ih2]))*dot8[ih1]+dot24*dot9[ih2]))*M+
			       (-2.*(1.+a12)*dot12*dot13+a22*dot16)*pow<3,1>(M)+
			       dot1*(-4.*dot11[ih2]*dot15+(-2.+a1)*dot17-8.*dot15*dot7[ih2]+6.*a1*dot15*dot7[ih2]+4.*a1*dot18*dot8[ih1]-
				     4.*dot18*(-(a2*dot4[ih1])+dot8[ih1])-2.*dot16*M+a1*dot16*M-2.*dot12*((-4.+a1)*dot14+2.*dot19+(-4.+a1)*dot13*M))+
			       (dot17-2.*dot12*(dot14+a12*dot14-2.*a1*dot19)+2.*dot15*dot7[ih2]+a12*(dot17-6.*dot15*dot7[ih2])+
				4.*dot18*dot8[ih1]+a1*(8.*dot11[ih2]*dot15-2.*dot17+4.*dot15*dot7[ih2]-4.*dot18*dot8[ih1]))*M2-
			       4.*dot10[ih1]*(-(dot7[ih2]*(dot14-dot19+dot13*M))-a12*dot7[ih2]*(dot14+2.*dot13*M)+
					      dot11[ih2]*(-(a2*dot14)+dot19+(-1.+2.*a1)*dot13*M)+
					      a1*(-dot26+2.*dot14*dot7[ih2]+dot23*M+3.*dot13*dot7[ih2]*M+dot18*(-dot1+dot2+dot3+dot5+dot6+M2))))
	    /(a1*a2*(dot5+a1*M2)*(-dot6-a2*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot5+a1*M2)));
	  diag[37]=-0.0625*M2*(4.*dot12*dot14*dot2-4.*dot11[ih2]*dot15*dot2-dot17*dot2+a1*dot17*dot2-4.*dot12*dot19*dot2+
			       2.*dot11[ih2]*dot27-4.*dot11[ih2]*dot15*dot3+a1*dot17*dot3-4.*dot12*dot19*dot3+
			       4.*dot11[ih2]*dot14*dot4[ih1]-4.*dot11[ih2]*dot19*dot4[ih1]+6.*dot18*dot2*dot4[ih1]-
			       6.*a1*dot18*dot2*dot4[ih1]-2.*dot26*dot4[ih1]+4.*a1*dot26*dot4[ih1]+4.*dot18*dot3*dot4[ih1]-
			       6.*a1*dot18*dot3*dot4[ih1]-4.*dot11[ih2]*dot15*dot5+a1*dot17*dot5-6.*a1*dot18*dot4[ih1]*dot5+
			       4.*dot12*dot14*dot6-4.*dot11[ih2]*dot15*dot6-4.*dot12*dot19*dot6+4.*dot18*dot4[ih1]*dot6-
			       4.*a1*dot18*dot4[ih1]*dot6-4.*dot15*dot2*dot7[ih2]+4.*a1*dot15*dot2*dot7[ih2]+2.*dot27*dot7[ih2]-
			       4.*a1*dot27*dot7[ih2]-4.*dot15*dot3*dot7[ih2]+4.*a1*dot15*dot3*dot7[ih2]+4.*dot14*dot4[ih1]*dot7[ih2]-
			       8.*a1*dot14*dot4[ih1]*dot7[ih2]+4.*a12*dot14*dot4[ih1]*dot7[ih2]-4.*dot15*dot5*dot7[ih2]+
			       8.*a1*dot15*dot5*dot7[ih2]-4.*dot15*dot6*dot7[ih2]+4.*dot11[ih2]*dot19*dot8[ih1]-
			       6.*dot18*dot2*dot8[ih1]+6.*a1*dot18*dot2*dot8[ih1]+4.*dot26*dot8[ih1]-4.*a1*dot26*dot8[ih1]-
			       4.*dot18*dot3*dot8[ih1]+6.*a1*dot18*dot3*dot8[ih1]-4.*dot18*dot5*dot8[ih1]+
			       6.*a1*dot18*dot5*dot8[ih1]-4.*dot18*dot6*dot8[ih1]+4.*a1*dot18*dot6*dot8[ih1]-
			       4.*a12*dot14*dot7[ih2]*dot8[ih1]-2.*dot27*dot9[ih2]+4.*a1*dot14*dot4[ih1]*dot9[ih2]+
			       4.*dot15*dot5*dot9[ih2]-4.*dot15*dot6*dot9[ih2]-4.*a1*dot14*dot8[ih1]*dot9[ih2]+
			       (-(a2*dot16*dot2)+2.*(2.*dot11[ih2]*dot13+dot23)*dot4[ih1]+a1*dot16*(dot3+dot5)+
				4.*dot12*dot13*(dot2+dot6)+4.*dot13*dot4[ih1]*dot7[ih2]-2.*dot24*(dot11[ih2]+dot7[ih2])+
				8.*a12*dot13*dot7[ih2]*(dot4[ih1]-dot8[ih1])-4.*dot23*dot8[ih1]+
				4.*a1*(dot23*(-dot4[ih1]+dot8[ih1])+dot7[ih2]*(dot24+dot13*(-3.*dot4[ih1]+dot8[ih1])))+
				2.*dot24*dot9[ih2]+4.*(-1.+2.*a1)*dot13*(dot4[ih1]-dot8[ih1])*dot9[ih2])*M+
			       (4.*a2*dot12*dot13+a12*dot16)*pow<3,1>(M)-
			       dot1*((1.+a1)*dot17-4.*dot15*(dot11[ih2]+dot7[ih2]-a1*dot7[ih2])+
				     2.*dot18*(dot10[ih1]+3.*a1*dot10[ih1]+dot4[ih1]-3.*a1*dot4[ih1]-dot8[ih1])+
				     6.*a1*dot18*dot8[ih1]+(1.+a1)*dot16*M+4.*dot12*(dot14-dot19+dot13*M))+
			       (4.*a2*dot12*(dot14-dot19)+4.*dot18*dot4[ih1]-4.*dot15*(dot11[ih2]+dot7[ih2])+
				4.*a1*dot18*(-2.*dot4[ih1]+dot8[ih1])+
				a12*(dot17+8.*dot15*dot7[ih2]+2.*dot18*(-dot4[ih1]+dot8[ih1]))+
				8.*a1*dot15*dot9[ih2]-4.*(dot18*dot8[ih1]+dot15*dot9[ih2]))*M2+
			       2.*dot10[ih1]*((-1.+3.*a1)*dot18*dot2+3.*a1*dot18*(dot3+dot5)+2.*a1*dot18*dot6-
					      2.*a12*dot14*dot7[ih2]+2.*dot19*(dot11[ih2]+dot7[ih2])+2.*dot14*dot9[ih2]-
					      2.*(-1.+2.*a1)*dot13*(a1*dot7[ih2]+dot9[ih2])*M-
					      2.*a1*(dot26+dot14*(-dot7[ih2]+dot9[ih2])-dot23*M)+a1*(2.+a1)*dot18*M2))
	    /(a1*a2*(dot5+a1*M2)*(-(a2*dot1)+dot2+dot6-a1*(dot2+dot3+dot6)+M2+(-2.+a1)*a1*M2)*(dot6+a2*M2));
	  // colour flows
	  Complex flow[3] = {diag[0] + diag[8] + diag[11] + diag[12] + diag[13] + diag[16] + diag[17] + diag[20] + diag[21] - diag[24] + diag[25] + diag[36] + diag[37],
			     diag[1] + diag[3] + diag[5] + diag[7] + diag[29] + diag[30] + diag[31]  + diag[32] + diag[33] + 2.*diag[37] +
			     (diag[4] - diag[5] + diag[6] - diag[7] + diag[9] - 2.*diag[11] - 2.*diag[13] + diag[15] + diag[18] + diag[23] + 2.*diag[24] - diag[29] - diag[30] - diag[32])/9.,
			     diag[2] - diag[3] - diag[5] - diag[7] + diag[26] + diag[27] + diag[28]  + diag[34] + diag[35] + 2.*diag[36] +
			     (-diag[4] + diag[5] - diag[6] + diag[7] - 2.*diag[8] + diag[10] + diag[14] - 2.*diag[16] + diag[19] - 2.*diag[20] + diag[22] - diag[26] - diag[27] - diag[34])/9.};
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
