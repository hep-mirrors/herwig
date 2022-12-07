// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGamma2BC1S0QQbarAmplitude class.
//

#include "GammaGamma2BC1S0QQbarAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/MatrixElement/TwoToThreePhaseSpace.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Handlers/EventHandler.h"
#include <numeric>

using namespace Herwig;

IBPtr GammaGamma2BC1S0QQbarAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGamma2BC1S0QQbarAmplitude::fullclone() const {
  return new_ptr(*this);
}

void GammaGamma2BC1S0QQbarAmplitude::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2) << n_;
}

void GammaGamma2BC1S0QQbarAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2) >> n_;
}

void GammaGamma2BC1S0QQbarAmplitude::doinit() {
  GammaGammaAmplitude::doinit();
  // get the non-perturbative ME
  O1_ = params_->singletMEProduction<0>(bcbar,n_,0,0);
}

//The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<GammaGamma2BC1S0QQbarAmplitude,GammaGammaAmplitude>
describeHerwigGammaGamma2BC1S0QQbarAmplitude("Herwig::GammaGamma2BC1S0QQbarAmplitude",
					     "HwOniumParameters.so HwMEGammaGamma.so HwMEGammaGammaOnium.so");

void GammaGamma2BC1S0QQbarAmplitude::Init() {

  static ClassDocumentation<GammaGamma2BC1S0QQbarAmplitude> documentation
    ("The GammaGamma2BC1S0QQbarAmplitude class implements the matrix elements for gamma gamma -> Bc+ b cbar +cc");

  static Parameter<GammaGamma2BC1S0QQbarAmplitude,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &GammaGamma2BC1S0QQbarAmplitude::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Reference<GammaGamma2BC1S0QQbarAmplitude,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &GammaGamma2BC1S0QQbarAmplitude::params_, false, false, true, false, false);
}


Energy GammaGamma2BC1S0QQbarAmplitude::generateW(double r, const tcPDVector & partons,
						 Energy Wmax,Energy2 & jacW, Energy2) {
  // twice B_c mass
  Energy Wmin = 2.*partons[0]->mass();
  Energy W = Wmin*pow(Wmax/Wmin,r);
  jacW = 2.*sqr(W)*log(Wmax/Wmin);
  return W;
}

vector<DiagPtr> GammaGamma2BC1S0QQbarAmplitude::getDiagrams(unsigned int iopt) const {
  vector<DiagPtr> output;
  output.reserve(8);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr b = getParticleData(ParticleID::b);
  tcPDPtr c = getParticleData(ParticleID::c);
  tcPDPtr state = getParticleData(long(541 + (n_-1)*100000));
  // gamma gamma process
  if(iopt==0) {
    output.push_back(new_ptr((Tree2toNDiagram(3), gamma, c->CC(), gamma, 1, c      , 4, state      , 2, c->CC(), 4, b, -1)));  
    output.push_back(new_ptr((Tree2toNDiagram(3), gamma, b->CC(), gamma, 2, b->CC(), 4, state      , 4, c->CC(), 1, b, -2)));
    output.push_back(new_ptr((Tree2toNDiagram(3), gamma, c      , gamma, 2, c      , 4, state      , 1, c->CC(), 4, b, -3))); 
    output.push_back(new_ptr((Tree2toNDiagram(3), gamma, b      , gamma, 1, b->CC(), 4, state      , 4, c->CC(), 2, b, -4)));
    output.push_back(new_ptr((Tree2toNDiagram(3), gamma, b->CC(), gamma, 1, b      , 4, state->CC(), 2, b->CC(), 4, c, -1)));  
    output.push_back(new_ptr((Tree2toNDiagram(3), gamma, c->CC(), gamma, 2, c->CC(), 4, state->CC(), 4, b->CC(), 1, c, -2)));
    output.push_back(new_ptr((Tree2toNDiagram(3), gamma, b      , gamma, 2, b      , 4, state->CC(), 1, b->CC(), 4, c, -3))); 
    output.push_back(new_ptr((Tree2toNDiagram(3), gamma, c      , gamma, 1, c->CC(), 4, state->CC(), 4, b->CC(), 2, c, -4)));
  }
  // e+e- initiated
  else {
    cPDPair in = generator()->eventHandler()->incoming();
    if(in.first->charged() && in.second->charged()) {
      output.push_back(new_ptr((Tree2toNDiagram(5), in.first, gamma, c->CC(), gamma, in.second, 1, in.first, 4, in.second,
				2, c      , 8, state      , 3, c->CC(), 8, b, -1)));
      output.push_back(new_ptr((Tree2toNDiagram(5), in.first, gamma, b->CC(), gamma, in.second, 1, in.first, 4, in.second,
				5, b->CC(), 8, state      , 8, c->CC(), 2, b, -2)));
      output.push_back(new_ptr((Tree2toNDiagram(5), in.first, gamma, c      , gamma, in.second, 1, in.first, 4, in.second,
				5, c      , 8, state      , 2, c->CC(), 8, b, -3)));
      output.push_back(new_ptr((Tree2toNDiagram(5), in.first, gamma, b      , gamma, in.second, 1, in.first, 4, in.second,
				2, b->CC(), 8, state      , 8, c->CC(), 3, b, -4)));
      output.push_back(new_ptr((Tree2toNDiagram(5), in.first, gamma, b->CC(), gamma, in.second, 1, in.first, 4, in.second,
				2, b      , 8, state->CC(), 3, b->CC(), 8, c, -1)));
      output.push_back(new_ptr((Tree2toNDiagram(5), in.first, gamma, c->CC(), gamma, in.second, 1, in.first, 4, in.second,
				5, c->CC(), 8, state->CC(), 8, b->CC(), 2, c, -2)));
      output.push_back(new_ptr((Tree2toNDiagram(5), in.first, gamma, b      , gamma, in.second, 1, in.first, 4, in.second,
				5, b      , 8, state->CC(), 2, b->CC(), 8, c, -3)));
      output.push_back(new_ptr((Tree2toNDiagram(5), in.first, gamma, c      , gamma, in.second, 1, in.first, 4, in.second,
				2, c->CC(), 8, state->CC(), 8, b->CC(), 3, c, -4)));
    }
  }
  return output;
}

double GammaGamma2BC1S0QQbarAmplitude::generateKinematics(const double * r,
							  const Energy2 & scale, 
							  vector<Lorentz5Momentum> & momenta,
							  const tcPDVector & partons) {
  double a1= partons[1]->mass()/(partons[2]->mass()+partons[1]->mass());
  Energy ecm=sqrt(scale);
  vector<Energy> mOut = {a1*partons[0]->mass(),(1.-a1)*partons[0]->mass(),partons[0]->mass()};
  double jac = TwoToThreePhaseSpace::twoToThreeFS(ecm,mOut,r,momenta[1],
  						  momenta[2],momenta[0],1.,r[4]);
  jac /= 0.5*pow(Constants::twopi,4);
  return jac;
}

ProductionMatrixElement
GammaGamma2BC1S0QQbarAmplitude::helicityAmplitude(const Energy2 & scale,
						  const vector<VectorWaveFunction> & v1,
						  const vector<VectorWaveFunction> & v2,
						  const vector<Lorentz5Momentum> & out,
						  const vector<SpinorWaveFunction> v4,
						  const vector<SpinorBarWaveFunction> ubar5,
						  double & output, DVector & dweights) const {
  // invariants
  Energy M  = out[0].mass();
  Energy2 M2 = sqr(M);
  double a1 = out[1].mass()/M, a2 = out[2].mass()/M;
  double a12(sqr(a1)),a22(sqr(a2));
  double eb(-1./3.),ec(2./3.);
  if(abs(v4[0].particle()->id())!=4) swap(eb,ec);
  // storage of the matrix element
  vector<PDT::Spin> spins = {PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half};
  vector<unsigned int> ihMax(4,0);
  ProductionMatrixElement me = bookME(ihMax,v1.size(),v2.size(),spins);
  dweights.resize(4,0.);
  Energy2 dot1 = v1[0].momentum()*v2[0].momentum();
  Energy2 dot2 = v1[0].momentum()*out[1];
  Energy2 dot3 = v1[0].momentum()*out[2];
  Energy2 dot6 = out[0]*out[1];
  Energy2 dot7 = out[0]*out[2];
  vector<Complex> diag(20,0.);
  for(unsigned int ih1A=0;ih1A<ihMax[0];++ih1A) {
    for(unsigned int ih1B=0;ih1B<ihMax[1];++ih1B) {
      unsigned int ih1 = 2*ih1A+ih1B;
      complex<Energy> dot5 = v2[0].momentum()*v1[ih1].wave();
      complex<Energy> dot8 = out[1]*v1[ih1].wave();
      complex<Energy> dot10 = out[2]*v1[ih1].wave();
      for(unsigned int ih2A=0;ih2A<ihMax[2];++ih2A) {
  	for(unsigned int ih2B=0;ih2B<ihMax[3];++ih2B) {
  	  unsigned int ih2 = 2*ih2A+ih2B;
	  complex<Energy> dot4 = v1[0].momentum()*v2[ih2].wave();
	  complex<Energy> dot9 = out[1]*v2[ih2].wave();
	  complex<Energy> dot11 = out[2]*v2[ih2].wave();
	  Complex dot12 = v1[ih1].wave()*v2[ih2].wave();
  	  for(unsigned int ih4=0;ih4<2;++ih4) {
  	    for(unsigned int ih5=0;ih5<2;++ih5) {
	      complex<Energy> dot13=v4[ih4].dimensionedWave().pseudoScalar(ubar5[ih5].dimensionedWave());
	      LorentzPolarizationVectorE vec1 = v4[ih4].dimensionedWave().generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	      LorentzPolarizationVectorE vec2 = v4[ih4].dimensionedWave().slash(v2[ih2].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	      LorentzPolarizationVectorE vec3 = v4[ih4].dimensionedWave().slash(v1[ih1].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	      LorentzVector<complex<Energy2> > vec4 = v4[ih4].dimensionedWave().slash(v1[0].momentum()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	      LorentzVector<complex<Energy2> > vec5 = v4[ih4].dimensionedWave().slash(v2[0].momentum()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	      LorentzPolarizationVectorE vec6 = v4[ih4].dimensionedWave().slash(v2[ih2].wave()).slash(v1[ih1].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	      LorentzVector<complex<Energy2> > vec7 = v4[ih4].dimensionedWave().slash(v2[ih2].wave()).slash(v1[0].momentum()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	      LorentzVector<complex<Energy2> > vec8 = v4[ih4].dimensionedWave().slash(v2[ih2].wave()).slash(v2[0].momentum()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	      LorentzVector<complex<Energy2> > vec9 = v4[ih4].dimensionedWave().slash(v1[0].momentum()).slash(v2[ih2].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	      LorentzVector<complex<Energy2> > vec10 = v4[ih4].dimensionedWave().slash(v2[0].momentum()).slash(v2[ih2].wave()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	      LorentzVector<complex<Energy2> > vec11 = v4[ih4].dimensionedWave().slash(v1[ih1].wave()).slash(v1[0].momentum()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	      LorentzVector<complex<Energy2> > vec12 = v4[ih4].dimensionedWave().slash(v1[ih1].wave()).slash(v2[0].momentum()).generalCurrent(ubar5[ih5].dimensionedWave(),-1,1);
	      LorentzVector<complex<Energy2> > vec13 = v4[ih4].dimensionedWave().slash(v2[ih2].wave()).slash(v1[ih1].wave()).slash(v2[0].momentum()).generalCurrent(ubar5[ih5].dimensionedWave(),1,-1);
	      complex<Energy2> dot14 = vec1*v1[0].momentum();
	      complex<Energy2> dot15 = vec1*v2[0].momentum();
	      complex<Energy> dot16 = vec1*v1[ih1].wave();
	      complex<Energy> dot17 = vec1*v2[ih2].wave();
	      complex<Energy2> dot18 = vec10*v1[ih1].wave();
	      complex<Energy2> dot19 = vec11*v2[ih2].wave();
	      complex<Energy3> dot20 = vec12*v1[0].momentum();
	      complex<Energy2> dot21 = vec12*v2[ih2].wave();
	      complex<Energy3> dot22 = vec13*v1[0].momentum();
	      complex<Energy2> dot23 = vec2*v1[0].momentum();
	      complex<Energy2> dot24 = vec2*v2[0].momentum();
	      complex<Energy> dot25 = vec2*v1[ih1].wave();
	      complex<Energy2> dot26 = vec3*v1[0].momentum();
	      complex<Energy2> dot27 = vec3*v2[0].momentum();
	      complex<Energy2> dot28 = vec4*v1[ih1].wave();
	      complex<Energy2> dot29 = vec4*v2[ih2].wave();
	      complex<Energy3> dot30 = vec5*v1[0].momentum();
	      complex<Energy2> dot31 = vec5*v1[ih1].wave();
	      complex<Energy2> dot32 = vec5*v2[ih2].wave();
	      complex<Energy2> dot33 = vec6*v1[0].momentum();
	      complex<Energy2> dot34 = vec6*v2[0].momentum();
	      complex<Energy2> dot35 = vec7*v1[ih1].wave();
	      complex<Energy3> dot36 = vec8*v1[0].momentum();
	      complex<Energy2> dot37 = vec8*v1[ih1].wave();
	      complex<Energy2> dot38 = vec9*v1[ih1].wave();
	      diag[0]=eb*ec*M2*(2.*dot3*dot34-2.*(dot1-dot2)*(dot34+dot35+dot37)+2.*a2*dot20*dot4-4.*dot16*dot3*dot4+4.*a1*dot16*dot3*dot4+
				4.*dot1*dot17*dot5-4.*a1*dot1*dot17*dot5-4.*dot17*dot2*dot5+4.*a1*dot17*dot2*dot5-4.*dot17*dot3*dot5+
				2.*a1*dot36*dot5-4.*dot14*dot4*dot5+8.*a1*dot14*dot4*dot5-4.*a12*dot14*dot4*dot5+4.*a1*dot15*dot4*dot5-
				4.*a12*dot15*dot4*dot5+2.*dot18*dot6-2.*a1*dot18*dot6+2.*dot35*dot6-2.*a1*dot35*dot6+2.*dot37*dot6-
				2.*a1*dot37*dot6+2.*dot38*dot6-2.*a1*dot38*dot6-4.*dot16*dot4*dot6+4.*a1*dot16*dot4*dot6+2.*dot34*dot7+
				2.*dot35*dot7+2.*dot37*dot7-4.*dot17*dot5*dot7+4.*a1*dot17*dot5*dot7+4.*a1*dot1*dot17*dot8-4.*a1*dot17*dot2*dot8-
				2.*a1*dot36*dot8-4.*a1*dot14*dot4*dot8+4.*a12*dot14*dot4*dot8-4.*a1*dot15*dot4*dot8+4.*a12*dot15*dot4*dot8-
				4.*a1*dot17*dot7*dot8+2.*a1*dot20*dot9+4.*dot16*dot3*dot9-4.*a1*dot16*dot3*dot9-8.*a1*dot14*dot5*dot9+
				4.*a12*dot14*dot5*dot9+4.*a12*dot15*dot5*dot9+4.*a1*dot14*dot8*dot9-4.*a12*dot14*dot8*dot9-4.*a12*dot15*dot8*dot9-3.*dot22*M+
				2.*(dot1*dot25+dot12*dot30-dot27*dot4+4.*a1*dot27*dot4-2.*a12*dot27*dot4-2.*dot28*dot4+4.*a1*dot28*dot4-
				    2.*a12*dot28*dot4-2.*dot31*dot4+4.*a1*dot31*dot4-2.*a12*dot31*dot4+2.*dot23*dot5-a1*dot24*dot5+
				    2.*dot13*dot4*dot5-6.*a1*dot13*dot4*dot5+2.*a12*dot13*dot4*dot5-dot23*dot8+a1*dot23*dot8-dot24*dot8+
				    2.*a1*dot24*dot8-dot29*dot8+a1*dot29*dot8-dot32*dot8+a1*dot32*dot8+2.*dot13*dot4*dot8-
				    4.*a1*dot13*dot4*dot8+2.*a12*dot13*dot4*dot8+(-1.+2.*a1)*a2*dot26*(dot4-dot9)-
				    2.*a2*(-(a2*dot27)-dot28-dot31+2.*dot13*dot5+a1*(dot28+dot31-dot13*(dot5+dot8)))*dot9)*M-
				2.*a2*dot11*(dot20-2.*dot16*dot3-2.*a2*dot14*dot5-2.*a1*dot14*dot8-(dot26+2.*(dot27+dot28+dot31-2.*dot13*dot5))*M+
					     2.*a1*(dot15*(dot5-dot8)+(dot26+dot27+dot28+dot31-dot13*(dot5+dot8))*M))-
				a2*(-2.*(-1.+2.*a1)*dot12*(dot14+dot15)-3.*a2*dot18+dot19+dot21-2.*dot34-4.*dot35-4.*dot37-
				    3.*dot38+4.*dot16*dot4+2.*dot17*dot5+
				    a1*(-2.*dot19-2.*dot21+dot33+dot34+2.*dot35+2.*dot37+3.*dot38-2.*dot16*dot4-2.*dot17*dot5+4.*dot17*dot8))*M2+
				2.*dot10*(dot36+2.*dot14*dot4-2.*dot14*dot9-(dot23+dot24+dot29+dot32-2.*dot13*dot4)*M+
					  2.*a2*dot11*(-(a2*dot14)+a1*(dot15+dot13*M))+
					  2.*a12*((dot14+dot15)*(dot4-dot9)+dot13*(dot4-dot9)*M+dot17*M2)+
					  a1*(2.*dot1*dot17-dot36-2.*(2.*dot14+dot15)*dot4+4.*dot14*dot9+
					      (dot23+2.*dot24+dot29+dot32-4.*dot13*dot4+2.*dot13*dot9)*M-
					      2.*dot17*(dot2+dot7+M2))))
		/(4.*a1*a2*(dot1-dot2-dot3)*(-dot1+dot2+dot3+dot6+dot7+M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	      diag[1]=-0.25*sqr(ec)*M2*(2.*dot10*dot36-2.*dot1*(dot34+dot35+dot37-2.*dot17*dot5)-4.*dot10*dot14*dot9-dot22*M+
					4.*a12*dot17*(dot10-dot5+dot8)*M2+
					2.*(dot2*(dot34+dot35+dot37-2.*dot17*dot5)+(dot34+dot35+dot37-2.*dot17*dot5)*dot7+
					    dot3*(dot34-2.*dot17*dot5+2.*dot16*dot9)+dot23*dot5*M+
					    (dot27+dot28+dot31-2.*dot13*dot5)*dot9*M+(dot34+dot35+dot37-2.*dot17*dot5)*M2)+
					2.*a1*((dot10-dot5+dot8)*(2.*dot1*dot17-dot36-2.*dot17*(dot2+dot7)+2.*dot14*dot9)-
					       (dot10-dot5+dot8)*(dot24-2.*dot13*dot9)*M-(dot34+dot35+dot37+2.*dot17*(dot10-2.*dot5+dot8))*M2))
		/(a1*a2*(dot1-dot2-dot3)*(dot1-dot2-dot7-a2*M2)*(dot7+a2*M2));
	      diag[2]=-0.25*eb*ec*M2*(2.*(-(dot18*dot3)-dot3*dot38-2.*dot17*dot3*dot5+2.*a1*dot17*dot3*dot5+dot36*dot5-a1*dot36*dot5+
					  2.*dot15*dot4*dot5-2.*a1*dot15*dot4*dot5-dot18*dot6+dot34*dot6-dot38*dot6-2.*dot17*dot5*dot6+
					  2.*a1*dot17*dot5*dot6+2.*dot12*(dot15*dot3+dot14*(dot3+dot6))+
					  a1*(dot36+2.*dot15*dot4-2.*dot17*(dot3+dot6))*dot8+2.*dot11*(dot16*dot3+dot14*(dot5-a1*dot5+a1*dot8)))+
				      (dot22-2.*dot12*dot30-2.*dot27*dot4+2.*a1*dot24*(-dot5+dot8)+
				       2.*dot11*(dot27+dot28+dot31-2.*dot13*(dot5+a1*dot5-a1*dot8)))*M+
				      2.*dot1*(-2.*dot11*dot16+dot18+dot38-dot25*M-2.*dot12*(dot14+dot15-dot13*M))+
				      2.*a1*(2.*dot12*dot14-dot18+dot34-dot38-2.*a2*dot17*dot5-2.*a1*dot17*dot8)*M2+
				      2.*dot10*(2.*dot1*dot17-dot36-2.*dot15*dot4+2.*dot11*(-(a2*dot14)+a1*dot13*M)+
						a1*(dot36+2.*dot15*dot4+dot24*M-2.*dot17*(dot3+dot6+a1*M2))))
		/(a1*(dot1-dot2-dot3)*(-dot1+dot3+dot6+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	      diag[3]=sqr(ec)*M2*(2.*(dot11*dot20+2.*dot12*(dot14+dot15)*dot3-dot3*(dot18+dot38)+
				      dot1*(-2.*dot12*(dot14+dot15)-2.*dot11*dot16+2.*dot10*dot17+dot18+dot38)-
				      dot20*dot4+2.*dot16*dot3*dot4-dot10*(dot36+2.*(dot14+dot15)*dot4)-2.*dot17*dot3*dot5+
				      dot36*dot5+2.*dot14*dot4*dot5+2.*dot15*dot4*dot5-
				      (2.*dot12*dot15+dot18+dot19+dot21-dot33-2.*dot34+dot38-2.*dot16*dot4+2.*dot17*dot5)*dot6)+
				  (3.*dot22+4.*dot1*(dot12*dot13-dot25)-4.*dot12*dot30+
				   2.*(dot11*(dot26+dot27+dot28+dot31)+(-dot27+dot28+dot31)*dot4+(dot26+dot27+dot28+dot31)*dot9)+
				   2.*dot5*(dot24+dot29+dot32-2.*dot13*(dot11+dot4+dot9)))*M+
				  (-2.*dot12*(dot14+dot15)+dot18-dot19-dot21+dot33+dot34+dot35+dot37+dot38)*M2-
				  a1*(-4.*dot1*dot12*dot14-4.*dot1*dot12*dot15-4.*dot1*dot11*dot16+2.*dot1*dot18+
				      4.*dot12*dot14*dot2+4.*dot12*dot15*dot2+4.*dot11*dot16*dot2-2.*dot18*dot2+
				      2.*dot11*dot20+4.*dot12*dot14*dot3+4.*dot12*dot15*dot3-2.*dot18*dot3+2.*dot1*dot38-
				      2.*dot2*dot38-2.*dot3*dot38-2.*dot20*dot4+4.*dot16*dot3*dot4-4.*dot11*dot14*dot5-
				      4.*dot11*dot15*dot5-4.*dot17*dot3*dot5+2.*dot36*dot5+8.*dot14*dot4*dot5+8.*dot15*dot4*dot5-
				      4.*dot12*dot14*dot6-4.*dot12*dot15*dot6-2.*dot19*dot6-2.*dot21*dot6+2.*dot33*dot6+
				      2.*dot34*dot6+4.*dot16*dot4*dot6-4.*dot17*dot5*dot6+4.*dot11*dot15*dot8+4.*dot17*dot3*dot8-
				      2.*dot36*dot8-4.*dot14*dot4*dot8-8.*dot15*dot4*dot8+4.*dot17*dot6*dot8+2.*dot20*dot9-
				      4.*dot16*dot3*dot9-4.*dot14*dot5*dot9+
				      2.*(dot25*(-dot1+dot2+dot3)+(3.*dot23+2.*dot24+3.*(dot29+dot32))*dot5+
					  dot11*(2.*dot26+dot27+dot28+dot31-2.*dot13*dot5)-(2.*dot23+dot24+2.*(dot29+dot32))*dot8+
					  dot4*(dot27+dot28+dot31-6.*dot13*dot5+2.*dot13*dot8)+2.*dot26*dot9+(dot27+dot28+dot31-2.*dot13*dot5)*dot9)*M+
				      2.*dot10*(2.*dot11*(dot14+dot15)-dot36-4.*(dot14+dot15)*dot4+2.*dot17*(dot1-dot2+dot6)+
						2.*dot14*dot9-(2.*dot23+dot24+2.*(dot29+dot32-dot13*dot4))*M)+
				      (-2.*dot12*(5.*dot14+3.*dot15)+6.*dot18-3.*dot19-3.*dot21+3.*dot33+
				       dot34+4.*dot35+4.*dot37+6.*dot38-2.*dot16*dot4+4.*dot17*dot5)*M2)+
				  a12*(4.*(dot14+dot15)*(dot10-dot5+dot8)*(dot11-dot4+dot9)-
				       4.*(dot10-dot5+dot8)*(dot23+dot24+dot29+dot32-dot13*(dot11+dot4+dot9))*M+
				       (-6.*dot12*(dot14+dot15)+4.*dot18-3.*dot19-3.*dot21+3.*dot33+3.*dot34+
					4.*(dot35+dot37+dot38)-2.*dot16*dot4-4.*dot17*(dot10-dot5+dot8))*M2))
		/(4.*a1*a2*(dot1-dot2-dot3)*(-dot7-a2*M2)*(-dot1+a1*(dot6+dot7+M2)));
	      diag[4]=eb*ec*M2*(2.*(2.*dot12*dot15*dot2-dot11*dot20+2.*dot12*dot15*dot3+dot19*dot3+dot21*dot3-dot2*dot34-dot3*dot34+
				    2.*dot11*dot15*dot5+2.*dot17*dot3*dot5-dot36*dot5-2.*dot15*dot4*dot5+2.*dot12*dot15*dot6+dot19*dot6+
				    dot21*dot6-dot34*dot6+2.*dot17*dot5*dot6+dot10*(-2.*dot11*dot15+dot36+2.*dot15*dot4-2.*dot17*(dot3+dot6))-
				    2.*(dot11*dot15+dot17*(dot3+dot6))*dot8-
				    dot1*(2.*dot12*dot15-2.*dot11*dot16+dot19+dot21-dot34-2.*dot17*dot8))+
				(-3.*dot22+4.*dot1*(-(dot12*dot13)+dot25)+4.*dot12*dot30+4.*dot27*dot4+
				 2.*dot23*dot5-4.*dot13*dot4*dot5-2.*dot24*(dot10-dot5+dot8))*M+
				a12*(-4.*(dot14+dot15)*(dot10-dot5+dot8)*(dot11-dot4+dot9)+
				     4.*(dot10-dot5+dot8)*(dot23+dot24+dot29+dot32-dot13*(dot11+dot4+dot9))*M+
				     (4.*dot10*dot17-dot18+dot33+dot34-dot38+2.*dot16*dot4-6.*dot17*dot5+4.*dot17*dot8)*M2)+
				a1*(-4.*dot1*dot16*(dot11-dot4+dot9)+
				    2.*(-(dot20*dot4)-2.*dot17*dot3*dot5+dot36*dot5+2.*dot14*dot4*dot5+4.*dot15*dot4*dot5-
					2.*dot17*dot5*dot6+2.*dot12*dot14*dot7+2.*dot12*dot15*dot7+dot19*dot7+dot21*dot7-
					dot33*dot7-dot34*dot7-2.*dot16*dot4*(dot2+dot7)+2.*dot17*dot3*dot8-dot36*dot8-
					4.*dot15*dot4*dot8+2.*dot17*dot6*dot8+
					dot11*(2.*dot16*dot2+dot20-2.*(dot14+2.*dot15)*dot5+4.*dot15*dot8)+
					2.*dot16*dot2*dot9+(dot20-2.*(dot14+dot15)*dot5+2.*dot15*dot8)*dot9)+
				    2.*((dot23+dot29+dot32)*dot5-(2.*dot23+dot24+2.*(dot29+dot32))*dot8+
					dot11*(dot26+2.*dot13*(-dot5+dot8))-2.*dot13*dot5*dot9+dot26*(-dot4+dot9)+
					2.*dot13*dot8*(dot4+dot9))*M+
				    (2.*dot12*(dot14+3.*dot15)+3.*dot19+3.*dot21-dot33-3.*dot34-2.*dot16*dot4+4.*dot17*dot5-4.*dot17*dot8)*M2+
				    dot10*(-2.*dot36+4.*(dot14+2.*dot15)*(dot11-dot4)+4.*dot17*(dot3+dot6)+
					   4.*(dot14+dot15)*dot9+2.*(-2.*dot23-dot24-2.*(dot29+dot32)+2.*dot13*(dot11+dot4+dot9))*M-4.*dot17*M2)))
		/(4.*a1*a2*(dot1-dot2-dot3)*(-dot1+dot2+dot3+dot6+dot7+M2)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2));
	      diag[5]=-0.25*sqr(ec)*M2*(2.*(dot2*dot34+dot11*(dot20+2.*dot15*dot8)-(dot3+dot6)*(dot19+dot21-dot34-2.*dot17*dot8)-
					    a1*((dot19+dot21-dot33-dot34)*dot7+dot11*(2.*dot16*dot2+dot20+2.*dot15*dot8)-
						dot4*(dot20+2.*dot16*(dot2+dot7)+2.*dot15*dot8)+(2.*dot16*dot2+dot20+2.*dot15*dot8)*dot9))+
					(dot22-2.*dot27*dot4+2.*dot24*dot8+2.*a1*(dot26-2.*dot13*dot8)*(dot11-dot4+dot9))*M+
					2.*dot1*(-2.*a2*dot11*dot16+dot19+dot21-dot34-2.*dot17*dot8+2.*a1*dot16*(-dot4+dot9)-
						 dot25*M+2.*dot12*(dot15+dot13*M))+
					a1*((-3.+a1)*dot19+(-3.+a1)*dot21+dot33+3.*dot34+2.*dot16*dot4-a1*(dot33+dot34+2.*dot16*dot4)+4.*dot17*dot8)*M2+
					2.*dot12*(-2.*dot15*(dot2+dot3+dot6)-2.*a1*(dot14+dot15)*dot7-dot30*M+a1*(-(a2*dot14)+(-3.+a1)*dot15)*M2))
		/(a1*a2*dot2*(-dot1+dot2+dot3+dot6+dot7+M2)*(-dot7-a2*M2));
	      diag[6]=eb*ec*M2*(2.*(-(dot3*dot38)-a2*dot20*dot4+2.*dot16*dot3*dot4-2.*a1*dot16*dot3*dot4+2.*dot14*dot4*dot5-
				    2.*a1*dot14*dot4*dot5+dot11*(dot20-a1*dot20+2.*a1*dot16*dot3-2.*a2*dot14*dot5)+dot34*dot6+
				    a1*dot35*dot6+a1*dot37*dot6-dot38*dot6+a1*dot38*dot6+2.*dot16*dot4*dot6-2.*a1*dot16*dot4*dot6-
				    2.*dot17*dot5*dot6-dot18*(dot3+dot6-a1*dot6)+a1*(-dot20+2.*dot16*dot3+2.*dot14*dot5)*dot9)+
				(dot22-2.*dot23*dot5-2.*a1*(dot27+dot28+dot31-2.*dot13*dot5)*(dot11-dot4+dot9))*M+
				a1*(-dot18+2.*dot34+dot35+dot37-dot38+2.*dot16*dot4-4.*dot17*dot5)*M2+
				dot10*(-4.*dot17*(dot3+dot6)+4.*dot15*(-(a2*(dot11-dot4))+a1*dot9)-
				       2.*dot24*M+4.*a1*dot13*(dot11-dot4+dot9)*M-4.*a1*dot17*M2))
		/(4.*a1*dot3*(-dot1+dot2+dot3+dot6+dot7+M2)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2));
	      diag[7]=-0.25*sqr(ec)*M2*(-2.*dot11*dot20+4.*dot11*(dot14+dot15)*dot5+
					2.*(dot10*dot36+dot20*dot4-(dot36+2.*(dot14+dot15)*dot4)*dot5+
					    dot3*(dot38-2.*dot16*dot4+2.*dot17*dot5)+
					    (-dot34+dot38-2.*dot16*dot4+4.*dot17*dot5)*dot6+dot18*(dot3+dot6))+
					(-3.*dot22+2.*(dot1*dot25+dot12*dot30+dot27*dot4+(2.*dot23+dot24-2.*dot13*dot4)*dot5))*M+
					a12*(-4.*(dot14+dot15)*(dot10-dot5+dot8)*(dot11-dot4+dot9)+
					     4.*(dot26+dot27+dot28+dot31-dot13*(dot10+dot5+dot8))*(dot11-dot4+dot9)*M+
					     (-4.*dot12*(dot14+dot15)+4.*dot10*dot17+2.*dot18-2.*dot19-
					      2.*dot21+dot33+dot34+dot35+dot37+2.*dot38-2.*dot17*dot5+4.*dot17*dot8)*M2)+
					a1*(-4.*dot1*dot11*dot16+2.*dot1*dot18+4.*dot11*dot16*dot2-2.*dot18*dot2+2.*dot11*dot20-2.*dot18*dot3+
					    2.*dot1*dot38-2.*dot2*dot38-2.*dot3*dot38-2.*dot20*dot4+4.*dot16*dot3*dot4-8.*dot11*dot14*dot5-
					    8.*dot11*dot15*dot5-4.*dot17*dot3*dot5+2.*dot36*dot5+8.*dot14*dot4*dot5+8.*dot15*dot4*dot5+
					    4.*dot11*dot16*dot6-4.*dot18*dot6-2.*dot35*dot6-2.*dot37*dot6-4.*dot38*dot6+4.*dot16*dot4*dot6-
					    4.*dot17*dot5*dot6+4.*dot11*dot16*dot7-2.*dot18*dot7-2.*dot38*dot7+4.*dot11*dot15*dot8+
					    4.*dot17*dot3*dot8-2.*dot36*dot8-4.*dot15*dot4*dot8+4.*dot17*dot6*dot8+2.*dot20*dot9-
					    4.*dot16*dot3*dot9-8.*dot14*dot5*dot9-4.*dot15*dot5*dot9+
					    2.*(dot1*(2.*dot12*dot13-dot25)+2.*dot11*dot13*dot5-dot24*dot5-2.*dot13*dot4*dot5-
						2.*dot12*dot13*(dot2+dot3+dot6+dot7)+dot25*(dot2+dot3+dot6+dot7)+dot23*dot8+
						2.*dot24*dot8+dot29*dot8+dot32*dot8-2.*dot13*dot4*dot8+2.*dot13*dot5*dot9-
						(dot26+2.*(dot27+dot28+dot31))*(dot11-dot4+dot9))*M+
					    2.*(-2.*dot12*dot13+dot25)*pow<3,1>(M)+
					    (2.*dot12*(dot14+dot15)+4.*dot11*dot16-2.*dot18+dot19+dot21-2.*dot34-
					     dot35-dot37-2.*dot38-2.*dot16*dot4+6.*dot17*dot5)*M2+
					    2.*dot10*(2.*dot11*(dot14+dot15)+2.*dot1*dot17-dot36-2.*(dot14+dot15)*dot4+
						      2.*dot14*dot9+(dot23+2.*dot24+dot29+dot32-2.*dot13*dot4)*M-2.*dot17*(dot2+dot7+M2))))
		/(a1*a2*(-dot1+dot2+dot3+dot6+dot7+M2)*(-dot7-a2*M2)*(-dot1+a1*(dot6+dot7+M2)));
	      diag[8]=eb*ec*M2*(-2.*dot1*(dot33+dot34+2.*dot17*(-(a2*dot10)-a1*dot5-a2*dot8)+2.*dot16*dot9)+
				2.*(dot33*dot7-a2*dot10*(dot36+2.*dot17*(dot2+dot7)-2.*dot14*dot9)+
				    dot2*(dot33+dot34-2.*dot17*(a1*dot5+dot8-a1*dot8)+2.*dot16*dot9)-
				    (dot5-dot8)*(-2.*dot17*dot7+a1*(dot36+2.*dot17*dot7-2.*dot14*dot9)))+
				(dot22-2.*dot23*dot5-2.*a2*dot24*(dot10-dot5+dot8)+
				 2.*(dot26+2.*a2*dot13*(dot10-dot5+dot8))*dot9)*M+
				2.*a2*(dot33-2.*a2*dot17*(dot10-dot5+dot8))*M2)
		/(4.*a2*(dot1-dot2-dot3)*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot7-a2*M2));
	      diag[9]=-0.25*sqr(eb)*M2*(-4.*dot12*dot15*(dot2+dot3)+2.*dot1*(2.*dot12*dot15+dot19+dot21-dot34)+
					2.*(-((dot19+dot21)*dot3)+dot2*dot34+dot3*dot34-2.*dot17*dot3*dot5+2.*a1*dot17*dot3*dot5+dot36*dot5-
					    a1*dot36*dot5+2.*dot15*dot4*dot5-2.*a1*dot15*dot4*dot5-2.*dot11*(dot16*dot2-a2*dot14*dot5)-
					    (2.*dot12*dot15+dot19+dot21-dot34+2.*a2*dot17*dot5)*dot6+
					    (2.*dot17*(-dot1+dot3+dot6)+a1*(2.*dot11*dot14+dot36+2.*dot15*dot4-2.*dot17*(dot3+dot6)))*dot8)+
					(dot22+dot1*(4.*dot12*dot13-2.*dot25)-2.*dot12*dot30-2.*(dot11*dot26+dot27*dot4)+2.*a2*(2.*dot11*dot13+dot24)*(dot5-dot8))*M-
					2.*a1*(2.*dot12*dot15+dot19+dot21-dot34+2.*a2*dot17*(dot5-dot8))*M2-
					2.*a2*dot10*(dot36+2.*dot15*dot4+dot24*M+2.*dot11*(dot14+dot13*M)-2.*dot17*(dot3+dot6+a1*M2)))
		/(a1*a2*(dot1-dot2-dot3)*(dot6+a1*M2)*(-dot1+dot3+dot6+a1*M2));
	      diag[10]=-0.25*sqr(eb)*M2*(4.*dot1*dot11*dot16-4.*a1*dot1*dot11*dot16+2.*a1*dot1*dot18-4.*dot11*dot16*dot2+4.*a1*dot11*dot16*dot2-
					 2.*a1*dot18*dot2-2.*dot11*dot20+2.*a1*dot11*dot20-2.*a1*dot18*dot3+2.*dot3*dot33+2.*dot3*dot34+
					 2.*a1*dot1*dot38-2.*a1*dot2*dot38-2.*a1*dot3*dot38-2.*a1*dot20*dot4+4.*a1*dot16*dot3*dot4+
					 4.*dot11*dot14*dot5-8.*a1*dot11*dot14*dot5+4.*a12*dot11*dot14*dot5+4.*dot11*dot15*dot5-
					 8.*a1*dot11*dot15*dot5+4.*a12*dot11*dot15*dot5-4.*a1*dot17*dot3*dot5+2.*a1*dot36*dot5+
					 4.*a1*dot14*dot4*dot5-4.*a12*dot14*dot4*dot5+4.*a1*dot15*dot4*dot5-4.*a12*dot15*dot4*dot5-
					 2.*a1*dot18*dot6+2.*a1*dot33*dot6+2.*a1*dot34*dot6-2.*a1*dot38*dot6+4.*a1*dot16*dot4*dot6-
					 4.*a1*dot17*dot5*dot6-2.*dot33*dot7+2.*a1*dot33*dot7+2.*a1*dot34*dot7+2.*a1*dot35*dot7+
					 2.*a1*dot37*dot7-4.*dot17*dot5*dot7+4.*a1*dot11*dot14*dot8-4.*a12*dot11*dot14*dot8-
					 4.*dot11*dot15*dot8+8.*a1*dot11*dot15*dot8-4.*a12*dot11*dot15*dot8-4.*dot17*dot3*dot8+
					 4.*a1*dot17*dot3*dot8-2.*a1*dot36*dot8+4.*a12*dot14*dot4*dot8-4.*a1*dot15*dot4*dot8+
					 4.*a12*dot15*dot4*dot8+4.*dot17*dot7*dot8-4.*a1*dot17*dot7*dot8+2.*a1*dot20*dot9+
					 4.*dot16*dot3*dot9-4.*a1*dot16*dot3*dot9-8.*a1*dot14*dot5*dot9+4.*a12*dot14*dot5*dot9-
					 4.*a1*dot15*dot5*dot9+4.*a12*dot15*dot5*dot9+4.*a1*dot14*dot8*dot9-4.*a12*dot14*dot8*dot9+
					 4.*a1*dot15*dot8*dot9-4.*a12*dot15*dot8*dot9+
					 (-4.*a2*dot1*dot12*dot13-3.*dot22-2.*(-2.+a1)*dot1*dot25+2.*dot12*(2.*a2*dot13*(dot2+dot3)+dot30)+
					  2.*(-(a2*dot2*dot25)+dot11*dot27+dot11*dot28-dot25*dot3+dot11*dot31+dot27*dot4+2.*dot23*dot5+
					      dot24*dot5-2.*dot13*dot4*dot5-2.*dot11*dot13*dot8-dot24*dot8+(dot27+dot28+dot31-2.*dot13*dot8)*dot9+
					      2.*a12*(dot26+dot27+dot28+dot31-dot13*(dot5+dot8))*(dot11-dot4+dot9)-
					      a1*(3.*dot11*dot26-dot25*dot3-dot26*dot4+(dot23+2.*dot24+dot29+dot32)*dot5-dot24*dot8-
						  2.*dot4*(dot27+dot28+dot31-dot13*dot8)+4.*dot11*(dot27+dot28+dot31-dot13*(dot5+dot8))+
						  3.*dot26*dot9+4.*(dot27+dot28+dot31-dot13*(dot5+dot8))*dot9)))*M-
					 (2.*(dot33+2.*dot17*dot5)+a1*(-2.*(-1.+2.*a1)*dot12*(dot14+dot15)+dot19+dot21-
								       5.*dot33-3.*dot34-2.*(dot35+dot37+dot16*dot4+2.*dot17*dot5)+
								       a1*(3.*dot18-2.*dot19-2.*dot21+2.*dot33+2.*dot34+3.*dot35+
									   3.*dot37+3.*dot38-2.*dot16*dot4+4.*dot17*dot5))-4.*a22*dot17*dot8)*M2+
					 2.*a2*dot10*(dot36-dot24*M-2.*dot9*(dot14+dot13*M)-2.*a2*dot11*(dot14+dot15+dot13*M)-
						      2.*a1*(dot4-dot9)*(dot14+dot15+dot13*M)+2.*dot17*(-dot1+dot2+dot7+M2-a1*M2)))
		/(a1*a2*(dot1-dot2-dot3)*(dot6+a1*M2)*(dot1-a2*(dot6+dot7+M2)));
	      diag[11]=-0.25*eb*ec*M2*(-2.*dot19*dot2-2.*dot2*dot21-2.*dot19*dot3-2.*dot21*dot3+2.*dot3*dot33+2.*dot3*dot34+
				       4.*dot16*dot2*dot4-4.*a1*dot16*dot2*dot4-2.*a1*dot20*dot4+4.*dot16*dot3*dot4-2.*dot19*dot6-
				       2.*dot21*dot6+2.*dot33*dot6+2.*dot34*dot6+4.*dot16*dot4*dot6-
				       4.*dot12*(dot14+dot15)*(dot3+dot6)-4.*dot12*dot14*dot7+4.*a1*dot12*dot14*dot7-8.*dot12*dot15*dot7+
				       4.*a1*dot12*dot15*dot7-4.*dot19*dot7+2.*a1*dot19*dot7-4.*dot21*dot7+2.*a1*dot21*dot7+
				       2.*dot33*dot7-2.*a1*dot33*dot7+4.*dot34*dot7-2.*a1*dot34*dot7+4.*dot16*dot4*dot7-
				       4.*a1*dot16*dot4*dot7+4.*dot17*dot2*dot8-4.*a1*dot15*dot4*dot8+4.*dot17*dot7*dot8-
				       4.*dot16*dot2*dot9+4.*a1*dot16*dot2*dot9+2.*a1*dot20*dot9+4.*a1*dot15*dot8*dot9-
				       dot22*M+2.*(dot12*dot30+dot27*dot4-dot8*(dot24-2.*a2*dot13*(dot4-dot9))-a2*dot26*(dot4-dot9))*M+
				       2.*dot1*(dot19+dot21-dot33-dot34-2.*dot16*(-(a2*(dot11-dot4))+a1*dot9)+dot25*M+2.*dot12*(dot14+dot15-dot13*M))-
				       2.*a2*dot11*(2.*dot16*dot2+dot20-dot26*M+2.*dot8*(dot15+dot13*M))+
				       (-2.*dot12*(dot14+a1*dot14+3.*dot15-a1*dot15)-3.*dot19-3.*dot21+dot33+3.*dot34+
					2.*dot16*dot4+4.*dot17*dot8+a1*(dot19+dot21+dot33-dot34+2.*dot16*dot4-4.*dot17*dot8))*M2)
		/(a2*dot2*(-dot1+dot2+dot3+dot6+dot7+M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	      diag[12]=-0.25*sqr(eb)*M2*(-2.*(dot1-dot2)*(dot34+dot35+dot37)+2.*dot18*dot6+
					 2.*(dot20*(dot4-a1*dot4)+2.*dot17*(dot1-dot2)*dot5-a1*dot18*dot6+dot35*dot6-a1*dot35*dot6+dot37*dot6-
					     a1*dot37*dot6+dot38*dot6-a1*dot38*dot6-2.*dot16*dot4*dot6+2.*a1*dot16*dot4*dot6+dot34*dot7+
					     dot35*dot7+dot37*dot7-2.*dot17*dot5*dot7+dot3*(dot34-2.*dot17*dot5-2.*a2*dot16*(dot4-dot9))+
					     a1*dot20*dot9+2.*dot14*dot5*(-(a2*dot4)-a1*dot9))-
					 dot22*M+2.*(dot23*dot5+a2*(dot27+dot28+dot31-2.*dot13*dot5)*(dot4-dot9))*M-
					 2.*a2*dot11*(dot20-2.*dot16*dot3+(dot27+dot28+dot31)*M-2.*dot5*(dot14+dot13*M))-
					 a2*((-2.+a1)*dot18+a1*(dot35+dot37+dot38-2.*dot16*dot4)-
					     2.*(dot34+2.*dot35+2.*dot37+dot38-2.*dot16*dot4-2.*dot17*dot5))*M2+
					 2.*dot10*(-2.*(-(a2*dot15*(dot11-dot4))+dot17*(-dot1+dot2+dot7)+a1*dot15*dot9)+
						   (dot24+2.*a2*dot13*(dot11-dot4+dot9))*M-2.*a2*dot17*M2))
		/(a1*a2*dot3*(-dot1+dot2+dot3+dot6+dot7+M2)*(dot6+a1*M2));
	      diag[13]=sqr(eb)*M2*(-4.*dot1*dot12*dot14+4.*a1*dot1*dot12*dot14-4.*dot1*dot12*dot15+4.*a1*dot1*dot12*dot15-4.*dot1*dot11*dot16+
				   4.*a1*dot1*dot11*dot16-2.*a1*dot1*dot18-4.*a1*dot12*dot14*dot2-4.*a1*dot12*dot15*dot2+4.*dot11*dot16*dot2-
				   4.*a1*dot11*dot16*dot2+2.*a1*dot18*dot2+2.*dot11*dot20-2.*a1*dot11*dot20+4.*dot12*dot14*dot3-
				   4.*a1*dot12*dot14*dot3+4.*dot12*dot15*dot3-4.*a1*dot12*dot15*dot3+2.*a1*dot18*dot3+2.*dot1*dot33-
				   2.*dot3*dot33+2.*dot1*dot34-2.*dot3*dot34-2.*a1*dot1*dot38+2.*a1*dot2*dot38+2.*a1*dot3*dot38+
				   2.*a1*dot20*dot4-4.*a1*dot16*dot3*dot4+4.*a1*dot11*dot14*dot5-4.*a12*dot11*dot14*dot5+
				   4.*a1*dot11*dot15*dot5-4.*a12*dot11*dot15*dot5+4.*a1*dot17*dot3*dot5-2.*a1*dot36*dot5-
				   4.*a1*dot14*dot4*dot5+4.*a12*dot14*dot4*dot5-4.*a1*dot15*dot4*dot5+4.*a12*dot15*dot4*dot5+
				   4.*dot12*dot14*dot6-4.*a1*dot12*dot14*dot6+4.*dot12*dot15*dot6-4.*a1*dot12*dot15*dot6+
				   2.*a1*dot18*dot6-2.*dot33*dot6-2.*dot34*dot6+2.*a1*dot38*dot6-4.*a1*dot16*dot4*dot6+
				   4.*a1*dot17*dot5*dot6+4.*dot12*dot14*dot7-8.*a1*dot12*dot14*dot7+8.*dot12*dot15*dot7-
				   8.*a1*dot12*dot15*dot7+2.*a1*dot18*dot7+2.*dot19*dot7-2.*a1*dot19*dot7+2.*dot21*dot7-
				   2.*a1*dot21*dot7-2.*dot33*dot7+2.*a1*dot33*dot7-4.*dot34*dot7+2.*a1*dot34*dot7+2.*a1*dot38*dot7-
				   4.*a1*dot11*dot14*dot8+4.*a12*dot11*dot14*dot8+4.*dot11*dot15*dot8-8.*a1*dot11*dot15*dot8+
				   4.*a12*dot11*dot15*dot8-4.*dot1*dot17*dot8+4.*dot17*dot3*dot8-4.*a1*dot17*dot3*dot8+
				   2.*a1*dot36*dot8+4.*a1*dot14*dot4*dot8-4.*a12*dot14*dot4*dot8+8.*a1*dot15*dot4*dot8-
				   4.*a12*dot15*dot4*dot8+4.*dot17*dot6*dot8-4.*a1*dot17*dot6*dot8+4.*dot1*dot16*dot9-
				   2.*a1*dot20*dot9-4.*dot16*dot3*dot9+4.*a1*dot16*dot3*dot9+4.*a1*dot14*dot5*dot9-
				   4.*a12*dot14*dot5*dot9-4.*a12*dot15*dot5*dot9-4.*dot16*dot6*dot9+4.*a1*dot16*dot6*dot9-
				   4.*dot16*dot7*dot9+4.*a1*dot16*dot7*dot9-4.*a1*dot14*dot8*dot9+4.*a12*dot14*dot8*dot9-
				   4.*a1*dot15*dot8*dot9+4.*a12*dot15*dot8*dot9+
				   (4.*dot1*dot12*dot13+3.*dot22+2.*(-3.+a1)*dot1*dot25+
				    2.*(dot2*(dot25-a1*dot25)+dot11*dot26-2.*dot12*dot30+(dot24+dot29+dot32)*dot5+dot25*(dot3+dot6+dot7)-
					(-2.*dot11*dot13+2.*dot23+dot24+2.*(dot29+dot32))*dot8-
					dot4*(dot26+2.*dot27-2.*dot13*dot8)+(dot26+2.*dot13*dot8)*dot9-
					a1*((3.*dot23+2.*dot24+3.*(dot29+dot32))*dot5+dot25*(dot3+dot6+dot7)-
					    (4.*dot23+3.*dot24+4.*(dot29+dot32))*dot8-dot4*(dot26+4.*dot13*dot5-4.*dot13*dot8)+
					    dot11*(dot26-2.*dot13*dot5+4.*dot13*dot8)+(dot26-2.*dot13*dot5+4.*dot13*dot8)*dot9)+
					2.*a12*(dot5-dot8)*(dot23+dot24+dot29+dot32-dot13*(dot11+dot4+dot9))))*M+2.*a2*dot25*pow<3,1>(M)+
				   (-2.*a2*dot12*(-(a2*dot14)+(-3.+a1)*dot15)+dot19+dot21-dot33-3.*dot34+2.*dot16*dot4+
				    a1*(2.*dot18+(-2.+a1)*dot19+(-2.+a1)*dot21+2.*(dot34+dot38-2.*dot16*dot4)-
					a1*(dot33+dot34+2.*dot16*dot4-4.*dot17*dot5)+4.*a2*dot17*dot8)-4.*a2*dot16*dot9)*M2-
				   2.*a2*dot10*(dot36-2.*a2*(dot14+dot15)*(dot11-dot4)+2.*dot17*(-dot1+dot2+dot7)+
						2.*(-(a2*dot14)+a1*dot15)*dot9+
						(-2.*a2*dot11*dot13+2.*dot23+dot24+2.*(dot29+dot32)-
						 2.*(dot13*(dot4+dot9)+a1*(dot23+dot24+dot29+dot32-dot13*(dot4+dot9))))*M+2.*a2*dot17*M2))
		/(4.*a1*a2*(-dot1+dot2+dot3+dot6+dot7+M2)*(dot6+a1*M2)*(dot1-a2*(dot6+dot7+M2)));
	      diag[14]=-0.25*eb*ec*M2*(4.*dot11*dot16*dot2-2.*dot2*dot34-2.*dot3*dot34-2.*dot34*dot6+
				       2.*(dot19+dot21)*(dot3+dot6)+4.*dot12*dot15*(dot2+dot3+dot6)-4.*dot17*(dot3+dot6)*dot8+
				       (dot22-2.*(dot11*dot26+dot12*dot30+dot27*dot4)+2.*(2.*dot11*dot13+dot24)*dot8)*M-
				       2.*dot1*(dot19+dot21-dot34-2.*dot17*dot8+dot25*M+2.*dot12*(dot15-dot13*M))+
				       2.*a1*(2.*dot12*dot15+dot19+dot21-dot34-2.*dot17*dot8)*M2)
		/(dot2*(-dot1+dot3+dot6+a1*M2)*(-dot2+a1*(-dot1+dot2+dot3+dot6+a1*M2)));
	      diag[15]=-0.25*sqr(ec)*M2*(-2.*dot7*(2.*dot12*dot15+dot19+dot21-dot34-2.*dot17*dot8)+
					 (dot22+dot1*(4.*dot12*dot13-2.*dot25)+2.*dot2*dot25-2.*dot12*(2.*dot13*dot2+dot30)-
					  2.*(dot26+dot27)*dot4+2.*dot24*dot8+4.*dot13*dot8*(dot4-dot9)+2.*dot26*dot9)*M+
					 2.*(-2.*dot12*dot15-dot19-dot21+dot34+a12*(2.*dot12*(dot14+dot15)+dot19+dot21-dot33-dot34-2.*dot16*dot4)+
					     2.*dot17*dot8+a1*(2.*dot12*dot15+dot19+dot21-dot34-2.*dot17*dot8))*M2)/(a2*dot2*sqr(dot7+a2*M2));
	      diag[16]=-0.25*eb*ec*M2*(2.*dot1*(dot34+dot35+dot37-2.*dot17*(dot10+dot5))+4.*dot10*dot17*(dot2+dot7)-
				       2.*(dot2*(dot34+dot35+dot37-2.*dot17*dot5)+(dot34+dot35+dot37-2.*dot17*dot5)*dot7+
					   dot3*(dot34-2.*dot17*dot5+2.*dot16*dot9))+
				       (-dot22+2.*(dot10*dot24+dot23*dot5+(dot27+dot28+dot31-2.*dot13*(dot10+dot5))*dot9))*M-
				       2.*a2*(dot34+dot35+dot37-2.*dot17*(dot10+dot5))*M2)
		/(dot3*(-(a2*dot1)+dot2+dot7-a1*(dot2+dot3+dot7)+M2+(-2.+a1)*a1*M2)*(dot1-dot2-dot7-a2*M2));
	      diag[17]=-0.25*sqr(ec)*M2*(2.*dot7*(dot34-2.*dot17*dot5+2.*dot16*dot9)+
					 (dot22-2.*dot23*dot5+2.*dot25*(-dot1+dot2+dot7)+2.*dot26*dot9-2.*(dot5-dot8)*(dot24-2.*dot13*dot9))*M+
					 2.*a2*dot25*pow<3,1>(M)+2.*a2*(dot34-2.*dot17*dot5+2.*dot16*dot9)*M2)/(a2*(dot1-dot2-dot7-a2*M2)*sqr(dot7+a2*M2));
	      diag[18]=-0.25*sqr(eb)*M2*(-2.*(dot34+dot35+dot37-2.*dot17*(dot10+dot5))*dot6+
					 (dot22-2.*dot25*dot3+2.*(dot27+dot28+dot31)*(dot11-dot4)-
					  2.*dot10*(2.*dot11*dot13+dot24-2.*dot13*dot4)-2.*(2.*dot11*dot13+dot23-2.*dot13*dot4)*dot5)*M+
					 2.*(dot18+dot35+dot37+dot38-2.*dot16*dot4+a12*(dot18+dot35+dot37+dot38-2.*dot16*dot4)-
					     a1*(2.*dot18+dot34+3.*dot35+3.*dot37+2.*dot38-4.*dot16*dot4-2.*dot17*(dot10+dot5)))*M2)/(a1*dot3*sqr(dot6+a1*M2));
	      diag[19]=sqr(eb)*M2*(dot22*M+2.*dot11*(-2.*dot10*dot13+dot27+dot28+dot31)*M-2.*(dot10*dot24+dot25*dot3+dot27*dot4-dot24*dot5)*M+
				   4.*dot11*dot16*(dot6+a1*M2)-2.*(dot34+dot25*M)*(dot6+a1*M2)+
				   dot12*(4.*dot13*dot3*M-2.*dot30*M+4.*(dot15+dot13*M)*(dot6+a1*M2)))/(4.*a1*sqr(dot6+a1*M2)*(-dot1+dot3+dot6+a1*M2));
   	      Complex amp = std::accumulate(diag.begin(),diag.end(),Complex(0.));
 	      dweights[0]+=norm(diag[17]);
 	      dweights[1]+=norm(diag[18]);
 	      dweights[2]+=norm(diag[15]);
 	      dweights[3]+=norm(diag[19]);
   	      output += norm(amp);
	      if(v1.size()==2 && v2.size()==2)  me(2*ih1B,2*ih2B,0,ih4,ih5) = amp;
	      else if(v1.size()==2)             me(2*ih1B,ih2A,ih2B,0,ih4,ih5) = amp;
	      else if(v2.size()==2)             me(ih1A,ih1B,2*ih2B,0,ih4,ih5) = amp;
	      else                              me(ih1A,ih1B,ih2A,ih2B,0,ih4,ih5) = amp;
  	    }
  	  }
  	}
      }
    }
  }
  // final factors
  output *= 512./27.*O1_*scale/M2/M2/M*
    sqr(Constants::pi*generator()->standardModel()->alphaEM())*
    sqr(Constants::pi*generator()->standardModel()->alphaS(scale));
  // output = O1_*scale/M2/M2/M;
  return me;
}

Selector<const ColourLines *>
GammaGamma2BC1S0QQbarAmplitude::colourGeometries(unsigned int iopt, const cPDVector &,
						 tcDiagPtr diag) const {
  static ColourLines c1[4]={ColourLines("-6 -2 4 7"),
			    ColourLines("7 -2 -4 -6"),
			    ColourLines("-6 2 4 7"),
			    ColourLines("7 2 -4 -6")};
  static ColourLines c2[4]={ColourLines("-10 -3 8 11"),
			    ColourLines("11 -3 -8 -10"),
			    ColourLines("-10 3 8 11"),
			    ColourLines("11 3 -8 -10")};
  Selector<const ColourLines *> sel;
  if(iopt==0) {
    sel.insert(1.0, &c1[abs(diag->id())-1]);
  }
  else {
    sel.insert(1.0, &c2[abs(diag->id())-1]);
  }
  return sel;
}

double GammaGamma2BC1S0QQbarAmplitude::me2(const vector<VectorWaveFunction> & v1,
					   const vector<VectorWaveFunction> & v2,
					   const Energy2 & , const Energy2 & ,
					   const Energy2 & scale, 
					   const vector<Lorentz5Momentum> & momenta,
					   const cPDVector & partons,
					   DVector & dweights ) const {
  double output(0.);
  SpinorWaveFunction      q2w(momenta[1],partons[1],outgoing);
  SpinorBarWaveFunction   q3w(momenta[2],partons[2],outgoing);
  vector<SpinorWaveFunction> v4;
  vector<SpinorBarWaveFunction> ubar5;
  for(unsigned int ix=0;ix<2;++ix) {
    q2w.reset(ix);
    v4.push_back(q2w);
    q3w.reset(ix);
    ubar5.push_back(q3w);
  }
  helicityAmplitude(scale,v1,v2,momenta,v4,ubar5,output,dweights);
  return output;
}
