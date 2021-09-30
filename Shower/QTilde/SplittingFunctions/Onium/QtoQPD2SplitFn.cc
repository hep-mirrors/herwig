// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQPD2SplitFn class.
//

#include "QtoQPD2SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQPD2SplitFn::pOver_ = 100.;

IBPtr QtoQPD2SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQPD2SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQPD2SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*sqr(GeV*GeV2)) << n_ << theta_ << sTheta_ << cTheta_ << fixedAlphaS_;
}

void QtoQPD2SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*sqr(GeV*GeV2)) >> n_ >> theta_ >> sTheta_ >> cTheta_ >> fixedAlphaS_;
}

void QtoQPD2SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  sTheta_ = sin(theta_/180.*Constants::pi);
  cTheta_ = cos(theta_/180.*Constants::pi);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQPD2SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQPD2SplitFn("Herwig::QtoQPD2SplitFn", "HwOniumShower.so HwOniumParameters.so");

void QtoQPD2SplitFn::Init() {

  static ClassDocumentation<QtoQPD2SplitFn> documentation
    ("The QtoQPD2SplitFn class implements the branching q-> q' D2");

  static Parameter<QtoQPD2SplitFn,Energy7> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQPD2SplitFn::O1_, GeV*GeV2*GeV2*GeV2, 0.131*GeV*GeV2*GeV2*GeV2, 0.0*GeV*GeV2*GeV2*GeV2, 10.0*GeV*GeV2*GeV2*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQPD2SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQPD2SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<QtoQPD2SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQPD2SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Parameter<QtoQPD2SplitFn,double> interfaceTheta
    ("Theta",
     "Mixing angle between the 1D2 and 3D2 states (in degrees)",
     &QtoQPD2SplitFn::theta_, 34.4, 0., 360.,
     false, false, Interface::limited);

}

void QtoQPD2SplitFn::guesstz(Energy2 t1,unsigned int iopt,
			     const IdList &ids,
			     double enhance,bool ident,
			     double detune, 
			     Energy2 &t_main, double &z_main) {
  unsigned int pdfopt = iopt!=1 ? 0 : pdfFactor();
  double lower = integOverP(zLimits().first ,ids,pdfopt);
  double upper = integOverP(zLimits().second,ids,pdfopt);
  Energy M = ids[0]->mass()+ids[1]->mass();
  double a2 = ids[1]->mass()/M;
  double aS2 = fixedAlphaS_ < 0 ? sqr(alpha()->overestimateValue()) : sqr(fixedAlphaS_);
  Energy2 pre =  4./27.*aS2*O1_/pow(a2,4)/M/sqr(M*M);
  Energy2 c = (upper - lower) * colourFactor() * pre * enhance * detune;
  double r = UseRandom::rnd();
  assert(iopt<=2);
  if(iopt==1) {
    c *= pdfMax();
    //symmetry of FS gluon splitting
    if(ident) c*= 2;
  }
  else if(iopt==2) c*=-1.;
  // guess t
  t_main = t1/(1.-t1/c*log(r));
  // guess z
  z_main = invIntegOverP(lower + UseRandom::rnd()*(upper - lower),ids,pdfopt);
}


double QtoQPD2SplitFn::ratioP(const double z, const Energy2 t,
			      const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M;
  double r = sqr(M)/t;
  // 1D2 coefficients
  double W1D2[5];
  W1D2[0] = z*(30 + 2*pow(a1,6)*pow(-1 + z,6) - 4*pow(a1,5)*pow(-1 + z,5)*(1 + 2*z) + 
	       (-3 + z)*z*(24 + 5*(-3 + z)*z) + 6*pow(a1,4)*pow(-1 + z,4)*(-1 - 2*z + 4*pow(z,2)) - 
	       4*pow(a1,3)*(-2 + z)*pow(-1 + z,3)*(7 + z*(-4 + 5*z)) + 2*a1*(-1 + z)*(54 + z*(-102 + z*(73 + (-20 + z)*z))) + 
	       pow(a1,2)*sqr(1.-z)*(134 + z*(-188 + z*(109 + 5*(-6 + z)*z))))/
    (12.*pow(-1 + a1,2)*pow(a1,2)*pow(1.-a1*(1.-z),6));
  W1D2[1] =(5 - 4*pow(a1,6)*pow(-1 + z,4)*(-15 + 11*z) + z*(-93 + 98*z - 30*pow(z,2)) + 
	    8*pow(a1,5)*pow(-1 + z,3)*(31 + z*(-46 + 11*z)) - pow(a1,4)*sqr(1.-z)*(-397 + z*(771 + z*(-319 + 41*z))) - 
	    4*a1*(7 + z*(-53 + z*(52 + 7*(-3 + z)*z))) + 4*pow(a1,3)*(-1 + z)*(77 + z*(-162 + z*(54 + 5*z*(2 + z)))) - 
	    2*pow(a1,2)*(-61 + z*(169 + z*(-58 + z*(-54 + z + 11*pow(z,2))))))/
    (12.*pow(-1 + a1,2)*pow(a1,2)*pow(1.-a1*(1.-z),4));
  W1D2[2] =((224 + a1*(-133 + 2*a1*(-155 + 16*a1*(5 + a1))))/(-1 + a1) + (32*pow(-1 + a1,2))/pow(1.-a1*(1.-z),3) + 
	    (12*(-1 + a1)*(23 + 8*a1))/sqr(1.-a1*(1.-z)) + (468 - 48*a1*(1 + 6*a1))/(1.-a1*(1.-z)) - 
	    (a1*(57 + 2*a1*(-81 + 8*a1*(5 + a1)))*z)/(-1 + a1))/(6.*pow(a1,3));
  W1D2[3] = (16*(1 + a1*(4 - 13*z + a1*(-2*a1*(-3 + z)*(-1 + z) + pow(1 + z,2)))))/(3.*a1*(1.-a1*(1.-z)));
  W1D2[4] = -128*(1.-a1)*a1/3.;
  // 3D2 coefficients
  double W3D2[5];
  W3D2[0]=z*(90 + 6*pow(a1,6)*pow(-1 + z,6) - 12*pow(a1,5)*pow(-1 + z,5)*(-5 + 2*z) + 
	     18*pow(a1,4)*pow(-1 + z,4)*(15 + 2*(-5 + z)*z) + 3*(-3 + z)*z*(24 + 5*(-3 + z)*z) - 
	     12*pow(a1,3)*pow(-1 + z,3)*(-50 + z*(51 + z*(-15 + 2*z))) + 6*a1*(-1 + z)*(66 + z*(-126 + z*(87 + (-22 + z)*z))) + 
	     pow(a1,2)*sqr(1.-z)*(690 + z*(-996 + z*(459 + z*(-66 + 7*z)))))/
    (24.*pow(-1 + a1,2)*pow(a1,2)*pow(1.-a1*(1.-z),6));
  W3D2[1] = (-16*pow(a1,6)*pow(-1 + z,5) + 20*pow(a1,5)*pow(-1 + z,4)*(1 + 3*z) - 
     3*pow(a1,4)*pow(-1 + z,3)*(-75 + z*(26 + 57*z)) + 3*(5 + z*(-93 + 98*z - 30*pow(z,2))) + 
     4*pow(a1,3)*sqr(1.-z)*(95 + z*(-75 + z*(-69 + 41*z))) - 12*a1*(-2 + z*(-52 + z*(90 + z*(-47 + 7*z)))) - 
     2*pow(a1,2)*(-1 + z)*(-115 + z*(-84 + z*(306 + z*(-172 + 21*z)))))/
    (24.*pow(-1 + a1,2)*pow(a1,2)*pow(1.-a1*(1.-z),4));
  W3D2[2]=  (33 - 135*z + a1*(51 + 32*pow(a1,5)*pow(-1 + z,3) - 90*z + 99*pow(z,2) - 
			      16*pow(a1,4)*sqr(1.-z)*(4 + (-11 + z)*z) + 6*pow(a1,3)*(-1 + z)*(-39 + z*(103 + (-49 + z)*z)) + 
			      a1*(-191 + 37*z*(17 + 7*(-3 + z)*z)) + pow(a1,2)*(-31 + 3*z*(16 + z*(98 - z*(80 + 13*z))))))/
    (12.*(-1 + a1)*pow(a1,2)*pow(1.-a1*(1.-z),3));
  W3D2[3] = 4.*(27.-54./a1+16.*a1+60.*(1.-a1)/(a1*(1.-a1*(1.-z)))-z)/3.;
  W3D2[4] = -128.*(1.-a1)*a1/3.;
  // mixing coefficients
  double Wmixed[5];
  Wmixed[0] =  (z*(-30 + 2*pow(a1,5)*pow(-1 + z,5) - 2*pow(a1,4)*pow(-1 + z,4)*(-3 + 4*z) - 
		   (-3 + z)*z*(24 + 5*(-3 + z)*z) + 4*pow(a1,3)*pow(-1 + z,3)*(-6 + z*(-1 + 3*z)) - 
		   2*pow(a1,2)*sqr(1.-z)*(44 + z*(-44 + 5*z*(1 + z))) + a1*(-1 + z)*(-90 + z*(156 + z*(-91 + 3*z*(4 + z))))))/
    (2.*pow(-1 + a1,2)*pow(a1,2)*pow(1.-a1*(1.-z),5));
  Wmixed[1] = (-5 - 8*pow(a1,5)*(-2 + z)*pow(-1 + z,3) + z*(93 - 98*z + 30*pow(z,2)) + 
	       2*pow(a1,4)*sqr(1.-z)*(29 + z*(-50 + 33*z)) - pow(a1,3)*(-1 + z)*(-73 + z*(265 + z*(-311 + 95*z))) + 
	       a1*(5 - z*(112 + z*(-93 + 2*z*(4 + z)))) + pow(a1,2)*(31 + z*(-159 + z*(329 + z*(-233 + 40*z)))))/
    (2.*pow(-1 + a1,2)*pow(a1,2)*pow(1.-a1*(1.-z),3));
  Wmixed[2] = (-11 + 45*z + a1*(8 - 8*pow(a1,4)*pow(-1 + z,3) + 20*z - 78*pow(z,2) + 2*pow(a1,3)*sqr(1.-z)*(7 + 5*z) + 
          6*pow(a1,2)*(-1 + z)*(11 + z*(-13 + 14*z)) + a1*(47 + z*(-167 + 5*(47 - 15*z)*z))))/
     ((-1 + a1)*pow(a1,2)*sqr(1.-a1*(1.-z)));
  Wmixed[3] = 16*(1 - 2/a1 - 2*a1*(-1 + z) + 5*z);
  Wmixed[4] = 0.;
  double ratio = 0., rr=1.;
  int itest = (abs(ids[2]->id())%100000)/10000;
  double mix1 = itest==1 ? sTheta_ :  cTheta_;
  double mix2 = itest==1 ? cTheta_ : -sTheta_;
  double ors=sqrt(1./6.);
  for(unsigned int ix=0;ix<5;++ix) {
    ratio += rr*(sqr(mix1)*W1D2[ix]+sqr(mix2)*W3D2[ix]+ors*mix1*mix2*Wmixed[ix]); 
    rr*=r;
  }
  ratio /= pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQPD2SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQPD2SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQPD2SplitFn::matrixElement(const double z, const Energy2 t, 
					 const IdList & ids, const double phi, bool) {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M, a2=1-a1;
  double r = sqr(M)/t;
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t+sqr(M)*(sqr(a1)*z*(1.-z)-sqr(a2)*(1.-z)-z));
  double rz = sqrt(z);
  double r6 = sqrt(6.);
  double r23= sqrt(2./3.);
  int itest = (abs(ids[2]->id())%100000)/10000;
  double mix1 = itest==1 ? sTheta_ :  cTheta_;
  double mix2 = itest==1 ? cTheta_ : -sTheta_;
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin2)));
  (*kernal)(0,0,0) =
    mix1*r*pow(phase,2)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(-1 + sqr(a1)*(4*r - 2*sqr(1.-z)) + 2*pow(a1,3)*r*(-1 + sqr(z)) + a1*(3 - 2*r - 5*z + 2*sqr(z)))/
    ((-1 + a1)*a1*(1.-a1*(1.-z))*sqr(1.-z)*rz)
    -mix2*r*pow(phase,2)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(-3 + 4*sqr(a1)*sqr(1.-z) + 8*pow(a1,4)*r*sqr(1.-z) + a1*(-1 + 4*r + 5*z - 4*sqr(z)) - 4*pow(a1,3)*r*(3 - 4*z + sqr(z)))/(r6*(-1 + a1)*a1*(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(0,0,1) =
    -0.5*(pT*mix1*r*phase*(1 + z - 3*sqr(z) + pow(z,3) + 8*pow(a1,5)*sqr(r)*pow(-1 + z,3)*(1 + z) - r*(3 + z) + 2*pow(a1,4)*r*sqr(1.-z)*(-5 + 6*z - sqr(z) + 8*r*(2 + z)) - a1*(8*sqr(r) + r*(-19 + 24*z - 5*sqr(z)) + sqr(1.-z)*(4 - z + sqr(z))) + pow(a1,3)*(-1 + z)*(48*sqr(r) - 2*pow(-1 + z,3) + r*(-33 + 57*z - 31*sqr(z) + 7*pow(z,3))) + sqr(a1)*(-16*sqr(r)*(-2 + z) + pow(-1 + z,3)*(-5 + 3*z) + r*(-39 + 83*z - 57*sqr(z) + 13*pow(z,3)))))/((-1 + a1)*a1*M*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz)+
    pT*mix2*r*phase*(16*pow(a1,6)*sqr(r)*pow(-1 + z,4) + 8*pow(a1,5)*r*pow(-1 + z,3)*(-1 + 6*r + z) + 4*pow(a1,4)*r*sqr(1.-z)*(-3 + 2*z + sqr(z) + 4*r*(2 + z)) + 3*(1 + z - 3*sqr(z) + pow(z,3) - r*(3 + z)) + a1*(16*sqr(r) - 3*r*(-5 + 2*z + 3*sqr(z)) + sqr(1.-z)*(-2 - 15*z + 5*sqr(z))) - sqr(a1)*(-1 + z)*(-48*sqr(r) + sqr(1.-z)*(-5 + 9*z) + r*(1 - 10*z + 9*sqr(z))) + pow(a1,3)*(-1 + z)*(4*pow(-1 + z,3) + 16*sqr(r)*(-2 + 3*z) + r*(3 - 17*z + 17*sqr(z) - 3*pow(z,3))))/(2.*r6*(-1 + a1)*a1*M*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(0,0,2) =
    -((mix1*(12*a1*pow(r,3)*pow(1.-a1*(1.-z),5)*(-1 + a1 + a1*z) + r*pow(1.-a1*(1.-z),3)*(-1 + z)*(-1 + 4*sqr(a1)*(-1 + z) - 6*z + sqr(z) + a1*(5 - 8*z + 3*sqr(z))) + sqr(1.-z)*z*(-3 + pow(a1,3)*pow(-1 + z,3) + 3*z - sqr(z) - sqr(a1)*sqr(1.-z)*(1 + 2*z) + a1*(5 - 6*z + pow(z,3))) + sqr(r)*pow(1.-a1*(1.-z),3)*(-3*(1 + z) + 2*pow(a1,3)*sqr(1.-z)*(7 + z) + 2*a1*(10 - 11*z + sqr(z)) + sqr(a1)*(-31 + 51*z - 33*sqr(z) + 13*pow(z,3)))))/(r6*(-1 + a1)*a1*pow(1.-a1*(1.-z),3)*sqr(1.-z)*rz))+
    mix2*(8*a1*pow(r,3)*pow(1.-a1*(1.-z),5) + r*sqr(1.-a1*(1.-z))*(-1 + z)*(-1 + a1*sqr(1.-z) - 6*z + sqr(z)) - sqr(1.-z)*z*(3 + sqr(a1)*sqr(1.-z) - 3*z + sqr(z) - 2*a1*(2 - 3*z + sqr(z))) + sqr(r)*pow(1.-a1*(1.-z),3)*(4*sqr(a1)*sqr(1.-z) - 3*(1 + z) + a1*(-1 - 4*z + 5*sqr(z))))/(2.*(-1 + a1)*a1*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(0,0,3) =
    pT*mix1*r*(8*a1*sqr(r)*pow(1.-a1*(1.-z),3)*(-1 + a1 + a1*z) + (-1 + z)*z*(-3 + a1 + 2*sqr(a1)*sqr(1.-z) + z - a1*sqr(z)) + r*sqr(1.-a1*(1.-z))*(-1 - 3*z + 2*sqr(a1)*(-3 + 2*z + sqr(z)) + a1*(7 - 12*z + 5*sqr(z))))/(2.*(-1 + a1)*a1*M*phase*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz)+
    pT*mix2*r*(16*a1*sqr(r)*pow(1.-a1*(1.-z),3)*(-1 - 2*a1*(-1 + z) + sqr(a1)*(-1 + z)) + r*sqr(1.-a1*(1.-z))*(3 - 12*sqr(a1)*sqr(1.-z) + 8*pow(a1,3)*sqr(1.-z) + 9*z + a1*(1 + 6*z - 7*sqr(z))) - (-1 + z)*z*(3*(-3 + z) + a1*(9 - 10*z + sqr(z))))/(2.*r6*(-1 + a1)*a1*M*phase*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(0,0,4) = mix1*r*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(-2*a1*r + 4*sqr(a1)*r - z + a1*(-1 + z)*z + 2*pow(a1,3)*r*(-1 + sqr(z)))/((-1 + a1)*a1*pow(phase,2)*(1.-a1*(1.-z))*sqr(1.-z)*rz) +
    mix2*r*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(-16*sqr(a1)*r*(-1 + z) + 8*pow(a1,4)*r*sqr(1.-z) + 3*z + a1*(-4*r + z - sqr(z)) - 4*pow(a1,3)*r*(5 - 8*z + 3*sqr(z)))/(r6*(-1 + a1)*a1*pow(phase,2)*(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(0,1,0) =
    -pT*mix1*r*phase*(2*a1*sqr(r)*pow(1.-a1*(1.-z),3) + sqr(1.-z)*z + r*(-1 + z)*(-1 + 4*a1 + 2*pow(a1,3)*sqr(1.-z) + sqr(a1)*(-5 + 4*z + sqr(z))))/((-1 + a1)*a1*M*(1.-a1*(1.-z))*sqr(1.-z)*rz)
    -pT*mix2*r*phase*(-((3 + a1*(-3 + z))*sqr(1.-z)*z) + r*sqr(1.-a1*(1.-z))*(-1 + z)*(3 + a1 + 4*sqr(a1)*(-1 + z) - 9*a1*z) + 4*a1*sqr(r)*pow(1.-a1*(1.-z),3)*(1 + a1 + 2*sqr(a1)*(-1 + z) - 3*a1*z))/(r6*(-1 + a1)*a1*M*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(0,1,1) =
    -0.5*mix1*(8*a1*pow(r,3)*pow(1.-a1*(1.-z),5) + (2 + 2*a1*(-1 + z) - z)*pow(-1 + z,3)*z + r*sqr(1.-a1*(1.-z))*sqr(1.-z)*(-1 + 2*sqr(a1)*(-1 + z) + z + a1*(3 + 5*z)) + sqr(r)*pow(1.-a1*(1.-z),3)*(-1 + z)*(-3 + 10*sqr(a1)*(-1 + z) + a1*(13 + 5*z)))/((-1 + a1)*a1*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz)
    -0.5*mix2*(r*pow(1.-a1*(1.-z),3)*sqr(1.-z)*(3 + a1 + 4*sqr(a1)*(-1 + z) - 3*z - 13*a1*z) + 16*a1*pow(r,3)*pow(1.-a1*(1.-z),5)*(1 + sqr(a1)*(-1 + z) - 2*a1*z) - pow(-1 + z,3)*z*(6 - 6*sqr(a1)*(-1 + z) - 3*z + a1*(-12 + 9*z + sqr(z))) + sqr(r)*pow(1.-a1*(1.-z),3)*(9*(-1 + z) - 4*pow(a1,3)*pow(-1 + z,3) + 8*pow(a1,4)*pow(-1 + z,3) - 6*a1*(-1 + sqr(z)) + sqr(a1)*(7 - 21*z + 45*sqr(z) - 31*pow(z,3))))/(r6*(-1 + a1)*a1*pow(1.-a1*(1.-z),3)*sqr(1.-z)*rz);
  (*kernal)(0,1,2) =
    pT*mix1*r*((-3*r + sqr(1.-z))*(-1 + z) + 12*pow(a1,4)*sqr(r)*pow(-1 + z,3) + 2*pow(a1,3)*r*sqr(1.-z)*(-7 + 18*r + 7*z) + a1*(12*sqr(r) - (-5 + z)*sqr(1.-z) - 2*r*(10 - 11*z + sqr(z))) + sqr(a1)*(-1 + z)*(36*sqr(r) + 4*sqr(1.-z) + r*(-31 + 30*z + sqr(z))))/(r6*(-1 + a1)*a1*M*phase*(1.-a1*(1.-z))*sqr(1.-z)*rz)
    -0.5*pT*mix2*r*((-3*r + sqr(1.-z))*(-1 + z) + 8*pow(a1,5)*sqr(r)*pow(-1 + z,3)*(1 + z) + pow(a1,3)*r*(-1 + z)*(-9 + 48*r + 7*z - 3*sqr(z) + 5*pow(z,3)) + 4*pow(a1,4)*r*sqr(1.-z)*(-1 + sqr(z) + 4*r*(2 + z)) - a1*(8*sqr(r) + r*(5 - 6*z + sqr(z)) + sqr(1.-z)*(-2 - z + sqr(z))) + sqr(a1)*(-16*sqr(r)*(-2 + z) + pow(-1 + z,3)*(1 + 3*z) + r*(-3 + 5*z - 9*sqr(z) + 7*pow(z,3))))/((-1 + a1)*a1*M*phase*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(0,1,3) =
    mix1*r*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(1 + 8*sqr(a1)*r*(-1 + z) - z + a1*(-6 + 8*r + 6*z))/(2.*(-1 + a1)*a1*pow(phase,2)*sqr(1.-z)*rz)-
    0.5*mix2*r*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(3 + 16*pow(a1,4)*r*sqr(1.-z) - 3*z - 4*sqr(a1)*(-3 + z)*(-1 + 4*r + z) + 8*pow(a1,3)*(-1 + z)*(-1 + 6*r + z) + a1*(1 - 16*r - 6*z + 5*sqr(z)))/(r6*(-1 + a1)*a1*pow(phase,2)*(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(0,1,4) =
    -2*pT*mix1*sqr(r)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)/((-1 + a1)*pow(phase,3)*M*sqr(1.-z)*rz)+
    2*r23*(-1 + 2*a1)*pT*mix2*sqr(r)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)/((-1 + a1)*pow(phase,3)*M*sqr(1.-z)*rz);
  (*kernal)(1,0,0) =
    -2.*phase*sqr(phase)*double(pT/M)*mix1*sqr(r)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)/((-1 + a1)*sqr(1.-z)*rz)+
    2*r23*(-1 + 2*a1)*pow(phase,3)*pT*mix2*sqr(r)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)/((-1 + a1)*M*sqr(1.-z)*rz);
  (*kernal)(1,0,1) =
    -0.5*(mix1*r*pow(phase,2)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(1 + 8*sqr(a1)*r*(-1 + z) - z + a1*(-6 + 8*r + 6*z)))/((-1 + a1)*a1*sqr(1.-z)*rz)+
    mix2*r*pow(phase,2)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(3 + 16*pow(a1,4)*r*sqr(1.-z) - 3*z - 4*sqr(a1)*(-3 + z)*(-1 + 4*r + z) + 8*pow(a1,3)*(-1 + z)*(-1 + 6*r + z) + a1*(1 - 16*r - 6*z + 5*sqr(z)))/(2.*r6*(-1 + a1)*a1*(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(1,0,2) =
    pT*mix1*r*phase*((-3*r + sqr(1.-z))*(-1 + z) + 12*pow(a1,4)*sqr(r)*pow(-1 + z,3) + 2*pow(a1,3)*r*sqr(1.-z)*(-7 + 18*r + 7*z) + a1*(12*sqr(r) - (-5 + z)*sqr(1.-z) - 2*r*(10 - 11*z + sqr(z))) + sqr(a1)*(-1 + z)*(36*sqr(r) + 4*sqr(1.-z) + r*(-31 + 30*z + sqr(z))))/(r6*(-1 + a1)*a1*M*(1.-a1*(1.-z))*sqr(1.-z)*rz)-
    0.5*pT*mix2*r*phase*((-3*r + sqr(1.-z))*(-1 + z) + 8*pow(a1,5)*sqr(r)*pow(-1 + z,3)*(1 + z) + pow(a1,3)*r*(-1 + z)*(-9 + 48*r + 7*z - 3*sqr(z) + 5*pow(z,3)) + 4*pow(a1,4)*r*sqr(1.-z)*(-1 + sqr(z) + 4*r*(2 + z)) - a1*(8*sqr(r) + r*(5 - 6*z + sqr(z)) + sqr(1.-z)*(-2 - z + sqr(z))) + sqr(a1)*(-16*sqr(r)*(-2 + z) + pow(-1 + z,3)*(1 + 3*z) + r*(-3 + 5*z - 9*sqr(z) + 7*pow(z,3))))/((-1 + a1)*a1*M*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz)
;
  (*kernal)(1,0,3) =
    mix1*(8*a1*pow(r,3)*pow(1.-a1*(1.-z),5) + (2 + 2*a1*(-1 + z) - z)*pow(-1 + z,3)*z + r*sqr(1.-a1*(1.-z))*sqr(1.-z)*(-1 + 2*sqr(a1)*(-1 + z) + z + a1*(3 + 5*z)) + sqr(r)*pow(1.-a1*(1.-z),3)*(-1 + z)*(-3 + 10*sqr(a1)*(-1 + z) + a1*(13 + 5*z)))/(2.*(-1 + a1)*a1*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz)+
    mix2*(r*pow(1.-a1*(1.-z),3)*sqr(1.-z)*(3 + a1 + 4*sqr(a1)*(-1 + z) - 3*z - 13*a1*z) + 16*a1*pow(r,3)*pow(1.-a1*(1.-z),5)*(1 + sqr(a1)*(-1 + z) - 2*a1*z) - pow(-1 + z,3)*z*(6 - 6*sqr(a1)*(-1 + z) - 3*z + a1*(-12 + 9*z + sqr(z))) + sqr(r)*pow(1.-a1*(1.-z),3)*(9*(-1 + z) - 4*pow(a1,3)*pow(-1 + z,3) + 8*pow(a1,4)*pow(-1 + z,3) - 6*a1*(-1 + sqr(z)) + sqr(a1)*(7 - 21*z + 45*sqr(z) - 31*pow(z,3))))/(2.*r6*(-1 + a1)*a1*pow(1.-a1*(1.-z),3)*sqr(1.-z)*rz);
  (*kernal)(1,0,4) =
    -pT*mix1*r*(2*a1*sqr(r)*pow(1.-a1*(1.-z),3) + sqr(1.-z)*z + r*(-1 + z)*(-1 + 4*a1 + 2*pow(a1,3)*sqr(1.-z) + sqr(a1)*(-5 + 4*z + sqr(z))))/((-1 + a1)*a1*M*phase*(1.-a1*(1.-z))*sqr(1.-z)*rz)
    -pT*mix2*r*(-((3 + a1*(-3 + z))*sqr(1.-z)*z) + r*sqr(1.-a1*(1.-z))*(-1 + z)*(3 + a1 + 4*sqr(a1)*(-1 + z) - 9*a1*z) + 4*a1*sqr(r)*pow(1.-a1*(1.-z),3)*(1 + a1 + 2*sqr(a1)*(-1 + z) - 3*a1*z))/
    (r6*(-1 + a1)*a1*M*phase*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(1,1,0) =
    -mix1*r*pow(phase,2)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(-2*a1*r + 4*sqr(a1)*r - z + a1*(-1 + z)*z + 2*pow(a1,3)*r*(-1 + sqr(z)))/((-1 + a1)*a1*(1.-a1*(1.-z))*sqr(1.-z)*rz)
    -mix2*r*pow(phase,2)*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(-16*sqr(a1)*r*(-1 + z) + 8*pow(a1,4)*r*sqr(1.-z) + 3*z + a1*(-4*r + z - sqr(z)) - 4*pow(a1,3)*r*(5 - 8*z + 3*sqr(z)))/(r6*(-1 + a1)*a1*(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(1,1,1) =
    pT*mix1*r*phase*(8*a1*sqr(r)*pow(1.-a1*(1.-z),3)*(-1 + a1 + a1*z) + (-1 + z)*z*(-3 + a1 + 2*sqr(a1)*sqr(1.-z) + z - a1*sqr(z)) + r*sqr(1.-a1*(1.-z))*(-1 - 3*z + 2*sqr(a1)*(-3 + 2*z + sqr(z)) + a1*(7 - 12*z + 5*sqr(z))))/(2.*(-1 + a1)*a1*M*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz)+
    pT*mix2*r*phase*(16*a1*sqr(r)*pow(1.-a1*(1.-z),3)*(-1 - 2*a1*(-1 + z) + sqr(a1)*(-1 + z)) + r*sqr(1.-a1*(1.-z))*(3 - 12*sqr(a1)*sqr(1.-z) + 8*pow(a1,3)*sqr(1.-z) + 9*z + a1*(1 + 6*z - 7*sqr(z))) - (-1 + z)*z*(3*(-3 + z) + a1*(9 - 10*z + sqr(z))))/(2.*r6*(-1 + a1)*a1*M*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(1,1,2) =
    mix1*(12*a1*pow(r,3)*pow(1.-a1*(1.-z),5)*(-1 + a1 + a1*z) + r*pow(1.-a1*(1.-z),3)*(-1 + z)*(-1 + 4*sqr(a1)*(-1 + z) - 6*z + sqr(z) + a1*(5 - 8*z + 3*sqr(z))) + sqr(1.-z)*z*(-3 + pow(a1,3)*pow(-1 + z,3) + 3*z - sqr(z) - sqr(a1)*sqr(1.-z)*(1 + 2*z) + a1*(5 - 6*z + pow(z,3))) + sqr(r)*pow(1.-a1*(1.-z),3)*(-3*(1 + z) + 2*pow(a1,3)*sqr(1.-z)*(7 + z) + 2*a1*(10 - 11*z + sqr(z)) + sqr(a1)*(-31 + 51*z - 33*sqr(z) + 13*pow(z,3))))/(r6*(-1 + a1)*a1*pow(1.-a1*(1.-z),3)*sqr(1.-z)*rz)
    -0.5*mix2*(8*a1*pow(r,3)*pow(1.-a1*(1.-z),5) + r*sqr(1.-a1*(1.-z))*(-1 + z)*(-1 + a1*sqr(1.-z) - 6*z + sqr(z)) - sqr(1.-z)*z*(3 + sqr(a1)*sqr(1.-z) - 3*z + sqr(z) - 2*a1*(2 - 3*z + sqr(z))) + sqr(r)*pow(1.-a1*(1.-z),3)*(4*sqr(a1)*sqr(1.-z) - 3*(1 + z) + a1*(-1 - 4*z + 5*sqr(z))))/((-1 + a1)*a1*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(1,1,3) =
    -0.5*pT*mix1*r*(1 + z - 3*sqr(z) + pow(z,3) + 8*pow(a1,5)*sqr(r)*pow(-1 + z,3)*(1 + z) - r*(3 + z) + 2*pow(a1,4)*r*sqr(1.-z)*(-5 + 6*z - sqr(z) + 8*r*(2 + z)) - a1*(8*sqr(r) + r*(-19 + 24*z - 5*sqr(z)) + sqr(1.-z)*(4 - z + sqr(z))) + pow(a1,3)*(-1 + z)*(48*sqr(r) - 2*pow(-1 + z,3) + r*(-33 + 57*z - 31*sqr(z) + 7*pow(z,3))) + sqr(a1)*(-16*sqr(r)*(-2 + z) + pow(-1 + z,3)*(-5 + 3*z) + r*(-39 + 83*z - 57*sqr(z) + 13*pow(z,3))))/((-1 + a1)*a1*M*phase*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz)
    +pT*mix2*r*(16*pow(a1,6)*sqr(r)*pow(-1 + z,4) + 8*pow(a1,5)*r*pow(-1 + z,3)*(-1 + 6*r + z) + 4*pow(a1,4)*r*sqr(1.-z)*(-3 + 2*z + sqr(z) + 4*r*(2 + z)) + 3*(1 + z - 3*sqr(z) + pow(z,3) - r*(3 + z)) + a1*(16*sqr(r) - 3*r*(-5 + 2*z + 3*sqr(z)) + sqr(1.-z)*(-2 - 15*z + 5*sqr(z))) - sqr(a1)*(-1 + z)*(-48*sqr(r) + sqr(1.-z)*(-5 + 9*z) + r*(1 - 10*z + 9*sqr(z))) + pow(a1,3)*(-1 + z)*(4*pow(-1 + z,3) + 16*sqr(r)*(-2 + 3*z) + r*(3 - 17*z + 17*sqr(z) - 3*pow(z,3))))/(2.*r6*(-1 + a1)*a1*M*phase*sqr(1.-a1*(1.-z))*sqr(1.-z)*rz);
  (*kernal)(1,1,4) =
    -mix1*r*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(-1 + sqr(a1)*(4*r - 2*sqr(1.-z)) + 2*pow(a1,3)*r*(-1 + sqr(z)) + a1*(3 - 2*r - 5*z + 2*sqr(z)))/
    ((-1 + a1)*a1*pow(phase,2)*(1.-a1*(1.-z))*sqr(1.-z)*rz)
    +mix2*r*(r*sqr(1.-a1*(1.-z)) + (-1 + z)*z)*(-3 + 4*sqr(a1)*sqr(1.-z) + 8*pow(a1,4)*r*sqr(1.-z) + a1*(-1 + 4*r + 5*z - 4*sqr(z)) - 4*pow(a1,3)*r*(3 - 4*z + sqr(z)))/
    (r6*(-1 + a1)*a1*pow(phase,2)*(1.-a1*(1.-z))*sqr(1.-z)*rz);
  return kernal;
}
