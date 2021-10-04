// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQPP1SplitFn class.
//

#include "QtoQPP1SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"

using namespace Herwig;

const double QtoQPP1SplitFn::pOver_ = 2.;

IBPtr QtoQPP1SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQPP1SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQPP1SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*sqr(GeV2)) << oenum(state_) << n_ << theta_ << sTheta_ << cTheta_ << fixedAlphaS_;
}

void QtoQPP1SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*sqr(GeV2)) >> ienum(state_) >> n_ >> theta_ >> sTheta_ >> cTheta_ >> fixedAlphaS_;
}

void QtoQPP1SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletMEProduction<1>(state_,n_,0,1);
  double theta = params_->singletTripletMixing(n_,1);
  sTheta_ = sin(theta/180.*Constants::pi);
  cTheta_ = cos(theta/180.*Constants::pi);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQPP1SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQPP1SplitFn("Herwig::QtoQPP1SplitFn", "HwOniumParameters.so HwOniumShower.so");

void QtoQPP1SplitFn::Init() {

  static ClassDocumentation<QtoQPP1SplitFn> documentation
    ("The QtoQPP1SplitFn class implements the branching q-> q' P1");

  static Reference<QtoQPP1SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &QtoQPP1SplitFn::params_, false, false, true, false, false);
  
  static Parameter<QtoQPP1SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQPP1SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<QtoQPP1SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &QtoQPP1SplitFn::state_, ccbar, false, false);
  static SwitchOption interfaceStateccbar
    (interfaceState,
     "ccbar",
     "Charmonium state",
     ccbar);
  static SwitchOption interfaceStatebbbar
    (interfaceState,
     "bbbar",
     "Bottomonium state",
     bbbar);
  
  static Parameter<QtoQPP1SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQPP1SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQPP1SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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
  Energy2 pre = 32./243.*aS2*O1_/pow(a2,4)/M/sqr(M);
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

double QtoQPP1SplitFn::ratioP(const double z, const Energy2 t,
			      const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M;
  double r = sqr(M)/t;
  // 1P1 coefficients
  double W1P1[4];
  W1P1[0] = z*(6.-8.*z+3.*sqr(z)
	       +a1*(+2.*(-4.+7.*z-4.*sqr(z)+pow(z,3))
		    +sqr(1.-z)*a1*((-1.+3*sqr(z))+a1*(1.-z)*(2.*(1.+z) +a1*(1.-z)))))/
    (16.*sqr(a1)*pow(1.-a1*(1.-z),4));
  W1P1[1] = (3.-9.*z
	     +a1*(+2.*(-8.+9.*z+5.*sqr(z))
		  +a1*(55 - 123*z +   73*sqr(z) -   13*pow(z,3)
		       +2.*a1*(-37 + 101*z -   83*sqr(z) +   19*pow(z,3)
			       -4.*a1*sqr(1.-z)*(-4.+3*z)))))/(16.*sqr(a1)*pow(1 + a1*(-1 + z),2));
  W1P1[2] = -0.5*((-1 + a1)*(-1 +   2*a1*(-2 + 5*z) +   2*pow(a1,3)*   (2 - 3*z +     sqr(z)) -   sqr(a1)*   (-1 + 2*z +     sqr(z))))/(a1*(1.-a1*(1.-z)));
  W1P1[3] = 4*sqr(1.-a1)*a1;
  // 3P1 coefficients
  double W3P1[4];
  W3P1[0]=z*(6 + sqr(a1)*(-5 + z)*(-3 + z)*sqr(1.-z) - 2*pow(a1,3)*(-3 + z)*pow(-1 + z,3) +  pow(a1,4)*pow(-1 + z,4) + 
	     z*(-8 + 3*z) + 2*a1*(-1 + z)*(8 + (-7 + z)*z))/(8.*sqr(a1)*pow(1 + a1*(-1 + z),4));
  W3P1[1]=(3 - 9*z + a1*(4 - 2*pow(a1,3)*pow(-1 + z,3) + 10*(-1 + z)*z + 2*sqr(a1)*sqr(1.-z)*(3 + z) - a1*(-1 + z)*(-15 + z*(20 + z))))/
    (8.*sqr(a1)*sqr(1.-a1*(1.-z)));
  W3P1[2]=((-1 + a1)*(2 + a1*(5 - 11*z + a1*(-5 + 2*a1*(-1 + z) - (-8 + z)*z))))/(2.*a1*(1 + a1*(-1 + z)));
  W3P1[3]=4.*sqr(1.-a1)*a1;
  // mixing coefficients
  double Wmixed[4];
  Wmixed[0] = z*(-6.+(8.-3.*z)*z
		 +a1*(6 + z*(-8 + z + sqr(z))
		      +a1*sqr(1.-z)*((1.-2.*z)-a1*(1.-z))))/
    (4.*sqr(a1)*pow(1.-a1*(1.-z),3));
  Wmixed[1] = (-3.+9.*z
	       +a1*(- 2.*(1.+z)*(-3.+5.*z)
		    + a1*(5.+z*(-25.+(35.-11.*z)*z)
			  -4.*a1*(1.-z)*(2*(2.+z*(-3.+2*z))
					 +a1*(1.-z)*(-2.+z)))))/(4.*sqr(a1)*sqr(1-a1*(1.-z)));
  Wmixed[2] = 2./a1*(1. - 2.*a1*(1.+z) + 3.*sqr(a1)*z + pow(a1,3)*(1.-z));


  Wmixed[3] = 0.;
  double ratio = 0., rr=1.;
  int itest = (abs(ids[2]->id())%100000)/10000;
  double mix1 = itest==1 ? sTheta_ :  cTheta_;
  double mix2 = itest==1 ? cTheta_ : -sTheta_;
  double ort=sqrt(0.5);
  for(unsigned int ix=0;ix<4;++ix) {
    ratio += rr*(sqr(mix1)*W1P1[ix]+sqr(mix2)*W3P1[ix]+ort*mix1*mix2*Wmixed[ix]); 
    rr*=r;
  }
  ratio /= pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQPP1SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQPP1SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQPP1SplitFn::matrixElement(const double z, const Energy2 t, 
					 const IdList & ids, const double phi, bool) {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M, a2=1-a1;
  double r = sqr(M)/t;
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t+sqr(M)*(sqr(a1)*z*(1.-z)-sqr(a2)*(1.-z)-z));
  double rz = sqrt(z);
  double r2 = sqrt(2.);
  int itest = (abs(ids[2]->id())%100000)/10000;
  double mix1 = itest==1 ? sTheta_ :  cTheta_;
  double mix2 = itest==1 ? cTheta_ : -sTheta_;
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  (*kernal)(0,0,0) = ii*phase*double(pT/M)*r/(1.-z)/rz/r2/a1*
    (-mix1   *(r*a1*(-1.+a1*(1.+z))
	       -0.5*(1.+2.*sqr(a1)*sqr(1.-z)+a1*(-3.+5.*z-2.*sqr(z)))/(1.-a1*(1.-z)))
     +mix2/r2*(r*a1*(1.+a1*(1.-z)-2.*sqr(a1)*(1.-z))
	       +(-1.+sqr(a1)*sqr(1.-z)+a1*(1.-z)*z)/(1.-a1*(1.-z))));
  (*kernal)(0,0,1) = ii/rz/a1*
    (mix1*(sqr(r)*a1*(1.-a1*(1.-z))*(1.-a1*(1.+z))/(1.-z)
	   + 0.25*z*(-2.+a1+sqr(a1)*sqr(1.-z)+z-a1*sqr(z))/sqr(1.-a1*(1.-z))
	   - 0.25*r*(-1.-z+a1*(1.-z)*(5.-4.*a1-3.*z))/(1.-z))
     +mix2/r2*(+ sqr(r)*a1*sqr(1.-a1*(1.-z))/(1.-z)
	       + 0.5*r*(-1.+a1*sqr(1.-z)-z)/(1.-z)
	       + 0.5*(2.-a1*(1.-z)-z)*z/(1.-a1*(1.-z))));
  (*kernal)(0,0,2) = ii*double(pT/M)/phase*r/rz/r2/(1.-z)/a1*
    (mix1   *(-0.5*(1.+a1*(1.-z))*z/(1.-a1*(1.-z))-r*a1*(1.-a1*(1.+z)))-
     mix2/r2*(r*a1*(1.+a1*(-3.+2*a1)*(1.-z))-z/(1.-a1*(1.-z))));
  (*kernal)(0,1,0) = ii/rz/a1/r2*
    (mix1*(- sqr(r)*a1*sqr(1.-a1*(1.-z))/(1.-z)
	   - 0.5*(1.-z)*z/(1.-a1*(1.-z))
	   + 0.5*r*(-1.+a1*(3.-2*a1*(1.-z)+z)))
     +mix2/r2*(r*(1.-sqr(a1)*(1.-z)-2.*a1*z) + (1.-a1)*(1.-z)*z/sqr(1.-a1*(1.-z))
	       -sqr(r)*a1*(1.-a1*(1.-z))*(1.+a1-2.*sqr(a1)*(1.-z)-3.*a1*z)/(1.-z)));
  (*kernal)(0,1,1) = ii*r*double(pT/M)/phase/rz/a1*
    (mix1   *(0.25*(1.-4.*a1)+r*a1*(1.-a1*(1.-z))/(1.-z))-
     mix2/r2*(1.-a1*(1.+z))*(0.5/(1.-a1*(1.-z))-r*a1/(1.-z)));
  (*kernal)(0,1,2) = ii/sqr(phase)*r/rz/r2*
    (mix1*(r*sqr(1.-a1*(1.-z))/(1.-z)-z)
    -mix2/r2*(1.-2.*a1)*(-r*sqr(1.-a1*(1.-z))/(1.-z)+z));
  (*kernal)(1,0,0) = -sqr(phase)*ii*r/rz/r2*
    ( mix1    *(r*sqr(1.-a1*(1.-z))/(1.-z)-z) +
      mix2/r2*(1.-2*a1)*(r*sqr(1.-a1*(1.-z))/(1.-z)-z));
  (*kernal)(1,0,1) = ii*phase*double(pT/M)*r/rz/a1*
    (mix1*((-0.25*(-1 + 4*a1)) - r*a1*(1.-a1*(1.-z))/((-1 + z)))-
     mix2*(1.-a1*(1.+z))/r2*(0.5/(1.-a1*(1.-z))- r*a1/(1.-z)));
  (*kernal)(1,0,2) = ii/rz/a1/r2*
    (mix1    *(+sqr(r)*a1*sqr(1.-a1*(1.-z))/(1.-z)
	       +0.5*(1.-z)*z/(1.-a1*(1.-z))
	       -0.5*r*(-1.+a1*(3.-2.*a1*(1.-z)+z)))
     -mix2/r2*(r*(1.+sqr(a1)*(-1 + z)-2.*a1*z)
	       +(1.-a1)*(1.-z)*z/sqr(1.-a1*(1.-z))
	       +sqr(r)*a1*(1.-a1*(1.-z))*(1 + a1 + 2*sqr(a1)*(-1 + z) - 3*a1*z)/((-1 + z))));
  (*kernal)(1,1,0) = ii*phase*pT/M*r/rz/(1.-z)/r2*
    (-mix1/a1*(0.5*(1.+a1*(1.-z))*z/(1.-a1*(1.-z))+r*a1*(1.-a1*(1.+z)))
     +mix2/r2*(r*(-1.+3.*a1*(1.-z)-2.*sqr(a1)*(1.-z))+z/(a1*(1.-a1*(1.-z)))));
  (*kernal)(1,1,1) = ii/rz/a1*
    (mix1*(-sqr(r)*a1*(1.-a1*(1.-z))*(1.-a1*(1.+z))/(1.-z)
	   -0.25*z*(-2.+a1+sqr(a1)*sqr(1.-z)+z-a1*sqr(z))/sqr(1.-a1*(1.-z))
	   +0.25*r*(-1.-z+a1*(1.-z)*(5.-4.*a1-3.*z))/(1.-z))
     -mix2/r2*(sqr(r)*a1*sqr(1.-a1*(1.-z))/(1.-z)
	       + 0.5*r*(-1.+a1*sqr(1.-z)-z)/(1.-z)
	       + 0.5*(2.-a1*(1.-z)-z)*z/(1.-a1*(1.-z))));
  (*kernal)(1,1,2) = ii*double(pT/M)/phase*r/rz/(1.-z)/r2/a1*
    (mix1*(0.5*(1.-a1*(3.-2.*a1*(1.-z)-2.*z)*(1.-z))/(1.-a1*(1.-z))
	   + a1*r*(1.-a1*(1.+z)))
     +mix2/r2*(a1*r*(1.+a1*(1.-2.*a1)*(1.-z))
	       +(-1 + a1*(a1*(-1 + z) - z)*(-1 + z))/(1.-a1*(1.-z))));
  // testing code
  // DecayMEPtr test1P1(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  // (*test1P1)(0,0,0) =(Complex(0,1)*phase*pT*pow(r,2)*(-1 + a1 + a1*z))/(r2*M*(-1 + z)*rz) - (Complex(0,0.5)*phase*pT*r*(1 + 2*pow(a1,2)*pow(-1 + z,2) + a1*(-3 + 5*z - 2*pow(z,2))))/(r2*a1*M*(1 + a1*(-1 + z))*(-1 + z)*rz);
  // (*test1P1)(0,0,1) =(Complex(0,1)*pow(r,2)*(1 + a1*(-1 + z))*(-1 + a1 + a1*z))/((-1 + z)*rz) + (Complex(0,0.25)*rz*(-2 + a1 + pow(a1,2)*pow(-1 + z,2) + z - a1*pow(z,2)))/(a1*pow(1 + a1*(-1 + z),2)) + (Complex(0,0.25)*r*(-1 - z + a1*(-1 + z)*(-5 + 4*a1 + 3*z)))/(a1*(-1 + z)*rz);
  // (*test1P1)(0,0,2) =(Complex(0,-0.5)*pT*r*(-1 + a1*(-1 + z))*rz)/(r2*a1*phase*M*(1 + a1*(-1 + z))*(-1 + z)) - (Complex(0,1)*pT*pow(r,2)*(-1 + a1 + a1*z))/(r2*phase*M*(-1 + z)*rz);
  // (*test1P1)(0,1,0) =(Complex(0,1)*pow(r,2)*pow(1 + a1*(-1 + z),2))/(r2*(-1 + z)*rz) + (Complex(0,0.5)*(-1 + z)*rz)/(r2*a1*(1 + a1*(-1 + z))) + (Complex(0,0.5)*r*(-1 + a1*(3 + 2*a1*(-1 + z) + z)))/(r2*a1*rz);
  // (*test1P1)(0,1,1) =(Complex(0,0.25)*(1 - 4*a1)*pT*r)/(a1*phase*M*rz) - (Complex(0,1)*pT*pow(r,2)*(1 + a1*(-1 + z)))/(phase*M*(-1 + z)*rz);
  // (*test1P1)(0,1,2) =(Complex(0,-1)*pow(r,2)*pow(1 + a1*(-1 + z),2))/(r2*sqr(phase)*(-1 + z)*rz) - (Complex(0,1)*r*rz)/(r2*sqr(phase));
  // (*test1P1)(1,0,0) =(Complex(0,1)*sqr(phase)*pow(r,2)*pow(1 + a1*(-1 + z),2))/(r2*(-1 + z)*rz) + (Complex(0,1)*sqr(phase)*r*rz)/r2;
  // (*test1P1)(1,0,1) =(Complex(0,-0.25)*(-1 + 4*a1)*phase*pT*r)/(a1*M*rz) - (Complex(0,1)*phase*pT*pow(r,2)*(1 + a1*(-1 + z)))/(M*(-1 + z)*rz);
  // (*test1P1)(1,0,2) =(Complex(0,-1)*pow(r,2)*pow(1 + a1*(-1 + z),2))/(r2*(-1 + z)*rz) - (Complex(0,0.5)*(-1 + z)*rz)/(r2*a1*(1 + a1*(-1 + z))) - (Complex(0,0.5)*r*(-1 + a1*(3 + 2*a1*(-1 + z) + z)))/(r2*a1*rz);
  // (*test1P1)(1,1,0) =(Complex(0,-0.5)*phase*pT*r*(-1 + a1*(-1 + z))*rz)/(r2*a1*M*(1 + a1*(-1 + z))*(-1 + z)) - (Complex(0,1)*phase*pT*pow(r,2)*(-1 + a1 + a1*z))/(r2*M*(-1 + z)*rz);
  // (*test1P1)(1,1,1) =(Complex(0,-1)*pow(r,2)*(1 + a1*(-1 + z))*(-1 + a1 + a1*z))/((-1 + z)*rz) - (Complex(0,0.25)*rz*(-2 + a1 + pow(a1,2)*pow(-1 + z,2) + z - a1*pow(z,2)))/(a1*pow(1 + a1*(-1 + z),2)) - (Complex(0,0.25)*r*(-1 - z + a1*(-1 + z)*(-5 + 4*a1 + 3*z)))/(a1*(-1 + z)*rz);
  // (*test1P1)(1,1,2) =(Complex(0,-0.5)*pT*r*(1 + a1*(3 + 2*a1*(-1 + z) - 2*z)*(-1 + z)))/(r2*a1*phase*M*(1 + a1*(-1 + z))*(-1 + z)*rz) + (Complex(0,1)*pT*pow(r,2)*(-1 + a1 + a1*z))/(r2*phase*M*(-1 + z)*rz);
  // DecayMEPtr test3P1(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  // (*test3P1)(0,0,0) = (Complex(0,-0.5)*phase*pT*pow(r,2)*(1 + a1 + 2*pow(a1,2)*(-1 + z) - a1*z))/(M*(-1 + z)*rz) - (Complex(0,0.5)*phase*pT*r*(-1 + pow(a1,2)*pow(-1 + z,2) - a1*(-1 + z)*z))/(a1*M*(1 + a1*(-1 + z))*(-1 + z)*rz);
  // (*test3P1)(0,0,1) = (Complex(0,-1)*pow(r,2)*pow(1 + a1*(-1 + z),2))/(r2*(-1 + z)*rz) - (Complex(0,0.5)*r*(-1 + a1*pow(-1 + z,2) - z))/(r2*a1*(-1 + z)*rz) + (Complex(0,0.5)*(2 + a1*(-1 + z) - z)*rz)/(r2*a1*(1 + a1*(-1 + z)));	
  // (*test3P1)(0,0,2) = (Complex(0,-0.5)*pT*pow(r,2)*(-1 + a1*(-3 + 2*a1)*(-1 + z)))/(phase*M*(-1 + z)*rz) - (Complex(0,0.5)*pT*r*rz)/(a1*phase*M*(1 + a1*(-1 + z))*(-1 + z));								
  // (*test3P1)(0,1,0) = (Complex(0,0.5)*r*(1/a1 + a1*(-1 + z) - 2*z))/rz + (Complex(0,0.5)*(-1 + a1)*(-1 + z)*rz)/(a1*pow(1 + a1*(-1 + z),2)) + (Complex(0,0.5)*pow(r,2)*(1 + a1*(-1 + z))*(1 + a1 + 2*pow(a1,2)*(-1 + z) - 3*a1*z))/((-1 + z)*rz);	
  // (*test3P1)(0,1,1) = (Complex(0,0.5)*pT*r*(-1 + a1 + a1*z))/(r2*a1*phase*M*(1 + a1*(-1 + z))*rz) + (Complex(0,1)*pT*pow(r,2)*(-1 + a1 + a1*z))/(r2*phase*M*(-1 + z)*rz);							
  // (*test3P1)(0,1,2) = (Complex(0,0.5)*(-1 + 2*a1)*pow(r,2)*pow(1 + a1*(-1 + z),2))/(sqr(phase)*(-1 + z)*rz) + (Complex(0,0.5)*(-1 + 2*a1)*r*rz)/sqr(phase);											
  // (*test3P1)(1,0,0) = (Complex(0,-0.5)*(-1 + 2*a1)*sqr(phase)*pow(r,2)*pow(1 + a1*(-1 + z),2))/((-1 + z)*rz) - Complex(0,0.5)*(-1 + 2*a1)*sqr(phase)*r*rz;											
  // (*test3P1)(1,0,1) = (Complex(0,0.5)*phase*pT*r*(-1 + a1 + a1*z))/(r2*a1*M*(1 + a1*(-1 + z))*rz) + (Complex(0,1)*phase*pT*pow(r,2)*(-1 + a1 + a1*z))/(r2*M*(-1 + z)*rz);							
  // (*test3P1)(1,0,2) = (Complex(0,-0.5)*r*(1/a1 + a1*(-1 + z) - 2*z))/rz - (Complex(0,0.5)*(-1 + a1)*(-1 + z)*rz)/(a1*pow(1 + a1*(-1 + z),2)) - (Complex(0,0.5)*pow(r,2)*(1 + a1*(-1 + z))*(1 + a1 + 2*pow(a1,2)*(-1 + z) - 3*a1*z))/((-1 + z)*rz);
  // (*test3P1)(1,1,0) = (Complex(0,-0.5)*phase*pT*pow(r,2)*(-1 - 3*a1*(-1 + z) + 2*pow(a1,2)*(-1 + z)))/(M*(-1 + z)*rz) - (Complex(0,0.5)*phase*pT*r*rz)/(a1*M*(1 + a1*(-1 + z))*(-1 + z));						
  // (*test3P1)(1,1,1) = (Complex(0,1)*pow(r,2)*pow(1 + a1*(-1 + z),2))/(r2*(-1 + z)*rz) + (Complex(0,0.5)*r*(-1 + a1*pow(-1 + z,2) - z))/(r2*a1*(-1 + z)*rz) - (Complex(0,0.5)*(2 + a1*(-1 + z) - z)*rz)/(r2*a1*(1 + a1*(-1 + z)));	
  // (*test3P1)(1,1,2) = (Complex(0,-0.5)*pT*pow(r,2)*(1 + a1*(-1 + 2*a1)*(-1 + z)))/(phase*M*(-1 + z)*rz) - (Complex(0,0.5)*pT*r*(-1 + a1*(a1*(-1 + z) - z)*(-1 + z)))/(a1*phase*M*(1 + a1*(-1 + z))*(-1 + z)*rz);				
  // cerr << "testing in kernal " << a1 << " " << z << " " << r << " " << pT/M << " " << phi << "\n";
  // for(unsigned int ih1=0;ih1<2;++ih1) {
  //   for(unsigned int ih2=0;ih2<2;++ih2) {
  //     for(unsigned int ih3=0;ih3<3;++ih3) {
  // 	cerr << ih1 << " " << ih2 << " " << ih3 << " " << (*kernal)(ih1,ih2,ih3)  << "\n";
  // 	cerr << "testing diff " << abs((*kernal)(ih1,ih2,ih3)-mix1*(*test1P1)(ih1,ih2,ih3)-mix2*(*test3P1)(ih1,ih2,ih3)) << "\n";
  //     }
  //   }
  // }
  return kernal;
}
