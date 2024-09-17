// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQP3D1SplitFn class.
//

#include "QtoQP3D1SplitFn.h"
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

const double QtoQP3D1SplitFn::pOver_ = 25.;

IBPtr QtoQP3D1SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQP3D1SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQP3D1SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*sqr(GeV*GeV2)) << oenum(state_) << n_ << fixedAlphaS_;
}

void QtoQP3D1SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*sqr(GeV*GeV2)) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

void QtoQP3D1SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletMEProduction<2>(state_,n_,1,1);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQP3D1SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQP3D1SplitFn("Herwig::QtoQP3D1SplitFn", "HwOniumParameters.so HwOniumShower.so");

void QtoQP3D1SplitFn::Init() {

  static ClassDocumentation<QtoQP3D1SplitFn> documentation
    ("The QtoQP3D1SplitFn class implements the branching q-> q' 3D1");

  static Reference<QtoQP3D1SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &QtoQP3D1SplitFn::params_, false, false, true, false, false);
  
  static Parameter<QtoQP3D1SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQP3D1SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<QtoQP3D1SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &QtoQP3D1SplitFn::state_, ccbar, false, false);
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
  static SwitchOption interfaceStatebcbar
    (interfaceState,
     "bcbar",
     "B_c state",
     bcbar);
  
  static Parameter<QtoQP3D1SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQP3D1SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQP3D1SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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
  Energy2 pre =  8./243.*aS2*O1_/pow(a2,4)/M/sqr(M*M);
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

double QtoQP3D1SplitFn::ratioP(const double z, const Energy2 t,
			       const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M;
  double r = sqr(M)/t;
  // coefficients
  double W0 = z*(32*pow(a1,8)*pow(-1 + z,6) - 16*pow(a1,7)*pow(-1 + z,5)*(-19 + 7*z) + 
		 2*pow(a1,6)*pow(-1 + z,4)*(641 + z*(-498 + 97*z)) - 4*pow(a1,5)*pow(-1 + z,3)*(-777 + z*(912 + z*(-337 + 42*z))) + 
		 2*pow(a1,4)*pow(-1 + z,2)*(2361 + z*(-3644 + z*(1915 - 422*z + 30*pow(z,2)))) - 
		 4*pow(a1,3)*(-2 + z)*(-1 + z)*(575 + z*(-807 + z*(343 + z*(-65 + 2*z)))) + 25*(6 + z*(-8 + z*(11 + z*(-10 + 3*z)))) + 
		 10*a1*(-98 + z*(172 + z*(-151 + z*(103 + z*(-37 + 3*z))))) + pow(a1,2)*(2806 + z*(-6432 + z*(5701 + z*(-2756 + z*(814 + z*(-104 + 3*z)))))))/
    (120.*pow(-1 + a1,2)*pow(a1,2)*pow(1 + a1*(-1 + z),6));
  double W1 = (-64*pow(a1,7)*(-2 + z)*pow(-1 + z,4) + 16*pow(a1,6)*pow(-1 + z,3)*(9 + z + 6*pow(z,2)) - 25*(-51 + z*(115 + 2*z*(-31 + 9*z))) + 
	       4*pow(a1,5)*pow(-1 + z,2)*(-531 + z*(1201 + z*(-633 + 59*z))) - 20*a1*(326 + z*(-836 + z*(634 + z*(-193 + 21*z)))) - 
	       pow(a1,4)*(-1 + z)*(8651 + z*(-23336 + z*(20078 + z*(-6224 + 575*z)))) + 
	       2*pow(a1,2)*(6849 + z*(-19877 + z*(19126 + z*(-7502 - 57*(-21 + z)*z)))) + 
	       4*pow(a1,3)*(-3741 + z*(12251 + z*(-14204 + z*(7076 + z*(-1471 + 105*z))))))/(120.*pow(-1 + a1,2)*pow(a1,2)*pow(1 + a1*(-1 + z),4));
  double W2 =(-741 - 45/(-1 + a1) + 2*a1*(675 + 4*a1*(-169 + 88*a1)) + (192*pow(-1 + a1,2))/pow(1 + a1*(-1 + z),3) - 
	      (72*(-1 + a1)*(-29 + 4*a1))/pow(1 + a1*(-1 + z),2) + (24*(108 + a1*(-121 + 28*a1)))/(1 + a1*(-1 + z)) + 
	      (a1*(333 + 2*a1*(-357 + 20*a1*(7 + 2*a1)))*z)/(-1 + a1))/(60.*pow(a1,3));
  double W3 =4*(30 + a1*(149 - 183*z + a1*(-235 + 8*a1*(-7 + z)*(-1 + z) + 11*(20 - 3*z)*z)))/(15.*a1*(1 + a1*(-1 + z)));
  double W4 =-128*(1.-a1)*a1/5.;
  double ratio = (W0+r*(W1+r*(W2+r*(W3+r*W4))))/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQP3D1SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQP3D1SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQP3D1SplitFn::matrixElement(const double z, const Energy2 t, 
					  const IdList & ids, const double phi, bool) {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M, a2=1-a1;
  double r = sqr(M)/t;
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t+sqr(M)*(sqr(a1)*z*(1.-z)-sqr(a2)*(1.-z)-z));
  double rz = sqrt(z);
  double r15 = sqrt(15.);
  double r30= sqrt(30.);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  (*kernal)(0,0,0)=(pT*phase*r*(-48*sqr(a1)*sqr(r)*(1 + a1*(-1 + z)) + (r*(-15 + a1*(-29 + 4*a1*(15 + 4*a1*(-1 + z) - 13*z) + 33*z)))/(-1 + a1) + (8*pow(a1,3)*(-2 + z)*sqr(1.-z) + 5*(5 + 3*(-2 + z)*z) - sqr(a1)*(-1 + z)*(57 + z*(-58 + 17*z)) + a1*(-66 + z*(105 + z*(-56 + 9*z))))/((-1 + a1)*pow(1 + a1*(-1 + z),2))))/(2.*r30*a1*M*(-1 + z)*rz);
  (*kernal)(0,0,1)=((-48*sqr(a1)*pow(r,3)*pow(1 + a1*(-1 + z),2))/(-1 + z) + (r*(15 - 29*a1 - 8*pow(a1,3)*sqr(1.-z) + a1*(38 - 17*z)*z + 5*(-2 + z)*z + 2*sqr(a1)*(-1 + z)*(-11 + 9*z)))/((-1 + a1)*(-1 + z)) + (z*(5 - 4*pow(a1,4)*pow(-1 + z,3) - 5*(-1 + z)*z + pow(a1,3)*sqr(1.-z)*(-19 + 7*z) - a1*(3 + z)*(7 + (-7 + z)*z) - sqr(a1)*(-1 + z)*(31 + z*(-21 + 2*z))))/((-1 + a1)*pow(1 + a1*(-1 + z),3)) + (sqr(r)*(-15 + a1*(-6 + 10*z + a1*(65 - 52*a1*sqr(1.-z) + 8*sqr(a1)*sqr(1.-z) + z*(-98 + 41*z)))))/((-1 + a1)*(-1 + z)))/(2.*r15*a1*rz);
  (*kernal)(0,0,2)=(pT*r*(48*sqr(a1)*sqr(r)*(1 + a1*(-1 + z)) + (z*(-5*(1 + z) + a1*(21 - 8*a1*(-3 + z)*(-1 + z) + 8*sqr(a1)*sqr(1.-z) + (-14 + z)*z)))/((-1 + a1)*pow(1 + a1*(-1 + z),2)) + (r*(15 + a1*(13 - 17*z + 4*a1*(-11 - 4*a1*(-1 + z) + 9*z))))/(-1 + a1)))/(2.*r30*a1*M*phase*(-1 + z)*rz);
  (*kernal)(0,1,0)=((48*sqr(a1)*pow(r,3)*pow(1 + a1*(-1 + z),2))/(-1 + z) + (r*(-5*(-5 + z) + a1*(-41 + 16*a1 + 9*z)))/(-1 + a1) + ((-1 + z)*z*(5*(-2 + z) + a1*(28 + 8*sqr(a1)*sqr(1.-z) + (-23 + z)*z - 2*a1*(-1 + z)*(-13 + 4*z))))/((-1 + a1)*pow(1 + a1*(-1 + z),3)) + (sqr(r)*(15 + a1*(14 - 10*z + a1*(-89 - 16*sqr(a1)*sqr(1.-z) + (114 - 41*z)*z + 4*a1*(-1 + z)*(-19 + 15*z)))))/((-1 + a1)*(-1 + z)))/(2.*r30*a1*rz);
  (*kernal)(0,1,1)=(pT*r*((-48*sqr(a1)*sqr(r)*(1 + a1*(-1 + z)))/(-1 + z) + (r*(-15 + a1*(-21 + 4*a1*(11 + 2*a1*(-1 + z) - 7*z) + 17*z)))/((-1 + a1)*(-1 + z)) + (5*(-3 + z) + a1*(44 - 8*pow(a1,3)*sqr(1.-z) - z*(27 + z) + 2*sqr(a1)*(-1 + z)*(-15 + 7*z) + a1*(-51 + (50 - 7*z)*z)))/((-1 + a1)*pow(1 + a1*(-1 + z),2))))/(2.*r15*a1*M*phase*rz);
  (*kernal)(0,1,2)=(r*((-48*sqr(a1)*sqr(r)*pow(1 + a1*(-1 + z),2))/(-1 + z) + (z*(-15 + a1*(-13 + 4*a1*(11 + 4*a1*(-1 + z) - 7*z) + 9*z)))/((-1 + a1)*(1 + a1*(-1 + z))) + (r*(-15 + a1*(2 - 6*z + a1*(57 + 16*sqr(a1)*sqr(1.-z) - 4*a1*(-1 + z)*(-15 + 19*z) + z*(-98 + 57*z)))))/((-1 + a1)*(-1 + z))))/(2.*r30*a1*pow(phase,2)*rz);
  (*kernal)(1,0,0)=(pow(phase,2)*r*((-48*sqr(a1)*sqr(r)*pow(1 + a1*(-1 + z),2))/(-1 + z) + (z*(-15 + a1*(-13 + 4*a1*(11 + 4*a1*(-1 + z) - 7*z) + 9*z)))/((-1 + a1)*(1 + a1*(-1 + z))) + (r*(-15 + a1*(2 - 6*z + a1*(57 + 16*sqr(a1)*sqr(1.-z) - 4*a1*(-1 + z)*(-15 + 19*z) + z*(-98 + 57*z)))))/((-1 + a1)*(-1 + z))))/(2.*r30*a1*rz);
  (*kernal)(1,0,1)=(pT*phase*r*((48*sqr(a1)*sqr(r)*(1 + a1*(-1 + z)))/(-1 + z) - (r*(-15 + a1*(-21 + 4*a1*(11 + 2*a1*(-1 + z) - 7*z) + 17*z)))/((-1 + a1)*(-1 + z)) + (-5*(-3 + z) + a1*(-44 + 8*pow(a1,3)*sqr(1.-z) + z*(27 + z) - 2*sqr(a1)*(-1 + z)*(-15 + 7*z) + a1*(51 + z*(-50 + 7*z))))/((-1 + a1)*pow(1 + a1*(-1 + z),2))))/(2.*r15*a1*M*rz);
  (*kernal)(1,0,2)=((48*sqr(a1)*pow(r,3)*pow(1 + a1*(-1 + z),2))/(-1 + z) + (r*(-5*(-5 + z) + a1*(-41 + 16*a1 + 9*z)))/(-1 + a1) + ((-1 + z)*z*(5*(-2 + z) + a1*(28 + 8*sqr(a1)*sqr(1.-z) + (-23 + z)*z - 2*a1*(-1 + z)*(-13 + 4*z))))/((-1 + a1)*pow(1 + a1*(-1 + z),3)) + (sqr(r)*(15 + a1*(14 - 10*z + a1*(-89 - 16*sqr(a1)*sqr(1.-z) + (114 - 41*z)*z + 4*a1*(-1 + z)*(-19 + 15*z)))))/((-1 + a1)*(-1 + z)))/(2.*r30*a1*rz);
  (*kernal)(1,1,0)=(pT*phase*r*(-48*sqr(a1)*sqr(r)*(1 + a1*(-1 + z)) + (r*(-15 + a1*(-13 + 4*a1*(11 + 4*a1*(-1 + z) - 9*z) + 17*z)))/(-1 + a1) - (z*(-5*(1 + z) + a1*(21 - 8*a1*(-3 + z)*(-1 + z) + 8*sqr(a1)*sqr(1.-z) + (-14 + z)*z)))/((-1 + a1)*pow(1 + a1*(-1 + z),2))))/(2.*r30*a1*M*(-1 + z)*rz);
  (*kernal)(1,1,1)=((-48*sqr(a1)*pow(r,3)*pow(1 + a1*(-1 + z),2))/(-1 + z) + (r*(15 - 29*a1 - 8*pow(a1,3)*sqr(1.-z) + a1*(38 - 17*z)*z + 5*(-2 + z)*z + 2*sqr(a1)*(-1 + z)*(-11 + 9*z)))/((-1 + a1)*(-1 + z)) + (z*(5 - 4*pow(a1,4)*pow(-1 + z,3) - 5*(-1 + z)*z + pow(a1,3)*sqr(1.-z)*(-19 + 7*z) - a1*(3 + z)*(7 + (-7 + z)*z) - sqr(a1)*(-1 + z)*(31 + z*(-21 + 2*z))))/((-1 + a1)*pow(1 + a1*(-1 + z),3)) + (sqr(r)*(-15 + a1*(-6 + 10*z + a1*(65 - 52*a1*sqr(1.-z) + 8*sqr(a1)*sqr(1.-z) + z*(-98 + 41*z)))))/((-1 + a1)*(-1 + z)))/(2.*r15*a1*rz);
  (*kernal)(1,1,2)=(pT*r*(48*sqr(a1)*sqr(r)*(1 + a1*(-1 + z)) + (r*(15 + a1*(29 - 33*z + 4*a1*(-15 - 4*a1*(-1 + z) + 13*z))))/(-1 + a1) + (-8*pow(a1,3)*(-2 + z)*sqr(1.-z) - 5*(5 + 3*(-2 + z)*z) + sqr(a1)*(-1 + z)*(57 + z*(-58 + 17*z)) + a1*(66 + z*(-105 + (56 - 9*z)*z)))/((-1 + a1)*pow(1 + a1*(-1 + z),2))))/(2.*r30*a1*M*phase*(-1 + z)*rz);
  return kernal;
}
