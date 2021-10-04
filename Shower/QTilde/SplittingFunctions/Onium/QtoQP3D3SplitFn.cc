// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQP3D3SplitFn class.
//

#include "QtoQP3D3SplitFn.h"
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

const double QtoQP3D3SplitFn::pOver_ = 100.;

IBPtr QtoQP3D3SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQP3D3SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQP3D3SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*sqr(GeV*GeV2)) << oenum(state_) << n_ << fixedAlphaS_;
}

void QtoQP3D3SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*sqr(GeV*GeV2)) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

void QtoQP3D3SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletMEProduction<2>(state_,n_,1,3);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQP3D3SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQP3D3SplitFn("Herwig::QtoQP3D3SplitFn", "HwOniumParameters.so HwOniumShower.so");

void QtoQP3D3SplitFn::Init() {

  static ClassDocumentation<QtoQP3D3SplitFn> documentation
    ("The QtoQP3D3SplitFn class implements the branching q-> q' 3D3");

  static Reference<QtoQP3D3SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &QtoQP3D3SplitFn::params_, false, false, true, false, false);
  
  static Parameter<QtoQP3D3SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQP3D3SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<QtoQP3D3SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &QtoQP3D3SplitFn::state_, ccbar, false, false);
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
  
  static Parameter<QtoQP3D3SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQP3D3SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQP3D3SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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
  Energy2 pre =  8./189.*aS2*O1_/pow(a2,4)/M/sqr(M*M);
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

double QtoQP3D3SplitFn::ratioP(const double z, const Energy2 t,
			       const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M;
  double r = sqr(M)/t;
  // coefficients
  double W0 = 2.*z*(84 - 18*pow(a1,5)*(-2 + z)*pow(-1 + z,5) + 3*pow(a1,6)*pow(-1 + z,6) - 288*z + 444*sqr(z) - 384*pow(z,3) + 
		   196*pow(z,4) - 56*pow(z,5) + 7*pow(z,6) + 9*pow(a1,4)*pow(-1 + z,4)*(22 - 26*z + 9*sqr(z)) - 
		   6*pow(a1,3)*pow(-1 + z,3)*(-82 + 147*z - 97*sqr(z) + 22*pow(z,3)) + 
		   3*sqr(a1)*pow(-1 + z,2)*(201 - 474*z + 455*sqr(z) - 202*pow(z,3) + 35*pow(z,4)) - 
		   6*a1*(60 - 234*z + 392*sqr(z) - 361*pow(z,3) + 192*pow(z,4) - 56*pow(z,5) + 7*pow(z,6)))/(15.*pow(-1 + a1,2)*pow(1 + a1*(-1 + z),6));
  double W1 =-2.*(-7 + 6*pow(a1,5)*(-7 + z)*pow(-1 + z,4) + 421*z - 813*sqr(z) + 611*pow(z,3) - 
		  191*pow(z,4) + 21*pow(z,5) + pow(a1,4)*pow(-1 + z,3)*(-161 + 16*z + sqr(z)) - 
		  4*pow(a1,3)*pow(-1 + z,2)*(56 + 104*z - 157*sqr(z) + 51*pow(z,3)) - 
		  2*a1*(7 + 643*z - 1547*sqr(z) + 1393*pow(z,3) - 563*pow(z,4) + 85*pow(z,5)) + 
		  3*sqr(a1)*(42 + 386*z - 1257*sqr(z) + 1351*pow(z,3) - 637*pow(z,4) + 115*pow(z,5)))/(15.*pow(-1 + a1,2)*pow(1 + a1*(-1 + z),4));
  double W2 =-2.*(-62 + 857*z - 750*sqr(z) + 175*pow(z,3) + 4*pow(a1,4)*pow(-1 + z,3)*(-58 + 15*z) - 
		  2*pow(a1,3)*pow(-1 + z,2)*(317 - 124*z + 65*sqr(z)) + sqr(a1)*(510 + 93*z - 1383*sqr(z) + 963*pow(z,3) - 183*pow(z,4)) + 
		  a1*(-46 - 1710*z + 2517*sqr(z) - 1164*pow(z,3) + 231*pow(z,4)))/(15.*(-1 + a1)*pow(1 + a1*(-1 + z),3));
  double W3 =-8.*(-33 + 161*z + a1*(-95 + 40*z - 49*sqr(z)) + 8*sqr(a1)*(16 - 19*z + 3*sqr(z)))/(15.*(1 + a1*(-1 + z)));
  double W4 = 896.*(-1 + a1)*a1/15.;
  double ratio = (W0+r*(W1+r*(W2+r*(W3+r*W4))))/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQP3D3SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQP3D3SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQP3D3SplitFn::matrixElement(const double z, const Energy2 t, 
					  const IdList & ids, const double phi, bool) {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M, a2=1-a1;
  double r = sqr(M)/t;
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t+sqr(M)*(sqr(a1)*z*(1.-z)-sqr(a2)*(1.-z)-z));
  double rz = sqrt(z), r2=sqrt(2.), r3=sqrt(3.);
  double r215 = sqrt(2./15.);
  double r410= sqrt(0.4);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin3)));
  (*kernal)(0,0,0)=2.*r2*pT*pow(phase,3)*sqr(r)*(r*pow(1 + a1*(-1 + z),2) + (-1 + z)*z)/((-1 + a1)*M*pow(-1 + z,3)*rz);
  (*kernal)(0,0,1)=2.*sqr(phase)*r*(pow(r + a1*r*(-1 + z),2)*(5 + 2*pow(a1,2)*pow(-1 + z,2) + z - a1*(7 - 8*z + sqr(z))) + (pow(-1 + z,2)*z*(2 + 2*z - sqr(z) + a1*(-2 + z + sqr(z))))/(1 + a1*(-1 + z)) + r*(-1 + z)*(2 + 7*z + pow(a1,2)*pow(-1 + z,2)*(2 + 3*z) - 2*a1*(2 + 3*z - 6*sqr(z) + pow(z,3))))/(r3*(-1 + a1)*pow(-1 + z,3)*rz);
  (*kernal)(0,0,2)=(2*r215*pT*phase*r*(-(r*(-1 + z)*(6 + 4*pow(a1,2)*pow(-1 + z,2) + 12*z - 3*sqr(z) + 2*a1*(-5 + 4*z + sqr(z)))) - sqr(r)*(8*pow(a1,3)*pow(-1 + z,3) + 5*(2 + z) - pow(a1,2)*pow(-1 + z,2)*(-26 + 3*z) + 2*a1*(-14 + 13*z + sqr(z))) - (pow(-1 + z,2)*(1 + 5*z - 4*sqr(z) + pow(z,3) + pow(a1,2)*pow(-1 + z,2)*(1 + 2*z) - a1*(2 + 5*z - 10*sqr(z) + 3*pow(z,3))))/pow(1 + a1*(-1 + z),2)))/((-1 + a1)*M*pow(-1 + z,3)*rz);
  (*kernal)(0,0,3)=(r410*(-((r*(1 + 2*pow(a1,2)*pow(-1 + z,2) + 13*z - 2*sqr(z) + 3*a1*(-1 + sqr(z))))/(-1 + z)) - (2*pow(r,3)*pow(1 + a1*(-1 + z),2)*(6*pow(a1,2)*pow(-1 + z,2) + 5*(1 + z) - a1*(11 - 12*z + sqr(z))))/pow(-1 + z,3) - (2*sqr(r)*(3 + 4*pow(a1,3)*pow(-1 + z,3) + 12*z + pow(a1,2)*pow(-1 + z,2)*(11 + 4*z) - 2*a1*(5 + 3*z - 9*sqr(z) + pow(z,3))))/pow(-1 + z,2) + (z*(-4 + 3*pow(a1,2)*(-2 + z)*pow(-1 + z,2) - pow(a1,3)*pow(-1 + z,3) + 6*z - 4*sqr(z) + pow(z,3) - 3*a1*(-3 + 6*z - 4*sqr(z) + pow(z,3))))/pow(1 + a1*(-1 + z),3)))/((-1 + a1)*rz);
  (*kernal)(0,0,4)=(2*r215*pT*r*(r*(-1 + z)*(2 + 4*pow(a1,2)*pow(-1 + z,2) + 15*z - 2*sqr(z) + 6*a1*(-1 + sqr(z))) + (pow(-1 + z,2)*z*(6 + 3*pow(a1,2)*pow(-1 + z,2) - 4*z + sqr(z) - 3*a1*(3 - 4*z + sqr(z))))/pow(1 + a1*(-1 + z),2) + sqr(r)*(5 + 8*pow(a1,3)*pow(-1 + z,3) + 10*z + pow(a1,2)*pow(-1 + z,2)*(21 + 2*z) + 6*a1*(-3 + z + 2*sqr(z)))))/((-1 + a1)*M*phase*pow(-1 + z,3)*rz);
  (*kernal)(0,0,5)=(2*r*(((4 + 3*a1*(-1 + z) - z)*pow(-1 + z,2)*sqr(z))/(1 + a1*(-1 + z)) + pow(r + a1*r*(-1 + z),2)*(1 + 2*pow(a1,2)*pow(-1 + z,2) + 5*z + 3*a1*(-1 + sqr(z))) + r*(-1 + z)*z*(5 + 5*pow(a1,2)*pow(-1 + z,2) + 4*z + 2*a1*(-5 + 4*z + sqr(z)))))/(r3*(-1 + a1)*sqr(phase)*pow(-1 + z,3)*rz);
  (*kernal)(0,0,6)=(-2*r2*pT*sqr(r)*rz*(r*pow(1 + a1*(-1 + z),2) + (-1 + z)*z))/((-1 + a1)*M*pow(phase,3)*pow(-1 + z,3));
  (*kernal)(0,1,0)=(-2*r2*sqr(phase)*r*(2*a1*r*(1 + a1*(-1 + z))*(-1 + z)*sqr(z) + pow(-1 + z,2)*sqr(z) + sqr(r)*pow(1 + a1*(-1 + z),3)*(-1 + a1 + a1*z)))/((-1 + a1)*(1 + a1*(-1 + z))*pow(-1 + z,2)*rz);
  (*kernal)(0,1,1)=(2*pT*phase*r*(((3 + 3*a1*(-1 + z) - z)*z)/pow(1 + a1*(-1 + z),2) + (r*(-2 + 3*z + a1*(2 + 4*z)))/(-1 + z) + (sqr(r)*(-5 - 2*a1*(-6 + z) + 2*pow(a1,3)*pow(-1 + z,2) + 3*pow(a1,2)*(-3 + 2*z + sqr(z))))/pow(-1 + z,2)))/(r3*(-1 + a1)*M*rz);
  (*kernal)(0,1,2)=(2*r215*(r*(-1 + a1 + 4*z + 5*a1*z) + (2*pow(r,3)*pow(1 + a1*(-1 + z),2)*(-5 + 4*pow(a1,2)*(-1 + z) + a1*(9 + z)))/pow(-1 + z,2) + ((-1 + z)*z*(3 + 3*pow(a1,2)*pow(-1 + z,2) - 3*z + sqr(z) - 3*a1*(2 - 3*z + sqr(z))))/pow(1 + a1*(-1 + z),3) + (2*sqr(r)*(-3 + 2*pow(a1,3)*pow(-1 + z,2) + 2*a1*(4 + z + sqr(z)) + pow(a1,2)*(-7 + 2*z + 5*sqr(z))))/(-1 + z)))/((-1 + a1)*rz);
  (*kernal)(0,1,3)=(r410*pT*r*((-2*r*(-3 + 4*pow(a1,2)*(-1 + z) + z + a1*(7 + z)))/(-1 + z) - (2*sqr(r)*(-5 + a1*(16 - 6*z) + 6*pow(a1,3)*pow(-1 + z,2) - pow(a1,2)*(17 - 18*z + sqr(z))))/pow(-1 + z,2) + (1 - 2*pow(a1,3)*pow(-1 + z,2) - 3*z + sqr(z) + a1*(-4 + 5*z - 3*sqr(z)) + pow(a1,2)*(5 - 6*z + sqr(z)))/pow(1 + a1*(-1 + z),2)))/((-1 + a1)*M*phase*rz);
  (*kernal)(0,1,4)=(2*r215*r*(-(((-5 + a1*(13 - 3*z) + 8*pow(a1,2)*(-1 + z))*pow(r + a1*r*(-1 + z),2))/pow(-1 + z,2)) - (z*(-2 - 2*a1*(-3 + z) + 4*pow(a1,2)*(-1 + z) + z))/(1 + a1*(-1 + z)) - (2*r*(-1 + 2*pow(a1,3)*pow(-1 + z,2) - 2*z + a1*(4 + 4*z - sqr(z)) + pow(a1,2)*(-5 + 2*z + 3*sqr(z))))/(-1 + z)))/((-1 + a1)*sqr(phase)*rz);
  (*kernal)(0,1,5)=(2*(-1 + 2*a1)*pT*sqr(r)*(r*pow(1 + a1*(-1 + z),2) + (-1 + z)*z))/(r3*(-1 + a1)*M*pow(phase,3)*pow(-1 + z,2)*rz);
  (*kernal)(0,1,6)=0;
  (*kernal)(1,0,0)=0;
  (*kernal)(1,0,1)=(-2*(-1 + 2*a1)*pT*pow(phase,3)*sqr(r)*(r*pow(1 + a1*(-1 + z),2) + (-1 + z)*z))/(r3*(-1 + a1)*M*pow(-1 + z,2)*rz);
  (*kernal)(1,0,2)=(2*r215*sqr(phase)*r*(-(((-5 + a1*(13 - 3*z) + 8*pow(a1,2)*(-1 + z))*pow(r + a1*r*(-1 + z),2))/pow(-1 + z,2)) - (z*(-2 - 2*a1*(-3 + z) + 4*pow(a1,2)*(-1 + z) + z))/(1 + a1*(-1 + z)) - (2*r*(-1 + 2*pow(a1,3)*pow(-1 + z,2) - 2*z + a1*(4 + 4*z - sqr(z)) + pow(a1,2)*(-5 + 2*z + 3*sqr(z))))/(-1 + z)))/((-1 + a1)*rz);
  (*kernal)(1,0,3)=(r410*pT*phase*r*((2*r*(-3 + 4*pow(a1,2)*(-1 + z) + z + a1*(7 + z)))/(-1 + z) + (2*sqr(r)*(-5 + a1*(16 - 6*z) + 6*pow(a1,3)*pow(-1 + z,2) - pow(a1,2)*(17 - 18*z + sqr(z))))/pow(-1 + z,2) + (-1 + 2*pow(a1,3)*pow(-1 + z,2) + 3*z - sqr(z) - pow(a1,2)*(5 - 6*z + sqr(z)) + a1*(4 - 5*z + 3*sqr(z)))/pow(1 + a1*(-1 + z),2)))/((-1 + a1)*M*rz);
  (*kernal)(1,0,4)=(2*r215*(r*(-1 + a1 + 4*z + 5*a1*z) + (2*pow(r,3)*pow(1 + a1*(-1 + z),2)*(-5 + 4*pow(a1,2)*(-1 + z) + a1*(9 + z)))/pow(-1 + z,2) + ((-1 + z)*z*(3 + 3*pow(a1,2)*pow(-1 + z,2) - 3*z + sqr(z) - 3*a1*(2 - 3*z + sqr(z))))/pow(1 + a1*(-1 + z),3) + (2*sqr(r)*(-3 + 2*pow(a1,3)*pow(-1 + z,2) + 2*a1*(4 + z + sqr(z)) + pow(a1,2)*(-7 + 2*z + 5*sqr(z))))/(-1 + z)))/((-1 + a1)*rz);
  (*kernal)(1,0,5)=(2*pT*r*(-(((3 + 3*a1*(-1 + z) - z)*z)/pow(1 + a1*(-1 + z),2)) - (r*(-2 + 3*z + a1*(2 + 4*z)))/(-1 + z) - (sqr(r)*(-5 - 2*a1*(-6 + z) + 2*pow(a1,3)*pow(-1 + z,2) + 3*pow(a1,2)*(-3 + 2*z + sqr(z))))/pow(-1 + z,2)))/(r3*(-1 + a1)*M*phase*rz);
  (*kernal)(1,0,6)=(-2*r2*r*(2*a1*r*(1 + a1*(-1 + z))*(-1 + z)*sqr(z) + pow(-1 + z,2)*sqr(z) + sqr(r)*pow(1 + a1*(-1 + z),3)*(-1 + a1 + a1*z)))/((-1 + a1)*sqr(phase)*(1 + a1*(-1 + z))*pow(-1 + z,2)*rz);
  (*kernal)(1,1,0)=(2*r2*pT*pow(phase,3)*sqr(r)*rz*(r*pow(1 + a1*(-1 + z),2) + (-1 + z)*z))/((-1 + a1)*M*pow(-1 + z,3));
  (*kernal)(1,1,1)=2.*sqr(phase)*r*(((4 + 3*a1*(-1 + z) - z)*pow(-1 + z,2)*sqr(z))/(1 + a1*(-1 + z)) + pow(r + a1*r*(-1 + z),2)*(1 + 2*pow(a1,2)*pow(-1 + z,2) + 5*z + 3*a1*(-1 + sqr(z))) + r*(-1 + z)*z*(5 + 5*pow(a1,2)*pow(-1 + z,2) + 4*z + 2*a1*(-5 + 4*z + sqr(z))))/(r3*(-1 + a1)*pow(-1 + z,3)*rz);
  (*kernal)(1,1,2)=(2*r215*pT*phase*r*(-(r*(-1 + z)*(2 + 4*pow(a1,2)*pow(-1 + z,2) + 15*z - 2*sqr(z) + 6*a1*(-1 + sqr(z)))) - (pow(-1 + z,2)*z*(6 + 3*pow(a1,2)*pow(-1 + z,2) - 4*z + sqr(z) - 3*a1*(3 - 4*z + sqr(z))))/pow(1 + a1*(-1 + z),2) - sqr(r)*(5 + 8*pow(a1,3)*pow(-1 + z,3) + 10*z + pow(a1,2)*pow(-1 + z,2)*(21 + 2*z) + 6*a1*(-3 + z + 2*sqr(z)))))/((-1 + a1)*M*pow(-1 + z,3)*rz);
  (*kernal)(1,1,3)=(r410*(-((r*(1 + 2*pow(a1,2)*pow(-1 + z,2) + 13*z - 2*sqr(z) + 3*a1*(-1 + sqr(z))))/(-1 + z)) - (2*pow(r,3)*pow(1 + a1*(-1 + z),2)*(6*pow(a1,2)*pow(-1 + z,2) + 5*(1 + z) - a1*(11 - 12*z + sqr(z))))/pow(-1 + z,3) - (2*sqr(r)*(3 + 4*pow(a1,3)*pow(-1 + z,3) + 12*z + pow(a1,2)*pow(-1 + z,2)*(11 + 4*z) - 2*a1*(5 + 3*z - 9*sqr(z) + pow(z,3))))/pow(-1 + z,2) + (z*(-4 + 3*pow(a1,2)*(-2 + z)*pow(-1 + z,2) - pow(a1,3)*pow(-1 + z,3) + 6*z - 4*sqr(z) + pow(z,3) - 3*a1*(-3 + 6*z - 4*sqr(z) + pow(z,3))))/pow(1 + a1*(-1 + z),3)))/((-1 + a1)*rz);
  (*kernal)(1,1,4)=(2*r215*pT*r*(r*(-1 + z)*(6 + 4*pow(a1,2)*pow(-1 + z,2) + 12*z - 3*sqr(z) + 2*a1*(-5 + 4*z + sqr(z))) + sqr(r)*(8*pow(a1,3)*pow(-1 + z,3) + 5*(2 + z) - pow(a1,2)*pow(-1 + z,2)*(-26 + 3*z) + 2*a1*(-14 + 13*z + sqr(z))) + (pow(-1 + z,2)*(1 + 5*z - 4*sqr(z) + pow(z,3) + pow(a1,2)*pow(-1 + z,2)*(1 + 2*z) - a1*(2 + 5*z - 10*sqr(z) + 3*pow(z,3))))/pow(1 + a1*(-1 + z),2)))/((-1 + a1)*M*phase*pow(-1 + z,3)*rz);
  (*kernal)(1,1,5)=(2*r*(pow(r + a1*r*(-1 + z),2)*(5 + 2*pow(a1,2)*pow(-1 + z,2) + z - a1*(7 - 8*z + sqr(z))) + (pow(-1 + z,2)*z*(2 + 2*z - sqr(z) + a1*(-2 + z + sqr(z))))/(1 + a1*(-1 + z)) + r*(-1 + z)*(2 + 7*z + pow(a1,2)*pow(-1 + z,2)*(2 + 3*z) - 2*a1*(2 + 3*z - 6*sqr(z) + pow(z,3)))))/(r3*(-1 + a1)*sqr(phase)*pow(-1 + z,3)*rz);
  (*kernal)(1,1,6)=(-2*r2*pT*sqr(r)*(r*pow(1 + a1*(-1 + z),2) + (-1 + z)*z))/((-1 + a1)*M*pow(phase,3)*pow(-1 + z,3)*rz);
  return kernal;
}
