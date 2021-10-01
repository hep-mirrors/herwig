// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQ3D3SplitFn class.
//

#include "QtoQ3D3SplitFn.h"
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

const double QtoQ3D3SplitFn::pOver_ = 10.;

IBPtr QtoQ3D3SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQ3D3SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQ3D3SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*sqr(GeV*GeV2)) << oenum(state_) << n_ << fixedAlphaS_;
}

void QtoQ3D3SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*sqr(GeV*GeV2)) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

void QtoQ3D3SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletME<2>(state_,n_,1,3);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQ3D3SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQ3D3SplitFn("Herwig::QtoQ3D3SplitFn", "HwOniumShower.so HwOniumParameters.so");

void QtoQ3D3SplitFn::Init() {

  static ClassDocumentation<QtoQ3D3SplitFn> documentation
    ("The QtoQ3D3SplitFn class implements the branching q-> q 3D3");

  static Reference<QtoQ3D3SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &QtoQ3D3SplitFn::params_, false, false, true, false, false);
  
  static Parameter<QtoQ3D3SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQ3D3SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<QtoQ3D3SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &QtoQ3D3SplitFn::state_, ccbar, false, false);
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
  
  static Parameter<QtoQ3D3SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQ3D3SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQ3D3SplitFn::guesstz(Energy2 t1,unsigned int iopt,
			     const IdList &ids,
			     double enhance,bool ident,
			     double detune, 
			     Energy2 &t_main, double &z_main) {
  unsigned int pdfopt = iopt!=1 ? 0 : pdfFactor();
  double lower = integOverP(zLimits().first ,ids,pdfopt);
  double upper = integOverP(zLimits().second,ids,pdfopt);
  Energy m = ids[0]->mass();
  double aS2 = fixedAlphaS_ < 0 ? sqr(alpha()->overestimateValue()) : sqr(fixedAlphaS_);
  Energy2 pre =  2./27.*aS2*O1_/m/sqr(m*m);
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

double QtoQ3D3SplitFn::ratioP(const double z, const Energy2 t,
			      const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double W0 = 8.*z*(291 - 414*z + 477*sqr(z) - 324*pow(z,3) + 205*pow(z,4) - 62*pow(z,5) + 19*pow(z,6))/(15.*pow(1 + z,6));
  double W1 = 128.*(7 - 183*z + 119*sqr(z) + 39*pow(z,3) - 90*pow(z,4) + 12*pow(z,5))/(15.*pow(1 + z,4));
  double W2 = 128.*(-89 + 670*z + 240*sqr(z) - 514*pow(z,3) + 229*pow(z,4))/(15.*pow(1 + z,3));
  double W3 = 512.*(97 - 286*z + 37*sqr(z))/(15.*(1 + z));
  double W4 =-57344./15.;
  double ratio = (W0+r*(W1+r*(W2+r*(W3+r*W4))))/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQ3D3SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQ3D3SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQ3D3SplitFn::matrixElement(const double z, const Energy2 t, 
					 const IdList & ids, const double phi, bool) {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  Complex ii(0.,1.); Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t-sqr(m*(1.+z)));
  double rz = sqrt(z), r2=sqrt(2.), r3=sqrt(3.);
  double r215 = sqrt(2./15.);
  double r410= sqrt(0.4);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin3)));
  (*kernal)(0,0,0)=-32.*r2*pT*pow(phase,3)*sqr(r)*((-1 + z)*z + r*sqr(1.+z))/(m*pow(-1 + z,3)*rz);
  (*kernal)(0,0,1)=-16.*sqr(phase)*r*(8*sqr(r)*pow(1 + z,3)*(1 + 2*z) - sqr(1.-z)*z*(-2 - 5*z + sqr(z)) - r*(2 + 15*z + 18*sqr(z) - 16*pow(z,3) - 20*pow(z,4) + pow(z,5)))/(r3*pow(-1 + z,3)*rz*(1 + z));
  (*kernal)(0,0,2)=(8*r215*pT*phase*r*(sqr(1.-z)*(1 + 10*z + sqr(z)) + 4*sqr(r)*pow(1 + z,3)*(6 + 23*z + sqr(z)) - 4*r*sqr(1.+z)*(2 + 12*z - 15*sqr(z) + pow(z,3))))/(m*pow(-1 + z,3)*rz*sqr(1.+z));
  (*kernal)(0,0,3)=(2*r410*(48*r*sqr(1.-z)*z*pow(1 + z,3) + 32*pow(r,3)*pow(1 + z,5)*(1 + 8*z + sqr(z)) - pow(-1 + z,3)*z*(-7 + 3*z - 5*sqr(z) + pow(z,3)) + 8*sqr(r)*pow(1 + z,3)*(-1 - 23*z - 9*sqr(z) + 31*pow(z,3) + 2*pow(z,4))))/(rz*pow(-1 + sqr(z),3));
  (*kernal)(0,0,4)=(-8*r215*pT*r*(sqr(1.-z)*z*(9 + 2*z + sqr(z)) + 4*r*z*sqr(1.+z)*(-13 + 11*z + 2*sqr(z)) + 4*sqr(r)*pow(1 + z,3)*(1 + 23*z + 6*sqr(z))))/(m*phase*pow(-1 + z,3)*rz*sqr(1.+z));
  (*kernal)(0,0,5)=(-16*r*rz*(8*sqr(r)*pow(1 + z,3)*(2 + z) + sqr(1.-z)*z*(5 + z) + r*(-5 - 22*z - 4*sqr(z) + 22*pow(z,3) + 9*pow(z,4))))/(r3*sqr(phase)*pow(-1 + z,3)*(1 + z));
  (*kernal)(0,0,6)=(32*r2*pT*sqr(r)*rz*((-1 + z)*z + r*sqr(1.+z)))/(m*pow(phase,3)*pow(-1 + z,3));
  (*kernal)(0,1,0)=(32*r2*sqr(phase)*r*((-1 + z)*sqr(z) + 2*r*sqr(z)*(1 + z) + sqr(r)*pow(1 + z,3)))/(rz*(-1 + sqr(z)));
  (*kernal)(0,1,1)=(-16*pT*phase*r*(8*sqr(r)*pow(1 + z,3) + 2*r*sqr(1.+z)*(-1 + 5*z) + z*(-3 + 2*z + sqr(z))))/(r3*m*(-1 + z)*rz*sqr(1.+z));
  (*kernal)(0,1,2)=(-8*r215*(24*pow(r,3)*pow(1 + z,5) + sqr(1.-z)*z*(3 + sqr(z)) + 8*sqr(r)*pow(1 + z,3)*(-1 + 2*z + 5*sqr(z)) + r*pow(1 + z,3)*(1 - 14*z + 13*sqr(z))))/((-1 + z)*rz*pow(1 + z,3));
  (*kernal)(0,1,3)=(8*r410*pT*r*(8*sqr(r)*pow(1 + z,3) + 2*r*sqr(1.+z)*(-1 + 5*z) + z*(-3 + 2*z + sqr(z))))/(m*phase*(-1 + z)*rz*sqr(1.+z));
  (*kernal)(0,1,4)=(32*r215*r*((-1 + z)*sqr(z) + 2*r*sqr(z)*(1 + z) + sqr(r)*pow(1 + z,3)))/(sqr(phase)*rz*(-1 + sqr(z)));
  (*kernal)(0,1,5)=0;
  (*kernal)(0,1,6)=0;
  (*kernal)(1,0,0)=0;
  (*kernal)(1,0,1)=0;
  (*kernal)(1,0,2)=(32*r215*sqr(phase)*r*((-1 + z)*sqr(z) + 2*r*sqr(z)*(1 + z) + sqr(r)*pow(1 + z,3)))/(rz*(-1 + sqr(z)));
  (*kernal)(1,0,3)=(-8*r410*pT*phase*r*(8*sqr(r)*pow(1 + z,3) + 2*r*sqr(1.+z)*(-1 + 5*z) + z*(-3 + 2*z + sqr(z))))/(m*(-1 + z)*rz*sqr(1.+z));
  (*kernal)(1,0,4)=(-8*r215*(24*pow(r,3)*pow(1 + z,5) + sqr(1.-z)*z*(3 + sqr(z)) + 8*sqr(r)*pow(1 + z,3)*(-1 + 2*z + 5*sqr(z)) + r*pow(1 + z,3)*(1 - 14*z + 13*sqr(z))))/((-1 + z)*rz*pow(1 + z,3));
  (*kernal)(1,0,5)=(16*pT*r*(8*sqr(r)*pow(1 + z,3) + 2*r*sqr(1.+z)*(-1 + 5*z) + z*(-3 + 2*z + sqr(z))))/(r3*m*phase*(-1 + z)*rz*sqr(1.+z));
  (*kernal)(1,0,6)=(32*r2*r*((-1 + z)*sqr(z) + 2*r*sqr(z)*(1 + z) + sqr(r)*pow(1 + z,3)))/(sqr(phase)*rz*(-1 + sqr(z)));
  (*kernal)(1,1,0)=-32.*r2*pT*pow(phase,3)*sqr(r)*rz*((-1 + z)*z + r*sqr(1.+z))/(m*pow(-1 + z,3));
  (*kernal)(1,1,1)=-16.*sqr(phase)*r*rz*(8*sqr(r)*pow(1 + z,3)*(2 + z) + sqr(1.-z)*z*(5 + z) + r*(-5 - 22*z - 4*sqr(z) + 22*pow(z,3) + 9*pow(z,4)))/(r3*pow(-1 + z,3)*(1 + z));
  (*kernal)(1,1,2)=(8*r215*pT*phase*r*(sqr(1.-z)*z*(9 + 2*z + sqr(z)) + 4*r*z*sqr(1.+z)*(-13 + 11*z + 2*sqr(z)) + 4*sqr(r)*pow(1 + z,3)*(1 + 23*z + 6*sqr(z))))/(m*pow(-1 + z,3)*rz*sqr(1.+z));
  (*kernal)(1,1,3)=(2*r410*(48*r*sqr(1.-z)*z*pow(1 + z,3) + 32*pow(r,3)*pow(1 + z,5)*(1 + 8*z + sqr(z)) - pow(-1 + z,3)*z*(-7 + 3*z - 5*sqr(z) + pow(z,3)) + 8*sqr(r)*pow(1 + z,3)*(-1 - 23*z - 9*sqr(z) + 31*pow(z,3) + 2*pow(z,4))))/(rz*pow(-1 + sqr(z),3));
  (*kernal)(1,1,4)=(-8*r215*pT*r*(sqr(1.-z)*(1 + 10*z + sqr(z)) + 4*sqr(r)*pow(1 + z,3)*(6 + 23*z + sqr(z)) - 4*r*sqr(1.+z)*(2 + 12*z - 15*sqr(z) + pow(z,3))))/(m*phase*pow(-1 + z,3)*rz*sqr(1.+z));
  (*kernal)(1,1,5)=(-16*r*(8*sqr(r)*pow(1 + z,3)*(1 + 2*z) - sqr(1.-z)*z*(-2 - 5*z + sqr(z)) - r*(2 + 15*z + 18*sqr(z) - 16*pow(z,3) - 20*pow(z,4) + pow(z,5))))/(r3*sqr(phase)*pow(-1 + z,3)*rz*(1 + z));
  (*kernal)(1,1,6)=(32*r2*pT*sqr(r)*((-1 + z)*z + r*sqr(1.+z)))/(m*pow(phase,3)*pow(-1 + z,3)*rz);
  return kernal;
}
