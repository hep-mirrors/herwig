// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQ3D2SplitFn class.
//

#include "QtoQ3D2SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQ3D2SplitFn::pOver_ = 15.;

IBPtr QtoQ3D2SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQ3D2SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQ3D2SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*sqr(GeV*GeV2)) << n_ << fixedAlphaS_;
}

void QtoQ3D2SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*sqr(GeV*GeV2)) >> n_ >> fixedAlphaS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQ3D2SplitFn,Sudakov1to2FormFactor>
  describeHerwigQtoQ3D2SplitFn("Herwig::QtoQ3D2SplitFn", "HwOniumShower.so");

void QtoQ3D2SplitFn::Init() {

  static ClassDocumentation<QtoQ3D2SplitFn> documentation
    ("The QtoQ3D2SplitFn class implements the branching q-> q 3D2");

  static Parameter<QtoQ3D2SplitFn,Energy7> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQ3D2SplitFn::O1_, GeV*GeV2*GeV2*GeV2, 0.131*GeV*GeV2*GeV2*GeV2, 0.0*GeV*GeV2*GeV2*GeV2, 10.0*GeV*GeV2*GeV2*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3D2SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQ3D2SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3D2SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQ3D2SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQ3D2SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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

double QtoQ3D2SplitFn::ratioP(const double z, const Energy2 t,
			      const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double W0 = 4.*z*(147 - 54*z + 117*sqr(z) - 84*pow(z,3) + 77*pow(z,4) - 22*pow(z,5) + 11*pow(z,6))/(3.*pow(1 + z,6));
  double W1 = 8.*(61 - 621*z - 166*sqr(z) + 126*pow(z,3) + 9*pow(z,4) + 15*pow(z,5))/(3.*pow(1 + z,4));
  double W2 = 256.*(-19 + 59*z + 66*sqr(z) - 23*pow(z,3) + 5*pow(z,4))/(3.*pow(1 + z,3));
  double W3 =-256.*(-47 + 74*z + sqr(z))/(3.*(1 + z));
  double W4 =-8192./3.;
  double ratio = (W0+r*(W1+r*(W2+r*(W3+r*W4))))/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQ3D2SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQ3D2SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQ3D2SplitFn::matrixElement(const double z, const Energy2 t, 
					 const IdList & ids, const double phi, bool) {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  Complex ii(0.,1.); Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t-sqr(m*(1.+z)));
  double rz = sqrt(z), r23 = sqrt(2./3.);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin2)));
  (*kernal)(0,0,0)=(8*r23*sqr(phase)*r*(8*sqr(r)*pow(1 + z,3) + z*(5 - 6*z + 3*sqr(z) - 2*pow(z,3)) - r*(5 + 17*z + 5*sqr(z) - 5*pow(z,3) + 2*pow(z,4))))/(sqr(1.-z)*rz*(1 + z));
  (*kernal)(0,0,1)=(-2*r23*pT*phase*r*(5 + 6*z - 14*pow(z,3) + 3*pow(z,4) + 16*sqr(r)*pow(1 + z,3)*(3 + z) + 2*r*sqr(1.+z)*(-17 - 8*z + sqr(z))))/(m*rz*pow(-1 + sqr(z),2));
  (*kernal)(0,0,2)=(-64*pow(r,3)*pow(1 + z,5) + 2*sqr(1.-z)*z*(5 - 2*z + sqr(z)) - 8*sqr(r)*pow(1 + z,3)*(-5 - 14*z + 7*sqr(z)) - 4*r*sqr(1.+z)*(1 + 13*z - 17*sqr(z) + 3*pow(z,3)))/(rz*pow(-1 + sqr(z),2));
  (*kernal)(0,0,3)=(4*r23*pT*r*(8*sqr(r)*pow(1 + z,3)*(1 + 3*z) + r*sqr(1.+z)*(-3 - 32*z + 11*sqr(z)) + z*(9 - 5*z - 5*sqr(z) + pow(z,3))))/(m*phase*rz*pow(-1 + sqr(z),2));
  (*kernal)(0,0,4)=(8*r23*r*rz*(8*sqr(r)*pow(1 + z,3) + z*(7 - 8*z + sqr(z)) + r*(-7 - 21*z - 5*sqr(z) + 9*pow(z,3))))/(sqr(phase)*sqr(1.-z)*(1 + z));
  (*kernal)(0,1,0)=(-8*r23*pT*phase*r*(8*sqr(r)*pow(1 + z,3) + r*sqr(1.+z)*(-5 + 7*z) + z*(-3 + 2*z + sqr(z))))/(m*(-1 + z)*rz*sqr(1.+z));
  (*kernal)(0,1,1)=(-2*r23*(48*pow(r,3)*pow(1 + z,5) + 2*sqr(1.-z)*z*(3 + sqr(z)) + r*pow(1 + z,3)*(5 - 22*z + 17*sqr(z)) + 2*sqr(r)*pow(1 + z,3)*(-17 - 2*z + 31*sqr(z))))/((-1 + z)*rz*pow(1 + z,3));
  (*kernal)(0,1,2)=(4*pT*r*(1 - 5*z + 3*sqr(z) + pow(z,3) + 16*sqr(r)*pow(1 + z,3) + 2*r*sqr(1.+z)*(-5 + 7*z)))/(m*phase*(-1 + z)*rz*sqr(1.+z));
  (*kernal)(0,1,3)=(4*r23*r*(8*sqr(r)*pow(1 + z,3) + z*(3 - 8*z + 5*sqr(z)) + r*(-3 - 9*z + 7*sqr(z) + 13*pow(z,3))))/(sqr(phase)*rz*(-1 + sqr(z)));
  (*kernal)(0,1,4)=0;
  (*kernal)(1,0,0)=0;
  (*kernal)(1,0,1)=(-4*r23*sqr(phase)*r*(8*sqr(r)*pow(1 + z,3) + z*(3 - 8*z + 5*sqr(z)) + r*(-3 - 9*z + 7*sqr(z) + 13*pow(z,3))))/(rz*(-1 + sqr(z)));
  (*kernal)(1,0,2)=(4*pT*phase*r*(1 - 5*z + 3*sqr(z) + pow(z,3) + 16*sqr(r)*pow(1 + z,3) + 2*r*sqr(1.+z)*(-5 + 7*z)))/(m*(-1 + z)*rz*sqr(1.+z));
  (*kernal)(1,0,3)=(2*r23*(48*pow(r,3)*pow(1 + z,5) + 2*sqr(1.-z)*z*(3 + sqr(z)) + r*pow(1 + z,3)*(5 - 22*z + 17*sqr(z)) + 2*sqr(r)*pow(1 + z,3)*(-17 - 2*z + 31*sqr(z))))/((-1 + z)*rz*pow(1 + z,3));
  (*kernal)(1,0,4)=(-8*r23*pT*r*(8*sqr(r)*pow(1 + z,3) + r*sqr(1.+z)*(-5 + 7*z) + z*(-3 + 2*z + sqr(z))))/(m*phase*(-1 + z)*rz*sqr(1.+z));
  (*kernal)(1,1,0)=(-8*r23*sqr(phase)*r*rz*(8*sqr(r)*pow(1 + z,3) + z*(7 - 8*z + sqr(z)) + r*(-7 - 21*z - 5*sqr(z) + 9*pow(z,3))))/(sqr(1.-z)*(1 + z));
  (*kernal)(1,1,1)=(4*r23*pT*phase*r*(8*sqr(r)*pow(1 + z,3)*(1 + 3*z) + r*sqr(1.+z)*(-3 - 32*z + 11*sqr(z)) + z*(9 - 5*z - 5*sqr(z) + pow(z,3))))/(m*rz*pow(-1 + sqr(z),2));
  (*kernal)(1,1,2)=(64*pow(r,3)*pow(1 + z,5) - 2*sqr(1.-z)*z*(5 - 2*z + sqr(z)) + 8*sqr(r)*pow(1 + z,3)*(-5 - 14*z + 7*sqr(z)) + 4*r*sqr(1.+z)*(1 + 13*z - 17*sqr(z) + 3*pow(z,3)))/(rz*pow(-1 + sqr(z),2));
  (*kernal)(1,1,3)=(-2*r23*pT*r*(5 + 6*z - 14*pow(z,3) + 3*pow(z,4) + 16*sqr(r)*pow(1 + z,3)*(3 + z) + 2*r*sqr(1.+z)*(-17 - 8*z + sqr(z))))/(m*phase*rz*pow(-1 + sqr(z),2));
  (*kernal)(1,1,4)=(-8*r23*r*(8*sqr(r)*pow(1 + z,3) + z*(5 - 6*z + 3*sqr(z) - 2*pow(z,3)) - r*(5 + 17*z + 5*sqr(z) - 5*pow(z,3) + 2*pow(z,4))))/(sqr(phase)*sqr(1.-z)*rz*(1 + z));
  return kernal;
}

