// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQ3D1SplitFn class.
//

#include "QtoQ3D1SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQ3D1SplitFn::pOver_ = 10.;

IBPtr QtoQ3D1SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQ3D1SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQ3D1SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*sqr(GeV*GeV2)) << n_ << fixedAlphaS_;
}

void QtoQ3D1SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*sqr(GeV*GeV2)) >> n_ >> fixedAlphaS_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQ3D1SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQ3D1SplitFn("Herwig::QtoQ3D1SplitFn", "HwOniumShower.so");

void QtoQ3D1SplitFn::Init() {

  static ClassDocumentation<QtoQ3D1SplitFn> documentation
    ("The QtoQ3D1SplitFn class implements the branching q-> q 3D1");

  static Parameter<QtoQ3D1SplitFn,Energy7> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQ3D1SplitFn::O1_, GeV*GeV2*GeV2*GeV2, 0.131*GeV*GeV2*GeV2*GeV2, 0.0*GeV*GeV2*GeV2*GeV2, 10.0*GeV*GeV2*GeV2*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3D1SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQ3D1SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3D1SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQ3D1SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQ3D1SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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

double QtoQ3D1SplitFn::ratioP(const double z, const Energy2 t,
			      const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double W0 = 4.*z*(73 - 82*z + 1071*sqr(z) - 412*pow(z,3) - 105*pow(z,4) + 14*pow(z,5) + 17*pow(z,6))/(15.*pow(1 + z,6));
  double W1 =-8.*(-673 + 677*z + 974*sqr(z) + 994*pow(z,3) + 275*pow(z,4) + 57*pow(z,5))/(15.*pow(1 + z,4));
  double W2 =-64.*(557 - 305*z - 885*sqr(z) + 197*pow(z,3) + 28*pow(z,4))/(15.*pow(1 + z,3));
  double W3 =-256.*(-211 + 178*z + 29*sqr(z))/(15.*(1 + z));
  double W4 =-8192./5.;
  double ratio = (W0+r*(W1+r*(W2+r*(W3+r*W4))))/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQ3D1SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQ3D1SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQ3D1SplitFn::matrixElement(const double z, const Energy2 t, 
					 const IdList & ids, const double phi, bool) {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  Complex ii(0.,1.); Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t-sqr(m*(1.+z)));
  double rz = sqrt(z), r215 = sqrt(2./15.), r15 = sqrt(15.);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  (*kernal)(0,0,0)=(-2*r215*pT*phase*r*(17 - 5*z + 7*sqr(z) + 5*pow(z,3) + 22*r*(-3 + z)*pow(1 + z,2) + 48*sqr(r)*pow(1 + z,3)))/(m*(-1 + z)*rz*pow(1 + z,2));
  (*kernal)(0,0,1)=(-2*(96*pow(r,3)*pow(1 + z,5) + 4*r*pow(1 + z,3)*(5 + z) + 4*sqr(r)*pow(1 + z,3)*(-31 - 30*z + 17*sqr(z)) - z*(1 + 30*z - 36*sqr(z) + 2*pow(z,3) + 3*pow(z,4))))/(r15*(-1 + z)*rz*pow(1 + z,3));
  (*kernal)(0,0,2)=(4*r215*pT*r*(24*sqr(r)*pow(1 + z,3) + r*pow(1 + z,2)*(-25 + 3*z) + z*(-1 + 12*z + sqr(z))))/(m*phase*(-1 + z)*rz*pow(1 + z,2));
  (*kernal)(0,1,0)=(2*r215*(48*pow(r,3)*pow(1 + z,5) + 2*pow(-1 + z,2)*z*(3 + sqr(z)) + r*pow(1 + z,3)*(17 - 18*z + sqr(z)) + 2*sqr(r)*pow(1 + z,3)*(-33 - 34*z + 15*sqr(z))))/((-1 + z)*rz*pow(1 + z,3));
  (*kernal)(0,1,1)=(-8*pT*r*(5 - 4*z + sqr(z) - 2*pow(z,3) + 24*sqr(r)*pow(1 + z,3) + r*pow(1 + z,2)*(-31 + 5*z)))/(r15*m*phase*(-1 + z)*rz*pow(1 + z,2));
  (*kernal)(0,1,2)=(-4*r215*r*(24*sqr(r)*pow(1 + z,3) - z*(-25 + 24*z + sqr(z)) + r*(-25 - 75*z - 27*sqr(z) + 23*pow(z,3))))/(sqr(phase)*rz*(-1 + sqr(z)));
  (*kernal)(1,0,0)=(-4*r215*sqr(phase)*r*(24*sqr(r)*pow(1 + z,3) - z*(-25 + 24*z + sqr(z)) + r*(-25 - 75*z - 27*sqr(z) + 23*pow(z,3))))/(rz*(-1 + sqr(z)));
  (*kernal)(1,0,1)=(8*pT*phase*r*(5 - 4*z + sqr(z) - 2*pow(z,3) + 24*sqr(r)*pow(1 + z,3) + r*pow(1 + z,2)*(-31 + 5*z)))/(r15*m*(-1 + z)*rz*pow(1 + z,2));
  (*kernal)(1,0,2)=(2*r215*(48*pow(r,3)*pow(1 + z,5) + 2*pow(-1 + z,2)*z*(3 + sqr(z)) + r*pow(1 + z,3)*(17 - 18*z + sqr(z)) + 2*sqr(r)*pow(1 + z,3)*(-33 - 34*z + 15*sqr(z))))/((-1 + z)*rz*pow(1 + z,3));
  (*kernal)(1,1,0)=(-4*r215*pT*phase*r*(24*sqr(r)*pow(1 + z,3) + r*pow(1 + z,2)*(-25 + 3*z) + z*(-1 + 12*z + sqr(z))))/(m*(-1 + z)*rz*pow(1 + z,2));
  (*kernal)(1,1,1)=(-2*(96*pow(r,3)*pow(1 + z,5) + 4*r*pow(1 + z,3)*(5 + z) + 4*sqr(r)*pow(1 + z,3)*(-31 - 30*z + 17*sqr(z)) - z*(1 + 30*z - 36*sqr(z) + 2*pow(z,3) + 3*pow(z,4))))/(r15*(-1 + z)*rz*pow(1 + z,3));
  (*kernal)(1,1,2)=(2*r215*pT*r*(17 - 5*z + 7*sqr(z) + 5*pow(z,3) + 22*r*(-3 + z)*pow(1 + z,2) + 48*sqr(r)*pow(1 + z,3)))/(m*phase*(-1 + z)*rz*pow(1 + z,2));
  return kernal;
}
