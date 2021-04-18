// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQ1D2SplitFn class.
//

#include "QtoQ1D2SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQ1D2SplitFn::pOver_ = 20.;

IBPtr QtoQ1D2SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQ1D2SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQ1D2SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*sqr(GeV*GeV2)) << n_ << fixedAlphaS_;
}

void QtoQ1D2SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*sqr(GeV*GeV2)) >> n_ >> fixedAlphaS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQ1D2SplitFn,Sudakov1to2FormFactor>
  describeHerwigQtoQ1D2SplitFn("Herwig::QtoQ1D2SplitFn", "HwOniumShower.so");

void QtoQ1D2SplitFn::Init() {

  static ClassDocumentation<QtoQ1D2SplitFn> documentation
    ("The QtoQ1D2SplitFn class implements the branching q-> q 1D2");

  static Parameter<QtoQ1D2SplitFn,Energy7> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQ1D2SplitFn::O1_, GeV*GeV2*GeV2*GeV2, 0.131*GeV*GeV2*GeV2*GeV2, 0.0*GeV*GeV2*GeV2*GeV2, 10.0*GeV*GeV2*GeV2*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQ1D2SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQ1D2SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<QtoQ1D2SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQ1D2SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQ1D2SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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

double QtoQ1D2SplitFn::ratioP(const double z, const Energy2 t,
			      const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double W0 = 8*z*(73 - 42*z + 31*pow(z,2) - 44*pow(z,3) + 55*pow(z,4) - 10*pow(z,5) + pow(z,6))/(3.*pow(1 + z,6));
  double W1 =-128.*(-2 + 39*z - 8*pow(z,2) - 26*pow(z,3) + 6*pow(z,4) + 7*pow(z,5))/(3.*pow(1 + z,4));
  double W2 =-128.*(22 - 113*z - 63*pow(z,2) + 105*pow(z,3) + pow(z,4))/(3.*pow(1 + z,3));
  double W3 = 10240.*(1.-2*z)/(3.*(1.+z));
  double W4 =-8192./3.;
  double ratio = (W0+r*(W1+r*(W2+r*(W3+r*W4))))/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQ1D2SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQ1D2SplitFn " << ratio << "\n";
  return ratio;
}


DecayMEPtr QtoQ1D2SplitFn::matrixElement(const double z, const Energy2 t, 
					 const IdList & ids, const double phi, bool) {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  Complex ii(0.,1.); Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t-sqr(m*(1.+z)));
  double rz = sqrt(z), r23 = sqrt(2./3.);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin2)));
  (*kernal)(0,0,0)=-16.*sqr(phase)*r*(2*sqr(r)*(-1 + z)*pow(1 + z,3) + sqr(z)*(3 - 4*z + sqr(z)) + r*z*(-1 - 7*z - 3*sqr(z) + 3*pow(z,3)))/(sqr(1.-z)*rz*(1 + z));
  (*kernal)(0,0,1)= 16.*pT*phase*r*(-2*(-1 + z)*z + 4*sqr(r)*(-1 + z)*pow(1 + z,3) + r*sqr(1.+z)*(1 - 8*z + 3*sqr(z)))/(m*rz*sqr(1.-sqr(z)));
  (*kernal)(0,0,2)= 2.*r23*(48*pow(r,3)*(-1 + z)*pow(1 + z,5) + sqr(1.-z)*z*(-7 + 3*z - 5*sqr(z) + pow(z,3)) + 2*r*pow(1 + z,3)*(-1 + 19*z - 23*sqr(z) + 5*pow(z,3)) + 8*sqr(r)*pow(1 + z,3)*(2 - 9*z - 12*sqr(z) + 7*pow(z,3)))/(sqr(1.-z)*rz*pow(1 + z,3));
  (*kernal)(0,0,3)=-16.*pT*r*(-2*(-1 + z)*z + 4*sqr(r)*(-1 + z)*pow(1 + z,3) + r*sqr(1.+z)*(1 - 8*z + 3*sqr(z)))/(m*phase*rz*sqr(1.-sqr(z)));
  (*kernal)(0,0,4)=-16.*r*(2*sqr(r)*(-1 + z)*pow(1 + z,3) + sqr(z)*(3 - 4*z + sqr(z)) + r*z*(-1 - 7*z - 3*sqr(z) + 3*pow(z,3)))/(sqr(phase)*sqr(1.-z)*rz*(1 + z));
  (*kernal)(0,1,0)= 16.*pT*phase*r*(sqr(1.-z)*z + 2*sqr(r)*pow(1 + z,3) + 2*r*z*(-1 + sqr(z)))/(m*sqr(1.-z)*rz*(1 + z));
  (*kernal)(0,1,1)=  8.*(pow(-1 + z,3)*z + 8*pow(r,3)*pow(1 + z,5) + 4*r*z*sqr(1.-sqr(z)) + 2*sqr(r)*pow(1 + z,3)*(-1 - 4*z + 5*sqr(z)))/(rz*sqr(1.-sqr(z)));
  (*kernal)(0,1,2)=(-4*r23*pT*r*(24*sqr(r)*pow(1 + z,3) + sqr(1.-z)*(1 + 3*z) + 8*r*(-1 - 2*z + sqr(z) + 2*pow(z,3))))/(m*phase*sqr(1.-z)*rz*(1 + z));
  (*kernal)(0,1,3)=(-16.*r*(sqr(1.-z)*z + 4*sqr(r)*pow(1 + z,3) + r*(-1 - 5*z + sqr(z) + 5*pow(z,3))))/(sqr(phase)*sqr(1.-z)*rz);
  (*kernal)(0,1,4)=(32*pT*sqr(r)*((-1 + z)*z + r*sqr(1.+z)))/(m*pow(phase,3)*sqr(1.-z)*rz);
  (*kernal)(1,0,0)=(32*pT*pow(phase,3)*sqr(r)*((-1 + z)*z + r*sqr(1.+z)))/(m*sqr(1.-z)*rz);
  (*kernal)(1,0,1)=(16.*sqr(phase)*r*(sqr(1.-z)*z + 4*sqr(r)*pow(1 + z,3) + r*(-1 - 5*z + sqr(z) + 5*pow(z,3))))/(sqr(1.-z)*rz);
  (*kernal)(1,0,2)=(-4*r23*pT*phase*r*(24*sqr(r)*pow(1 + z,3) + sqr(1.-z)*(1 + 3*z) + 8*r*(-1 - 2*z + sqr(z) + 2*pow(z,3))))/(m*sqr(1.-z)*rz*(1 + z));
  (*kernal)(1,0,3)=(-8*(pow(-1 + z,3)*z + 8*pow(r,3)*pow(1 + z,5) + 4*r*z*sqr(1.-sqr(z)) + 2*sqr(r)*pow(1 + z,3)*(-1 - 4*z + 5*sqr(z))))/(rz*sqr(1.-sqr(z)));
  (*kernal)(1,0,4)=(16.*pT*r*(sqr(1.-z)*z + 2*sqr(r)*pow(1 + z,3) + 2*r*z*(-1 + sqr(z))))/(m*phase*sqr(1.-z)*rz*(1 + z));
  (*kernal)(1,1,0)=(16.*sqr(phase)*r*(2*sqr(r)*(-1 + z)*pow(1 + z,3) + sqr(z)*(3 - 4*z + sqr(z)) + r*z*(-1 - 7*z - 3*sqr(z) + 3*pow(z,3))))/(sqr(1.-z)*rz*(1 + z));
  (*kernal)(1,1,1)=(-16.*pT*phase*r*(-2*(-1 + z)*z + 4*sqr(r)*(-1 + z)*pow(1 + z,3) + r*sqr(1.+z)*(1 - 8*z + 3*sqr(z))))/(m*rz*sqr(1.-sqr(z)));
  (*kernal)(1,1,2)=(-2*r23*(48*pow(r,3)*(-1 + z)*pow(1 + z,5) + sqr(1.-z)*z*(-7 + 3*z - 5*sqr(z) + pow(z,3)) + 2*r*pow(1 + z,3)*(-1 + 19*z - 23*sqr(z) + 5*pow(z,3)) + 8*sqr(r)*pow(1 + z,3)*(2 - 9*z - 12*sqr(z) + 7*pow(z,3))))/(sqr(1.-z)*rz*pow(1 + z,3));
  (*kernal)(1,1,3)=(16.*pT*r*(-2*(-1 + z)*z + 4*sqr(r)*(-1 + z)*pow(1 + z,3) + r*sqr(1.+z)*(1 - 8*z + 3*sqr(z))))/(m*phase*rz*sqr(1.-sqr(z)));
  (*kernal)(1,1,4)=(16.*r*(2*sqr(r)*(-1 + z)*pow(1 + z,3) + sqr(z)*(3 - 4*z + sqr(z)) + r*z*(-1 - 7*z - 3*sqr(z) + 3*pow(z,3))))/(sqr(phase)*sqr(1.-z)*rz*(1 + z));
  return kernal;
}
