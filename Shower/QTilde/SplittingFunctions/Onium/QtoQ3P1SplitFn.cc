// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQ3P1SplitFn class.
//

#include "QtoQ3P1SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQ3P1SplitFn::pOver_ = 1.5;

IBPtr QtoQ3P1SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQ3P1SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQ3P1SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*sqr(GeV2)) << n_ << fixedAlphaS_;
}

void QtoQ3P1SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*sqr(GeV2)) >> n_ >> fixedAlphaS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQ3P1SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQ3P1SplitFn("Herwig::QtoQ3P1SplitFn", "HwOniumShower.so HwOniumParameters.so");

void QtoQ3P1SplitFn::Init() {

  static ClassDocumentation<QtoQ3P1SplitFn> documentation
    ("The QtoQ3P1SplitFn class implements the branching q-> q' P1");

  static Parameter<QtoQ3P1SplitFn,Energy5> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQ3P1SplitFn::O1_, GeV*GeV2*GeV2, 0.794*GeV*GeV2*GeV2, 0.0*GeV*GeV2*GeV2, 10.0*GeV*GeV2*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3P1SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQ3P1SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3P1SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQ3P1SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);
}

void QtoQ3P1SplitFn::guesstz(Energy2 t1,unsigned int iopt,
			     const IdList &ids,
			     double enhance,bool ident,
			     double detune, 
			     Energy2 &t_main, double &z_main) {
  unsigned int pdfopt = iopt!=1 ? 0 : pdfFactor();
  double lower = integOverP(zLimits().first ,ids,pdfopt);
  double upper = integOverP(zLimits().second,ids,pdfopt);
  Energy m = ids[0]->mass();
  double aS2 = fixedAlphaS_ < 0 ? sqr(alpha()->overestimateValue()) : sqr(fixedAlphaS_);
  Energy2 pre = 64./81.*aS2*O1_/m/sqr(m);
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

double QtoQ3P1SplitFn::ratioP(const double z, const Energy2 t,
			      const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double W0= 0.5*z*(17.-4.*z+6.*sqr(z)-4.*pow(z,3)+pow(z,4))/pow(1.+z,4);
  double W1 = (17.-55.*z+7.*sqr(z)-pow(z,3))/sqr(1.+z);
  double W2 =-4./(1.+z)*(12.-13.*z-sqr(z));
  double W3 = 32.;
  double ratio = (W0+r*(W1+r*(W2+r*W3)))/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQ3P1SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQ3P1SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQ3P1SplitFn::matrixElement(const double z, const Energy2 t, 
					 const IdList & ids, const double phi, bool) {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t-sqr(m*(1.+z)));
  double rz = sqrt(z), r2 = sqrt(2.);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  
  (*kernal)(0,0,0) = ii*phase*pT/m/rz*r/(1.-z)*(4.*r-(3.+sqr(z))/(1.+z));
  (*kernal)(0,0,1) =-ii/rz/r2*((-3.+z)*z/(1.+z) - 4.*sqr(r)*sqr(1.+z)/(1.-z) + 2.*r*(1.+4.*z-sqr(z))/(1.-z));
  (*kernal)(0,0,2) =-4.*ii*pT/m*r*rz/phase/(1.-z)*(r-1./(1.+z));
  (*kernal)(0,1,0) = ii/rz*(3.*r*(1.-z)+2.*(1.-z)*z/sqr(1.+z)-4.*sqr(r)*(1.+z));
  (*kernal)(0,1,1) = ii*r2*pT/m/rz/phase*r*(2.*r-(1.-z)/(1.+z));
  (*kernal)(0,1,2) = 0;								   
  (*kernal)(1,0,0) = 0;								   
  (*kernal)(1,0,1) = ii*r2*pT/m/rz*phase*r*(2.*r-(1.-z)/(1.+z));
  (*kernal)(1,0,2) =-ii/rz*(3.*r*(1.-z)+2.*(1.-z)*z/sqr(1.+z)-4.*sqr(r)*(1.+z));
  (*kernal)(1,1,0) =-4.*ii*phase*pT/m*r*rz/(1.-z)*(r-1./(1.+z));
  (*kernal)(1,1,1) = ii/rz/r2*((-3.+z)*z/(1.+z) -4.*sqr(r)*sqr(1.+z)/(1.-z) + 2.*r*(1.+4*z-sqr(z))/(1.-z));
  (*kernal)(1,1,2) = ii/phase*pT/m/rz*r/(1.-z)*(4.*r-(3.+sqr(z))/(1.+z));


  
  // // testing code
  // DecayMEPtr test(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  // (*test)(0,0,0) =(Complex(0,4)*phase*pT*pow(r,2))/(sqrt(z)*(m - m*z)) + (Complex(0,1)*phase*pT*r*(3 + pow(z,2)))/(m*(-1 + z)*sqrt(z)*(1 + z));
  // (*test)(0,0,1) =(Complex(0,-1)*(-3 + z)*sqrt(z))/(sqrt(2)*(1 + z)) - (Complex(0,2)*sqrt(2)*pow(r,2)*pow(1 + z,2))/((-1 + z)*sqrt(z)) - (Complex(0,1)*sqrt(2)*r*(-1 - 4*z + pow(z,2)))/((-1 + z)*sqrt(z));
  // (*test)(0,0,2) =(Complex(0,4)*pT*pow(r,2)*sqrt(z))/(phase*m*(-1 + z)) + (Complex(0,4)*pT*r*sqrt(z))/(phase*(m - m*pow(z,2)));
  // (*test)(0,1,0) =(Complex(0,-3)*r*(-1 + z))/sqrt(z) - (Complex(0,2)*(-1 + z)*sqrt(z))/pow(1 + z,2) - (Complex(0,4)*pow(r,2)*(1 + z))/sqrt(z);
  // (*test)(0,1,1) =(Complex(0,2)*sqrt(2)*pT*pow(r,2))/(phase*m*sqrt(z)) + (Complex(0,1)*sqrt(2)*pT*r*(-1 + z))/(phase*m*sqrt(z)*(1 + z));
  // (*test)(0,1,2) =0;
  // (*test)(1,0,0) =0;
  // (*test)(1,0,1) =(Complex(0,2)*sqrt(2)*phase*pT*pow(r,2))/(m*sqrt(z)) + (Complex(0,1)*sqrt(2)*phase*pT*r*(-1 + z))/(m*sqrt(z)*(1 + z));
  // (*test)(1,0,2) =(Complex(0,3)*r*(-1 + z))/sqrt(z) + (Complex(0,2)*(-1 + z)*sqrt(z))/pow(1 + z,2) + (Complex(0,4)*pow(r,2)*(1 + z))/sqrt(z);
  // (*test)(1,1,0) =(Complex(0,4)*phase*pT*pow(r,2)*sqrt(z))/(m*(-1 + z)) + (Complex(0,4)*phase*pT*r*sqrt(z))/(m - m*pow(z,2));
  // (*test)(1,1,1) =(Complex(0,1)*(-3 + z)*sqrt(z))/(sqrt(2)*(1 + z)) + (Complex(0,2)*sqrt(2)*pow(r,2)*pow(1 + z,2))/((-1 + z)*sqrt(z)) + (Complex(0,1)*sqrt(2)*r*(-1 - 4*z + pow(z,2)))/((-1 + z)*sqrt(z));
  // (*test)(1,1,2) =(Complex(0,4)*pT*pow(r,2))/(phase*sqrt(z)*(m - m*z)) + (Complex(0,1)*pT*r*(3 + pow(z,2)))/(phase*m*(-1 + z)*sqrt(z)*(1 + z));

  
  // cerr << "testing in kernal " << " " << z << " " << r << " " << pT/m << " " << phi << "\n";
  // for(unsigned int ih1=0;ih1<2;++ih1) {
  //   for(unsigned int ih2=0;ih2<2;++ih2) {
  //     for(unsigned int ih3=0;ih3<3;++ih3) {
  // 	cerr << ih1 << " " << ih2 << " " << ih3 << " " << (*kernal)(ih1,ih2,ih3)  << "\n";
  // 	cerr << "testing diff " << abs((*kernal)(ih1,ih2,ih3)-(*test)(ih1,ih2,ih3)) << "\n";
  //     }
  //   }
  // }
  return kernal;
}
