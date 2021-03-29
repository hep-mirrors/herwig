// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQ3P2SplitFn class.
//

#include "QtoQ3P2SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQ3P2SplitFn::pOver_ = 1.;

IBPtr QtoQ3P2SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQ3P2SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQ3P2SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*sqr(GeV2)) << n_ << fixedAlphaS_;
}

void QtoQ3P2SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*sqr(GeV2)) >> n_ >> fixedAlphaS_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQ3P2SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQ3P2SplitFn("Herwig::QtoQ3P2SplitFn", "HwOniumShower.so");

void QtoQ3P2SplitFn::Init() {

  static ClassDocumentation<QtoQ3P2SplitFn> documentation
    ("The QtoQ3P2SplitFn class implements the branching q-> q' 3P2");

  static Parameter<QtoQ3P2SplitFn,Energy5> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQ3P2SplitFn::O1_, GeV*GeV2*GeV2, 0.794*GeV*GeV2*GeV2, 0.0*GeV*GeV2*GeV2, 10.0*GeV*GeV2*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3P2SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQ3P2SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3P2SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQ3P2SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQ3P2SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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

double QtoQ3P2SplitFn::ratioP(const double z, const Energy2 t,
			       const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double W0 = z*(49.-68.*z+38.*sqr(z)-4.*z*sqr(z)+pow(z,4))/(6.*pow(1.+z,4));
  double W1 = (15.-165.*z+161.*sqr(z)-59.*z*sqr(z))/(3.*sqr(1.+z));
  double W2 = -4.*(28.-81.*z+11.*sqr(z))/(3.*(1.+z));
  double W3 = 160./3.;
  double ratio = (W0+r*(W1+r*(W2+r*W3)))/pOver_;
  // cerr << "testing values " << " " << z << " " << W0 << " " << W1 << " " << W2 << " " << W3 << "\n";
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQ3P2SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQ3P2SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQ3P2SplitFn::matrixElement(const double z, const Energy2 t, 
					  const IdList & ids, const double phi, bool) {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double rz=sqrt(z);
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t-sqr(m*(1.+z)));
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin2)));
  (*kernal)(0,0,0) = -ii*sqr(phase)*4.*r/rz/(1.-z)*(r*sqr(1.+z)/(1.-z)-z);
  (*kernal)(0,0,1) =  ii*phase*double(pT/m)*r/rz/(1.-z)*(-(1.+4.*z-sqr(z))/(1.+z)
							 + 4.*r*(1+3.*z)/(1.-z));
  (*kernal)(0,0,2) =  ii/sqrt(6.)/rz*(z*(5.-2.*z+sqr(z))/sqr(1.+z)-24.*z*r/(1.-z)
				      + sqr(r)*4.*(1+11.*z+11.*sqr(z)+pow(z,3))/sqr(1.-z));
  (*kernal)(0,0,3) = ii/phase*double(pT/m)*r*rz*4./(1.-z)*(1./(1.+z)-r*(3.+z)/(1.-z));
  (*kernal)(0,0,4) =-ii/sqr(phase)*4.*r/(1.-z)*rz*(r*sqr(1.+z)/(1.-z)-z);
  (*kernal)(0,1,0) =-ii*phase*double(pT/m)/rz*4.*r*(z/(1.+z)+r);
  (*kernal)(0,1,1) = ii/rz*(2.*(1.-z)*z/sqr(1.+z)+r*(1-5.*z)-4.*sqr(r)*(1.+z));
  (*kernal)(0,1,2) = ii*double(pT/m)/sqrt(6.)/phase*4.*r/rz*(r+z/(1.+z));
  (*kernal)(0,1,3) = 0.;
  (*kernal)(0,1,4) = 0.;
  (*kernal)(1,0,0) = 0.;  
  (*kernal)(1,0,1) = 0.;
  (*kernal)(1,0,2) =-ii*double(pT/m)/sqrt(6.)*phase*4.*r/rz*(r+z/(1.+z));
  (*kernal)(1,0,3) =+ii/rz*(2.*(1.-z)*z/sqr(1.+z)+r*(1.-5.*z)-4.*sqr(r)*(1.+z));
  (*kernal)(1,0,4) = ii/phase*double(pT/m)/rz*4.*r*(z/(1.+z)+r);
  (*kernal)(1,1,0) =-ii*sqr(phase)*4.*r/(1.-z)*rz*(r*sqr(1.+z)/(1.-z)-z);
  (*kernal)(1,1,1) = ii*phase*double(pT/m)*4.*r*rz/(1.-z)*(-1./(1.+z) + r*(3.+z)/(1.-z));
  (*kernal)(1,1,2) = ii/sqrt(6.)/rz* (z*(5.-2.*z+sqr(z))/sqr(1.+z) -24.*r*z/(1.-z)
				      +4.*sqr(r)*(1.+11.*z+11.*sqr(z)+pow(z,3))/sqr(1.-z));
  (*kernal)(1,1,3) = -ii/phase*double(pT/m)*r/rz/(1.-z)*(-(1.+4.*z-sqr(z))/(1.+z)+4.*r*(1.+3*z)/(1.-z));
  (*kernal)(1,1,4) =-ii/sqr(phase)*4.*r/rz/(1.-z)*(r*sqr(1.+z)/(1.-z)-z);



  // DecayMEPtr test(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin2)));
  // (*test)(0,0,0) = (Complex(0,-4)*sqr(phase)*r*sqrt(z))/(-1 + z) - (Complex(0,4)*sqr(phase)*pow(r,2)*pow(1 + z,2))/(pow(-1 + z,2)*sqrt(z))												  ;
  // (*test)(0,0,1) = (Complex(0,4)*phase*pT*pow(r,2)*(1 + 3*z))/(m*pow(-1 + z,2)*sqrt(z)) - (Complex(0,1)*phase*pT*r*(-1 - 4*z + pow(z,2)))/(m*(-1 + z)*sqrt(z)*(1 + z))								  ;
  // (*test)(0,0,2) = (Complex(0,4)*sqrt(6)*r*sqrt(z))/(-1 + z) + (Complex(0,1)*sqrt(z)*(5 - 2*z + pow(z,2)))/(sqrt(6)*pow(1 + z,2)) + (Complex(0,2)*sqrt(0.6666666666666666)*pow(r,2)*(1 + 11*z + 11*pow(z,2) + pow(z,3)))/(pow(-1 + z,2)*sqrt(z)) ;
  // (*test)(0,0,3) = (Complex(0,-4)*pT*pow(r,2)*sqrt(z)*(3 + z))/(phase*m*pow(-1 + z,2)) + (Complex(0,4)*pT*r*sqrt(z))/(phase*(m - m*pow(z,2)))											  ;
  // (*test)(0,0,4) = (Complex(0,-4)*r*pow(z,1.5))/(sqr(phase)*(-1 + z)) - (Complex(0,4)*pow(r,2)*sqrt(z)*pow(1 + z,2))/(sqr(phase)*pow(-1 + z,2))											  ;
  // (*test)(0,1,0) = (Complex(0,-4)*phase*pT*pow(r,2))/(m*sqrt(z)) - (Complex(0,4)*phase*pT*r*sqrt(z))/(m + m*z)																  ;
  // (*test)(0,1,1) = (Complex(0,-2)*(-1 + z)*sqrt(z))/pow(1 + z,2) - (Complex(0,4)*pow(r,2)*(1 + z))/sqrt(z) - (Complex(0,1)*r*(-1 + 5*z))/sqrt(z)													  ;
  // (*test)(0,1,2) = (Complex(0,2)*sqrt(0.6666666666666666)*pT*pow(r,2))/(phase*m*sqrt(z)) + (Complex(0,2)*sqrt(0.6666666666666666)*pT*r*sqrt(z))/(phase*(m + m*z))									  ;
  // (*test)(0,1,3) = 0																													  ;
  // (*test)(0,1,4) = 0																													  ;
  // (*test)(1,0,0) = 0																													  ;
  // (*test)(1,0,1) = 0																													  ;
  // (*test)(1,0,2) = (Complex(0,-2)*sqrt(0.6666666666666666)*phase*pT*pow(r,2))/(m*sqrt(z)) - (Complex(0,2)*sqrt(0.6666666666666666)*phase*pT*r*sqrt(z))/(m + m*z)									  ;
  // (*test)(1,0,3) = (Complex(0,-2)*(-1 + z)*sqrt(z))/pow(1 + z,2) - (Complex(0,4)*pow(r,2)*(1 + z))/sqrt(z) - (Complex(0,1)*r*(-1 + 5*z))/sqrt(z)													  ;
  // (*test)(1,0,4) = (Complex(0,4)*pT*pow(r,2))/(phase*m*sqrt(z)) + (Complex(0,4)*pT*r*sqrt(z))/(phase*(m + m*z))															  ;
  // (*test)(1,1,0) = (Complex(0,-4)*sqr(phase)*r*pow(z,1.5))/(-1 + z) - (Complex(0,4)*sqr(phase)*pow(r,2)*sqrt(z)*pow(1 + z,2))/pow(-1 + z,2)												  ;
  // (*test)(1,1,1) = (Complex(0,4)*phase*pT*pow(r,2)*sqrt(z)*(3 + z))/(m*pow(-1 + z,2)) + (Complex(0,4)*phase*pT*r*sqrt(z))/(m*(-1 + pow(z,2)))											  ;
  // (*test)(1,1,2) = (Complex(0,4)*sqrt(6)*r*sqrt(z))/(-1 + z) + (Complex(0,1)*sqrt(z)*(5 - 2*z + pow(z,2)))/(sqrt(6)*pow(1 + z,2)) + (Complex(0,2)*sqrt(0.6666666666666666)*pow(r,2)*(1 + 11*z + 11*pow(z,2) + pow(z,3)))/(pow(-1 + z,2)*sqrt(z)) ;
  // (*test)(1,1,3) = (Complex(0,-4)*pT*pow(r,2)*(1 + 3*z))/(phase*m*pow(-1 + z,2)*sqrt(z)) + (Complex(0,1)*pT*r*(-1 - 4*z + pow(z,2)))/(phase*m*(-1 + z)*sqrt(z)*(1 + z))								  ;
  // (*test)(1,1,4) = (Complex(0,-4)*r*sqrt(z))/(sqr(phase)*(-1 + z)) - (Complex(0,4)*pow(r,2)*pow(1 + z,2))/(sqr(phase)*pow(-1 + z,2)*sqrt(z))                                                                                               ;
  // cerr << "testing in kernal " << " " << z << " " << r << " " << pT/m << " " << phi << "\n";
  // for(unsigned int ih1=0;ih1<2;++ih1) {
  //   for(unsigned int ih2=0;ih2<2;++ih2) {
  //     for(unsigned int ih3=0;ih3<5;++ih3) {
  // 	cerr << ih1 << " " << ih2 << " " << ih3 << " " << (*kernal)(ih1,ih2,ih3)  << "\n";
  // 	cerr << "testing diff " << abs( (*kernal)(ih1,ih2,ih3) - (*test)(ih1,ih2,ih3) ) << "\n";
  //     }
  //   }
  // }
  return kernal;
}
