// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQ3P0SplitFn class.
//

#include "QtoQ3P0SplitFn.h"
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

const double QtoQ3P0SplitFn::pOver_ = 2.;

IBPtr QtoQ3P0SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQ3P0SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQ3P0SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*sqr(GeV2)) << oenum(state_) << n_ << fixedAlphaS_;
}

void QtoQ3P0SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*sqr(GeV2)) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

void QtoQ3P0SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletMEProduction<1>(state_,n_,1,0);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQ3P0SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQ3P0SplitFn("Herwig::QtoQ3P0SplitFn", "HwOniumShower.so HwOniumParameters.so");

void QtoQ3P0SplitFn::Init() {

  static ClassDocumentation<QtoQ3P0SplitFn> documentation
    ("The QtoQ3P0SplitFn class implements the branching q-> q 3P0");

  static Reference<QtoQ3P0SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &QtoQ3P0SplitFn::params_, false, false, true, false, false);
  
  static Parameter<QtoQ3P0SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQ3P0SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<QtoQ3P0SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &QtoQ3P0SplitFn::state_, ccbar, false, false);
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
  
  static Parameter<QtoQ3P0SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQ3P0SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQ3P0SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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

double QtoQ3P0SplitFn::ratioP(const double z, const Energy2 t,
			       const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double W0 = z*sqr(5.-2.*z+sqr(z))/(12.*pow(1+z,4));
  double W1 = 2.*(33.-12.*z-7.*sqr(z)-2.*pow(z,3))/(3.*sqr(1+z));
  double W2 = 4.*(-23.+3.*z+2.*sqr(z))/(3.*(1.+z));
  double W3 = 32./3.;
  double ratio = (W0+r*(W1+r*(W2+r*W3)))/pOver_;
  // cerr << "testing values " << " " << z << " " << W0 << " " << W1 << " " << W2 << " " << W3 << "\n";
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQ3P0SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQ3P0SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQ3P0SplitFn::matrixElement(const double z, const Energy2 t, 
					  const IdList & ids, const double phi, bool) {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t-sqr(m*(1.+z)));
  // calculate the kernal N.B. prefactor 1./sqrt(3.)/sqrt(z) removed
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  (*kernal)(0,0,0) = 0.5*z*(z*(z-2.)+5.)/sqr(1.+z)+6.*r-4.*sqr(r)*(1.+z);
  (*kernal)(1,1,0) = (*kernal)(0,0,0);
  (*kernal)(0,1,0) = 2.*double(pT/m)*r/phase*((z+3.)/(1.+z)-2.*r);
  (*kernal)(1,0,0) = -conj((*kernal)(0,1,0));
  return kernal;
}
