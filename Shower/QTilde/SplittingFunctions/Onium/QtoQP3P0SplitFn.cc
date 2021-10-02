// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQP3P0SplitFn class.
//

#include "QtoQP3P0SplitFn.h"
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

const double QtoQP3P0SplitFn::pOver_ = 1.;

IBPtr QtoQP3P0SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQP3P0SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQP3P0SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*sqr(GeV2)) << oenum(state_) << n_ << fixedAlphaS_;
}

void QtoQP3P0SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*sqr(GeV2)) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

void QtoQP3P0SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletMEProduction<1>(state_,n_,1,0);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQP3P0SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQP3P0SplitFn("Herwig::QtoQP3P0SplitFn", "HwOniumParameters.so HwOniumShower.so");

void QtoQP3P0SplitFn::Init() {

  static ClassDocumentation<QtoQP3P0SplitFn> documentation
    ("The QtoQP3P0SplitFn class implements the branching q-> q' 3P0");

  static Reference<QtoQP3P0SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &QtoQP3P0SplitFn::params_, false, false, true, false, false);
  
  static Parameter<QtoQP3P0SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQP3P0SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<QtoQP3P0SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &QtoQP3P0SplitFn::state_, ccbar, false, false);
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
  
  static Parameter<QtoQP3P0SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQP3P0SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQP3P0SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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
  Energy2 pre = 32./81.*aS2*O1_/pow(a2,4)/M/sqr(M);
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

double QtoQP3P0SplitFn::ratioP(const double z, const Energy2 t,
			       const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M;
  double r = sqr(M)/t;
  double W0 = z*sqr( 3.*(-2.+z) +
		     a1*(15.-12.*z+sqr(z) +
			 a1*((-13.+18.*z-5.*sqr(z)) +
			     4.*a1*sqr(1.-z))))/(48.*sqr(a1*sqr(1.-a1*(1.-z))));
  double W1 = (45.-27*z + a1*(6.*(6.-13.*z+5.*sqr(z)) +
			      a1*(-407.+591.*z-233.*sqr(z)+9.*z*sqr(z) +
				  a1*(590.-998.*z+498.*sqr(z)-58.*z*sqr(z) +
				      8.*a1*(-41.+84.*z-53.*sqr(z)+10.*z*sqr(z) +
					     4.*a1*sqr(1.-z)*(2.-z))))))/(48.*sqr(a1)*sqr(1 - a1*(1 - z)));
  double W2 = (1.-a1)*(-3.+a1*((-22.+12.*z) +
			       a1*((41.-30.*z+3.*sqr(z)) +
				   2.*a1*(-8.+9.*z-sqr(z)))))/(6.*a1*(1 - a1*(1 - z)));
  double W3 = 4./3.*sqr(1.-a1)*a1;
  double ratio = (W0+r*(W1+r*(W2+r*W3)))/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQP3P0SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQP3P0SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQP3P0SplitFn::matrixElement(const double z, const Energy2 t, 
					  const IdList & ids, const double phi, bool) {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M, a2=1-a1;
  double r = sqr(M)/t;
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t+sqr(M)*(sqr(a1)*z*(1.-z)-sqr(a2)*(1.-z)-z));
  // calculate the kernal N.B. prefactor 1./4./sqrt(3.)/sqrt(z) removed
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  (*kernal)(0,0,0) = z*(6.-3.*z-a1*(15.-a1*(13.-5.*z)*(1.-z)+4.*sqr(a1)*sqr(1.-z)-(12.-z)*z))/(a1*sqr(1.-a1*(1.-z)))
    +r*(3.+(1.-2.*a1)*a1*(7.-4.*a1*(1.-z)-5.*z))/a1 - 8.*(1.-a1)*a1*sqr(r)*(1.-a1*(1.-z));
  (*kernal)(1,1,0) = (*kernal)(0,0,0);
  (*kernal)(0,1,0) = double(pT/M)*r/phase*(-8.*(1.-a1)*a1*r+(3.+a1*(7.-2.*a1*(9.-4.*a1*(1.-z)-5.*z)-z))/(a1*(1.-a1*(1.-z))));
  (*kernal)(1,0,0) = -conj((*kernal)(0,1,0));
  return kernal;
}
