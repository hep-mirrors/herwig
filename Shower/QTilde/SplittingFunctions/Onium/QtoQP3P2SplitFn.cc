// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQP3P2SplitFn class.
//

#include "QtoQP3P2SplitFn.h"
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

const double QtoQP3P2SplitFn::pOver_ = 2.;

IBPtr QtoQP3P2SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQP3P2SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQP3P2SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*sqr(GeV2)) << oenum(state_) << n_ << fixedAlphaS_;
}

void QtoQP3P2SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*sqr(GeV2)) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

void QtoQP3P2SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletMEProduction<1>(state_,n_,1,2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQP3P2SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQP3P2SplitFn("Herwig::QtoQP3P2SplitFn", "HwOniumParameters.so HwOniumShower.so");

void QtoQP3P2SplitFn::Init() {

  static ClassDocumentation<QtoQP3P2SplitFn> documentation
    ("The QtoQP3P2SplitFn class implements the branching q-> q' 3P2");

  static Reference<QtoQP3P2SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &QtoQP3P2SplitFn::params_, false, false, true, false, false);
  
  static Parameter<QtoQP3P2SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQP3P2SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<QtoQP3P2SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &QtoQP3P2SplitFn::state_, ccbar, false, false);
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
  
  static Parameter<QtoQP3P2SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQP3P2SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQP3P2SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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

double QtoQP3P2SplitFn::ratioP(const double z, const Energy2 t,
			       const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M;
  double r = sqr(M)/t;
  double W0 = z*(30.-72.*z+69*sqr(z)-30.*z*sqr(z)+5.*sqr(sqr(z)) 
		 +a1*( -4.*(18.-51.*z+55.*sqr(z)-27*z*sqr(z)+5.*sqr(sqr(z)))
		       +a1*( +4.*sqr(1.-z)*(14.-17.*z+6.*sqr(z))
			     +a1*(-8.*pow(1.-z,3)*(2.-z) + 2.*a1*pow(1.-z,4) ))))/(12.*pow(1.-a1*(1.-z),4));
  double W1 = (5.-z*(93.-98*z+30.*sqr(z))
	       +a1*(+ 2.*(5.+46.*z-57.*sqr(z)+14.*z*sqr(z))
		    + a1*( -35.+45.*z-17.*sqr(z)+7.*z*sqr(z)
			   +4.*a1*sqr(1.-z)*(5.-z))))/(12.*sqr(1.-a1*(1.-z)));
  double W2 = (1.-a1)*(-11.+45.*z + a1*(3.*z*(4.-5.*z)-23. + 2.*a1*(17.-4.*z)*(1.-z)))/(6.*(1.-a1*(1.-z)));
  double W3 = 20.*sqr(1.-a1)*a1/3.;
  double ratio = (W0+r*(W1+r*(W2+r*W3)))/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQP3P2SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQP3P2SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQP3P2SplitFn::matrixElement(const double z, const Energy2 t, 
					  const IdList & ids, const double phi, bool) {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M, a2=1-a1;
  double r = sqr(M)/t;
  double rz=sqrt(z);
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t+sqr(M)*(sqr(a1)*z*(1.-z)-sqr(a2)*(1.-z)-z));
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin2)));
  (*kernal)(0,0,0) = -ii*sqr(phase)*r/rz/(1.-z)*(r*sqr(1.-a1*(1.-z))/(1.-z)-z);  
  (*kernal)(0,0,1) =  0.5*ii*phase*double(pT/M)*r/rz/(1.-z)*(-(1.+2.*z-sqr(z)+a1*(-1.+sqr(z)))/(1.-a1*(1.-z))
							     + r*(3.+2.*sqr(a1)*sqr(1.-z)+z-a1*(5.-6.*z+sqr(z)))/(1.-z));
  (*kernal)(0,0,2) =  ii/sqrt(6.)/rz*(z*(3.+sqr(a1)*sqr(1.-z)-3.*z+sqr(z)-2.*a1*(2.-3.*z+sqr(z)))/sqr(1.-a1*(1.-z))
				      -r*(1 + 2*sqr(a1)*sqr(1.-z) + 6*z - sqr(z) + a1*(-3 + 2*z + sqr(z)))/(1.-z)
				      + sqr(r)*(-sqr(a1)*(-11.+z)*sqr(1.-z)-4.*pow(a1,3)*pow(1.-z,3)+3.*(1.+z)+2.*a1*(-5.+4.*z+sqr(z)))/sqr(1.-z));
  (*kernal)(0,0,3) = -0.5*ii/phase*double(pT/M)*r/rz/(1.-z)*( -(3.-2*a1*(1.-z)-z)*z/((1.-a1*(1.-z)))
							      + r*(1.+2.*sqr(a1)*sqr(1.-z)+3.*z+a1*(-3.+2.*z+sqr(z)))/(1.-z));
  (*kernal)(0,0,4) =-ii/sqr(phase)*r/(1.-z)*rz*( r*sqr(1.-a1*(1.-z))/(1.-z) - z);
  (*kernal)(0,1,0) =-ii*phase*double(pT/M)/rz*r*( z/(1.-a1*(1.-z)) + r*(1.-a1-a1*z)/(1.-z));  
  (*kernal)(0,1,1) = 0.5*ii/rz*((2.-2.*a1*(1.-z)-z)*(1.-z)*z/sqr(1.-a1*(1.-z))
				+r*(1.-a1-z-3.*a1*z)
				+sqr(r)*(-3.-2.*a1*(-4.+z)+2.*pow(a1,3)*sqr(1.-z)+sqr(a1)*(-7.+6.*z+sqr(z)))/(1.-z));
  (*kernal)(0,1,2) = ii*double(pT/M)/sqrt(6.)/phase*r/rz*(r*(3.-a1*(7.-z)+4.*sqr(a1)*(1.-z))/(1.-z)
							  +(-1.+a1*(3.-z)-2.*sqr(a1)*(1.-z)+z)/(1.-a1*(1.-z)));
  (*kernal)(0,1,3) = 0.5*ii*(1.-2.*a1)/sqr(phase)*r/rz*(r*sqr(1.-a1*(1.-z))/(1.-z)-z);
  (*kernal)(0,1,4) = 0.;
  (*kernal)(1,0,0) = 0.;  
  (*kernal)(1,0,1) = 0.5*ii*(1.-2.*a1)*sqr(phase)*r/rz*(r*sqr(1.-a1*(1.-z))/(1.-z)-z);
  (*kernal)(1,0,2) =-ii*double(pT/M)/sqrt(6.)*phase*r/rz*(r*(3.-a1*(7.-z)+4.*sqr(a1)*(1.-z))/(1.-z)
							  +(-1.+a1*(3.-z)-2.*sqr(a1)*(1.-z)+z)/(1.-a1*(1.-z)));
  (*kernal)(1,0,3) =+0.5*ii/rz*((2.-2.*a1*(1.-z)-z)*(1.-z)*z/sqr(1.-a1*(1.-z))
				+r*(1.-a1-z-3.*a1*z)
				+sqr(r)*(-3.-2.*a1*(-4.+z)+2.*pow(a1,3)*sqr(1.-z)+sqr(a1)*(-7.+6.*z+sqr(z)))/(1.-z));
  (*kernal)(1,0,4) = ii/phase*double(pT/M)/rz*r*( z/(1.-a1*(1.-z)) + r*(1.-a1-a1*z)/(1.-z));
  (*kernal)(1,1,0) =-ii*sqr(phase)*r/(1.-z)*rz*( r*sqr(1.-a1*(1.-z))/(1.-z) - z);
  (*kernal)(1,1,1) = 0.5*ii*phase*double(pT/M)*r/rz/(1.-z)*(-(3.-2.*a1*(1.-z)-z)*z/(1.-a1*(1.-z))
							    + r*(1.+2.*sqr(a1)*sqr(1.-z)+3.*z+a1*(-3.+2.*z+sqr(z)))/(1.-z));
  (*kernal)(1,1,2) = ii/sqrt(6.)/rz* (z*(3.+sqr(a1)*sqr(1.-z)-3.*z+sqr(z)-2.*a1*(2.-3.*z+sqr(z)))/sqr(1.-a1*(1.-z))
				      -r*(1.+2.*sqr(a1)*sqr(1.-z)+6.*z-sqr(z)+a1*(-3.+2.*z+sqr(z)))/(1.-z)
				      +sqr(r)*(-(sqr(a1)*(-11.+z)*sqr(1.-z))-4.*pow(a1,3)*pow(1.-z,3)+3.*(1.+z)+2.*a1*(-5.+4.*z+sqr(z)))/(sqr(1.-z)));
  (*kernal)(1,1,3) = -0.5*ii/phase*double(pT/M)*r/rz/(1.-z)*(-(1.+2.*z-sqr(z)+a1*(-1.+sqr(z)))/(1.-a1*(1.-z))
							     +r*(3.+2.*sqr(a1)*sqr(1.-z)+z-a1*(5.-6.*z+sqr(z)))/(1.-z));
  (*kernal)(1,1,4) =-ii/sqr(phase)*r/rz/(1.-z)*(r*sqr(1.-a1*(1.-z))/(1.-z)-z);
  return kernal;
}
