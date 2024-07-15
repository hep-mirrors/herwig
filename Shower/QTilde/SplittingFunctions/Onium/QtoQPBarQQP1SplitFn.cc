// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQPBarQQP1SplitFn class.
//

#include "QtoQPBarQQP1SplitFn.h"
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

const double QtoQPBarQQP1SplitFn::pOver_ = 6.;

IBPtr QtoQPBarQQP1SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQPBarQQP1SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQPBarQQP1SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(R02_,GeV*GeV2) << oenum(state_) << n_ << fixedAlphaS_;
}

void QtoQPBarQQP1SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(R02_,GeV*GeV2) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

void QtoQPBarQQP1SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  R02_ = params_->radialWaveFunctionSquared(state_,n_);
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<QtoQPBarQQP1SplitFn,Sudakov1to2FormFactor>
  describeHerwigQtoQPBarQQP1SplitFn("Herwig::QtoQPBarQQP1SplitFn", "HwOniumParameters.so HwOniumShower.so");

void QtoQPBarQQP1SplitFn::Init() {

  static ClassDocumentation<QtoQPBarQQP1SplitFn> documentation
    ("The QtoQPBarQQP1SplitFn class implements the splitting function for q -> qbar' (qq')_1");

  static Reference<QtoQPBarQQP1SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &QtoQPBarQQP1SplitFn::params_, false, false, true, false, false);
  
  static Parameter<QtoQPBarQQP1SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQPBarQQP1SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<QtoQPBarQQP1SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &QtoQPBarQQP1SplitFn::state_, cc, false, false);
  static SwitchOption interfaceStateccbar
    (interfaceState,
     "cc",
     "cc diquark",
     cc);
  static SwitchOption interfaceStatebb
    (interfaceState,
     "bb",
     "bb diquark",
     bb);
  static SwitchOption interfaceStatebc
    (interfaceState,
     "bc",
     "bc diquark",
     bc);

  static Parameter<QtoQPBarQQP1SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQPBarQQP1SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQPBarQQP1SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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
  Energy2 pre = 1./9.*aS2*R02_/sqr(a2)/M/Constants::pi;
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

double QtoQPBarQQP1SplitFn::ratioP(const double z, const Energy2 t,
				   const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M;
  double r = sqr(M)/t;
  double W0 = z*(6.+sqr(a1*(1.-z))+2.*a1*(1.-z)*(z-2.)+ z*(3.*z-8.))/sqr(1.-a1*(1.-z));
  double W1 = (3.+2.*sqr(a1)*(1.-z)*(z-3.)-9.*z+a1*(3.-2.*z+3.*sqr(z)))/(1.-a1*(1.-z));
  double W2 = -12.*(1.-a1)*a1;
  double ratio =(W0+r*W1+sqr(r)*W2)/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQPBarQQP1SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQPBarQQP1SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQPBarQQP1SplitFn::matrixElement(const double z, const Energy2 t, 
					      const IdList & ids, const double phi, bool) {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M, a2=1-a1;
  double r = sqr(M)/t;
  double rz=sqrt(z);
  double r2=sqrt(2.);
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t+sqr(M)*(sqr(a1)*z*(1.-z)-sqr(a2)*(1.-z)-z));
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  (*kernal)(0,1,0) = -r2*phase*double(pT/M)*r/((1 - z)*rz);
  (*kernal)(0,1,1) = (1.+(1.-a1)*(1.-z))*rz/(1.-a1*(1.-z)) - r*(1.+a1*(-2.+(-1.+2.*a1)*(1.-z))*(1.-z)+z)/((1 - z)*rz);
  (*kernal)(0,1,2) = (r2*double(pT/M)*r*rz)/(phase*(1 - z));
  (*kernal)(0,0,0) = r2*( (1.-z)*rz/(1.-a1*(1.-z)) + r*(1.-a1*(1.+z))/rz);
  (*kernal)(0,0,1) = (-1.+2.*a1)*double(pT/M)*r/(phase*rz);
  (*kernal)(0,0,2) = 0.;
  (*kernal)(1,1,0) = 0.;
  (*kernal)(1,1,1) = ((1.-2.*a1)*phase*double(pT/M)*r)/rz;
  (*kernal)(1,1,2) = r2*((1.-z)*rz/(1.-a1*(1.-z)) + r*(1.-a1*(1.+z))/rz); 
  (*kernal)(1,0,0) = -r2*r*phase*double(pT/M)*rz/(1 - z);
  (*kernal)(1,0,1) = r*(-1-a1*(-2.+(-1.+2.*a1)*(1.-z))*(1.-z) - z)/((1 - z)*rz) + (1.+(1.-a1)*(1.-z))*rz/(1.-a1*(1.-z));
  (*kernal)(1,0,2) = r*r2*double(pT/M)/(phase*(1.-z)*rz);
  return kernal;
}


