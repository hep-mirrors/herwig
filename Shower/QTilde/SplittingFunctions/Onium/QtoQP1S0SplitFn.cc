// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQP1S0SplitFn class.
//

#include "QtoQP1S0SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQP1S0SplitFn::pOver_ = 3.;

IBPtr QtoQP1S0SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQP1S0SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQP1S0SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2) << n_ << fixedAlphaS_;
}

void QtoQP1S0SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2) >> n_ >> fixedAlphaS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQP1S0SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQP1S0SplitFn("Herwig::QtoQP1S0SplitFn", "HwOniumShower.so");

void QtoQP1S0SplitFn::Init() {

  static ClassDocumentation<QtoQP1S0SplitFn> documentation
    ("The QtoQP1S0SplitFn class implements the branching q-> q' 1S0");

  static Parameter<QtoQP1S0SplitFn,Energy3> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQP1S0SplitFn::O1_, GeV*GeV2, 0.573*GeV*GeV2, 0.0*GeV*GeV2, 10.0*GeV*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQP1S0SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQP1S0SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<QtoQP1S0SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQP1S0SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQP1S0SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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
  Energy2 pre = 8./81.*aS2 * O1_ / M / sqr(a2); 
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

double QtoQP1S0SplitFn::ratioP(const double z, const Energy2 t,
			       const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m1 = ids[0]->mass();
  Energy M  = m1 + ids[1]->mass();
  double a1 = m1/M;
  double r = sqr(M)/t;
  double W0 = sqr(1.+(1.-a1)*(1.-z))*z/sqr(1.-a1*(1.-z));
  double W1 = (-2.*sqr(a1*z)+2.*sqr(a1)+a1*sqr(z)+2.*a1*z-7.*a1-3.*z+5)/(1.-a1*(1.-z));
  double W2 = -4.*a1*(1.-a1);
  double ratio =(W0+r*W1+sqr(r)*W2)/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQP1S0SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQP1S0SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQP1S0SplitFn::matrixElement(const double z, const Energy2 t, 
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
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  (*kernal)(0,0,0) = -(1.+(1.-a1)*(1.-z))*rz/(1.-a1*(1.-z)) - r*(1.-a1-a1*z)/rz;
  (*kernal)(1,1,0) = -(*kernal)(0,0,0);
  (*kernal)(0,1,0) = -double(pT/M)*r/rz/phase;
  (*kernal)(1,0,0) = conj((*kernal)(0,1,0));
  return kernal;
}
