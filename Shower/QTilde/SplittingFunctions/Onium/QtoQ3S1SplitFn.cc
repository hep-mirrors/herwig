// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQ3S1SplitFn class.
//

#include "QtoQ3S1SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQ3S1SplitFn::pOver_ = 3.;

IBPtr QtoQ3S1SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQ3S1SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQ3S1SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2) << n_ << fixedAlphaS_;

}

void QtoQ3S1SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2) >> n_ >> fixedAlphaS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQ3S1SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQ3S1SplitFn("Herwig::QtoQ3S1SplitFn", "HwOniumShower.so");

void QtoQ3S1SplitFn::Init() {

  static ClassDocumentation<QtoQ3S1SplitFn> documentation
    ("The QtoQ3S1SplitFn class implements the branching q-> q 3S1");

  static Parameter<QtoQ3S1SplitFn,Energy3> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQ3S1SplitFn::O1_, GeV*GeV2, 0.573*GeV*GeV2, 0.0*GeV*GeV2, 10.0*GeV*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3S1SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQ3S1SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<QtoQ3S1SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQ3S1SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void QtoQ3S1SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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

double QtoQ3S1SplitFn::ratioP(const double z, const Energy2 t,
			       const IdList & ids, const bool, const RhoDMatrix &) const {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double W0 = z*(17.-22.*z+9.*sqr(z))/sqr(1.+z);
  double W1 = 8*(3.-8.*z+sqr(z))/(1+z);
  double W2 = -48.;
  double ratio =(W0+r*W1+sqr(r)*W2)/pOver_;
  if(ratio>1.) cerr << "ratio greater than 1 in QtoQ3S1SplitFn " << ratio << "\n";
  if(ratio<0.) cerr << "ratio negative       in QtoQ3S1SplitFn " << ratio << "\n";
  return ratio;
}

DecayMEPtr QtoQ3S1SplitFn::matrixElement(const double z, const Energy2 t, 
			 const IdList & ids, const double phi, bool) {
  Energy m = ids[0]->mass();
  double r = sqr(m)/t;
  double rz=sqrt(z);
  double r2=sqrt(2.);
  Complex ii(0.,1.);
  Complex phase = exp(ii*phi);
  Energy pT = sqrt(z*(1.-z)*t-sqr(m*(1.+z)));
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  (*kernal)(0,0,0) = -2.*r2*phase*double(pT/m)*r/((1.-z)*rz);
  (*kernal)(0,0,1) = rz*((3.-z)/(1.+z)-8.*r/(1.-z));
  (*kernal)(0,0,2) = 2.*r2*double(pT/m)*r*rz/(phase*(1.-z));
  (*kernal)(0,1,0) = 2.*r2*(1.-z)/(1.+z)*(z + r*(1.+z))/rz;
  (*kernal)(0,1,1) = 0.;
  (*kernal)(0,1,2) = 0.;
  (*kernal)(1,0,0) = 0.;
  (*kernal)(1,0,1) = 0.;
  (*kernal)(1,0,2) = 2.*r2*(1.-z)/(1.+z)*(z + r*(1.+z))/rz;
  (*kernal)(1,1,0) = -2.*r2*r*phase*double(pT/m)*rz/(1.-z);  
  (*kernal)(1,1,1) = rz*((3.-z)/(1.+z)-8.*r/(1.-z));
  (*kernal)(1,1,2) = 4.*r*r2*0.5*double(pT/m)/(phase*(1.-z)*rz);
  return kernal;
}

