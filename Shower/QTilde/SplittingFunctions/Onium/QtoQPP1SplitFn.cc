// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQPP1SplitFn class.
//

#include "QtoQPP1SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQPP1SplitFn::pOver_ = 1.;

IBPtr QtoQPP1SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQPP1SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQPP1SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*sqr(GeV2)) << n_ << theta_ << sTheta_ << cTheta_ << fixedAlphaS_;
}

void QtoQPP1SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*sqr(GeV2)) >> n_ >> theta_ >> sTheta_ >> cTheta_ >> fixedAlphaS_;
}

void QtoQPP1SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  sTheta_ = sin(theta_);
  cTheta_ = cos(theta_);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQPP1SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQPP1SplitFn("Herwig::QtoQPP1SplitFn", "HwOniumShower.so");

void QtoQPP1SplitFn::Init() {

  static ClassDocumentation<QtoQPP1SplitFn> documentation
    ("The QtoQPP1SplitFn class implements the branching q-> q' P1");

  static Parameter<QtoQPP1SplitFn,Energy5> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQPP1SplitFn::O1_, GeV*GeV2*GeV2, 0.794*GeV*GeV2*GeV2, 0.0*GeV*GeV2*GeV2, 10.0*GeV*GeV2*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQPP1SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQPP1SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<QtoQPP1SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &QtoQPP1SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Parameter<QtoQPP1SplitFn,double> interfaceTheta
    ("Theta",
     "Mixing angle between the 1P1 and 3P1 states (in degrees)",
     &QtoQPP1SplitFn::theta_, 25., 0., 360.,
     false, false, Interface::limited);

}

void QtoQPP1SplitFn::guesstz(Energy2 t1,unsigned int iopt,
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
