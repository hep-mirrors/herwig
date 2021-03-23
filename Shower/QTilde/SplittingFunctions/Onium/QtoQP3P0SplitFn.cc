// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQP3P0SplitFn class.
//

#include "QtoQP3P0SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const double QtoQP3P0SplitFn::pOver_ = 1.;

IBPtr QtoQP3P0SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQP3P0SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQP3P0SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*sqr(GeV2)) << n_ << fixedAlphaS_;
}

void QtoQP3P0SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*sqr(GeV2)) >> n_ >> fixedAlphaS_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQP3P0SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQP3P0SplitFn("Herwig::QtoQP3P0SplitFn", "HwOniumShower.so");

void QtoQP3P0SplitFn::Init() {

  static ClassDocumentation<QtoQP3P0SplitFn> documentation
    ("The QtoQP3P0SplitFn class implements the branching q-> q' 3P0");

  static Parameter<QtoQP3P0SplitFn,Energy5> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQP3P0SplitFn::O1_, GeV*GeV2*GeV2, 0.794*GeV*GeV2*GeV2, 0.0*GeV*GeV2*GeV2, 10.0*GeV*GeV2*GeV2,
     false, false, Interface::limited);
  
  static Parameter<QtoQP3P0SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &QtoQP3P0SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
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
