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

const double QtoQP1S0SplitFn::pOver_ = 1.;

IBPtr QtoQP1S0SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQP1S0SplitFn::fullclone() const {
  return new_ptr(*this);
}

void QtoQP1S0SplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2);
}

void QtoQP1S0SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QtoQP1S0SplitFn,Sudakov1to2FormFactor>
describeHerwigQtoQP1S0SplitFn("Herwig::QtoQP1S0SplitFn", "HwOniumShower.so");

void QtoQP1S0SplitFn::Init() {

  static ClassDocumentation<QtoQP1S0SplitFn> documentation
    ("The QtoQP1S0SplitFn class implements the branching q-> q' M");

  static Parameter<QtoQP1S0SplitFn,Energy3> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &QtoQP1S0SplitFn::O1_, GeV*GeV2, 0.573*GeV*GeV2, 0.0*GeV*GeV2, 10.0*GeV*GeV2,
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
  Energy2 pre = 8./81.*sqr(alpha()->overestimateValue()) * O1_ / M / sqr(a2); 
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
