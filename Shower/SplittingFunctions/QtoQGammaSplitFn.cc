// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGammaSplitFn class.
//

#include "QtoQGammaSplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtoQGammaSplitFn.tcc"
#endif


using namespace Herwig;

QtoQGammaSplitFn::~QtoQGammaSplitFn() {}

NoPIOClassDescription<QtoQGammaSplitFn> QtoQGammaSplitFn::initQtoQGammaSplitFn;
// Definition of the static class description member.

void QtoQGammaSplitFn::Init() {

  static ClassDocumentation<QtoQGammaSplitFn> documentation
    ("There is no documentation for the QtoQGammaSplitFn class");

}

double QtoQGammaSplitFn::P(const double z, const Energy2 qtilde2,
			   const IdList & ids) {
  Energy m = CurrentGenerator::current().getParticleData(ids[0])->mass();
  Energy2 m2 = sqr(m); 
  return (1. + sqr(z)- 2.*m2/(qtilde2*z))/(1.-z); 
}

double QtoQGammaSplitFn::overestimateP(const double z, const IdList &) {
  return 2./(1.-z); 
}


double QtoQGammaSplitFn::integOverP(const double z) {
  return -2.*log(1.-z); 
}


double QtoQGammaSplitFn::invIntegOverP(const double r) {
  return 1. - exp(-r/2.); 
}



void QtoQGammaSplitFn::colourConnection(const ShoColinePair &parent,
					ShoColinePair &first,
					ShoColinePair &second) {

  // Return immediately if the input is inconsistent.
  if((!parent.first && !parent.second) || (parent.first && parent.second))
    return;
  
  // second should be Gamma, doesn't get colour, of course. 
  first = parent;
  second = ShoColinePair();
}

