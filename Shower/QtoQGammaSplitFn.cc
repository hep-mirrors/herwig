// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGammaSplitFun class.
//

#include "QtoQGammaSplitFn.h"

using namespace Herwig;

ClassDescription<QtoQGammaSplitFn> QtoQGammaSplitFn::initQtoQGammaSplitFn;

QtoQGammaSplitFn::~QtoQGammaSplitFn() {}

double QtoQGammaSplitFn::P(const double z, const Energy2 qtilde2,
			   const IdList &) {
  return (1. + sqr(z))/(1.-z); 
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
