// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGSplitFn class.
//

#include "QtoQGSplitFn.h"

using namespace Herwig;

ClassDescription<QtoQGSplitFn> QtoQGSplitFn::initQtoQGSplitFn;

QtoQGSplitFn::~QtoQGSplitFn() {}

double QtoQGSplitFn::P(const double z, const Energy2 qtilde2,
			const IdList &ids) {
  Energy m = CurrentGenerator::current().getParticleData(ids[0])->mass();
  Energy2 m2 = sqr(m); 
  double val = 4./3.*(1. + sqr(z) - 2.*m2/(qtilde2*z))/(1.-z);
  return val;

}

double QtoQGSplitFn::overestimateP(const double z, const IdList &) { 
  return 8./3./(1.-z); 
}

double QtoQGSplitFn::integOverP(const double z) { return -8./3.*log(1.-z); }

double QtoQGSplitFn::invIntegOverP(const double r) {
  return 1. - exp(- 3.*r/8.); 
}


void QtoQGSplitFn::colourConnection(const ShoColinePair &parent,
				     ShoColinePair &first,
				     ShoColinePair &second) {

  // Return immediately if the input is inconsistent.
  if ((!parent.first && !parent.second) || (parent.first && parent.second)) 
    return;
  
  // Initialize
  first = second = ShoColinePair();

  // The first branching product is considered to be the quark 
  // and the second the gluon. The colour line of the parent
  // is one of the two colour lines of the gluon, whereas the
  // other one of the latter is a new colour line which is
  // also share by the first product (the quark).
  if(parent.first) { // the parent is a quark
    second.first = parent.first;
    first.first = second.second = new_ptr(ColourLine());
  } else if(parent.second) { // the parent is an antiquark
    second.second = parent.second;
    first.second = second.first = new_ptr(ColourLine());
  }
}

