// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoGQSplitFn class.
//

#include "QtoGQSplitFn.h"

using namespace Herwig;

ClassDescription<QtoGQSplitFn> QtoGQSplitFn::initQtoGQSplitFn;

QtoGQSplitFn::~QtoGQSplitFn() {}

double QtoGQSplitFn::P(const double z, const Energy2 qtilde2,
			const IdList &ids) {
  Energy m = CurrentGenerator::current().getParticleData(ids[0])->mass();
  Energy2 m2 = sqr(m); 
  // Lets switch everything from z->1-z
  //double val = 4./3.*(1. + sqr(z) - 2.*m2/(qtilde2*z))/(1.-z);
  double val = 4./3.*(2.*(1.-z)+sqr(z) - 2.*m2/(qtilde2*(1.-z)))/z;
  return val;
}

double QtoGQSplitFn::overestimateP(const double z, const IdList &) { 
  return 8./3./z; 
}

double QtoGQSplitFn::integOverP(const double z) { return 8./3.*log(z); }

double QtoGQSplitFn::invIntegOverP(const double r) {
  return exp(3.*r/8.); 
}


void QtoGQSplitFn::colourConnection(const ShoColinePair &parent,
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

