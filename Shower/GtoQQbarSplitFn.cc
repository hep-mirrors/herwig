// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoQQbarSplitFun class.
//

#include "GtoQQbarSplitFn.h"

using namespace Herwig;

ClassDescription<GtoQQbarSplitFn> GtoQQbarSplitFn::initGtoQQbarSplitFn; 

GtoQQbarSplitFn::~GtoQQbarSplitFn() {}

double GtoQQbarSplitFn::P(const double z, const Energy2 qtilde2, 
			   const IdList &ids) {
  Energy m = CurrentGenerator::current().getParticleData(ids[1])->mass();
  double zz = z*(1.-z);
  double term = 2.*sqr(m)/zz/qtilde2; 
  double val = 1./2.*(1.-2.*zz+term);
  return val;
}

double GtoQQbarSplitFn::overestimateP(const double, const IdList &) {
  return 1./2.; 
}

double GtoQQbarSplitFn::integOverP(const double z) { 
  return z/2.; 
}


double GtoQQbarSplitFn::invIntegOverP(const double r) {
  return 2.*r; 
}


void GtoQQbarSplitFn::colourConnection(const ShoColinePair &parent,
				       ShoColinePair &first,
				       ShoColinePair &second) {

  // Return immediately if the input is inconsistent.
  if ( ! parent.first  ||  ! parent.second ) return;
  
  // Initialize
  first = second = ShoColinePair();

  // The first branching product is considered to be the quark 
  // and the second the anti-quark. 
  first.first = parent.first;
  second.second = parent.second;

}
