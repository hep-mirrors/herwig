// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGSplitFn class.
//

#include "QtoQGSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtoQGSplitFn.tcc"
#endif


using namespace Herwig;

QtoQGSplitFn::~QtoQGSplitFn() {}

NoPIOClassDescription<QtoQGSplitFn> QtoQGSplitFn::initQtoQGSplitFn;
// Definition of the static class description member.

void QtoQGSplitFn::Init() {

  static ClassDocumentation<QtoQGSplitFn> documentation
    ("The QtoQGSplitFn class implements the q -> qg splitting function");

}

double QtoQGSplitFn::P(const double z, const Energy2 qtilde2,
			const IdList &ids) {
  Energy m = getParticleData(ids[0])->mass();
  Energy2 m2 = sqr(m); 
  double val = 4./3.*(1. + sqr(z) - 2.*m2/(qtilde2*z))/(1.-z);
  return val;

}

double QtoQGSplitFn::overestimateP(const double z, const IdList &) { 
  return 8./3./(1.-z); 
}

double QtoQGSplitFn::ratioP(const double z, const Energy2 qtilde2,
			    const IdList &ids) {
  Energy m = getParticleData(ids[0])->mass();
  Energy2 m2 = sqr(m); 
  double val = 0.5*(1. + sqr(z) - 2.*m2/(qtilde2*z));
  return val;
} 

double QtoQGSplitFn::integOverP(const double z) { return -8./3.*log(1.-z); }

double QtoQGSplitFn::invIntegOverP(const double r) {
  return 1. - exp(- 3.*r/8.); 
}

void QtoQGSplitFn::colourConnection(const ColinePair &parent,
				    ColinePair &first,
				    ColinePair &second) {

  // Return immediately if the input is inconsistent.
  if ((!parent.first && !parent.second) || (parent.first && parent.second)) 
    return;
  
  // Initialize
  first = second = ColinePair();

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

bool QtoQGSplitFn::accept(const IdList &ids)
{
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[1]||ids[2]!=ParticleID::g) return false;
  tcPDPtr q=getParticleData(ids[0]);
  return q->iSpin()==PDT::Spin1Half&&(q->iColour()==PDT::Colour3||
				      q->iColour()==PDT::Colour3bar);
}
