// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoGQSplitFn class.
//

#include "QtoGQSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtoGQSplitFn.tcc"
#endif


using namespace Herwig;

NoPIOClassDescription<QtoGQSplitFn> QtoGQSplitFn::initQtoGQSplitFn;
// Definition of the static class description member.

void QtoGQSplitFn::Init() {

  static ClassDocumentation<QtoGQSplitFn> documentation
    ("The QtoGQSplitFn class implements the splitting function for q -> g q");

}

double QtoGQSplitFn::P(const double z, const Energy2 t,
		       const IdList &ids, const bool mass) const {
  double val=(2.*(1.-z)+sqr(z))/z;
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=2.*sqr(m)/t;
  }
  return 4./3.*val;
}

double QtoGQSplitFn::overestimateP(const double z, const IdList &) const { 
  return 8./3./z; 
}

double QtoGQSplitFn::ratioP(const double z, const Energy2 t,
			    const IdList &ids,const bool mass) const {
  double val=2.*(1.-z)+sqr(z);
  if(mass) {
    Energy m=getParticleData(ids[0])->mass();
    val -=2.*sqr(m)*z/t;
  }
  return 0.5*val;
}

double QtoGQSplitFn::integOverP(const double z) const { return 8./3.*log(z); }

double QtoGQSplitFn::invIntegOverP(const double r) const {
  return exp(3.*r/8.); 
}

void QtoGQSplitFn::colourConnection(const ColinePair &parent,
				     ColinePair &first,
				     ColinePair &second) const {

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

bool QtoGQSplitFn::accept(const IdList &ids)
const {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[2]||ids[1]!=ParticleID::g) return false;
  tcPDPtr q=getParticleData(ids[0]);
  return q->iSpin()==PDT::Spin1Half&&(q->iColour()==PDT::Colour3||
				      q->iColour()==PDT::Colour3bar);
}


