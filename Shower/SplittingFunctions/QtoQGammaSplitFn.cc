// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGammaSplitFn class.
//

#include "QtoQGammaSplitFn.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtoQGammaSplitFn.tcc"
#endif


using namespace Herwig;

NoPIOClassDescription<QtoQGammaSplitFn> QtoQGammaSplitFn::initQtoQGammaSplitFn;
// Definition of the static class description member.

void QtoQGammaSplitFn::Init() {

  static ClassDocumentation<QtoQGammaSplitFn> documentation
    ("The QtoQGammaSplitFn class implements the splitting function for f->f gamma");

}

double QtoQGammaSplitFn::P(const double z, const Energy2 t,
			   const IdList & ids, const bool mass) const {
  double val=(1. + sqr(z))/(1.-z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=2.*sqr(m)/t;
  }
  double charge=getParticleData(ids[0])->iCharge()*3.;
  return sqr(charge)*val;; 
}

double QtoQGammaSplitFn::overestimateP(const double z, const IdList & ids) const {
  double charge=getParticleData(ids[0])->iCharge()*3.;
  return 2.*sqr(charge)/(1.-z); 
}

double QtoQGammaSplitFn::ratioP(const double z, const Energy2 t,
				const IdList & ids, const bool mass) const {
  double val=1. + sqr(z); 
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val -=2.*sqr(m)*(1.-z)/t;
  }
  return 0.5*val; 
}

double QtoQGammaSplitFn::integOverP(const double z) const {
  return -2.*log(1.-z); 
}

double QtoQGammaSplitFn::invIntegOverP(const double r) const {
  return 1. - exp(-r/2.); 
}

void QtoQGammaSplitFn::colourConnection(const ColinePair &parent,
					ColinePair &first,
					ColinePair &second) const {
  
  // Return immediately if the input is inconsistent.
  if((!parent.first && !parent.second) || (parent.first && parent.second))
    return;
  
  // second should be Gamma, doesn't get colour, of course. 
  first = parent;
  second = ColinePair();
}

bool QtoQGammaSplitFn::accept(const IdList & ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[1]||ids[2]!=ParticleID::gamma) return false;
  return getParticleData(ids[0])->charged();
}
  
  
