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

QtoQGammaSplitFn::~QtoQGammaSplitFn() {}

NoPIOClassDescription<QtoQGammaSplitFn> QtoQGammaSplitFn::initQtoQGammaSplitFn;
// Definition of the static class description member.

void QtoQGammaSplitFn::Init() {

  static ClassDocumentation<QtoQGammaSplitFn> documentation
    ("The QtoQGammaSplitFn class implements the splitting function for f->f gamma");

}

double QtoQGammaSplitFn::P(const double z, const Energy2 qtilde2,
			   const IdList & ids) {
  Energy m = getParticleData(ids[0])->mass();
  double charge=getParticleData(ids[0])->iCharge()*3.;
  Energy2 m2 = sqr(m); 
  return sqr(charge)*(1. + sqr(z)- 2.*m2/(qtilde2*z))/(1.-z); 
}

double QtoQGammaSplitFn::overestimateP(const double z, const IdList & ids) {
  double charge=getParticleData(ids[0])->iCharge()*3.;
  return 2.*sqr(charge)/(1.-z); 
}

double QtoQGammaSplitFn::ratioP(const double z, const Energy2 qtilde2,
				const IdList & ids) {
  Energy m = getParticleData(ids[0])->mass();
  Energy2 m2 = sqr(m); 
  cerr << "testing QtoQGammaRatio " << P(z,qtilde2,ids)/overestimateP(z,ids) 
       << " " << 0.5*(1. + sqr(z)- 2.*m2/(qtilde2*z)) << endl;
  return 0.5*(1. + sqr(z)- 2.*m2/(qtilde2*z)); 
}

double QtoQGammaSplitFn::integOverP(const double z) {
  return -2.*log(1.-z); 
}

double QtoQGammaSplitFn::invIntegOverP(const double r) {
  return 1. - exp(-r/2.); 
}

void QtoQGammaSplitFn::colourConnection(const ColinePair &parent,
					ColinePair &first,
					ColinePair &second) {
  
  // Return immediately if the input is inconsistent.
  if((!parent.first && !parent.second) || (parent.first && parent.second))
    return;
  
  // second should be Gamma, doesn't get colour, of course. 
  first = parent;
  second = ColinePair();
}

bool QtoQGammaSplitFn::accept(const IdList & ids) {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[1]||ids[2]!=ParticleID::gamma) return false;
  return getParticleData(ids[0])->charged();
}
  
  
