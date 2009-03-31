// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGammaSplitFn class.
//

#include "QtoQGammaSplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

IBPtr QtoQGammaSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr QtoQGammaSplitFn::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<QtoQGammaSplitFn> QtoQGammaSplitFn::initQtoQGammaSplitFn;
// Definition of the static class description member.

void QtoQGammaSplitFn::Init() {

  static ClassDocumentation<QtoQGammaSplitFn> documentation
    ("The QtoQGammaSplitFn class implements the splitting"
     " function for q -> q gamma");

}

double QtoQGammaSplitFn::P(const double z, const Energy2 t,
		       const IdList &ids, const bool mass) const {
  tPDPtr part=getParticleData(ids[0]);
  double charge2 = sqr(double(part->iCharge())/3.);
  double val = (1. + sqr(z))/(1.-z);
  if(mass) {
    Energy m = part->mass();  
    val -= 2.*sqr(m)/t;
  }
  return charge2*val;
}

double QtoQGammaSplitFn::overestimateP(const double z, const IdList & ids) const { 
  double charge2 = sqr(double(getParticleData(ids[0])->iCharge())/3.);
  return 2.*charge2/(1.-z); 
}

double QtoQGammaSplitFn::ratioP(const double z, const Energy2 t,
			    const IdList &ids, const bool mass) const {
  double val = 1. + sqr(z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val -= 2.*sqr(m)*(1.-z)/t;
  } 
  return 0.5*val;
}

double QtoQGammaSplitFn::integOverP(const double z, const IdList & ids,
				unsigned int PDFfactor) const {
  double charge2 = sqr(double(getParticleData(ids[0])->iCharge())/3.);
  switch (PDFfactor) {
  case 0:
    return -2.*charge2*log(1.-z);
  case 1:
    return 2.*charge2*log(z/(1.-z));
  case 2:
    return 2.*charge2/(1.-z);
  case 3:
  default:
    throw Exception() << "QtoQGammaSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

double QtoQGammaSplitFn::invIntegOverP(const double r, const IdList & ids,
				       unsigned int PDFfactor) const {
  double charge2 = sqr(double(getParticleData(ids[0])->iCharge())/3.);
  switch (PDFfactor) {
  case 0:
    return 1. - exp(-0.5*r/charge2); 
  case 1:
    return 1./(1.-exp(-0.5*r/charge2));
  case 2:
    return 1.-0.5/r/charge2;
  case 3:
  default:
    throw Exception() << "QtoQGammaSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

void QtoQGammaSplitFn::colourConnection(tShowerParticlePtr parent,
				    tShowerParticlePtr first,
				    tShowerParticlePtr ,
				    const bool back) const {
  if(!parent->data().coloured()) return;
  if(!back) {
    ColinePair cparent = ColinePair(parent->colourLine(), 
				    parent->antiColourLine());
    // ensure input consistency
    assert((!cparent.first &&  cparent.second) || 
	   ( cparent.first && !cparent.second));
    // q -> q g
    if(cparent.first) {
      cparent.first->addColoured(first);
    }
    // qbar -> qbar g
    else {
      cparent.second->addAntiColoured(first);
    }
  }
  else {
    ColinePair cfirst = ColinePair(first->colourLine(), 
				   first->antiColourLine());
    // ensure input consistency
    assert(( cfirst.first && !cfirst.second) ||
	   (!cfirst.first &&  cfirst.second)); 
    // q -> q g
    if(cfirst.first) {
      cfirst.first->addColoured(parent);
    }
    // qbar -> qbar g
    else {
      cfirst.second->addAntiColoured(parent);
    }
  }
}

bool QtoQGammaSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[1]||ids[2]!=ParticleID::gamma) return false;
  tcPDPtr q=getParticleData(ids[0]);
  return q->iSpin() == PDT::Spin1Half  &&  q->iCharge() != 0;
}
