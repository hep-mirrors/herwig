// -*- C++ -*-
//
// QtoQGammaSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGammaSplitFn class.
//

#include "QtoQGammaSplitFn.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

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

double QtoQGammaSplitFn::integOverP(const double z, unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return -2.*log(1.-z); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "QtoQGammaSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double QtoQGammaSplitFn::invIntegOverP(const double r, unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 1. - exp(-0.5*r); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "QtoQGammaSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

void QtoQGammaSplitFn::colourConnection(tShowerParticlePtr parent,
					tShowerParticlePtr first,
					tShowerParticlePtr ,
					const bool back) const {
  if(!back) {
    ColinePair cparent = ColinePair(parent->colourLine(), 
				    parent->antiColourLine());
    // ensure input consistency
    assert(( cparent.first && !cparent.second) || 
	   (!cparent.first &&  cparent.second));
    // q -> q gamma
    if(cparent.first) cparent.first->addColoured(first);
    // qbar -> qbar gamma
    else              cparent.second->addAntiColoured(first);
  }
  else {
    ColinePair cfirst = ColinePair(first->colourLine(), 
				   first->antiColourLine());
    // ensure input consistency
    assert(( cfirst.first && !cfirst.second) ||
	   (!cfirst.first &&  cfirst.second));
    // q -> q gamma
    if(cfirst.first) cfirst.first->addColoured(parent);
    // qbar -> qbar gamma
    else             cfirst.second->addAntiColoured(parent);
  }
}

bool QtoQGammaSplitFn::accept(const IdList & ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[1]||ids[2]!=ParticleID::gamma) return false;
  return getParticleData(ids[0])->charged();
}
  
  
