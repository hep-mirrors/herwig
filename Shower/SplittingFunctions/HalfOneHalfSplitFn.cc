// -*- C++ -*-
//
// HalfOneHalfSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfOneHalfSplitFn class.
//

#include "HalfOneHalfSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include <cassert>

using namespace Herwig;

NoPIOClassDescription<HalfOneHalfSplitFn> HalfOneHalfSplitFn::initHalfOneHalfSplitFn;
// Definition of the static class description member.

void HalfOneHalfSplitFn::Init() {

  static ClassDocumentation<HalfOneHalfSplitFn> documentation
    ("The HalfOneHalfSplitFn class implements the splitting function for q -> g q");

}

double HalfOneHalfSplitFn::P(const double z, const Energy2 t,
		       const IdList &ids, const bool mass) const {
  double val=(2.*(1.-z)+sqr(z))/z;
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=2.*sqr(m)/t;
  }
  return colourFactor()*val;
}

double HalfOneHalfSplitFn::overestimateP(const double z, const IdList &) const { 
  return 2.*colourFactor()/z; 
}

double HalfOneHalfSplitFn::ratioP(const double z, const Energy2 t,
			    const IdList &ids,const bool mass) const {
  double val=2.*(1.-z)+sqr(z);
  if(mass) {
    Energy m=getParticleData(ids[0])->mass();
    val -=2.*sqr(m)*z/t;
  }
  return 0.5*val;
}

double HalfOneHalfSplitFn::integOverP(const double z, const IdList & ,
				unsigned int PDFfactor) const { 
  switch(PDFfactor) {
  case 0:
    return 2.*colourFactor()*log(z); 
  case 1:
    return -2.*colourFactor()/z;
  case 2:
    return 2.*colourFactor()*log(z/(1.-z));
  case 3:
  default:
    throw Exception() << "HalfOneHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double HalfOneHalfSplitFn::invIntegOverP(const double r, const IdList & ,
				   unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return exp(0.5*r/colourFactor()); 
  case 1:
    return -2.*colourFactor()/r;
  case 2:
    return 1./(1.+exp(-0.5*r/colourFactor()));
  case 3:
  default:
    throw Exception() << "HalfOneHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool HalfOneHalfSplitFn::accept(const IdList &ids) const {
  // 3 particles and in and out fermion same
  if(ids.size()!=3 || ids[0]!=ids[2]) return false;
  tcPDPtr q=getParticleData(ids[0]);
  tcPDPtr g=getParticleData(ids[1]);
  if(q->iSpin()!=PDT::Spin1Half ||
     g->iSpin()!=PDT::Spin1) return false;
  return checkColours(ids);
}

