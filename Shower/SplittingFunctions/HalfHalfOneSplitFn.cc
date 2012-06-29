// -*- C++ -*-
//
// HalfHalfOneSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfHalfOneSplitFn class.
//

#include "HalfHalfOneSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

NoPIOClassDescription<HalfHalfOneSplitFn> HalfHalfOneSplitFn::initHalfHalfOneSplitFn;
// Definition of the static class description member.

void HalfHalfOneSplitFn::Init() {

  static ClassDocumentation<HalfHalfOneSplitFn> documentation
    ("The HalfHalfOneSplitFn class implements the q -> qg splitting function");

}

double HalfHalfOneSplitFn::P(const double z, const Energy2 t,
			     const IdList &ids, const bool mass) const {
  double val = (1. + sqr(z))/(1.-z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();  
    val -= 2.*sqr(m)/t;
  }
  return colourFactor()*val;
}

double HalfHalfOneSplitFn::overestimateP(const double z, const IdList &) const { 
  return 2.*colourFactor()/(1.-z); 
}

double HalfHalfOneSplitFn::ratioP(const double z, const Energy2 t,
				  const IdList &ids, const bool mass) const {
  double val = 1. + sqr(z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val -= 2.*sqr(m)*(1.-z)/t;
  } 
  return 0.5*val;
} 

double HalfHalfOneSplitFn::integOverP(const double z, const IdList & ,
				      unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return -2.*colourFactor()*log(1.-z);
  case 1:
    return  2.*colourFactor()*log(z/(1.-z));
  case 2:
    return  2.*colourFactor()/(1.-z);
  case 3:
  default:
    throw Exception() << "HalfHalfOneSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

double HalfHalfOneSplitFn::invIntegOverP(const double r, const IdList & ,
				   unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return 1. - exp(- 0.5*r/colourFactor()); 
  case 1:
    return 1./(1.-exp(-0.5*r/colourFactor()));
  case 2:
    return 1.-2.*colourFactor()/r;
  case 3:
  default:
    throw Exception() << "HalfHalfOneSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

bool HalfHalfOneSplitFn::accept(const IdList &ids) const {
  // 3 particles and in and out fermion same
  if(ids.size()!=3 || ids[0]!=ids[1]) return false;
  tcPDPtr q=getParticleData(ids[0]);
  tcPDPtr g=getParticleData(ids[2]);
  if(q->iSpin()!=PDT::Spin1Half ||
     g->iSpin()!=PDT::Spin1) return false;
  return checkColours(ids);
}
