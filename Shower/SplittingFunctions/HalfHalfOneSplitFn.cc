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
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeNoPIOClass<HalfHalfOneSplitFn,Herwig::SplittingFunction>
describeHalfHalfOneSplitFn ("Herwig::HalfHalfOneSplitFn","HwShower.so");

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
  return colourFactor(ids)*val;
}

double HalfHalfOneSplitFn::overestimateP(const double z,
					 const IdList & ids) const { 
  return 2.*colourFactor(ids)/(1.-z); 
}

double HalfHalfOneSplitFn::ratioP(const double z, const Energy2 t,
				  const IdList & ids, const bool mass) const {
  double val = 1. + sqr(z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val -= 2.*sqr(m)*(1.-z)/t;
  } 
  return 0.5*val;
} 

double HalfHalfOneSplitFn::integOverP(const double z,
				      const IdList & ids,
				      unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return -2.*colourFactor(ids)*Math::log1m(z);
  case 1:
    return  2.*colourFactor(ids)*log(z/(1.-z));
  case 2:
    return  2.*colourFactor(ids)/(1.-z);
  case 3:
  default:
    throw Exception() << "HalfHalfOneSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

double HalfHalfOneSplitFn::invIntegOverP(const double r, const IdList & ids,
					 unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return 1. - exp(- 0.5*r/colourFactor(ids)); 
  case 1:
    return 1./(1.-exp(-0.5*r/colourFactor(ids)));
  case 2:
    return 1.-2.*colourFactor(ids)/r;
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
