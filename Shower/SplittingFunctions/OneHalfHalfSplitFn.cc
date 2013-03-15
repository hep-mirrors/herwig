// -*- C++ -*-
//
// OneHalfHalfSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneHalfHalfSplitFn class.
//

#include "OneHalfHalfSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeNoPIOClass<OneHalfHalfSplitFn,Herwig::SplittingFunction>
describeOneHalfHalfSplitFn ("Herwig::OneHalfHalfSplitFn","HwShower.so");

void OneHalfHalfSplitFn::Init() {

  static ClassDocumentation<OneHalfHalfSplitFn> documentation
    ("The OneHalfHalfSplitFn class implements the splitting function for g->q qbar");

}

double OneHalfHalfSplitFn::P(const double z, const Energy2 t, 
			  const IdList &ids, const bool mass) const {
  double zz = z*(1.-z);
  double val=1.-2.*zz;
  if(mass) {
    Energy m = getParticleData(ids[1])->mass();
    val +=2.*sqr(m)/t;
  }
  return colourFactor()*val;
}

double OneHalfHalfSplitFn::overestimateP(const double, const IdList &) const {
  return colourFactor(); 
}

double OneHalfHalfSplitFn::ratioP(const double z, const Energy2 t, 
			       const IdList &ids, const bool mass) const {
  double zz = z*(1.-z);
  double val = 1.-2.*zz;
  if(mass) {
    Energy m = getParticleData(ids[1])->mass();
    val+= 2.*sqr(m)/t;
  }
  return val;
}

double OneHalfHalfSplitFn::integOverP(const double z, const IdList & ,
				   unsigned int PDFfactor) const { 
  switch(PDFfactor) {
  case 0:
    return colourFactor()*z; 
  case 1:
    return colourFactor()*log(z);
  case 2:
    return -colourFactor()*log(1.-z);
  case 3:
    return colourFactor()*log(z/(1.-z));
  default:
    throw Exception() << "OneHalfHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double OneHalfHalfSplitFn::invIntegOverP(const double r, const IdList & ,
				      unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return r/colourFactor(); 
  case 1:
    return exp(r/colourFactor());
  case 2:
    return 1.-exp(-r/colourFactor());
  case 3:
    return 1./(1.+exp(-r/colourFactor()));
  default:
    throw Exception() << "OneHalfHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool OneHalfHalfSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[1]!=-ids[2]) return false;
  tcPDPtr q=getParticleData(ids[1]);
  if(q->iSpin()!=PDT::Spin1Half) return false;
  tcPDPtr g=getParticleData(ids[0]);
  if(g->iSpin()!=PDT::Spin1) return false;
  return checkColours(ids);
}
