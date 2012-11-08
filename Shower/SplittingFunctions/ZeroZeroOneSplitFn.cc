// -*- C++ -*-
//
// PhitoPhiGSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZeroZeroOneSplitFn class.
//

#include "ZeroZeroOneSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeNoPIOClass<ZeroZeroOneSplitFn,Herwig::SplittingFunction>
describeZeroZeroOneSplitFn ("Herwig::ZeroZeroOneSplitFn","HwShower.so");

void ZeroZeroOneSplitFn::Init() {

  static ClassDocumentation<ZeroZeroOneSplitFn> documentation
    ("The ZeroZeroOneSplitFn class implements the splitting function for the "
     "radiation of a gluon by a scalar coloured particle");

}

double ZeroZeroOneSplitFn::P(const double z, const Energy2 t,
			   const IdList &ids, const bool mass) const {
  double val = z/(1.-z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=  sqr(m)/t;
  }
  return 2.*colourFactor()*val;
}

double ZeroZeroOneSplitFn::overestimateP(const double z, const IdList &) const { 
  return 2.*colourFactor()/(1.-z); 
}

double ZeroZeroOneSplitFn::ratioP(const double z, const Energy2 t,
				const IdList &ids,const bool mass) const { 
  double val = z;
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=sqr(m)*(1.-z)/t;
  }
  return val;
} 

double ZeroZeroOneSplitFn::integOverP(const double z, const IdList & ,
				    unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return -2.*colourFactor()*log(1.-z); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "ZeroZeroOneSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double ZeroZeroOneSplitFn::invIntegOverP(const double r, const IdList & ,
				       unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 1. - exp(- 0.5*r/colourFactor());
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "ZeroZeroOneSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

bool ZeroZeroOneSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[1]) return false;
  tcPDPtr q=getParticleData(ids[0]);
  tcPDPtr g=getParticleData(ids[2]);
  if(q->iSpin()!=PDT::Spin0 ||
     g->iSpin()!=PDT::Spin1) return false;
  return checkColours(ids);
}
