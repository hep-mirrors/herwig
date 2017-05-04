// -*- C++ -*-
//
// PhitoPhiGSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;

DescribeNoPIOClass<ZeroZeroOneSplitFn,Herwig::SplittingFunction>
describeZeroZeroOneSplitFn ("Herwig::ZeroZeroOneSplitFn","HwShower.so");

void ZeroZeroOneSplitFn::Init() {

  static ClassDocumentation<ZeroZeroOneSplitFn> documentation
    ("The ZeroZeroOneSplitFn class implements the splitting function for the "
     "radiation of a gluon by a scalar coloured particle");

}

double ZeroZeroOneSplitFn::P(const double z, const Energy2 t,
			     const IdList &ids, const bool mass, const RhoDMatrix &) const {
  double val = z/(1.-z);
  if(mass) {
    Energy m = ids[0]->mass();
    val-=  sqr(m)/t;
  }
  return 2.*colourFactor(ids)*val;
}

double ZeroZeroOneSplitFn::overestimateP(const double z,
					 const IdList &ids) const { 
  return 2.*colourFactor(ids)/(1.-z); 
}

double ZeroZeroOneSplitFn::ratioP(const double z, const Energy2 t,
				  const IdList &ids,const bool mass, const RhoDMatrix &) const { 
  double val = z;
  if(mass) {
    Energy m = ids[0]->mass();
    val-=sqr(m)*(1.-z)/t;
  }
  return val;
} 

double ZeroZeroOneSplitFn::integOverP(const double z, const IdList & ids,
				    unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return -2.*colourFactor(ids)*log(1.-z); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "ZeroZeroOneSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double ZeroZeroOneSplitFn::invIntegOverP(const double r, const IdList & ids,
				       unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 1. - exp(- 0.5*r/colourFactor(ids));
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
  if(ids[0]->iSpin()!=PDT::Spin0 ||
     ids[2]->iSpin()!=PDT::Spin1) return false;
  return checkColours(ids);
}

vector<pair<int, Complex> > 
ZeroZeroOneSplitFn::generatePhiForward(const double, const Energy2, const IdList &,
				const RhoDMatrix &) {
  // scalar so no dependence
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}

vector<pair<int, Complex> > 
ZeroZeroOneSplitFn::generatePhiBackward(const double, const Energy2, const IdList &,
					const RhoDMatrix &) {
  // scalar so no dependence
  assert(false);
  return vector<pair<int, Complex> >(1,make_pair(0,1.));
}

DecayMEPtr ZeroZeroOneSplitFn::matrixElement(const double z, const Energy2 t, 
					     const IdList & ids, const double phi,
                                             bool timeLike) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  Energy m = timeLike ? ids[0]->mass() : ZERO;
  (*kernal)(0,0,0) = -exp(Complex(0.,1.)*phi)*sqrt(1.-(1.-z)*sqr(m)/z/t)*sqrt(z/(1.-z));
  (*kernal)(0,0,2) = -conj((*kernal)(0,0,0));
  return kernal;
}
