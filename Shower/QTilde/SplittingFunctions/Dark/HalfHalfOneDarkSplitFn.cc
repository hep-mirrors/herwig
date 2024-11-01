// -*- C++ -*-
//
// HalfHalfOneDarkSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfHalfOneDarkSplitFn class.
//

#include "HalfHalfOneDarkSplitFn.h"
#include "Herwig/Models/HiddenValley/HiddenValleyModel.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <cassert>

using namespace Herwig;

DescribeNoPIOClass<HalfHalfOneDarkSplitFn,Herwig::Sudakov1to2FormFactor>
describeHalfHalfOneDarkSplitFn ("Herwig::HalfHalfOneDarkSplitFn","HwDarkShower.so");

void HalfHalfOneDarkSplitFn::Init() {
  // documentation
  static ClassDocumentation<HalfHalfOneDarkSplitFn> documentation
    ("The HalfHalfOneDarkSplitFn class implements the dark coloured q -> qg splitting function");
}

bool HalfHalfOneDarkSplitFn::checkColours(const IdList & ids) const {
  if(colourStructure()==TripletTripletOctet) {
    if(ids[0]!=ids[1]) return false;
    if((ids[0]->iColour()==PDT::DarkColourFundamental||
	      ids[0]->iColour()==PDT::DarkColourAntiFundamental) &&
        ids[2]->iColour()==PDT::DarkColourAdjoint) return true;
    return false;
  }
  else if(colourStructure()==OctetOctetOctet) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(ids[ix]->iColour()!=PDT::DarkColourAdjoint) return false;
    }
    return true;
  }
  else if(colourStructure()==OctetTripletTriplet) {
    if(ids[0]->iColour()!=PDT::DarkColourAdjoint) return false;
    if(ids[1]->iColour()==PDT::DarkColourFundamental&&
       ids[2]->iColour()==PDT::DarkColourAntiFundamental)
      return true;
    if(ids[1]->iColour()==PDT::DarkColourAntiFundamental&&
       ids[2]->iColour()==PDT::DarkColourFundamental)
      return true;
    return false;
  }
  else if(colourStructure()==TripletOctetTriplet) {
    if(ids[0]!=ids[2]) return false;
    if((ids[0]->iColour()==PDT::DarkColourFundamental||
	ids[0]->iColour()==PDT::DarkColourAntiFundamental) &&
       ids[1]->iColour()==PDT::DarkColourAdjoint) return true;
    return false;
  }
  else {
    assert(false);
  }
  return false;
}

double HalfHalfOneDarkSplitFn::P(const double z, const Energy2 t,
			     const IdList &ids, const bool mass, const RhoDMatrix & ) const {
  double val = (1. + sqr(z))/(1.-z);
  if(mass) {
    Energy m = ids[0]->mass();
    val -= 2.*sqr(m)/t;
  }
  return val;
}

double HalfHalfOneDarkSplitFn::overestimateP(const double z,
                                             const IdList & ) const {
  return 2./(1.-z);
}

double HalfHalfOneDarkSplitFn::ratioP(const double z, const Energy2 t,
                                      const IdList & ids, const bool mass, const RhoDMatrix & ) const {
  double val = 1. + sqr(z);
  if(mass) {
    Energy m = ids[0]->mass();
    val -= 2.*sqr(m)*(1.-z)/t;
  }
  return 0.5*val;
}

double HalfHalfOneDarkSplitFn::integOverP(const double z,
                                          const IdList & ,
                                          unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return -2.*Math::log1m(z);
  case 1:
    return  2.*log(z/(1.-z));
  case 2:
    return  2./(1.-z);
  case 3:
  default:
    throw Exception() << "HalfHalfOneDarkSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double HalfHalfOneDarkSplitFn::invIntegOverP(const double r, const IdList &,
                                             unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return 1. - exp(- 0.5*r);
  case 1:
    return 1./(1.-exp(-0.5*r));
  case 2:
    return 1.-2./r;
  case 3:
  default:
    throw Exception() << "HalfHalfOneDarkSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool HalfHalfOneDarkSplitFn::accept(const IdList &ids) const {
  // 3 particles and in and out fermion same
  if(ids.size()!=3 || ids[0]!=ids[1]) return false;
  if(ids[0]->iSpin()!=PDT::Spin1Half ||
      ids[2]->iSpin()!=PDT::Spin1) return false;
  return checkColours(ids);
}

vector<pair<int, Complex> >
HalfHalfOneDarkSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return {{ {0, 1.} }};
}

vector<pair<int, Complex> >
HalfHalfOneDarkSplitFn::generatePhiBackward(const double, const Energy2, const IdList & ,
					const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return {{ {0, 1.} }};
}

DecayMEPtr HalfHalfOneDarkSplitFn::matrixElement(const double z, const Energy2 t,
                                             const IdList & ids, const double phi,
                                             bool timeLike) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  Energy m = !timeLike ? ZERO : ids[0]->mass();
  double mt = m/sqrt(t);
  double root = sqrt(1.-(1.-z)*sqr(m)/z/t);
  double romz = sqrt(1.-z);
  double rz   = sqrt(z);
  Complex phase = exp(Complex(0.,1.)*phi);
  (*kernal)(0,0,0) = -root/romz*phase;
  (*kernal)(1,1,2) =  -conj((*kernal)(0,0,0));
  (*kernal)(0,0,2) =  root/romz*z/phase;
  (*kernal)(1,1,0) = -conj((*kernal)(0,0,2));
  (*kernal)(1,0,2) =  mt*(1.-z)/rz;
  (*kernal)(0,1,0) =  conj((*kernal)(1,0,2));
  (*kernal)(0,1,2) =  0.;
  (*kernal)(1,0,0) =  0.;
  return kernal;
}
