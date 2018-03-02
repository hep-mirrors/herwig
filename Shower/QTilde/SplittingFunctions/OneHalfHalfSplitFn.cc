// -*- C++ -*-
//
// OneHalfHalfSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;

DescribeNoPIOClass<OneHalfHalfSplitFn,Herwig::SplittingFunction>
describeOneHalfHalfSplitFn ("Herwig::OneHalfHalfSplitFn","HwShower.so");

void OneHalfHalfSplitFn::Init() {

  static ClassDocumentation<OneHalfHalfSplitFn> documentation
    ("The OneHalfHalfSplitFn class implements the splitting function for g->q qbar");

}

double OneHalfHalfSplitFn::P(const double z, const Energy2 t, 
			     const IdList &ids, const bool mass, const RhoDMatrix &) const {
  double zz = z*(1.-z);
  double val=1.-2.*zz;
  if(mass) {
    Energy m = ids[1]->mass();
    val +=2.*sqr(m)/t;
  }
  return colourFactor(ids)*val;
}

double OneHalfHalfSplitFn::overestimateP(const double,
					 const IdList &ids) const {
  return colourFactor(ids); 
}

double OneHalfHalfSplitFn::ratioP(const double z, const Energy2 t, 
			       const IdList &ids, const bool mass, const RhoDMatrix &) const {
  double zz = z*(1.-z);
  double val = 1.-2.*zz;
  if(mass) {
    Energy m = ids[1]->mass();
    val+= 2.*sqr(m)/t;
  }
  return val;
}

double OneHalfHalfSplitFn::integOverP(const double z, const IdList & ids,
				      unsigned int PDFfactor) const { 
  switch(PDFfactor) {
  case 0:
    return colourFactor(ids)*z; 
  case 1:
    return colourFactor(ids)*log(z);
  case 2:
    return -colourFactor(ids)*log(1.-z);
  case 3:
    return colourFactor(ids)*log(z/(1.-z));
  default:
    throw Exception() << "OneHalfHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double OneHalfHalfSplitFn::invIntegOverP(const double r,
					 const IdList & ids,
					 unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return r/colourFactor(ids); 
  case 1:
    return exp(r/colourFactor(ids));
  case 2:
    return 1.-exp(-r/colourFactor(ids));
  case 3:
    return 1./(1.+exp(-r/colourFactor(ids)));
  default:
    throw Exception() << "OneHalfHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool OneHalfHalfSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[1]!=ids[2]->CC()) return false;
  if(ids[1]->iSpin()!=PDT::Spin1Half) return false;
  if(ids[0]->iSpin()!=PDT::Spin1) return false;
  return checkColours(ids);
}

vector<pair<int, Complex> > 
OneHalfHalfSplitFn::generatePhiForward(const double z, const Energy2 t, const IdList & ids,
				       const RhoDMatrix & rho) { 
  assert(rho.iSpin()==PDT::Spin1);
  double modRho = abs(rho(0,2));
  Energy mq = ids[1]->mass();
  Energy2 mq2 = sqr(mq);
  double fact = z*(1.-z)-mq2/t;
  double max = 1.+2.*fact*(-1.+2.*modRho);
  vector<pair<int, Complex> > output;
  output.push_back(make_pair( 0,(rho(0,0)+rho(2,2))*(1.-2.*fact)/max));
  output.push_back(make_pair(-2,2.*fact*rho(0,2)/max));
  output.push_back(make_pair( 2,2.*fact*rho(2,0)/max));
  return output;
}

vector<pair<int, Complex> > 
OneHalfHalfSplitFn::generatePhiBackward(const double, const Energy2, const IdList &,
					const RhoDMatrix & ) { 
  // no dependance
  return {{ {0, 1.} }};
}

DecayMEPtr OneHalfHalfSplitFn::matrixElement(const double z, const Energy2 t, 
					     const IdList & ids, const double phi,
                                             bool timeLike) {
  static const Complex ii(0.,1.);
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  double mt = !timeLike ? ZERO : ids[1]->mass()/sqrt(t);
  double root =sqrt(1.-sqr(mt)/z/(1.-z));
  (*kernal)(0,0,0) = mt/sqrt(z*(1.-z));
  (*kernal)(2,1,1) = (*kernal)(0,0,0);
  (*kernal)(0,0,1) = -z*root*exp(-ii*phi);
  (*kernal)(2,1,0) = -conj((*kernal)(0,0,1));
  (*kernal)(0,1,0) = (1.-z)*exp(-ii*phi)*root;
  (*kernal)(2,0,1) = -conj((*kernal)(0,1,0));
  (*kernal)(0,1,1) = 0.;
  (*kernal)(2,0,0) = 0.;
  return kernal;
}
