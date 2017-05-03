// -*- C++ -*-
//
// OneOneOneSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneOneOneSplitFn class.
//

#include "OneOneOneSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;

DescribeNoPIOClass<OneOneOneSplitFn,Herwig::SplittingFunction>
describeOneOneOneSplitFn ("Herwig::OneOneOneSplitFn","HwShower.so");

void OneOneOneSplitFn::Init() {

  static ClassDocumentation<OneOneOneSplitFn> documentation
    ("The OneOneOneSplitFn class implements the g -> gg splitting function");

}

double OneOneOneSplitFn::P(const double z, const Energy2,
			   const IdList & ids, const bool, const RhoDMatrix &)const {
  // (this is historically important! the first physics - two years
  // after the birth of the project - in the Herwig shower! Alberto
  // & Stefan, 25/04/2002).
  return colourFactor(ids)*sqr(1.-z*(1.-z))/(z*(1.-z));
}

double OneOneOneSplitFn::overestimateP(const double z,
				       const IdList & ids) const {
  return colourFactor(ids)*(1/z + 1/(1.-z)); 
}


double OneOneOneSplitFn::ratioP(const double z, const Energy2,
				const IdList & , const bool, const RhoDMatrix &) const {
  return sqr(1.-z*(1.-z));
}

double OneOneOneSplitFn::invIntegOverP(const double r,
				       const IdList & ids,
				       unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 1./(1.+exp(-r/colourFactor(ids))); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneOneOneSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
} 

double OneOneOneSplitFn::integOverP(const double z, const IdList & ids,
				    unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    assert(z>0.&&z<1.);
    return colourFactor(ids)*log(z/(1.-z)); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneOneOneSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool OneOneOneSplitFn::accept(const IdList & ids) const {
  if(ids.size()!=3) return false;
  for(unsigned int ix=0;ix<ids.size();++ix) {
    if(ids[0]->iSpin()!=PDT::Spin1) return false;
  }
  return checkColours(ids);
}

vector<pair<int, Complex> > 
OneOneOneSplitFn::generatePhiForward(const double z, const Energy2, const IdList &,
				     const RhoDMatrix & rho) {
  assert(rho.iSpin()==PDT::Spin1);
  double modRho = abs(rho(0,2));
  double max = 2.*z*modRho*(1.-z)+sqr(1.-(1.-z)*z)/(z*(1.-z));
  vector<pair<int, Complex> > output;
  output.push_back(make_pair( 0,(rho(0,0)+rho(2,2))*sqr(1.-(1.-z)*z)/(z*(1.-z))/max));
  output.push_back(make_pair(-2,-rho(0,2)*z*(1.-z)/max));
  output.push_back(make_pair( 2,-rho(2,0)*z*(1.-z)/max));
  return output;
}

vector<pair<int, Complex> > 
OneOneOneSplitFn::generatePhiBackward(const double z, const Energy2, const IdList &,
			      const RhoDMatrix & rho) {
  assert(rho.iSpin()==PDT::Spin1);
  double diag = sqr(1 - (1 - z)*z)/(1 - z)/z;
  double off  = (1.-z)/z;
  double max  = 2.*abs(rho(0,2))*off+diag;
  vector<pair<int, Complex> > output;
  output.push_back(make_pair( 0, (rho(0,0)+rho(2,2))*diag/max));
  output.push_back(make_pair( 2,-rho(0,2)           * off/max));
  output.push_back(make_pair(-2,-rho(2,0)           * off/max));
  return output;
}

DecayMEPtr OneOneOneSplitFn::matrixElement(const double z, const Energy2, 
					   const IdList &, const double phi,
                                           bool) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  double omz = 1.-z;
  double root = sqrt(z*omz);
  Complex phase = exp(Complex(0.,1.)*phi);
  (*kernal)(0,0,0) =  phase/root;
  (*kernal)(2,2,2) = -conj((*kernal)(0,0,0));
  (*kernal)(0,0,2) = -sqr(z)/root/phase;
  (*kernal)(2,2,0) = -conj((*kernal)(0,0,2));
  (*kernal)(0,2,0) = -sqr(omz)/root/phase;
  (*kernal)(2,0,2) = -conj((*kernal)(0,2,0));
  (*kernal)(0,2,2) = 0.;
  (*kernal)(2,0,0) = 0.;
  return kernal;
}
