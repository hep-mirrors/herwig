// -*- C++ -*-
//
// OneOneOneSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
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

DescribeNoPIOClass<OneOneOneSplitFn,Herwig::SudakovFormFactor>
describeOneOneOneSplitFn ("Herwig::OneOneOneSplitFn","HwShower.so");

void OneOneOneSplitFn::Init() {

  static ClassDocumentation<OneOneOneSplitFn> documentation
    ("The OneOneOneSplitFn class implements the g -> gg splitting function");

}

vector<pair<int, Complex> > 
OneOneOneSplitFn::generatePhiForward(const double z, const Energy2, const IdList &,
				     const RhoDMatrix & rho) {
  assert(rho.iSpin()==PDT::Spin1);
  double modRho = abs(rho(0,2));
  double max = 2.*z*modRho*(1.-z)+sqr(1.-(1.-z)*z)/(z*(1.-z));
  vector<pair<int, Complex> > output;
  output.reserve(3);
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
  output.reserve(3);
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
