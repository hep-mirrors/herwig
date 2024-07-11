// -*- C++ -*-
//
// HalfOneHalfSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfOneHalfSplitFn class.
//

#include "HalfOneHalfSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;

DescribeNoPIOClass<HalfOneHalfSplitFn,Herwig::Sudakov1to2FormFactor>
describeHalfOneHalfSplitFn ("Herwig::HalfOneHalfSplitFn","HwShower.so");

void HalfOneHalfSplitFn::Init() {

  static ClassDocumentation<HalfOneHalfSplitFn> documentation
    ("The HalfOneHalfSplitFn class implements the splitting "
     "function for q -> g q");

}
double HalfOneHalfSplitFn::integOverP(const double z, const IdList &,
				      unsigned int PDFfactor) const { 
  switch(PDFfactor) {
  case 0:
    return 2.*log(z); 
  case 1:
    return -2./z;
  case 2:
    return 2.*log(z/(1.-z));
  case 3:
  default:
    throw Exception() << "HalfOneHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double HalfOneHalfSplitFn::invIntegOverP(const double r, 
					 const IdList & ,
					 unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return exp(0.5*r); 
  case 1:
    return -2./r;
  case 2:
    return 1./(1.+exp(-0.5*r));
  case 3:
  default:
    throw Exception() << "HalfOneHalfSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool HalfOneHalfSplitFn::accept(const IdList &ids) const {
  // 3 particles and in and out fermion same
  if(ids.size()!=3 || ids[0]!=ids[2]) return false;
  if(ids[0]->iSpin()!=PDT::Spin1Half ||
     ids[1]->iSpin()!=PDT::Spin1) return false;
  return checkColours(ids);
}


vector<pair<int, Complex> > 
HalfOneHalfSplitFn::generatePhiForward(const double, const Energy2, const IdList & ,
				       const RhoDMatrix &) {
  // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
  // and rest = splitting function for Tr(rho)=1 as required by defn
  return {{ {0, 1.} }};
}

vector<pair<int, Complex> > 
HalfOneHalfSplitFn::generatePhiBackward(const double z, const Energy2 t, const IdList & ids,
					const RhoDMatrix & rho) {
  assert(rho.iSpin()==PDT::Spin1);
  double mt = sqr(ids[0]->mass())/t;
  double diag = (1.+sqr(1.-z))/z - 2.*mt;
  double off  = 2.*(1.-z)/z*(1.-mt*z/(1.-z));
  double max = diag+2.*abs(rho(0,2))*off;
  vector<pair<int, Complex> > output;
  output.push_back(make_pair( 0, (rho(0,0)+rho(2,2))*diag/max));
  output.push_back(make_pair( 2, -rho(0,2)          * off/max)); 
  output.push_back(make_pair(-2, -rho(2,0)          * off/max));
  return output;
}

DecayMEPtr HalfOneHalfSplitFn::matrixElement(const double z, const Energy2 t, 
					     const IdList & ids, const double phi,
                                             bool timeLike) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1,PDT::Spin1Half)));
  Energy m = !timeLike ? ZERO : ids[0]->mass();
  double mt = m/sqrt(t);
  double root = sqrt(1.-z*sqr(m)/(1.-z)/t);
  double romz = sqrt(1.-z); 
  double rz   = sqrt(z);
  Complex phase = exp(-Complex(0.,1.)*phi);
  (*kernal)(0,0,0) = -root/rz/phase;
  (*kernal)(1,2,1) = -conj((*kernal)(0,0,0));
  (*kernal)(0,2,0) =  root/rz*(1.-z)*phase;
  (*kernal)(1,0,1) = -conj((*kernal)(0,2,0));
  (*kernal)(1,2,0) =  mt*z/romz;
  (*kernal)(0,0,1) =  conj((*kernal)(1,2,0));
  (*kernal)(0,2,1) =  0.;
  (*kernal)(1,0,0) =  0.;
  return kernal;
}
