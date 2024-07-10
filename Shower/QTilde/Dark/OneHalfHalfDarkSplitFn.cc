// -*- C++ -*-
//
// OneHalfHalfDarkSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneHalfHalfDarkSplitFn class.
//

#include "OneHalfHalfDarkSplitFn.h"
#include "HiddenValleyModel.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeNoPIOClass<OneHalfHalfDarkSplitFn,Herwig::Sudakov1to2FormFactor>
describeOneHalfHalfDarkSplitFn ("Herwig::OneHalfHalfDarkSplitFn","HwShower.so");

void OneHalfHalfDarkSplitFn::Init() {

  static ClassDocumentation<OneHalfHalfDarkSplitFn> documentation
    ("The OneHalfHalfDarkSplitFn class implements the splitting function for g->q qbar");

}

bool OneHalfHalfDarkSplitFn::checkColours(const IdList & ids) const {
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

double OneHalfHalfDarkSplitFn::P(const double z, const Energy2 t,
			     const IdList &ids, const bool mass, const RhoDMatrix &) const {
  double zz = z*(1.-z);
  double val=1.-2.*zz;
  if(mass) {
    Energy m = ids[1]->mass();
    val +=2.*sqr(m)/t;
  }
  return val;
}

double OneHalfHalfDarkSplitFn::overestimateP(const double,
                                             const IdList &) const {
  return 1.;
}

double OneHalfHalfDarkSplitFn::ratioP(const double z, const Energy2 t,
                                      const IdList &ids, const bool mass, const RhoDMatrix &) const {
  double zz = z*(1.-z);
  double val = 1.-2.*zz;
  if(mass) {
    Energy m = ids[1]->mass();
    val+= 2.*sqr(m)/t;
  }
  return val;
}

double OneHalfHalfDarkSplitFn::integOverP(const double z, const IdList & ,
                                          unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return z;
  case 1:
    return log(z);
  case 2:
    return -log(1.-z);
  case 3:
    return log(z/(1.-z));
  case 4:
    return 2.*sqrt(z);
  case 5:
    return (2./3.)*z*sqrt(z);
  default:
    throw Exception() << "OneHalfHalfDarkSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double OneHalfHalfDarkSplitFn::invIntegOverP(const double r,
                                             const IdList & ,
                                             unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return r;
  case 1:
    return exp(r);
  case 2:
    return 1.-exp(-r);
  case 3:
    return 1./(1.+exp(-r));
  case 4:
    return 0.25*sqr(r);
  case 5:
    return pow(1.5*r,2./3.);
  default:
    throw Exception() << "OneHalfHalfDarkSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool OneHalfHalfDarkSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[1]!=ids[2]->CC()) return false;
  if(ids[1]->iSpin()!=PDT::Spin1Half) return false;
  if(ids[0]->iSpin()!=PDT::Spin1) return false;
  return OneHalfHalfDarkSplitFn::checkColours(ids);
}

vector<pair<int, Complex> >
OneHalfHalfDarkSplitFn::generatePhiForward(const double z, const Energy2 t, const IdList & ids,
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
OneHalfHalfDarkSplitFn::generatePhiBackward(const double, const Energy2, const IdList &,
					const RhoDMatrix & ) {
  // no dependance
  return {{ {0, 1.} }};
}

DecayMEPtr OneHalfHalfDarkSplitFn::matrixElement(const double z, const Energy2 t,
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
