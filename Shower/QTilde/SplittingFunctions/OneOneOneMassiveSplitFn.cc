// -*- C++ -*-
//
// OneOneOneMassiveSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneOneOneMassiveSplitFn class.
//

#include "OneOneOneMassiveSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

DescribeNoPIOClass<OneOneOneMassiveSplitFn,Herwig::SplittingFunction>
describeOneOneOneMassiveSplitFn ("Herwig::OneOneOneMassiveSplitFn","HwShower.so");

void OneOneOneMassiveSplitFn::Init() {

  static ClassDocumentation<OneOneOneMassiveSplitFn> documentation
    ("The OneOneOneMassiveSplitFn class implements the g -> gg splitting function");

}

double OneOneOneMassiveSplitFn::P(const double z, const Energy2 t,
				  const IdList & ids, const bool, const RhoDMatrix & rho )const {
  Energy2 m2 = sqr(ids[0]->mass());
  double rho00 = rho(1,1).real();
  return 2.*colourFactor(ids)*(z/(1.-z)-m2/t
			       + (1.-rho00)*( (1.-z)/z + z*(1.-z) -sqr(1.-z)*m2/t )
			       + 2.*rho00*sqr(1.-z)*m2/t);
}

double OneOneOneMassiveSplitFn::overestimateP(const double z,
					      const IdList & ids) const {
  return 2.*colourFactor(ids)*(1/z + 1/(1.-z)); 
}


double OneOneOneMassiveSplitFn::ratioP(const double z, const Energy2 t,
				       const IdList & ids , const bool, const RhoDMatrix & rho) const {
  Energy2 m2 = sqr(ids[0]->mass());
  double rho00 = rho(1,1).real();
  return (sqr(z)-z*(1.-z)*m2/t
	  + (1.-rho00)*sqr(1.-z)*( 1. + sqr(z) -z*(1.-z)*m2/t )
	  + 2.*rho00*z*(1.-z)*sqr(1.-z)*m2/t);
}

double OneOneOneMassiveSplitFn::invIntegOverP(const double r,
				       const IdList & ids,
				       unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 1./(1.+exp(-0.5*r/colourFactor(ids))); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneOneOneMassiveSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
} 

double OneOneOneMassiveSplitFn::integOverP(const double z, const IdList & ids,
				    unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    assert(z>0.&&z<1.);
    return 2.*colourFactor(ids)*log(z/(1.-z)); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "OneOneOneMassiveSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool OneOneOneMassiveSplitFn::accept(const IdList & ids) const {
  if(ids.size()!=3) return false;
  for(unsigned int ix=0;ix<ids.size();++ix) {
    if(ids[0]->iSpin()!=PDT::Spin1) return false;
  }
  return checkColours(ids);
}

vector<pair<int, Complex> > 
OneOneOneMassiveSplitFn::generatePhiForward(const double z, const Energy2, const IdList &,
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
OneOneOneMassiveSplitFn::generatePhiBackward(const double z, const Energy2, const IdList &,
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

DecayMEPtr OneOneOneMassiveSplitFn::matrixElement(const double z, const Energy2 t, 
						  const IdList & ids, const double phi,
						  bool) {
  // calculate the kernal
  DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  Energy2 m2 = sqr(ids[0]->mass());
  double omz = 1.-z;
  double root = sqrt(z*omz);
  Complex phase = exp(Complex(0.,1.)*phi);
  double mroot = sqrt(1.-m2*(1.-z)/z/t);
  (*kernal)(0,0,0) =  phase/root*mroot;
  (*kernal)(2,2,2) = -conj((*kernal)(0,0,0));
  (*kernal)(0,0,2) = -sqr(z)/root/phase*mroot;
  (*kernal)(2,2,0) = -conj((*kernal)(0,0,2));
  (*kernal)(0,2,0) = -sqr(omz)/root/phase*mroot;
  (*kernal)(2,0,2) = -conj((*kernal)(0,2,0));
  (*kernal)(0,2,2) = 0.;
  (*kernal)(2,0,0) = 0.;
  (*kernal)(1,0,0) = 0.;
  (*kernal)(1,2,2) = 0.;
  (*kernal)(1,0,2) = sqrt(2.*m2/t)*(1.-z);
  (*kernal)(1,2,0) = (*kernal)(1,0,2);  
  return kernal;
}
