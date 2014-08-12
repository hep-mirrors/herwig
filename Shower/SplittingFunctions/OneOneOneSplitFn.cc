// -*- C++ -*-
//
// OneOneOneSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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

using namespace Herwig;

DescribeNoPIOClass<OneOneOneSplitFn,Herwig::SplittingFunction>
describeOneOneOneSplitFn ("Herwig::OneOneOneSplitFn","HwShower.so");

void OneOneOneSplitFn::Init() {

  static ClassDocumentation<OneOneOneSplitFn> documentation
    ("The OneOneOneSplitFn class implements the g -> gg splitting function");

}

double OneOneOneSplitFn::P(const double z, const Energy2,
			   const IdList & ids, const bool)const {
  // (this is historically important! the first physics - two years
  // after the birth of the project - in the Herwig++ shower! Alberto
  // & Stefan, 25/04/2002).
  return colourFactor(ids)*sqr(1.-z*(1.-z))/(z*(1.-z));
}

double OneOneOneSplitFn::overestimateP(const double z,
				       const IdList & ids) const {
  return colourFactor(ids)*(1/z + 1/(1.-z)); 
}


double OneOneOneSplitFn::ratioP(const double z, const Energy2,
				const IdList & , const bool) const {
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
    tcPDPtr part = getParticleData(ids[ix]);
    if(part->iSpin()!=PDT::Spin1) return false;
  }
  return checkColours(ids);
}

vector<pair<int, Complex> > 
OneOneOneSplitFn::generatePhi(ShowerParticle & ,ShoKinPtr ,
				const double z, const Energy2, const IdList &,
				const RhoDMatrix & rho) { 
  assert(rho.iSpin()==PDT::Spin1);
  double modRho = abs(rho(0,2));
  double max = 2.*z*modRho*(1.-z)+sqr(+1.-(1.-z)*z)/(z*(1.-z));
  vector<pair<int, Complex> > output;
  output.push_back(make_pair( 0,(rho(0,0)+rho(2,2))*sqr(1.-(1.-z)*z)/(z*(1.-z))/max));
  output.push_back(make_pair(-2,-rho(0,2)*z*(1.-z)));
  output.push_back(make_pair( 2,-rho(2,0)*z*(1.-z)));
  return output;
}

DecayMatrixElement OneOneOneSplitFn::matrixElement(ShowerParticle &,ShoKinPtr,
						     const double z, const Energy2, 
						     const IdList &, const double phi) {
  // calculate the kernal
  DecayMatrixElement kernal(PDT::Spin1,PDT::Spin1,PDT::Spin1);
  double root = sqrt(z*(1.-z));
  double omz = 1.-z;
  double fact = 1./z/(1.-z)-1.;
  for(int lamA=-1;lamA<2;lamA+=2) {
    Complex factA = exp(Complex(0.,1.)*double(lamA)*phi);
    for(int lamB=-1;lamB<2;lamB+=2) {
      Complex factB = exp(-Complex(0.,1.)*double(lamB)*phi);
      for(int lamC=-1;lamC<2;lamC+=2) {
	Complex factC = exp(-Complex(0.,1.)*double(lamC)*phi);
	kernal(1+lamA,1+lamB,1+lamC) = -0.5*root*factA*factB*factC*
	  (lamA*lamB*lamC*fact+lamA+lamB/z+lamC/omz);
      }
    }
  }
  return kernal;
}
