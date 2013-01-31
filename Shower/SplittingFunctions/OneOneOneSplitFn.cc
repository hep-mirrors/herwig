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
