// -*- C++ -*-
//
// GtoGGSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoGGSplitFn class.
//

#include "GtoGGSplitFn.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

NoPIOClassDescription<GtoGGSplitFn> GtoGGSplitFn::initGtoGGSplitFn;
// Definition of the static class description member.

void GtoGGSplitFn::Init() {

  static ClassDocumentation<GtoGGSplitFn> documentation
    ("The GtoGGSplitFn class implements the g -> gg splitting function");

}

double GtoGGSplitFn::P(const double z, const Energy2, const IdList &,
		       const bool)const {
  double val = 3.*sqr(1.-z*(1.-z))/(z*(1.-z));
  // (this is historically important! the first physics - two years
  // after the birth of the project - in the Herwig++ shower! Alberto
  // & Stefan, 25/04/2002).
  return val;
}

double GtoGGSplitFn::overestimateP(const double z, const IdList &) const {
  return 3.*(1/z + 1/(1.-z)); 
}

double GtoGGSplitFn::ratioP(const double z, const Energy2, const IdList &, 
			    const bool) const {
  return sqr(1.-z*(1.-z));
}

double GtoGGSplitFn::invIntegOverP(const double r, unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 1./(1.+exp(-r/3.)); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "GtoGGSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
} 

double GtoGGSplitFn::integOverP(const double z, unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 3.*log(z/(1.-z)); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "GtoGGSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

void GtoGGSplitFn::colourConnection(tShowerParticlePtr parent,
				    tShowerParticlePtr first,
				    tShowerParticlePtr second,
				    const bool back) const {
  if(!back) {
    ColinePair cparent = ColinePair(parent->colourLine(), 
				    parent->antiColourLine());
    // ensure input consistency
    assert(cparent.first&&cparent.second);
    // Randomly decide which of the two gluon products take the
    // colour line passing for the colour of the parent gluon
    // (the other will take the one passing for the anticolour of
    //  the parent gluon).
    if(UseRandom::rndbool()) {
      ColinePtr newline=new_ptr(ColourLine());
      cparent.first->addColoured(first);
      cparent.second->addAntiColoured(second);
      newline->addColoured(second);
      newline->addAntiColoured(first);
    }
    else {
      ColinePtr newline=new_ptr(ColourLine());
      cparent.first->addColoured(second);
      cparent.second->addAntiColoured(first);
      newline->addColoured(first);
      newline->addAntiColoured(second);
    }
  }
  else {
    ColinePair cfirst = ColinePair(first->colourLine(), 
				   first->antiColourLine());
    // ensure input consistency
    assert(cfirst.first&&cfirst.second);
    // Randomly decide which of the two gluon products take the
    // colour line passing for the colour of the parent gluon
    // (the other will take the one passing for the anticolour of
    //  the parent gluon).
    if (UseRandom::rndbool()) {
      ColinePtr newline=new_ptr(ColourLine());
      cfirst.first->addColoured(parent);
      cfirst.second->addColoured(second);
      newline->addAntiColoured(second);
      newline->addAntiColoured(parent);
    }
    else {
      ColinePtr newline=new_ptr(ColourLine());
      cfirst.first->addAntiColoured(second);
      cfirst.second->addAntiColoured(parent);
      newline->addColoured(parent);
      newline->addColoured(second);
    }
  }
}

bool GtoGGSplitFn::accept(const IdList & ids)
const {
  if(ids.size()!=3) return false;
  for(unsigned int ix=0;ix<ids.size();++ix)
    {if(ids[ix]!=ParticleID::g) return false;}
  return true;
}
