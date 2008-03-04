// -*- C++ -*-
//
// GluinotoGluinoGSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GluinotoGluinoGSplitFn class.
//

#include "GluinotoGluinoGSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Repository/UseRandom.h"
#include <cassert>

using namespace Herwig;

NoPIOClassDescription<GluinotoGluinoGSplitFn> 
GluinotoGluinoGSplitFn::initGluinotoGluinoGSplitFn;
// Definition of the static class description member.

void GluinotoGluinoGSplitFn::Init() {

  static ClassDocumentation<GluinotoGluinoGSplitFn> documentation
    ("The GluinotoGluinoGSplitFn class implements the splitting function for "
     "gluino -> gluino g");

}


double GluinotoGluinoGSplitFn::P(const double z, const Energy2 t,
		       const IdList &ids, const bool mass) const {
  double val = (1. + sqr(z))/(1.-z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();  
    val -= 2.*sqr(m)/t;
  }
  return 3.*val;

}

double GluinotoGluinoGSplitFn::overestimateP(const double z, const IdList &) const {
  return 6./(1.-z); 
}

double GluinotoGluinoGSplitFn::ratioP(const double z, const Energy2 t,
			    const IdList &ids, const bool mass) const {
  double val = 1. + sqr(z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val -= 2.*sqr(m)*(1.-z)/t;
  } 
  return 0.5*val;
} 

double GluinotoGluinoGSplitFn::integOverP(const double z, const IdList & ,
					  unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return -6.*log(1.-z); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "GluinotoGluinoGSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double GluinotoGluinoGSplitFn::invIntegOverP(const double r, const IdList & ,
					     unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 1. - exp(-r/6.); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "GluinotoGluinoGSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

void GluinotoGluinoGSplitFn::colourConnection(tShowerParticlePtr parent,
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

bool GluinotoGluinoGSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[1]||ids[2]!=ParticleID::g) return false;
  tcPDPtr q=getParticleData(ids[0]);
  return q->iSpin()==PDT::Spin1Half&&q->iColour()==PDT::Colour8;
}
