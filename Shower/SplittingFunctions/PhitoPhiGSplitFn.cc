// -*- C++ -*-
//
// PhitoPhiGSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhitoPhiGSplitFn class.
//

#include "PhitoPhiGSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

NoPIOClassDescription<PhitoPhiGSplitFn> PhitoPhiGSplitFn::initPhitoPhiGSplitFn;
// Definition of the static class description member.

void PhitoPhiGSplitFn::Init() {

  static ClassDocumentation<PhitoPhiGSplitFn> documentation
    ("The PhitoPhiGSplitFn class implements the splitting function for the "
     "radiation of a gluon by a scalar coloured particle");

}

double PhitoPhiGSplitFn::P(const double z, const Energy2 t,
			   const IdList &ids, const bool mass) const {
  double val = z/(1.-z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=  sqr(m)/t;
  }
  return 8./3.*val;
}

double PhitoPhiGSplitFn::overestimateP(const double z, const IdList &) const { 
  return 8./3./(1.-z); 
}

double PhitoPhiGSplitFn::ratioP(const double z, const Energy2 t,
				const IdList &ids,const bool mass) const { 
  double val = z;
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=sqr(m)*(1.-z)/t;
  }
  return val;
} 

double PhitoPhiGSplitFn::integOverP(const double z, unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return -8./3.*log(1.-z); 
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "PhitoPhiGSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double PhitoPhiGSplitFn::invIntegOverP(const double r, unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 1. - exp(- 3.*r/8.);
  case 1:
  case 2:
  case 3:
  default:
    throw Exception() << "PhitoPhiGSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

void PhitoPhiGSplitFn::colourConnection(tShowerParticlePtr parent,
					tShowerParticlePtr first,
					tShowerParticlePtr second,
					const bool back) const {
  if(!back) {
    ColinePair cparent = ColinePair(parent->colourLine(), 
				    parent->antiColourLine());
    // ensure input consistency
    assert(( cparent.first && !cparent.second) || 
	   (!cparent.first &&  cparent.second));
    // q~ -> q~ g
    if(cparent.first) {
      ColinePtr newline=new_ptr(ColourLine());
      cparent.first->addColoured(second);
      newline->addColoured     ( first);
      newline->addAntiColoured (second);
    }
    // q~bar -> q~bar g
    else {
      ColinePtr newline=new_ptr(ColourLine());
      cparent.second->addAntiColoured(second);
      newline->addColoured(second);
      newline->addAntiColoured(first);
    }
  }
  else {
    ColinePair cfirst = ColinePair(first->colourLine(), 
				   first->antiColourLine());
    // ensure input consistency
    assert(( cfirst.first && !cfirst.second) ||
	   (!cfirst.first &&  cfirst.second)); 
    // q~ -> q~ g
    if(cfirst.first) {
      ColinePtr newline=new_ptr(ColourLine());
      cfirst.first->addAntiColoured(second);
      newline->addColoured(second);
      newline->addColoured(parent);
    }
    // q~bar -> q~bar g
    else {
      ColinePtr newline=new_ptr(ColourLine());
      cfirst.second->addColoured(second);
      newline->addAntiColoured(second);
      newline->addAntiColoured(parent);
    }
  }
}

bool PhitoPhiGSplitFn::accept(const IdList &ids) const
{
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[1]||ids[2]!=ParticleID::g) return false;
  tcPDPtr q=getParticleData(ids[0]);
  return q->iSpin()==PDT::Spin0&&(q->iColour()==PDT::Colour3||
				  q->iColour()==PDT::Colour3bar);
}
