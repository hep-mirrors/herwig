// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoQQbarSplitFn class.
//

#include "GtoQQbarSplitFn.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

NoPIOClassDescription<GtoQQbarSplitFn> GtoQQbarSplitFn::initGtoQQbarSplitFn;
// Definition of the static class description member.

void GtoQQbarSplitFn::Init() {

  static ClassDocumentation<GtoQQbarSplitFn> documentation
    ("The GtoQQbarSplitFn class implements the splitting function for g->q qbar");

}

double GtoQQbarSplitFn::P(const double z, const Energy2 t, 
			  const IdList &ids, const bool mass) const {
  double zz = z*(1.-z);
  double val=1.-2.*zz;
  if(mass) {
    Energy m = getParticleData(ids[1])->mass();
    val +=2.*sqr(m)/t;
  }
  return 0.5*val;
}

double GtoQQbarSplitFn::overestimateP(const double, const IdList &) const {
  return 0.5; 
}

double GtoQQbarSplitFn::ratioP(const double z, const Energy2 t, 
			       const IdList &ids, const bool mass) const {
  double zz = z*(1.-z);
  double val = 1.-2.*zz;
  if(mass) {
    Energy m = getParticleData(ids[1])->mass();
    val+= 2.*sqr(m)/t;
  }
  return val;
}

double GtoQQbarSplitFn::integOverP(const double z, unsigned int PDFfactor) const { 
  switch(PDFfactor) {
  case 0:
    return z/2.; 
  case 1:
    return 0.5*log(z);
  case 2:
    return -0.5*log(1.-z);
  case 3:
    return 0.5*log(z/(1.-z));
  default:
    throw Exception() << "GtoQQbarSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double GtoQQbarSplitFn::invIntegOverP(const double r, unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return 2.*r; 
  case 1:
    return exp(2.*r);
  case 2:
    return 1.-exp(-2.*r);
  case 3:
    return 1./(1.+exp(-2.*r));
  default:
    throw Exception() << "GtoQQbarSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

void GtoQQbarSplitFn::colourConnection(tShowerParticlePtr parent,
				       tShowerParticlePtr first,
				       tShowerParticlePtr second,
				       const bool back) const {
  if(!back) {
    ColinePair cparent = ColinePair(parent->colourLine(), 
				    parent->antiColourLine());
    // ensure input consistency
    assert(cparent.first&&cparent.second);
    cparent.first ->addColoured    ( first);
    cparent.second->addAntiColoured(second);
  }
  else {
    ColinePair cfirst = ColinePair(first->colourLine(), 
				   first->antiColourLine());
    // ensure input consistency
    assert(( cfirst.first && !cfirst.second) ||
	   (!cfirst.first &&  cfirst.second));
    // g -> q qbar
    if(cfirst.first) {
      ColinePtr newline=new_ptr(ColourLine());
      cfirst.first->addColoured(parent);
      newline->addAntiColoured(second);
      newline->addAntiColoured(parent);
    }
    // g -> qbar q
    else {
      ColinePtr newline=new_ptr(ColourLine());
      cfirst.second->addAntiColoured(parent);
      newline->addColoured(second);
      newline->addColoured(parent);
    }
  }
}

bool GtoQQbarSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[1]!=-ids[2]||ids[0]!=ParticleID::g) return false;
  tcPDPtr q=getParticleData(ids[1]);
  return q->iSpin()==PDT::Spin1Half&&(q->iColour()==PDT::Colour3||
				      q->iColour()==PDT::Colour3bar);
}
