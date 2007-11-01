// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoGQSplitFn class.
//

#include "QtoGQSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

NoPIOClassDescription<QtoGQSplitFn> QtoGQSplitFn::initQtoGQSplitFn;
// Definition of the static class description member.

void QtoGQSplitFn::Init() {

  static ClassDocumentation<QtoGQSplitFn> documentation
    ("The QtoGQSplitFn class implements the splitting function for q -> g q");

}

double QtoGQSplitFn::P(const double z, const Energy2 t,
		       const IdList &ids, const bool mass) const {
  double val=(2.*(1.-z)+sqr(z))/z;
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=2.*sqr(m)/t;
  }
  return 4./3.*val;
}

double QtoGQSplitFn::overestimateP(const double z, const IdList &) const { 
  return 8./3./z; 
}

double QtoGQSplitFn::ratioP(const double z, const Energy2 t,
			    const IdList &ids,const bool mass) const {
  double val=2.*(1.-z)+sqr(z);
  if(mass) {
    Energy m=getParticleData(ids[0])->mass();
    val -=2.*sqr(m)*z/t;
  }
  return 0.5*val;
}

double QtoGQSplitFn::integOverP(const double z, unsigned int PDFfactor) const { 
  switch(PDFfactor) {
  case 0:
    return 8./3.*log(z); 
  case 1:
    return -8./3./z;
  case 2:
    return 8./3.*log(z/(1.-z));
  case 3:
  default:
    throw Exception() << "QtoGQSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double QtoGQSplitFn::invIntegOverP(const double r, unsigned int PDFfactor) const {
  switch(PDFfactor) {
  case 0:
    return exp(3.*r/8.); 
  case 1:
    return -8./3./r;
  case 2:
    return 1./(1.+exp(-3.*r/8.));
  case 3:
  default:
    throw Exception() << "QtoGQSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

void QtoGQSplitFn::colourConnection(tShowerParticlePtr parent,
				    tShowerParticlePtr first,
				    tShowerParticlePtr second,
				    const bool back) const {
  if(!back) {
    ColinePair cparent = ColinePair(parent->colourLine(), 
				    parent->antiColourLine());
    // ensure input consistency
    assert(( cparent.first && !cparent.second) || 
	   (!cparent.first &&  cparent.second));
    // q -> g q
    if(cparent.first) {
      ColinePtr newline=new_ptr(ColourLine());
      cparent.first->addColoured(first);
      newline->addColoured    (second);
      newline->addAntiColoured( first);
    }
    // qbar -> g qbar
    else {
      ColinePtr newline=new_ptr(ColourLine());
      cparent.second->addAntiColoured(first);
      newline->addColoured    ( first);
      newline->addAntiColoured(second);
    }
  }
  else {
    ColinePair cfirst = ColinePair(first->colourLine(), 
				   first->antiColourLine());
    // ensure input consistency
    assert(cfirst.first&&cfirst.second);
    // q -> g q
    if(parent->id()>0) {
      cfirst.first ->addColoured(parent);
      cfirst.second->addColoured(second);
    }
    else {
      cfirst.first ->addAntiColoured(second);
      cfirst.second->addAntiColoured(parent);
    }
  }
}

bool QtoGQSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[2]||ids[1]!=ParticleID::g) return false;
  tcPDPtr q=getParticleData(ids[0]);
  return q->iSpin()==PDT::Spin1Half&&(q->iColour()==PDT::Colour3||
				      q->iColour()==PDT::Colour3bar);
}

