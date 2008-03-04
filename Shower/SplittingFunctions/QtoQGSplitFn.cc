// -*- C++ -*-
//
// QtoQGSplitFn.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGSplitFn class.
//

#include "QtoQGSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

NoPIOClassDescription<QtoQGSplitFn> QtoQGSplitFn::initQtoQGSplitFn;
// Definition of the static class description member.

void QtoQGSplitFn::Init() {

  static ClassDocumentation<QtoQGSplitFn> documentation
    ("The QtoQGSplitFn class implements the q -> qg splitting function");

}

double QtoQGSplitFn::P(const double z, const Energy2 t,
		       const IdList &ids, const bool mass) const {
  double val = (1. + sqr(z))/(1.-z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();  
    val -= 2.*sqr(m)/t;
  }
  return 4./3.*val;

}

double QtoQGSplitFn::overestimateP(const double z, const IdList &) const
{ 
  return 8./3./(1.-z); 
}

double QtoQGSplitFn::ratioP(const double z, const Energy2 t,
			    const IdList &ids, const bool mass) const
{
  double val = 1. + sqr(z);
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val -= 2.*sqr(m)*(1.-z)/t;
  } 
  return 0.5*val;
} 

double QtoQGSplitFn::integOverP(const double z, const IdList & ,
				unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return -8./3.*log(1.-z);
  case 1:
    return 8./3.*log(z/(1.-z));
  case 2:
    return 8./3./(1.-z);
  case 3:
  default:
    throw Exception() << "QtoQGSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

double QtoQGSplitFn::invIntegOverP(const double r, const IdList & ,
				   unsigned int PDFfactor) const {
  switch (PDFfactor) {
  case 0:
    return 1. - exp(- 3.*r/8.); 
  case 1:
    return 1./(1.-exp(-3.*r/8.));
  case 2:
    return 1.-8./3./r;
  case 3:
  default:
    throw Exception() << "QtoQGSplitFn::invIntegOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  } 
}

void QtoQGSplitFn::colourConnection(tShowerParticlePtr parent,
				    tShowerParticlePtr first,
				    tShowerParticlePtr second,
				    const bool back) const {
  if(!back) {
    ColinePair cparent = ColinePair(parent->colourLine(), 
				    parent->antiColourLine());
    // ensure input consistency
    assert((!cparent.first &&  cparent.second) || 
	   ( cparent.first && !cparent.second));
    // q -> q g
    if(cparent.first) {
      ColinePtr newline=new_ptr(ColourLine());
      cparent.first->addColoured(second);
      newline->addColoured     ( first);
      newline->addAntiColoured (second);
    }
    // qbar -> qbar g
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
    // q -> q g
    if(cfirst.first) {
      ColinePtr newline=new_ptr(ColourLine());
      cfirst.first->addAntiColoured(second);
      newline->addColoured(second);
      newline->addColoured(parent);
    }
    // qbar -> qbar g
    else {
      ColinePtr newline=new_ptr(ColourLine());
      cfirst.second->addColoured(second);
      newline->addAntiColoured(second);
      newline->addAntiColoured(parent);
    }
  }
}

bool QtoQGSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[1]||ids[2]!=ParticleID::g) return false;
  tcPDPtr q=getParticleData(ids[0]);
  return q->iSpin()==PDT::Spin1Half&&(q->iColour()==PDT::Colour3||
				      q->iColour()==PDT::Colour3bar);
}

double QtoQGSplitFn::generatePhi(ShowerParticle & ,ShoKinPtr ,
				 const double , const Energy ,
				 const IdList & , const RhoDMatrix &, 
				 const double ) {
  return 1.;
}

DecayMatrixElement QtoQGSplitFn::
matrixElement(ShowerParticle &, ShoKinPtr,
	      const double z, const Energy qtilde, 
	      const IdList & ids,
	      const RhoDMatrix & mapping, const double phi) {
  // compute the splitting matrix 
  Energy m(getParticleData(ids[0])->mass());
  // temporarily as an array and store in DecayMatrixElement class permanently
  double terma(m/qtilde*sqrt(1.-z)/z);
  double termb(sqrt((1.-sqr(m/(qtilde*z)))/(1.-z)));
  Complex ephi(cos(phi)+Complex(0.,1.)*sin(phi));
  Complex me[2][2][2];
  // +++
  me[1][1][1]=termb*conj(ephi);
  // -++
  me[0][1][1]=0.;
  // ++-
  me[1][1][0]=-z*termb*ephi;
  // -+-
  me[0][1][0]=terma;
  // +-+
  me[1][0][1]=terma;
  // --+
  me[0][0][1]=z*termb*conj(ephi);
  // +--
  me[1][0][0]=0.;
  // ---
  me[0][0][0]=-termb*ephi;
  // compute the decay matrix element
  DecayMatrixElement bme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1);
  for(unsigned int in1=0;in1<2;++in1) {
    for(unsigned int in2=0;in2<2;++in2) {
      for(unsigned int out1=0;out1<2;++out1) {
	for(unsigned int out2=0;out2<2;++out2) {
	  bme(in1,out1,2*out2)+=mapping(in1,in2)*me[in2][out1][out2];
	}
      }
    }
  }
  return bme;
}
