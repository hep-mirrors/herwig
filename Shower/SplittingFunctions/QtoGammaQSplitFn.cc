// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoGammaQSplitFn class.
//

#include "QtoGammaQSplitFn.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include <cassert>

using namespace Herwig;

NoPIOClassDescription<QtoGammaQSplitFn> QtoGammaQSplitFn::initQtoGammaQSplitFn;
// Definition of the static class description member.

void QtoGammaQSplitFn::Init() {

  static ClassDocumentation<QtoGammaQSplitFn> documentation
    ("The QtoGammaQSplitFn class implements the splitting function for q -> gamma q");

}

double QtoGammaQSplitFn::P(const double z, const Energy2 t,
		       const IdList &ids, const bool mass) const {
  double val=(2.*(1.-z)+sqr(z))/z;
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=2.*sqr(m)/t;
  }
  double charge=getParticleData(ids[0])->iCharge()/3.;
  return sqr(charge)*val;
}

double QtoGammaQSplitFn::overestimateP(const double z, const IdList & ids) const { 
  double charge=getParticleData(ids[0])->iCharge()/3.;
  return 2.*sqr(charge)/z; 
}

double QtoGammaQSplitFn::ratioP(const double z, const Energy2 t,
			    const IdList &ids,const bool mass) const {
  double val=2.*(1.-z)+sqr(z);
  if(mass) {
    Energy m=getParticleData(ids[0])->mass();
    val -=2.*sqr(m)*z/t;
  }
  return 0.5*val;
}

double QtoGammaQSplitFn::integOverP(const double z, const IdList & ids,
				    unsigned int PDFfactor) const { 
  double charge=getParticleData(ids[0])->iCharge()/3.;
  switch(PDFfactor) {
  case 0:
    return 2.*sqr(charge)*log(z); 
  case 1:
    return -2.*sqr(charge)/z;
  case 2:
    return 2.*sqr(charge)*log(z/(1.-z));
  case 3:
  default:
    throw Exception() << "QtoGammaQSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

double QtoGammaQSplitFn::invIntegOverP(const double r, const IdList & ids,
				       unsigned int PDFfactor) const {
  double charge=getParticleData(ids[0])->iCharge()/3.;
  switch(PDFfactor) {
  case 0:
    return exp(0.5*r/sqr(charge)); 
  case 1:
    return -2.*sqr(charge)/r;
  case 2:
    return 1./(1.+exp(-0.5*r/sqr(charge)));
  case 3:
  default:
    throw Exception() << "QtoGammaQSplitFn::integOverP() invalid PDFfactor = "
		      << PDFfactor << Exception::runerror;
  }
}

bool QtoGammaQSplitFn::accept(const IdList &ids) const {
  if(ids.size()!=3) return false;
  if(ids[0]!=ids[2]||ids[1]!=ParticleID::gamma) return false;
  return getParticleData(ids[0])->charged();
}

void QtoGammaQSplitFn::colourConnection(tShowerParticlePtr parent,
					tShowerParticlePtr first,
					tShowerParticlePtr second,
					const bool back) const {
  if(!back) {
    ColinePair cparent = ColinePair(parent->colourLine(), 
				    parent->antiColourLine());
    // ensure input consistency
    assert(( cparent.first && !cparent.second) || 
	   (!cparent.first &&  cparent.second));
    // q -> gamma q
    if(cparent.first) cparent.first->addColoured(second);
    // qbar -> gamma qbar 
    else              cparent.second->addAntiColoured(second);
  }
  else {
    ColinePtr newline=new_ptr(ColourLine());
    // q -> gamma q
    if(parent->id()>0) {
      newline->addColoured(parent);
      newline->addColoured(second);
    }
    else {
      newline->addAntiColoured(parent);
      newline->addAntiColoured(second);
    }
    cerr << "testing did backward \n" 
	 << *parent << "\n" 
	 << *first << "\n" << *second << "\n";
  }
}
