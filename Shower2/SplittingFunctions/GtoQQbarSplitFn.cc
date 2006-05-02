// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoQQbarSplitFn class.
//

#include "GtoQQbarSplitFn.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GtoQQbarSplitFn.tcc"
#endif


using namespace Herwig;

GtoQQbarSplitFn::~GtoQQbarSplitFn() {}

NoPIOClassDescription<GtoQQbarSplitFn> GtoQQbarSplitFn::initGtoQQbarSplitFn;
// Definition of the static class description member.

void GtoQQbarSplitFn::Init() {

  static ClassDocumentation<GtoQQbarSplitFn> documentation
    ("The GtoQQbarSplitFn class implements the splitting function for g->q qbar");

}

double GtoQQbarSplitFn::P(const double z, const Energy2 qtilde2, 
			   const IdList &ids) {
  Energy m = getParticleData(ids[1])->mass();
  double zz = z*(1.-z);
  double term = 2.*sqr(m)/zz/qtilde2; 
  double val = 1./2.*(1.-2.*zz+term);
  return val;
}

double GtoQQbarSplitFn::overestimateP(const double, const IdList &) {
  return 1./2.; 
}

double GtoQQbarSplitFn::ratioP(const double z, const Energy2 qtilde2, 
			   const IdList &ids) {
  Energy m = getParticleData(ids[1])->mass();
  double zz = z*(1.-z);
  double term = 2.*sqr(m)/zz/qtilde2; 
  double val = (1.-2.*zz+term);
  return val;
}

double GtoQQbarSplitFn::integOverP(const double z) { 
  return z/2.; 
}


double GtoQQbarSplitFn::invIntegOverP(const double r) {
  return 2.*r; 
}

void GtoQQbarSplitFn::colourConnection(const ColinePair &parent,
				       ColinePair &first,
				       ColinePair &second) {

  // Return immediately if the input is inconsistent.
  if ( ! parent.first  ||  ! parent.second ) return;
  
  // Initialize
  first = second = ColinePair();

  // The first branching product is considered to be the quark 
  // and the second the anti-quark. 
  first.first = parent.first;
  second.second = parent.second;

}

bool GtoQQbarSplitFn::accept(const IdList &ids) {
  if(ids.size()!=3) return false;
  if(ids[1]!=-ids[2]||ids[0]!=ParticleID::g) return false;
  tcPDPtr q=getParticleData(ids[1]);
  return q->iSpin()==PDT::Spin1Half&&(q->iColour()==PDT::Colour3||
				      q->iColour()==PDT::Colour3bar);
}
