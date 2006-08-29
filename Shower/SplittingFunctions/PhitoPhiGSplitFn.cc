// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhitoPhiGSplitFn class.
//

#include "PhitoPhiGSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PhitoPhiGSplitFn.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PhitoPhiGSplitFn::~PhitoPhiGSplitFn() {}

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

double PhitoPhiGSplitFn::overestimateP(const double z, const IdList &) const
{ 
  return 8./3./(1.-z); 
}

double PhitoPhiGSplitFn::ratioP(const double z, const Energy2 t,
				const IdList &ids,const bool mass) const
{ 
  double val = z;
  if(mass) {
    Energy m = getParticleData(ids[0])->mass();
    val-=sqr(m)*(1.-z)/t;
  }
  return val;
} 

double PhitoPhiGSplitFn::integOverP(const double z) const {
  return -8./3.*log(1.-z); 
}

double PhitoPhiGSplitFn::invIntegOverP(const double r) const {
  return 1. - exp(- 3.*r/8.); 
}

void PhitoPhiGSplitFn::colourConnection(const ColinePair &parent,
					ColinePair &first,
					ColinePair &second) const {

  // Return immediately if the input is inconsistent.
  if ((!parent.first && !parent.second) || (parent.first && parent.second)) 
    return;
  
  // Initialize
  first = second = ColinePair();

  // The first branching product is considered to be the squark 
  // and the second the gluon. The colour line of the parent
  // is one of the two colour lines of the gluon, whereas the
  // other one of the latter is a new colour line which is
  // also share by the first product (the squark).
  if(parent.first) { // the parent is a squark
    second.first = parent.first;
    first.first = second.second = new_ptr(ColourLine());
  } else if(parent.second) { // the parent is an antisquark
    second.second = parent.second;
    first.second = second.first = new_ptr(ColourLine());
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
