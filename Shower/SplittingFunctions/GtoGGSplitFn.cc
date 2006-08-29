// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoGGSplitFn class.
//

#include "GtoGGSplitFn.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GtoGGSplitFn.tcc"
#endif


using namespace Herwig;

NoPIOClassDescription<GtoGGSplitFn> GtoGGSplitFn::initGtoGGSplitFn;
// Definition of the static class description member.

void GtoGGSplitFn::Init() {

  static ClassDocumentation<GtoGGSplitFn> documentation
    ("The GtoGGSplitFn class implements the g -> gg splitting function");

}

double GtoGGSplitFn::P(const double z, const Energy2, const IdList &,
		       const bool )const {
  double val = 2.*3.*sqr(1.-z*(1.-z))/(z*(1.-z));
  // (this is historically important! the first physics - two years
  // after the birth of the project - in the Herwig++ shower! Alberto
  // & Stefan, 25/04/2002).
  return val;
}

double GtoGGSplitFn::overestimateP(const double z, const IdList &) const {
  return 2.*3.*(1/z + 1/(1.-z)); 
}

double GtoGGSplitFn::ratioP(const double z, const Energy2,
			    const IdList &, const bool ) const {
  return sqr(1.-z*(1.-z));
}

double GtoGGSplitFn::invIntegOverP(const double r) const {
  return exp(r/3./2.)/(1.+exp(r/3./2.)); 
} 

double GtoGGSplitFn::integOverP(const double z) const {
  return 2.*3.*log(z/(1.-z)); 
}

void GtoGGSplitFn::colourConnection(const ColinePair &parent,
				    ColinePair &first,
				    ColinePair &second) const {

  // Return immediately if the input is inconsistent.
  if(!parent.first || !parent.second) return;
  
  // Randomly decide which of the two gluon products take the
  // colour line passing for the colour of the parent gluon
  // (the other will take the one passing for the anticolour of
  //  the parent gluon).
  if(UseRandom::rndbool()) {
    first.first = parent.first;
    second.second = parent.second;
    first.second = second.first = new_ptr(ColourLine());    
  } else {
    second.first = parent.first;
    first.second = parent.second;
    first.first = second.second = new_ptr(ColourLine());    
  }

}

bool GtoGGSplitFn::accept(const IdList & ids)
const {
  if(ids.size()!=3) return false;
  for(unsigned int ix=0;ix<ids.size();++ix)
    {if(ids[ix]!=ParticleID::g) return false;}
  return true;
}
