// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoGGSplitFn class.
//

#include "GtoGGSplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GtoGGSplitFn.tcc"
#endif


using namespace Herwig;

GtoGGSplitFn::~GtoGGSplitFn() {}

NoPIOClassDescription<GtoGGSplitFn> GtoGGSplitFn::initGtoGGSplitFn;
// Definition of the static class description member.

void GtoGGSplitFn::Init() {

  static ClassDocumentation<GtoGGSplitFn> documentation
    ("There is no documentation for the GtoGGSplitFn class");

}


double GtoGGSplitFn::P(const double z, const Energy2 qtilde2, const IdList &){
  // ACHTUNG! factor two included by hand! 
  double val = 2.*3.*sqr(1.-z*(1.-z))/(z*(1.-z));
  //double val = 3.*sqr(1.-z*(1.-z))/(z*(1.-z));
  // Here we write the LO splitting function P(z) for g -> gg
  // splittings that is well-known from the text books 

  // (this is historically important! the first physics - two years
  // after the birth of the project - in the Herwig++ shower! Alberto
  // & Stefan, 25/04/2002).
  return val;
}

double GtoGGSplitFn::overestimateP(const double z, const IdList &) {
  // ACHTUNG! factor two included by hand! 
  return 2.*3.*(1/z + 1/(1.-z)); 
  //return 3.*(1/z + 1/(1.-z)); 
}

double GtoGGSplitFn::integOverP(const double z) {
  return 2.*3.*log(z/(1.-z)); 
  //return 3.*log(z/(1.-z)); 
}


double GtoGGSplitFn::invIntegOverP(const double r) {
  return exp(r/3./2.)/(1.+exp(r/3./2.)); 
  //return exp(r/3.)/(1.+exp(r/3.)); 
} 


void GtoGGSplitFn::colourConnection(const ShoColinePair &parent,
				     ShoColinePair &first,
				     ShoColinePair &second) {

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

