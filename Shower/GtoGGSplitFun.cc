// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoGGSplitFun class.
//

#include "GtoGGSplitFun.h"
#include "Pythia7/Repository/UseRandom.h"

using namespace Herwig;


GtoGGSplitFun::~GtoGGSplitFun() {}


Complex GtoGGSplitFun::fullFun( const double z, const Energy2 qtilde2, const double phi ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER G->GG

  return val;

}


Complex GtoGGSplitFun::integratedFun( const double z, const Energy2 qtilde2 ) {

  double val = 3.*sqr(1.-z*(1.-z))/(z*(1.-z));

  // Here we write the LO splitting function P(z) for g -> gg
  // splittings that is well-known from the text books 

  // (this is historically important! the first physics - two years
  // after the birth of the project - in the Herwig++ shower! Alberto
  // & Stefan, 25/04/2002).
  
  return val;

}


Complex GtoGGSplitFun::fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER G->GG
  
  return val;

}


Complex GtoGGSplitFun::integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER G->GG

  return val;

}


Complex GtoGGSplitFun::overestimateIntegratedFun( const double z ) {
  return 3.*(1/z + 1/(1.-z)); 
}


Complex GtoGGSplitFun::integOverIntegratedFun(const double z) {
  return 3.*log(z/(1.-z)); 
}


Complex GtoGGSplitFun::invIntegOverIntegratedFun(const double r) {
  return exp(r/3.)/(1.+exp(r/3.)); 
} 


void GtoGGSplitFun::colourConnection( const ShoColinePair & parentShoColinePair,
				      ShoColinePair & firstProductShoColinePair,
				      ShoColinePair & secondProductShoColinePair ) {

  // Return immediately if the input is inconsistent.
  if ( ! parentShoColinePair.first  ||  ! parentShoColinePair.second ) return;
  
  // Randomly decide which of the two gluon products take the
  // colour line passing for the colour of the parent gluon
  // (the other will take the one passing for the anticolour of
  //  the parent gluon).
  if ( UseRandom::rndbool() ) {
    firstProductShoColinePair.first = parentShoColinePair.first;
    secondProductShoColinePair.second = parentShoColinePair.second;
    firstProductShoColinePair.second = secondProductShoColinePair.first 
      = new_ptr( ShowerColourLine() );    
  } else {
    firstProductShoColinePair.second = parentShoColinePair.first;
    secondProductShoColinePair.first = parentShoColinePair.second;
    firstProductShoColinePair.first = secondProductShoColinePair.second 
      = new_ptr( ShowerColourLine() );    
  }

}

