// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoGGSplitFun class.
//

#include "GtoGGSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


GtoGGSplitFun::~GtoGGSplitFun() {}


AbstractClassDescription<GtoGGSplitFun> GtoGGSplitFun::initGtoGGSplitFun;
// Definition of the static class description member.


void GtoGGSplitFun::Init() {

  static ClassDocumentation<GtoGGSplitFun> documentation
    ("This abstract class defines the exact LO splitting function G->GG.");

}


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


