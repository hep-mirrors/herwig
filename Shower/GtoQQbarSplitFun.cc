// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoQQbarSplitFun class.
//

#include "GtoQQbarSplitFun.h"

using namespace Herwig;


GtoQQbarSplitFun::~GtoQQbarSplitFun() {}


double GtoQQbarSplitFun::fullFun( const double z, const Energy2 qtilde2, const double phi ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER G->QQbar

  return val;

}


double GtoQQbarSplitFun::integratedFun( const double z, const Energy2 qtilde2 ) {

  double val = 1./2.*(sqr(z) + sqr(1.-z));

  // ***ACHTUNG*** generally this will depend on the q, qbar masses
  // and the scale qtilde2!  may be dangerous in the case of heavy
  // quark production that wouldn't be kinematically allowed!

  return val;

}


double GtoQQbarSplitFun::fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi,
						 const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER G->QQbar
  
  return val;

}


double GtoQQbarSplitFun::integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER G->QQbar

  return val;

}


double GtoQQbarSplitFun::overestimateIntegratedFun( const double z ) {
  return 1./2.; 
}


double GtoQQbarSplitFun::integOverIntegratedFun(const double z) {
  return z/2.; 
}


double GtoQQbarSplitFun::invIntegOverIntegratedFun(const double r) {
  return 2.*r; 
}


void GtoQQbarSplitFun::colourConnection( const ShoColinePair & parentShoColinePair,
					 ShoColinePair & firstProductShoColinePair,
					 ShoColinePair & secondProductShoColinePair ) {

  // Return immediately if the input is inconsistent.
  if ( ! parentShoColinePair.first  ||  ! parentShoColinePair.second ) return;
  
  // Initialize
  firstProductShoColinePair = secondProductShoColinePair = ShoColinePair();

  // The first branching product is considered to be the quark 
  // and the second the anti-quark. 
  firstProductShoColinePair.first = parentShoColinePair.first;
  secondProductShoColinePair.second = parentShoColinePair.second;

}
