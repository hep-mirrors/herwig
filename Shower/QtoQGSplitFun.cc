// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGSplitFun class.
//

#include "QtoQGSplitFun.h"

using namespace Herwig;


QtoQGSplitFun::~QtoQGSplitFun() {}


double QtoQGSplitFun::fullFun( const double z, const Energy2 qtilde2, const double phi ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QG

  return val;

}


double QtoQGSplitFun::integratedFun( const double z, const Energy2 qtilde2 ) {

  Energy2 m2 = sqr(massEmitter()); 

  m2 = 0;  // ***ACHTUNG*** generally neq 0!

  double val = 4./3.*(1. + sqr(z) - 2.*m2/(qtilde2*z))/(1.-z);

  return val;

}


double QtoQGSplitFun::fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QG
  
  return val;

}


double QtoQGSplitFun::integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QG

  return val;

}


double QtoQGSplitFun::overestimateIntegratedFun( const double z ) {
  return 2.*4./3./(1.-z); 
}


double QtoQGSplitFun::integOverIntegratedFun(const double z) {
  return -8./3.*log(1.-z); 
}


double QtoQGSplitFun::invIntegOverIntegratedFun(const double r) {
  return 1. - exp(- 3.*r/8.); 
}


void QtoQGSplitFun::colourConnection( const ShoColinePair & parentShoColinePair,
				      ShoColinePair & firstProductShoColinePair,
				      ShoColinePair & secondProductShoColinePair ) {

  // Return immediately if the input is inconsistent.
  if ( ( ! parentShoColinePair.first  &&  ! parentShoColinePair.second ) ||
       ( parentShoColinePair.first  &&  parentShoColinePair.second ) ) {
    return;
  }
  
  // Initialize
  firstProductShoColinePair = secondProductShoColinePair = ShoColinePair();

  // The first branching product is considered to be the quark 
  // and the second the gluon. The colour line of the parent
  // is one of the two colour lines of the gluon, whereas the
  // other one of the latter is a new colour line which is
  // also share by the first product (the quark).
  if ( parentShoColinePair.first ) { // the parent is a quark
    secondProductShoColinePair.first = parentShoColinePair.first;
    firstProductShoColinePair.first = secondProductShoColinePair.second 
      = new_ptr( ColourLine() );
  } else if ( parentShoColinePair.second ) { // the parent is an antiquark
    secondProductShoColinePair.second = parentShoColinePair.second;
    firstProductShoColinePair.second = secondProductShoColinePair.first 
      = new_ptr( ColourLine() );
  }

}

