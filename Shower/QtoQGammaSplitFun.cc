// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGammaSplitFun class.
//

#include "QtoQGammaSplitFun.h"

using namespace Herwig;


QtoQGammaSplitFun::~QtoQGammaSplitFun() {}


double QtoQGammaSplitFun::fullFun( const double z, const Energy2 qtilde2, const double phi ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QGamma

  return val;

}


double QtoQGammaSplitFun::integratedFun( const double z, const Energy2 qtilde2) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QGamma
  
  return val;

}


double QtoQGammaSplitFun::fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QGamma
  
  return val;

}


double QtoQGammaSplitFun::integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QGamma

  return val;

}


void QtoQGammaSplitFun::colourConnection( const ShoColinePair & parentShoColinePair,
					  ShoColinePair & firstProductShoColinePair,
					  ShoColinePair & secondProductShoColinePair ) {

  // Return immediately if the input is inconsistent.
  if ( ( ! parentShoColinePair.first  &&  ! parentShoColinePair.second ) ||
       ( parentShoColinePair.first  &&  parentShoColinePair.second ) ) {
    return;
  }
  
  firstProductShoColinePair = parentShoColinePair;
  secondProductShoColinePair = ShoColinePair();

}
