// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGSplitFun class.
//

#include "QtoQGSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


QtoQGSplitFun::~QtoQGSplitFun() {}


AbstractClassDescription<QtoQGSplitFun> QtoQGSplitFun::initQtoQGSplitFun;
// Definition of the static class description member.


void QtoQGSplitFun::Init() {

  static ClassDocumentation<QtoQGSplitFun> documentation
    ("This abstract class defines the exact LO splitting function Q->QG.");

}


Complex QtoQGSplitFun::fullFun( const double z, const Energy2 qtilde2, const double phi ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QG

  return val;

}


Complex QtoQGSplitFun::integratedFun( const double z, const Energy2 qtilde2 ) {

  Energy2 m2 = sqr(massEmitter()); 

  m2 = 0;  // ***ACHTUNG*** generally neq 0!

  double val = 4./3.*(1. + sqr(z) - 2.*m2/(qtilde2*z))/(1.-z);

  return val;

}


Complex QtoQGSplitFun::fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QG
  
  return val;

}


Complex QtoQGSplitFun::integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QG

  return val;

}


Complex QtoQGSplitFun::overestimateIntegratedFun( const double z ) {
  return 2.*4./3./(1.-z); 
}


Complex QtoQGSplitFun::integOverIntegratedFun(const double z) {
  return -8./3.*log(1.-z); 
}


Complex QtoQGSplitFun::invIntegOverIntegratedFun(const double r) {
  return 1. - exp(- 3.*r/8.); 
}


