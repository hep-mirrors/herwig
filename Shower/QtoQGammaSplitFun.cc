// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtoQGammaSplitFun class.
//

#include "QtoQGammaSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


QtoQGammaSplitFun::~QtoQGammaSplitFun() {}


AbstractClassDescription<QtoQGammaSplitFun> QtoQGammaSplitFun::initQtoQGammaSplitFun;
// Definition of the static class description member.


void QtoQGammaSplitFun::Init() {

  static ClassDocumentation<QtoQGammaSplitFun> documentation
    ("This abstract class defines the exact LO splitting function Q->QGamma.");

}


Complex QtoQGammaSplitFun::fullFun( const double z, const Energy2 qtilde2, const double phi ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QGamma

  return val;

}


Complex QtoQGammaSplitFun::integratedFun( const double z, const Energy2 qtilde2) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QGamma
  
  return val;

}


Complex QtoQGammaSplitFun::fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QGamma
  
  return val;

}


Complex QtoQGammaSplitFun::integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER Q->QGamma

  return val;

}

