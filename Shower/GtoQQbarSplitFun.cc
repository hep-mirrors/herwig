// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoQQbarSplitFun class.
//

#include "GtoQQbarSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


GtoQQbarSplitFun::~GtoQQbarSplitFun() {}


AbstractClassDescription<GtoQQbarSplitFun> GtoQQbarSplitFun::initGtoQQbarSplitFun;
// Definition of the static class description member.


void GtoQQbarSplitFun::Init() {

  static ClassDocumentation<GtoQQbarSplitFun> documentation
    ("This abstract class defines the exact LO splitting function G->QQbar.");

}


Complex GtoQQbarSplitFun::fullFun( const double z, const Energy2 qtilde2, const double phi ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER G->QQbar

  return val;

}


Complex GtoQQbarSplitFun::integratedFun( const double z, const Energy2 qtilde2 ) {

  double val = 1./2.*(sqr(z) + sqr(1.-z));

  // ***ACHTUNG*** generally this will depend on the q, qbar masses
  // and the scale qtilde2!  may be dangerous in the case of heavy
  // quark production that wouldn't be kinematically allowed!

  return val;

}


Complex GtoQQbarSplitFun::fullFunWithHelicities( const double z, const Energy2 qtilde2, const double phi,
						 const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER G->QQbar
  
  return val;

}


Complex GtoQQbarSplitFun::integratedFunWithHelicities( const double z, const Energy2 qtilde2, const int h0, const int h1, const int h2 ) {

  double val = 0.0;

  //***LOOKHERE*** WRITE THE CODE FOR THE LEADING ORDER G->QQbar

  return val;

}


Complex GtoQQbarSplitFun::overestimateIntegratedFun( const double z ) {
  return 1./2.; 
}


Complex GtoQQbarSplitFun::integOverIntegratedFun(const double z) {
  return z/2.; 
}


Complex GtoQQbarSplitFun::invIntegOverIntegratedFun(const double r) {
  return 2.*r; 
}

