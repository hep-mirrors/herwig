// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplitFun1to2 class.
//

#include "SplitFun1to2.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


SplitFun1to2::~SplitFun1to2() {}


AbstractClassDescription<SplitFun1to2> SplitFun1to2::initSplitFun1to2;
// Definition of the static class description member.


void SplitFun1to2::Init() {

  static ClassDocumentation<SplitFun1to2> documentation
    ("This is the abstract class from which all 1->2 splitting function ",
     "classes inherit from.");

}


Complex SplitFun1to2::overestimateFullFun( const double z, const double phi ) {
  return fullFun( z, Energy2(), phi );
}


Complex SplitFun1to2::overestimateIntegratedFun( const double z ) {
  return integratedFun( z, Energy2() );
}


Complex SplitFun1to2::
overestimateFullFunWithHelicities( const double z, const double phi,
				  const int h0, const int h1, const int h2 ) {
  return fullFunWithHelicities( z, Energy2(), phi, h0, h1, h2 );
}


Complex SplitFun1to2::
overestimateIntegratedFunWithHelicities( const double z,
					const int h0, const int h1, const int h2 ) {
  return integratedFunWithHelicities( z, Energy2(), h0, h1, h2 );
}


Complex SplitFun1to2::integOverIntegratedFun(const double z) {
  return 0.;
} 


Complex SplitFun1to2::invIntegOverIntegratedFun(const double r) {
  return 0.;
}

