// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeBodyOnShellME class.
//
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceChannel.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThreeBodyOnShellME.h"

// the functions for the integrands of the width
namespace Herwig {
using namespace Genfun;

FUNCTION_OBJECT_IMP(ThreeBodyOnShellME)

ThreeBodyOnShellME::
 ThreeBodyOnShellME(DecayIntegratorPtr in,int inmode){_decayer=in;_mode=inmode;}
  
  
ThreeBodyOnShellME::~ThreeBodyOnShellME() {
}
  
ThreeBodyOnShellME::
ThreeBodyOnShellME(const ThreeBodyOnShellME & right) {  }
  
unsigned int ThreeBodyOnShellME::dimensionality() const {return 7;}

double ThreeBodyOnShellME::operator ()(const Argument & a) const 
 {return _decayer->threeBodyMatrixElement(_mode,a[0],a[1],a[2],a[3],a[4],a[5],a[6]);}

}



