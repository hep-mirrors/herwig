// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeBodyDGammaDs class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceChannel.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThreeBodyDGammaDs.h"

namespace Herwig {
using namespace ThePEG;
using namespace Genfun;

FUNCTION_OBJECT_IMP(ThreeBodyDGammaDs)

ThreeBodyDGammaDs::
 ThreeBodyDGammaDs(DecayIntegratorPtr in,int inmode){_decayer=in;_mode=inmode;}

ThreeBodyDGammaDs::~ThreeBodyDGammaDs() {}


ThreeBodyDGammaDs::
ThreeBodyDGammaDs(const ThreeBodyDGammaDs & right) {  }
  
unsigned int ThreeBodyDGammaDs::dimensionality() const {return 5;}

double ThreeBodyDGammaDs::operator ()(const Argument & a) const 
 {return _decayer->threeBodydGammads(_mode,a[0],a[1],a[2],a[3],a[4]);}

}

