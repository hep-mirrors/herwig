// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlpha class.
//

#include "ShowerAlpha.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/Interface/Parameter.h"
#include "Pythia7/Interface/Reference.h" 

using namespace Herwig;


ShowerAlpha::~ShowerAlpha() {}


void ShowerAlpha::persistentOutput(PersistentOStream & os) const {
  os << _pointerShowerConstrainer
     << _scaleFactor;
}


void ShowerAlpha::persistentInput(PersistentIStream & is, int) {
  is >> _pointerShowerConstrainer
     >> _scaleFactor;
}


AbstractClassDescription<ShowerAlpha> ShowerAlpha::initShowerAlpha;
// Definition of the static class description member.


void ShowerAlpha::Init() {

  static ClassDocumentation<ShowerAlpha> documentation
    ("This is the abstract class from which the various types of running alphas.",
     "inherit from.");

  static Parameter<ShowerAlpha,double> interfaceShowerAlpha 
    ("ScaleFactor", "Factor that multiplies the scale argument, mu, of the running alpha.",
     &ShowerAlpha::_scaleFactor, 0, 1.0 , 0.0 , 10.0);

  static Reference<ShowerAlpha,ShowerConstrainer> 
    interfaceShowerConstrainer("ShowerConstrainer", 
			       "A reference to the ShowerConstrainer object", 
			       &Herwig::ShowerAlpha::_pointerShowerConstrainer,
			       false, false, true, false);
}

