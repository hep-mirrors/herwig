// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlpha class.
//

#include "ShowerAlpha.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h" 

using namespace Herwig;


ShowerAlpha::~ShowerAlpha() {}


void ShowerAlpha::persistentOutput(PersistentOStream & os) const {
  os << _scaleFactor;
}


void ShowerAlpha::persistentInput(PersistentIStream & is, int) {
  is >> _scaleFactor;
}


AbstractClassDescription<ShowerAlpha> ShowerAlpha::initShowerAlpha;
// Definition of the static class description member.

void ShowerAlpha::setSV(ShowerVarsPtr scp) {
  _pointerShowerVariables = scp; 
}


void ShowerAlpha::Init() {

  static ClassDocumentation<ShowerAlpha> documentation
    ("This is the abstract class from which the various types of running alphas.",
     "inherit from.");

  static Parameter<ShowerAlpha,double> interfaceShowerAlpha 
    ("ScaleFactor", "Factor that multiplies the scale argument, mu, of the running alpha.",
     &ShowerAlpha::_scaleFactor, 0, 1.0 , 0.0 , 10.0,false,false,false);


//   static Reference<ShowerAlpha,ShowerVariables> 
//     interfaceShowerVariables("ShowerVariables", 
// 			       "A reference to the ShowerVariables object", 
// 			       &Herwig::ShowerAlpha::_pointerShowerVariables,
// 			       false, false, true, false);
}

