// -*- C++ -*-
//
// ShowerAlpha.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlpha class.
//

#include "ShowerAlpha.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerAlpha.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void ShowerAlpha::persistentOutput(PersistentOStream & os) const {
  os << _scaleFactor;
}

void ShowerAlpha::persistentInput(PersistentIStream & is, int) {
  is >> _scaleFactor;
}

AbstractClassDescription<ShowerAlpha> ShowerAlpha::initShowerAlpha;
// Definition of the static class description member.

void ShowerAlpha::Init() {

  static ClassDocumentation<ShowerAlpha> documentation
    ("This is the abstract class from which the various types of running alphas.",
     "inherit from.");

  static Parameter<ShowerAlpha,double> interfaceShowerAlpha 
    ("ScaleFactor", "Factor that multiplies the scale argument, mu, "
     "of the running alpha.",
     &ShowerAlpha::_scaleFactor, 0, 1.0 , 0.0 , 10.0,false,false,false);

}

