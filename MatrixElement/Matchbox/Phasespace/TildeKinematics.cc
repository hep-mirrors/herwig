// -*- C++ -*-
//
// TildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TildeKinematics class.
//

#include "TildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Rebinder.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

TildeKinematics::TildeKinematics() 
 : HandlerBase() {}

TildeKinematics::~TildeKinematics() {}

Energy TildeKinematics::lastScale() const {
  if ( ( theDipole->bornEmitter() < 2 && theDipole->bornSpectator() > 1 ) ||
       ( theDipole->bornEmitter() > 1 && theDipole->bornSpectator() < 2 ) ) {
    return -(bornEmitterMomentum()-bornSpectatorMomentum()).m();
  }
  return (bornEmitterMomentum()+bornSpectatorMomentum()).m();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void TildeKinematics::rebind(const TranslationMap & trans) {
  theDipole = trans.translate(theDipole);
  HandlerBase::rebind(trans);
}

IVector TildeKinematics::getReferences() {
  IVector ret = HandlerBase::getReferences();
  ret.push_back(theDipole);
  return ret;
}



void TildeKinematics::persistentOutput(PersistentOStream & os) const {
  os << theDipole << theRealXComb << theBornXComb
     << ounit(theBornEmitterMomentum,GeV) << ounit(theBornSpectatorMomentum,GeV);
}

void TildeKinematics::persistentInput(PersistentIStream & is, int) {
  is >> theDipole >> theRealXComb >> theBornXComb
     >> iunit(theBornEmitterMomentum,GeV) >> iunit(theBornSpectatorMomentum,GeV);
}

void TildeKinematics::Init() {

  static ClassDocumentation<TildeKinematics> documentation
    ("TildeKinematics is the base class for the 'tilde' "
     "kinematics being used for subtraction terms in the "
     "formalism of Catani and Seymour.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<TildeKinematics,HandlerBase>
describeTildeKinematics("Herwig::TildeKinematics", "Herwig.so");
