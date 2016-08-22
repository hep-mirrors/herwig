// -*- C++ -*-
//
// DipoleSplittingReweight.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleSplittingReweight class.
//

#include <config.h>
#include "DipoleSplittingReweight.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/DipoleShower/DipoleShowerHandler.h"

using namespace Herwig;

DipoleSplittingReweight::DipoleSplittingReweight() 
  : HandlerBase() {}

DipoleSplittingReweight::~DipoleSplittingReweight() {}

void DipoleSplittingReweight::updateCurrentHandler() {
  if ( ShowerHandler::currentHandler() != theCurrentHandler ) {
    Ptr<ShowerHandler>::tptr sptr = ShowerHandler::currentHandler();
    theCurrentHandler = 
      dynamic_ptr_cast<Ptr<DipoleShowerHandler>::tptr>(sptr);
  }
}

Ptr<DipoleShowerHandler>::tptr DipoleSplittingReweight::currentHandler() const {
  return theCurrentHandler;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleSplittingReweight::persistentOutput(PersistentOStream & ) const {
}

void DipoleSplittingReweight::persistentInput(PersistentIStream &, int) {
}

AbstractClassDescription<DipoleSplittingReweight> DipoleSplittingReweight::initDipoleSplittingReweight;
// Definition of the static class description member.

void DipoleSplittingReweight::Init() {

  static ClassDocumentation<DipoleSplittingReweight> documentation
    ("DipoleSplittingReweight is used by the dipole shower "
     "to reweight splittings from a given dipole splitting kernel.");

}

