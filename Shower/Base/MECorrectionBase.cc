// -*- C++ -*-
//
// MECorrectionBase.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MECorrectionBase class.
//

#include "MECorrectionBase.h"
#include "Evolver.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "KinematicsReconstructor.h"
#include "PartnerFinder.h"

using namespace Herwig;

AbstractClassDescription<MECorrectionBase> 
MECorrectionBase::initMECorrectionBase;
// Definition of the static class description member.

void MECorrectionBase::Init() {

  static ClassDocumentation<MECorrectionBase> documentation
    ("The MECorrectionBase class is base class for"
     " implementation of the classic"
     " matrix element correction in Herwig++");

  static Reference<MECorrectionBase,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &MECorrectionBase::_alpha, false, false, true, false, false);

}

void MECorrectionBase::persistentOutput(PersistentOStream & os) const {
  os << _alpha;
}

void MECorrectionBase::persistentInput(PersistentIStream & is, int) {
  is >> _alpha;
}

tEvolverPtr MECorrectionBase::evolver() const {
  return _evolver;
}

void MECorrectionBase::evolver(tEvolverPtr in) {
  _evolver=in;
}
