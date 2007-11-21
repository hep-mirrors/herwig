// -*- C++ -*-
//
// SSGSGSGVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGSGSGVertex class.
//

#include "SSGSGSGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

void SSGSGSGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS;
}

void SSGSGSGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS;
  _couplast = 0.;
  _q2last = 0.*GeV2;
}

ClassDescription<SSGSGSGVertex> SSGSGSGVertex::initSSGSGSGVertex;
// Definition of the static class description member.

void SSGSGSGVertex::Init() {

  static ClassDocumentation<SSGSGSGVertex> documentation
    ("This class implements the gluon-gluino-gluino vertex");

}

void SSGSGSGVertex::setCoupling(Energy2 q2,tcPDPtr part1,
				tcPDPtr part2,tcPDPtr part3) {
  if((part1->id() == ParticleID::g && part2->id() == ParticleID::SUSY_g &&
      part3->id() == ParticleID::SUSY_g) || 
     (part2->id() == ParticleID::g && part1->id() == ParticleID::SUSY_g &&
      part3->id() == ParticleID::SUSY_g) ||
     (part3->id() == ParticleID::g && part1->id() == ParticleID::SUSY_g &&
      part2->id() == ParticleID::SUSY_g)) {
    if(q2 != _q2last) {
      _couplast = sqrt(4.*Constants::pi*_theSS->alphaS(q2));
      _q2last = q2;
    }
    setNorm(_couplast);
    setLeft(1.);setRight(1.);
  }
  else {
    throw HelicityConsistencyError() 
      << "SSGSGSGVertex::setCoupling() - Incorrect particle found. "
      << part1->id() << "  " << part2->id()
      << "  " << part3->id()
      << Exception::warning;
    setNorm(0.);
    setLeft(0.); setRight(0);
  }
}
