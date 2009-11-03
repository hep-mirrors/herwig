// -*- C++ -*-
//
// SSGSSVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGSSVertex class.
//

#include "SSGSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGSSVertex::SSGSSVertex() : _couplast(0.),_q2last(ZERO) {
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(21,ix,-ix);
  }
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(21,ix,-ix);
  }
}

void SSGSSVertex::doinit() {
  orderInGs(1);
  orderInGem(0);
  VSSVertex::doinit();
}

NoPIOClassDescription<SSGSSVertex> SSGSSVertex::initSSGSSVertex;
// Definition of the static class description member.

void SSGSSVertex::Init() {

  static ClassDocumentation<SSGSSVertex> documentation
    ("There is no documentation for the SSGSSVertex class");

}

void SSGSSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2, tcPDPtr) {
  long isf(0);
  if(part1->id() == ParticleID::g) {
    isf = abs(part2->id());
  }
  else if(part2->id() == ParticleID::g) {
    isf = abs(part1->id());
  }
  else {
    isf = abs(part1->id());
  }
  if((isf >= 1000001 && isf <= 1000006) || 
     (isf>=2000001 && isf <= 2000006) ) {
    if(q2 != _q2last || _couplast == 0.) {
      _couplast = strongCoupling(q2);
      _q2last = q2;
    }
    norm(_couplast);
}
  else {
    throw  HelicityConsistencyError() 
      << "SSGSSVertex::setCoupling() - Incorrect particle(s) in vertex. "
      << part1->id() << " " << part2->id()
      << Exception::warning;
  }
}
