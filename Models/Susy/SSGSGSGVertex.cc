// -*- C++ -*-
//
// SSGSGSGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGSGSGVertex class.
//

#include "SSGSGSGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGSGSGVertex::SSGSGSGVertex() : _couplast(0.),_q2last(ZERO) {
  orderInGs(1);
  orderInGem(0);
}

// Static variable needed for the type description system in ThePEG.
DescribeNoPIOClass<SSGSGSGVertex,FFVVertex>
describeHerwigSSGSGSGVertex("Herwig::SSGSGSGVertex", "HwSusy.so");

void SSGSGSGVertex::Init() {

  static ClassDocumentation<SSGSGSGVertex> documentation
    ("This class implements the gluon-gluino-gluino vertex");

}

void SSGSGSGVertex::setCoupling(Energy2 q2,tcPDPtr part1,
				tcPDPtr part2,tcPDPtr part3) {
  assert(part1->id()==ParticleID::SUSY_g &&
	 part2->id()==ParticleID::SUSY_g &&
	 part3->id() == ParticleID::g);
  if(q2 != _q2last || _couplast==0.) {
    _couplast = strongCoupling(q2);
    _q2last = q2;
  }
  norm(_couplast);
  left(1.);
  right(1.);
}

void SSGSGSGVertex::doinit() {
  addToList(1000021, 1000021, 21);
  FFVVertex::doinit();
}
