// -*- C++ -*-
//
// SSGSSVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
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
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGSSVertex::SSGSSVertex() : _couplast(0.),_q2last(ZERO) {
  orderInGs(1);
  orderInGem(0);
}

void SSGSSVertex::doinit() {
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(21,ix,-ix);
  }
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(21,ix,-ix);
  }
  VSSVertex::doinit();
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeNoPIOClass<SSGSSVertex,Helicity::VSSVertex>
describeSSGSSVertex("Herwig::SSGSSVertex", "HwSusy.so");

void SSGSSVertex::Init() {

  static ClassDocumentation<SSGSSVertex> documentation
    ("The SSGSSVertex class implements the coupling"
     " of the gluon to the squarks");

}

void SSGSSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2, tcPDPtr) {
  assert(part1->id()==ParticleID::g);
  long isf = abs(part2->id());
  assert( (isf >= 1000001 && isf <= 1000006) || 
	  (isf >= 2000001 && isf <= 2000006) );
  if(q2 != _q2last || _couplast == 0.) {
    _couplast = strongCoupling(q2);
    _q2last = q2;
  }
  norm(_couplast);
}
