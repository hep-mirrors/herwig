// -*- C++ -*-
//
// HEJMEPP2Jets.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HEJMEPP2Jets class.
//

#include "HEJMEPP2Jets.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/DiagramDrawer.h"

using namespace Herwig;

HEJMEPP2Jets::HEJMEPP2Jets() {}

HEJMEPP2Jets::~HEJMEPP2Jets() {}

IBPtr HEJMEPP2Jets::clone() const {
  return new_ptr(*this);
}

IBPtr HEJMEPP2Jets::fullclone() const {
  return new_ptr(*this);
}

void HEJMEPP2Jets::getDiagrams() const {

  tcPDPtr g = getParticleData(ParticleID::g);

  Tree2toNDiagram diag(nGluons()+3); // n + 1 t-channel props and two incoming

  diag.operator,(firstIncoming());
  for ( unsigned int i = 0; i < nGluons()+1; ++i ) {
    diag.operator,(g);
  }
  diag.operator,(secondIncoming());

  diag.operator,(1); diag.operator,(firstIncoming());

  for ( unsigned int i = 2; i < nGluons()+2; ++i ) {
    diag.operator,(i); diag.operator,(g);
  }

  diag.operator,(nGluons()+2); diag.operator,(secondIncoming());

  diag.operator,(-1);

  add(new_ptr(diag));

  /*
  cerr << name() << " generated diagram:\n";
  DiagramTester::drawDiag(cerr,diag);
  cerr << "\n" << flush;
  */

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HEJMEPP2Jets::persistentOutput(PersistentOStream &) const {
}

void HEJMEPP2Jets::persistentInput(PersistentIStream &, int) {
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<HEJMEPP2Jets,Herwig::HEJMEBase>
  describeHerwigHEJMEPP2Jets("Herwig::HEJMEPP2Jets", "HwHEJ.so");

void HEJMEPP2Jets::Init() {

  static ClassDocumentation<HEJMEPP2Jets> documentation
    ("HEJMEPP2Jets interfaces to HEJ jet production.");

}

