// -*- C++ -*-
//
// PhasespaceCouplings.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhasespaceCouplings class.
//

#include "PhasespaceCouplings.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PhasespaceCouplings::PhasespaceCouplings() {}

PhasespaceCouplings::~PhasespaceCouplings() {}

IBPtr PhasespaceCouplings::clone() const {
  return new_ptr(*this);
}

IBPtr PhasespaceCouplings::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void PhasespaceCouplings::persistentOutput(PersistentOStream & os) const {
  os << theCouplings.size();
  for ( map<LTriple,double>::const_iterator cit = 
	  theCouplings.begin(); cit != theCouplings.end(); ++cit )
    os << cit->first << cit->second;
}

void PhasespaceCouplings::persistentInput(PersistentIStream & is, int) {
  theCouplings.clear();
  size_t size;
  LTriple k;
  is >> size;
  while ( size-- && is ) {
    is >> k;
    is >> theCouplings[k];
  }
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PhasespaceCouplings,HandlerBase>
  describeHerwigPhasespaceCouplings("Herwig::PhasespaceCouplings", "Herwig.so");

void PhasespaceCouplings::Init() {

  static ClassDocumentation<PhasespaceCouplings> documentation
    ("Store couplings for the phasespace generator.");

}

