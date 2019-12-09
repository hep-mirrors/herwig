// -*- C++ -*-
//
// MatchboxFactoryMatcher.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxFactoryMatcher class.
//

#include "MatchboxFactoryMatcher.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxFactoryMatcher::MatchboxFactoryMatcher()
  : theGroup("") {}

MatchboxFactoryMatcher::~MatchboxFactoryMatcher() {}

IBPtr MatchboxFactoryMatcher::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxFactoryMatcher::fullclone() const {
  return new_ptr(*this);
}

PMPtr MatchboxFactoryMatcher::pmclone() const {
  return new_ptr(*this);
}

bool MatchboxFactoryMatcher::check(const ParticleData & data) const {
  return theIds.find(data.id()) != theIds.end();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxFactoryMatcher::persistentOutput(PersistentOStream & os) const {
  os << theGroup << theIds;
}

void MatchboxFactoryMatcher::persistentInput(PersistentIStream & is, int) {
  is >> theGroup >> theIds;
}

void MatchboxFactoryMatcher::doinit() {
  if ( !MatchboxFactory::isMatchboxRun() )
    return;
  if ( !MatchboxFactory::currentFactory() )
    throw Exception()
      << "MatchboxFactoryMatcher::doinit(): No factory object is available in the matcher '"
      << name() << "'" << Exception::runerror;
  map<string,PDVector>::const_iterator grp
    = MatchboxFactory::currentFactory()->particleGroups().find(theGroup);
  if ( grp == MatchboxFactory::currentFactory()->particleGroups().end() )
    throw Exception()
      << "MatchboxFactoryMatcher::doinit(): Particle group '" << theGroup << "' not defined in factory object '"
      << MatchboxFactory::currentFactory()->name() << "'" << Exception::runerror;
  theIds.clear();
  for ( PDVector::const_iterator p = grp->second.begin();
	p != grp->second.end(); ++p )
    theIds.insert((**p).id());
  MatcherBase::doinit();
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxFactoryMatcher,ThePEG::MatcherBase>
  describeHerwigMatchboxFactoryMatcher("Herwig::MatchboxFactoryMatcher", "Herwig.so");

void MatchboxFactoryMatcher::Init() {

  static ClassDocumentation<MatchboxFactoryMatcher> documentation
    ("MatchboxFactoryMatcher matches particles according to MatchboxFactory particle groups");

  static Parameter<MatchboxFactoryMatcher,string> interfaceGroup
    ("Group",
     "Set the group name to match.",
     &MatchboxFactoryMatcher::theGroup, "",
     false, false);

}

