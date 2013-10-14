// -*- C++ -*-
//
// MatchboxHybridAmplitude.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxHybridAmplitude class.
//

#include "MatchboxHybridAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

using namespace Herwig;

MatchboxHybridAmplitude::MatchboxHybridAmplitude() {}

MatchboxHybridAmplitude::~MatchboxHybridAmplitude() {}

IBPtr MatchboxHybridAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxHybridAmplitude::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxHybridAmplitude::isConsistent() const {
  return
    treeLevelAmplitude()->sortOutgoing() == 
    oneLoopAmplitude()->sortOutgoing() &&
    treeLevelAmplitude()->hasInitialAverage() == 
    oneLoopAmplitude()->hasInitialAverage() &&
    treeLevelAmplitude()->hasFinalStateSymmetry() == 
    oneLoopAmplitude()->hasFinalStateSymmetry() &&
    !treeLevelAmplitude()->isOLPTree() &&
    !treeLevelAmplitude()->isOLPLoop() &&
    oneLoopAmplitude()->haveOneLoop() &&
    treeLevelAmplitude()->orderInGs() ==
    oneLoopAmplitude()->orderInGs() &&
    treeLevelAmplitude()->orderInGem() ==
    oneLoopAmplitude()->orderInGem() &&
    !(treeLevelAmplitude()->nDimAdditional() != 0 &&
      oneLoopAmplitude()->nDimAdditional() != 0);
}

void MatchboxHybridAmplitude::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  treeLevelAmplitude()->prepareAmplitudes(me);
}

void MatchboxHybridAmplitude::prepareOneLoopAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  oneLoopAmplitude()->prepareOneLoopAmplitudes(me);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxHybridAmplitude::persistentOutput(PersistentOStream & os) const {
  os << theTreeLevelAmplitude << theOneLoopAmplitude;
}

void MatchboxHybridAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> theTreeLevelAmplitude >> theOneLoopAmplitude;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxHybridAmplitude,Herwig::MatchboxAmplitude>
  describeHerwigMatchboxHybridAmplitude("Herwig::MatchboxHybridAmplitude", "HwMatchbox.so");

void MatchboxHybridAmplitude::Init() {

  static ClassDocumentation<MatchboxHybridAmplitude> documentation
    ("MatchboxHybridAmplitude unifies two amplitude objects to "
     "provide tree and one-loop matrix elements.");


  static Reference<MatchboxHybridAmplitude,MatchboxAmplitude> interfaceTreeLevelAmplitude
    ("TreeLevelAmplitude",
     "Set the tree level amplitude to be used.",
     &MatchboxHybridAmplitude::theTreeLevelAmplitude, false, false, true, false, false);

  static Reference<MatchboxHybridAmplitude,MatchboxAmplitude> interfaceOneLoopAmplitude
    ("OneLoopAmplitude",
     "Set the tree level amplitude to be used.",
     &MatchboxHybridAmplitude::theOneLoopAmplitude, false, false, true, false, false);

}

