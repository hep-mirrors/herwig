// -*- C++ -*-
//
// MPIXSecReweighter.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MPIXSecReweighter class.
//

#include "MPIXSecReweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "MPIHandler.h"

using namespace Herwig;

MPIXSecReweighter::MPIXSecReweighter() 
  : sumWeights(0.), xSec(ZERO) {}

MPIXSecReweighter::~MPIXSecReweighter() {}

void MPIXSecReweighter::
handle(EventHandler & eh, const tPVector &,
       const Hint &) {

  if ( MPIHandler::currentHandler() &&
       dynamic_cast<StandardEventHandler*>(&eh) ) {

    double corr = 1.;
    StandardEventHandler& seh = 
      dynamic_cast<StandardEventHandler&>(eh);
    CrossSection mpiXSec = MPIHandler::currentHandler()->inelasticXSec();
    double weight = seh.currentEvent()->weight();

    if ( weight == 0. )
      return;

    CrossSection next = seh.integratedXSecNoReweight();
    assert(next != ZERO);

    if ( xSec != ZERO ) {      
      corr = (mpiXSec/next)*( 1. +(sumWeights/weight)*(1.-next/xSec) );
      sumWeights += weight;
      xSec = next;
    } else {
      xSec = next;
      sumWeights = weight;
      corr = mpiXSec/xSec;
    }
    seh.reweight(corr);
  }

}

IBPtr MPIXSecReweighter::clone() const {
  return new_ptr(*this);
}

IBPtr MPIXSecReweighter::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MPIXSecReweighter::persistentOutput(PersistentOStream &) const {}

void MPIXSecReweighter::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MPIXSecReweighter,StepHandler>
  describeHerwigMPIXSecReweighter("Herwig::MPIXSecReweighter", "JetCuts.so SimpleKTCut.so HwMPI.so");

void MPIXSecReweighter::Init() {

  static ClassDocumentation<MPIXSecReweighter> documentation
    ("MPIXSecReweighter sets up the proper minimum bias cross section.");

}

