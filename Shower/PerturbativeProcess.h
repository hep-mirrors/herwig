// -*- C++ -*-
//
// PerturbativeProcess.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PerturbativeProcess_H
#define HERWIG_PerturbativeProcess_H

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Config/Pointers.h"
#include "PerturbativeProcess.fh"

namespace Herwig {
using namespace ThePEG;

/**
 *  Simple struct for hard processes and decays
 */
class PerturbativeProcess : public Base {

public:

  /**
   *  The incoming particles
   */
  vector<pair<PPtr,tPerturbativeProcessPtr> > & incoming() {
    return incoming_;
  }

  /**
   *  The outgoing particles
   */
  vector<pair<PPtr, PerturbativeProcessPtr> > & outgoing() {
    return outgoing_;
  }

protected:

  /**
   *  The incoming particles
   */
  vector<pair<PPtr,tPerturbativeProcessPtr> > incoming_;

  /*
   *   The outgoing particles
   */
  vector<pair<PPtr,PerturbativeProcessPtr> > outgoing_;

  /**
   *   The subprocess
   */
  tSubProPtr subprocess_;

  /**
   *  The XComb
   */
  XCombPtr xcomb_;
};

/**
 *  Typedef for map of PerturbativeProcess for decays
 */
typedef multimap<Energy,PerturbativeProcessPtr,std::greater<Energy> > DecayProcessMap;

}

#endif
