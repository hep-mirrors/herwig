// -*- C++ -*-
//
// StandardSelectors.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_StandardSelectors_H
#define HERWIG_StandardSelectors_H

/**
 * This file contains declarations of
 * standard selector classes for Herwig in addition to those of ThePEG.
 * The classes contain only static
 * functions and are assumed to be used as template arguments to the
 * ParticleSelector class.
 */

#include "ThePEG/EventRecord/SelectorBase.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/ParticleTraits.h"

namespace Herwig {
using namespace ThePEG;

/**
 *  Selector to select weakly decaying B hadrons
 */
struct WeakBHadronSelector: public SelectorBase {

  /**
   * Return true if the particle should be extracted.
   */
  virtual bool check(const Particle &  p) const;

};
}

#endif /* HERWIG_StandardSelectors_H */
