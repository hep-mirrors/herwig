// -*- C++ -*-
#ifndef HERWIG_StandardSelectors_H
#define HERWIG_StandardSelectors_H

/**
 * @file Utilities/StandardSelectors.h This file contains declarations of
 * standard selector classes for Herwig++ in addition to those of ThePEG.
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
  inline virtual bool check(const Particle &  p) const;

};

#include "StandardSelectors.icc"
}

#endif /* HERWIG_StandardSelectors_H */
