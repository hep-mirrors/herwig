// -*- C++ -*-
//
// DecayRadiationGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DecayRadiationGenerator_H
#define HERWIG_DecayRadiationGenerator_H
//
// This is the declaration of the DecayRadiationGenerator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Decay/DecayIntegrator.fh"
#include "DecayRadiationGenerator.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The DecayRadiationGenerator class is the base class for classes generating
 * QED radiation in particle decays in Herwig. Classes implementing specific
 * algorithms must inherit from this class and implement the virtual generatePhotons
 * member.
 * 
 * @see \ref DecayRadiationGeneratorInterfaces "The interfaces"
 * defined for DecayRadiationGenerator.
 */
class DecayRadiationGenerator: public Interfaced {

public:

  /**
   *  Member to generate the photons in the decay. This must be implemented
   *  in classes inheriting from this one to produce the radiation.
   * @param p The decaying particle
   * @param children The decay products
   * @param decayer The decayer which would normally generate this decay
   * @return The decay products with additional radiation
   */
  virtual ParticleVector generatePhotons(const Particle & p,
					 ParticleVector children,
					 tDecayIntegratorPtr decayer)=0;

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DecayRadiationGenerator & operator=(const DecayRadiationGenerator &) = delete;

};

}

#endif /* HERWIG_DecayRadiationGenerator_H */
