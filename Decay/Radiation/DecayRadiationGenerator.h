// -*- C++ -*-
//
// DecayRadiationGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
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
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with no persistent data.
   */
  static AbstractNoPIOClassDescription<DecayRadiationGenerator> initDecayRadiationGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DecayRadiationGenerator & operator=(const DecayRadiationGenerator &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DecayRadiationGenerator. */
template <>
struct BaseClassTrait<Herwig::DecayRadiationGenerator,1> {
  /** Typedef of the first base class of DecayRadiationGenerator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DecayRadiationGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DecayRadiationGenerator>
  : public ClassTraitsBase<Herwig::DecayRadiationGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DecayRadiationGenerator"; }
};

/** @endcond */

}

#endif /* HERWIG_DecayRadiationGenerator_H */
