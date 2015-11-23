// -*- C++ -*-
//
// HwDecayHandler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HwDecayHandler_H
#define HERWIG_HwDecayHandler_H
//
// This is the declaration of the HwDecayHandler class.
//
#include "ThePEG/Handlers/DecayHandler.h"
#include "ThePEG/EventRecord/Particle.h"


namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 * The <code>HwDecayHandler</code> is the Herwig decay handler which 
 *  administers the decays of unstable particles in Herwig. It
 * is derived from ThePEG::DecayHandler and includes a different handle
 * method in order to simulate decays including spin correlations.
 *
 * The handle method decays all particles in the current step, including
 * spin correlations. Another feature of the DecayHandler is that it correctly
 * handles mutlistep decays where a Decayer supplys intermediate decay products
 * in addition to the outgoing particles.
 *
 * @see ThePEG::StepHandler
 * @see ThePEG::CollisionHandler
 * @see ThePEG::SubProcessHandler
 * @see ThePEG::DecayHandler
 * 
 */

class HwDecayHandler: public DecayHandler {

public:

  /**
   * Default constructor
   */
  HwDecayHandler()  : DecayHandler(), _newstep(true) 
  {}

public:

  /**
   * Look through all \a tagged particled and decay all unstable ones.
   * @param eh the EventHandler in charge of the generation.
   * @param tagged the vector of particles to consider. If empty, all
   * final state particles in the current Step is considered.
   * @param hint a possible Hint which is ignored in this implementation.
   */
  virtual void handle(EventHandler & eh, const tPVector & tagged,
		      const Hint & hint)
   ;

  /**
   * Perform the decay of one unstable particle.
   * @param parent the particle to be decayed.
   * @param s the Step where decay products are inserted.
   * @throws Veto if the Handler requires the current step to be discarded.
   * @throws Exception if something goes wrong.
   */
  virtual void performDecay(tPPtr parent, Step & s) const
   ;
  
  /**
   * add the decay products of in intermediate particle produced in a decay
   * @param parent the particle which has been decayed.
   * @param s the Step where decay products are inserted.
   * @throws Veto if the Handler requires the current step to be discarded.
   * @throws Exception if something goes wrong.
   */
  void addDecayedParticle(tPPtr parent, Step & s) const
   ;

  /**
   * Standard Init function
   */
  static void Init();

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /**
   *  Develop a stable particle
   */
  void develop(tPPtr particle) const {
    tcSpinPtr hwspin = particle->spinInfo();
    if ( hwspin ) 
      hwspin->develop();
  }

private:

  /**
   * Describe a concrete class with persistent date/
   */
  static ClassDescription<HwDecayHandler> initHwDecayHandler;

  /**
   *  Private and non-existent assignment operator.
   */
  HwDecayHandler & operator=(const HwDecayHandler &);

private:

  /**
   *  Option for adding particles in a new Step
   */
  bool _newstep;

  /**
   *  Particles which should not be decayed
   */
  set<tcPDPtr> _excluded;

  /**
   *  Vector to fill the set as an interface
   */
  vector<PDPtr> _excludedVector;

};
}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * Hw64Decayer.
 */
template <>
struct BaseClassTrait<Herwig::HwDecayHandler,1> {
  /** Typedef of the base class of Hw64Decayer. */
  typedef DecayHandler NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Hw64Decayer class.
 */
template <>
struct ClassTraits<Herwig::HwDecayHandler>: public ClassTraitsBase<Herwig::HwDecayHandler> {
  /** Return the class name. */
  static string className() { return "Herwig::HwDecayHandler"; }
};

/** @endcond */

}

#endif /* HERWIG_HwDecayHandler_H */
