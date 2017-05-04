// -*- C++ -*-
#ifndef Herwig_EvtGenDecayer_H
#define Herwig_EvtGenDecayer_H
//
// This is the declaration of the EvtGenDecayer class.
//

#include "ThePEG/PDT/Decayer.h"
#include "EvtGenInterface.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The EvtGenDecayer class is designed to allow the EvtGen decay package to be used
 * as a Decayer in the Herwig structure.
 *
 * It is a simple wrapper which uses members of the Herwig EvtGen class to perform
 * the decay
 *
 * @see EvtGenInterface
 * @see \ref EvtGenDecayerInterfaces "The interfaces"
 * defined for EvtGenDecayer.
 */
class EvtGenDecayer: public Decayer {

public:

  /**
   * Standard constructors
   */
  EvtGenDecayer() : check_(0) , evtOpt_(0)
  {}

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
  //@}

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

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  Method to check conservation of charge and momentum in the decay
   *  for testing only
   * @param parent The decaying particle
   */
  void checkDecay(PPtr parent) const;

  /**
   *  Method to rescale the momenta of the decay products if required to
   *  conserve 4-momentum
   */
  bool rescale(const Particle & parent,
	       const ParticleVector & children) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EvtGenDecayer & operator=(const EvtGenDecayer &);

private:

  /**
   *  Pointer to the EvtGen interface object
   */
  EvtGenInterfacePtr evtgen_;

  /**
   *  Perform checks ?
   */
  unsigned int check_;

  /**
   *  Option for how EvtGen is used
   */
  unsigned int evtOpt_;

};

}

#endif /* Herwig_EvtGenDecayer_H */
