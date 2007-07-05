// -*- C++ -*-
#ifndef HERWIG_EvtGenDecayer_H
#define HERWIG_EvtGenDecayer_H
//
// This is the declaration of the EvtGenDecayer class.
//

#include "ThePEG/PDT/Decayer.h"
#include "EvtGen.h"
#include "EvtGenDecayer.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The EvtGenDecayer class is designed to allow the EvtGen decay package to be used
 * as a Decayer in the Herwig++ structure.
 *
 * It is a simple wrapper which uses members of the Herwig++ EvtGen class to perform
 * the decay
 *
 * @see EvtGen
 * @see \ref EvtGenDecayerInterfaces "The interfaces"
 * defined for EvtGenDecayer.
 */
class EvtGenDecayer: public Decayer {

public:

  /**
   * The default constructor.
   */
  inline EvtGenDecayer();

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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<EvtGenDecayer> initEvtGenDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EvtGenDecayer & operator=(const EvtGenDecayer &);

private:

  /**
   *  Pointer to the EvtGen interface object
   */
  EvtGenPtr _evtgen;

  /**
   *  Option for how EvtGen is used
   */
  unsigned int _evtopt;

  /**
   *  Perform checks ?
   */
  bool _check;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EvtGenDecayer. */
template <>
struct BaseClassTrait<Herwig::EvtGenDecayer,1> {
  /** Typedef of the first base class of EvtGenDecayer. */
  typedef Decayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EvtGenDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EvtGenDecayer>
  : public ClassTraitsBase<Herwig::EvtGenDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EvtGenDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the EvtGenDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwEvtGen.so"; }
};

/** @endcond */

}

#include "EvtGenDecayer.icc"

#endif /* HERWIG_EvtGenDecayer_H */
