// -*- C++ -*-
//
// BtoSGammaDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_BtoSGammaDecayer_H
#define HERWIG_BtoSGammaDecayer_H
//
// This is the declaration of the BtoSGammaDecayer class.
//

#include "PartonicDecayerBase.h"
#include "Herwig/Decay/FormFactors/BtoSGammaHadronicMass.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the BtoSGammaDecayer class.
 *
 */
class BtoSGammaDecayer: public PartonicDecayerBase {

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;

  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<BtoSGammaDecayer> initBtoSGammaDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BtoSGammaDecayer & operator=(const BtoSGammaDecayer &) = delete;

private:

  /**
   *  Pointer to the object which generates the hadronic mass spectrum
   */
  BtoSGammaHadronicMassPtr _hadronicmass;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BtoSGammaDecayer. */
template <>
struct BaseClassTrait<Herwig::BtoSGammaDecayer,1> {
  /** Typedef of the first base class of BtoSGammaDecayer. */
  typedef Herwig::PartonicDecayerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BtoSGammaDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BtoSGammaDecayer>
  : public ClassTraitsBase<Herwig::BtoSGammaDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BtoSGammaDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the BtoSGammaDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwPartonicDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_BtoSGammaDecayer_H */
