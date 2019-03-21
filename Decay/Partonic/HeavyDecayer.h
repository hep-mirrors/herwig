// -*- C++ -*-
//
// HeavyDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HeavyDecayer_H
#define HERWIG_HeavyDecayer_H
// This is the declaration of the HeavyDecayer class.

#include "PartonicDecayerBase.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 *  This class is designed for the partonic decay of a bottom or charm mesons
 *  and baryons and is intended to be the same as that in HERWIG6.4. 
 *  Only four body partonic decays are supported.
 *
 *  Two types of matrix element are supported for this decay
 *
 *  - MECode=0   flat phase space
 *  - MECode=100 V-A matrix element for the heavy quark decay in the spectator model.
 *
 * The class decays a particle based on the DecayMode given.
 * The basic idea is that the heavy parton in the heavy meson will decay
 * weakly while the other parton will be more or less unchanged. This produces
 * 4 partons, the spectator plus the result of the weak decay. The W then
 * decays again into two more partons,
 * e.g. a decay of \f$B^0\to d,\bar{c},\bar{d}, u \f$.
 *
 *    \f$\bar{b}\to\bar{c}\f$ (colour connected to spectator)
 *   and \f$W^+\to\bar{d}u\f$  (the W decay products are colour connected)
 *
 * The resulting partons then need to be hadronized and decayed again.
 *
 * @see QuarkoniumDecayer
 * @see Hw64Decayer
 * @see Decayer
 * 
 */
class HeavyDecayer: public PartonicDecayerBase {

public:

  /**
   * Default constructor
   */
  HeavyDecayer();

  /**
   * Check if this decayer can perfom the decay for a particular mode
   * @param parent The decaying particle
   * @param children The decay products
   * @return true If this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;
  
  /**
   *  Perform the decay of the particle to the specified decay products
   * @param parent The decaying particle
   * @param children The decay products
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

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
   * Weighting of phase space for V-A matrix elements
   */
  static double VAWt(Energy2, Energy2, Energy2, InvEnergy4);

private:

  /**
   *  Describe a concrete class with persistent data.
   */
  static ClassDescription<HeavyDecayer> initHeavyDecayer;

  /**
   *  Private and non-existent assignment operator.
   */
  const HeavyDecayer & operator=(const HeavyDecayer &) = delete;

private:

  /**
   *  The code for the matrix element being used.
   */
  int MECode;
};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * HeavyDecayer.
 */
template <>
struct BaseClassTrait<Herwig::HeavyDecayer,1> {
  /** Typedef of the base class of HeavyDecayer. */
  typedef Herwig::PartonicDecayerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * HeavyDecayer class.
 */
template <>
struct ClassTraits<Herwig::HeavyDecayer>: public ClassTraitsBase<Herwig::HeavyDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig::HeavyDecayer"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwPartonicDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_HeavyDecayer_H */
