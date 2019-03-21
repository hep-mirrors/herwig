// -*- C++ -*-
//
// HerwigStrategy.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2008-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_HerwigStrategy_H
#define Herwig_HerwigStrategy_H
// This is the declaration of the HerwigStrategy class.

#include "ThePEG/Repository/Strategy.h"
#include <string>

namespace Herwig {

using namespace ThePEG;

/**
 * The HerwigStrategy class is a sub-class of the Strategy class,
 * simply implementing the correct citation for Herwig in the
 * ClassDocumentation interface.
 *
 * @see Strategy
 * 
 */
class HerwigStrategy: public ThePEG::Strategy {

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

  /**
   * Freeform version string
   */
  static const std::string version;

  /**
   * Version string
   */
  virtual const std::string versionstring() const;

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
   * Describe concrete class without persistent data.
   */
  static NoPIOClassDescription<HerwigStrategy> initHerwigStrategy;

  /**
   *  Private and non-existent assignment operator.
   */
  HerwigStrategy & operator=(const HerwigStrategy &) = delete;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of HerwigStrategy. */
template <>
struct BaseClassTrait<Herwig::HerwigStrategy,1>: public ClassTraitsType {
  /** Typedef of the first base class of HerwigStrategy. */
  typedef Strategy NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  HerwigStrategy class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<Herwig::HerwigStrategy>: public ClassTraitsBase<Herwig::HerwigStrategy> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HerwigStrategy"; }
};

/** @endcond */

}

#endif /* Herwig_HerwigStrategy_H */
