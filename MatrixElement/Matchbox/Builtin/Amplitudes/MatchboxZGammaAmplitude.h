// -*- C++ -*-
//
// MatchboxZGammaAmplitude.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxZGammaAmplitude_H
#define Herwig_MatchboxZGammaAmplitude_H
//
// This is the declaration of the MatchboxZGammaAmplitude class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MatchboxZGammaAmplitude class.
 *
 * @see \ref MatchboxZGammaAmplitudeInterfaces "The interfaces"
 * defined for MatchboxZGammaAmplitude.
 */
class MatchboxZGammaAmplitude: public MatchboxAmplitude {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxZGammaAmplitude();

  /**
   * The destructor.
   */
  virtual ~MatchboxZGammaAmplitude();
  //@}

public:

  /**
   * Return true, if the Z contribution should be taken into account
   */
  bool includeZ() const { return theIncludeZ; }

  /**
   * Return true, if the gamma contribution should be taken into account
   */
  bool includeGamma() const { return theIncludeGamma; }

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxZGammaAmplitude & operator=(const MatchboxZGammaAmplitude &) = delete;

  /**
   * True, if the Z contribution should be taken into account
   */
  bool theIncludeZ;

  /**
   * True, if the gamma contribution should be taken into account
   */
  bool theIncludeGamma;

};

}

#endif /* Herwig_MatchboxZGammaAmplitude_H */
