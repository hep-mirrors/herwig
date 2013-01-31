// -*- C++ -*-
//
// MEMatching.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MEMatching_H
#define Herwig_MEMatching_H
//
// This is the declaration of the MEMatching class.
//

#include "Herwig++/MatrixElement/Matchbox/Matching/ShowerApproximation.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MEMatching implements NLO matching with matrix element correction (aka Powheg).
 *
 */
class MEMatching: public Herwig::ShowerApproximation {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MEMatching();

  /**
   * The destructor.
   */
  virtual ~MEMatching();
  //@}

public:

  /**
   * Return true, if this shower approximation will require a
   * splitting generator
   */
  virtual bool needsSplittingGenerator() const { return true; }

public:

  /**
   * Return the shower approximation to the real emission cross
   * section for the given pair of Born and real emission
   * configurations.
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Return the shower approximation splitting kernel for the given
   * pair of Born and real emission configurations in units of the
   * Born center of mass energy squared, and including a weight to
   * project onto the splitting given by the dipole used.
   */
  virtual double me2() const;

  /**
   * Generate a weight for the given dipole channel
   */
  double channelWeight(int emitter, int emission, int spectator) const;

  /**
   * Generate a normalized weight taking into account all channels
   */
  double channelWeight() const;

  /**
   * Return true, if 'Born screening' is taken into account
   */
  bool bornScreening() const { return theBornScreening; }

  /**
   * Return the power of pt used in the screening term
   */
  double screeningPower() const { return theScreeningPower; }

  /**
   * Return the screening `matrix element squared'
   */
  double screeningME2() const;

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * True, if 'Born screening' should be done
   */
  bool theBornScreening;

  /**
   * The power of pt used in the screening term
   */
  double theScreeningPower;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEMatching & operator=(const MEMatching &);

};

}

#endif /* Herwig_MEMatching_H */
