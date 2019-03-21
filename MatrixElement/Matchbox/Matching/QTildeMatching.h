// -*- C++ -*-
//
// QTildeMatching.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_QTildeMatching_H
#define Herwig_QTildeMatching_H
//
// This is the declaration of the QTildeMatching class.
//

#include "Herwig/MatrixElement/Matchbox/Matching/ShowerApproximation.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Shower/QTilde/Default/QTildeFinder.h"
#include "Herwig/Shower/QTilde/Default/QTildeSudakov.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief QTildeMatching implements NLO matching with the default shower.
 *
 */
class QTildeMatching: public Herwig::ShowerApproximation {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  QTildeMatching();

  /**
   * The destructor.
   */
  virtual ~QTildeMatching();
  //@}

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
   * Determine if the configuration is below or above the cutoff.
   */
  virtual void checkCutoff();

  /**
   * Determine all kinematic variables which are not provided by the
   * dipole kinematics; store all shower variables in the respective
   * dipole object for later use.
   */
  virtual void getShowerVariables();

protected:

  /**
   * Return true, if the shower was able to generate an emission
   * leading from the given Born to the given real emission process.
   */
  virtual bool isInShowerPhasespace() const;

  /**
   * Return true, if the shower emission leading from the given Born
   * to the given real emission process would have been generated
   * above the shower's infrared cutoff.
   */
  virtual bool isAboveCutoff() const;

  /**
   * Calculate qtilde^2 and z for the splitting considered
   */
  void calculateShowerVariables() const;

  /**
   * Return the splitting function as a function of the kinematic
   * variables
   */
  double splitFn(const pair<Energy2,double>&) const; 

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeMatching & operator=(const QTildeMatching &) = delete;

  /**
   * The shower handler to be used
   */
  Ptr<ShowerHandler>::ptr theShowerHandler;

  /**
   * The qtilde partner finder for calculating the hard scales
   */
  Ptr<QTildeFinder>::ptr theQTildeFinder;

  /**
   * The qtilde Sudakov to access the cutoff
   */
  Ptr<QTildeSudakov>::ptr theQTildeSudakov;

  /**
   * True, if PDF weight should be corrected for z/x mismatch at the
   * hard phase space boundary
   */
  bool theCorrectForXZMismatch;

};

}

#endif /* Herwig_QTildeMatching_H */
