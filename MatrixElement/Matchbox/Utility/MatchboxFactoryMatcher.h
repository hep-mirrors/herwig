// -*- C++ -*-
//
// MatchboxFactoryMatcher.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxFactoryMatcher_H
#define Herwig_MatchboxFactoryMatcher_H
//
// This is the declaration of the MatchboxFactoryMatcher class.
//

#include "ThePEG/PDT/MatcherBase.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxFactoryMatcher matches particles according to
 * MatchboxFatory particle groups
 *
 * @see \ref MatchboxFactoryMatcherInterfaces "The interfaces"
 * defined for MatchboxFactoryMatcher.
 */
class MatchboxFactoryMatcher: public ThePEG::MatcherBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxFactoryMatcher();

  /**
   * The destructor.
   */
  virtual ~MatchboxFactoryMatcher();
  //@}

public:

  /**
   * Check if a particle type meets the criteria.
   */
  virtual bool check(const ParticleData &) const;

  /**
   * Specialized clone method for MatcherBase used by the
   * Repository. A sub class must make sure that also the MatcherBase
   * object corresponding to the complex conjugate of this is cloned.
   */
  virtual PMPtr pmclone() const;

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
  //@}

private:

  /**
   * A pointer to the factory to be used
   */
  Ptr<MatchboxFactory>::ptr theFactory;

  /**
   * The particle group to be matched
   */
  string theGroup;

  /**
   * The set of particle ids to be matched
   */
  set<long> theIds;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxFactoryMatcher & operator=(const MatchboxFactoryMatcher &);

};

}

#endif /* Herwig_MatchboxFactoryMatcher_H */
