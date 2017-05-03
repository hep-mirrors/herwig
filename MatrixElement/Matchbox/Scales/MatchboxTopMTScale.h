// -*- C++ -*-
//
// MatchboxTopMTScale.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxTopMTScale_H
#define Herwig_MatchboxTopMTScale_H
//
// This is the declaration of the MatchboxTopMTScale class.
//

#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxScaleChoice.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Daniel Rauch
 *
 * \brief MatchboxTopMTScale implements a scale choice related to the average
 *        of the transverse masses of a top and antitop quark.
 *
 */
class MatchboxTopMTScale: public MatchboxScaleChoice {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxTopMTScale();

  /**
   * The destructor.
   */
  virtual ~MatchboxTopMTScale();
  //@}

public:

  /**
   * Return the renormalization scale.
   */
  virtual Energy2 renormalizationScale() const;

  /**
   * Return the factorization scale.
   */
  virtual Energy2 factorizationScale() const;

  /**
   * Return the shower hard scale.
   */
  virtual Energy2 showerScale() const;


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
   * Switch to choose the definition of the shower hard scale.
   */
  unsigned int theShowerScaleMode;

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxTopMTScale & operator=(const MatchboxTopMTScale &);

};

}

#endif /* Herwig_MatchboxTopMTScale_H */
