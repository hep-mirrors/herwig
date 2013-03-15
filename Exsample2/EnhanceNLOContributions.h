// -*- C++ -*-
//
// EnhanceNLOContributions.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_EnhanceNLOContributions_H
#define Herwig_EnhanceNLOContributions_H
//
// This is the declaration of the EnhanceNLOContributions class.
//

#include "SamplingBias.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Enhance NLO virtual or real emission contributions
 */
class EnhanceNLOContributions: public SamplingBias {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  EnhanceNLOContributions();

  /**
   * The destructor.
   */
  virtual ~EnhanceNLOContributions();
  //@}

public:

  /**
   * Return an enhancement factor for the number of initial points to
   * presample the given XComb.
   */
  virtual double enhanceInitialPoints(const StandardXComb&) const;

  /**
   * Return an oversampling factor for the the given XComb.
   */
  virtual double oversamplingFactor(const StandardXComb&) const;

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
   * The enhancement factor for the number of initial points to
   * presample the given XComb.
   */
  double theEnhanceInitialPointsVirtual;

  /**
   * The oversampling factor for the the given XComb.
   */
  double theOversamplingFactorVirtual;

  /**
   * The enhancement factor for the number of initial points to
   * presample the given XComb.
   */
  double theEnhanceInitialPointsReal;

  /**
   * The oversampling factor for the the given XComb.
   */
  double theOversamplingFactorReal;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EnhanceNLOContributions & operator=(const EnhanceNLOContributions &);

};

}

#endif /* Herwig_EnhanceNLOContributions_H */
