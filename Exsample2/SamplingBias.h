// -*- C++ -*-
//
// SamplingBias.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_SamplingBias_H
#define Herwig_SamplingBias_H
//
// This is the declaration of the SamplingBias class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Handlers/StandardXComb.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Bias sampling by certain criteria
 */
class SamplingBias: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SamplingBias();

  /**
   * The destructor.
   */
  virtual ~SamplingBias();
  //@}

public:

  /**
   * Return an enhancement factor for the number of initial points to
   * presample the given XComb.
   */
  virtual double enhanceInitialPoints(const StandardXComb&) const = 0;

  /**
   * Return an oversampling factor for the the given XComb.
   */
  virtual double oversamplingFactor(const StandardXComb&) const = 0;

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
  SamplingBias & operator=(const SamplingBias &);

};

}

#endif /* Herwig_SamplingBias_H */
