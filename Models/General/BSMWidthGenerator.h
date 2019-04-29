// -*- C++ -*-
//
// BSMWidthGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_BSMWidthGenerator_H
#define HERWIG_BSMWidthGenerator_H
//
// This is the declaration of the BSMWidthGenerator class.
//

#include "Herwig/PDT/GenericWidthGenerator.h"
#include "Herwig/Decay/General/GeneralTwoBodyDecayer.fh"
#include "BSMWidthGenerator.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class is designed to be able to calculate the running
 * width for a BSM particle given the decay mode and partial width
 * calculated from the decayer.
 *
 * @see \ref BSMWidthGeneratorInterfaces "The interfaces"
 * defined for BSMWidthGenerator.
 */
class BSMWidthGenerator: public GenericWidthGenerator {

public:

  /**
   * The default constructor.
   */
  BSMWidthGenerator() : theModes(0) {}

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

  /**
   * Perform the set up for a mode in classes inheriting from this one
   * @param mode The decay mode
   * @param decayer The decayer for the mode.
   * @param imode The number of the mode.
   */
  virtual void setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer, 
			 unsigned int imode);

  /**
   * The \f$1\to2\f$ width for outgoing particles which can be off-shell.
   * @param iloc The location of the mode in the list.
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the first outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @return The partial width.
   */
  virtual Energy partial2BodyWidth(int iloc,Energy m0,Energy m1,Energy m2) const;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BSMWidthGenerator & operator=(const BSMWidthGenerator &) = delete;

private:
  
  /**
   * A vector decay modes with their associated decayer
   * cast to GeneralTwoBodyDecayer.
   */
  vector<pair<tcDMPtr, tcGeneralTwoBodyDecayerPtr> > theModes;
};

  /** 
   * An exception class to report an error.
   */
  class BSMWidthException : public Exception {};
  
}

#endif /* HERWIG_BSMWidthGenerator_H */
