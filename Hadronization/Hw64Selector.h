// -*- C++ -*-
//
// Hw64Selector.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Hw64Selector_H
#define HERWIG_Hw64Selector_H
//
// This is the declaration of the Hw64Selector class.
//

#include "StandardModelHadronSpectrum.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup hadronization
 * The Hw64Selector class selects the hadrons produced in cluster decay using
 * the FORTRAN HERWIG variant of the cluster model.
 *
 * @see \ref Hw64SelectorInterfaces "The interfaces"
 * defined for Hw64Selector.
 */
class Hw64Selector: public StandardModelHadronSpectrum {

public:

  /**
   * The default constructor.
   */
  Hw64Selector() : StandardModelHadronSpectrum(0),
		   _pwtDIquarkS0( 1.0 ),_pwtDIquarkS1( 1.0 )
  {}

  /**
   * Method to return a pair of hadrons given the PDG codes of
   * two or three constituents
   * @param cluMass The mass of the cluster
   * @param par1 The particle pointer of the first constituent
   * @param par2 The particle pointer of the second constituent
   * @param par3 The particle pointer of the third constituent
   */
  virtual pair<tcPDPtr,tcPDPtr> chooseHadronPair(const Energy cluMass,
						 tcPDPtr par1, tcPDPtr par2) const;

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

protected :
  
  /**
   *  Weights for baryons
   */
  virtual double baryonWeight(long id) const;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Hw64Selector & operator=(const Hw64Selector &) = delete;

private:

  /**
   *  The weights for the diquarks
   */
  //@{
  /**
   * The probability of producting a spin-0 diquark.
   */
  double _pwtDIquarkS0;

  /**
   * The probability of producting a spin-1 diquark.
   */
  double _pwtDIquarkS1;
  //@}

};

}

#endif /* HERWIG_Hw64Selector_H */
