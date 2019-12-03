// -*- C++ -*-
//
// AlphaEM.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_AlphaEM_H
#define HERWIG_AlphaEM_H
//
// This is the declaration of the AlphaEM class.
//

#include "ThePEG/StandardModel/AlphaEMBase.h"

namespace Herwig {

using namespace ThePEG;


/**
 * The AlphaEM class is an exact reimplementation of the electromagentic
 * coupling in FORTRAN HERWIG and is mainly intended for testing.
 * It uses that same hadronic parameterisation as in the ThePEG::SimpleEM
 * but differs in the treatment of the top and leptonic contribution.
 *
 * @see \ref AlphaEMInterfaces "The interfaces"
 * defined for AlphaEM.
 */
class AlphaEM: public AlphaEMBase {

public:

  /**
   * The default constructor.
   */
  AlphaEM() : _me(),_mmu(),_mtau(), _mtop() {}

  /**
   * The \f$\alpha_{EM}\f$. Return the value of the coupling at a
   * given \a scale using the given standard model object, \a sm.
   */
  virtual double value(Energy2 scale, const StandardModelBase &) const;

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

  /**
   *  The real part of the photon self-energy
   * @param ratio The ratio of the mass squared of the fermion to the scale squared,
   * \f$m^2/Q^2\f$.
   */
  double realPi(double ratio) const;

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
  AlphaEM & operator=(const AlphaEM &) = delete;


private:

  /**
   *  Masses of the Standard Model fermions we need for the
   *  self-energies
   */
  //@{
  /**
   *  Electron mass squared
   */
  Energy2 _me;

  /**
   *  Muon mass squared
   */
  Energy2 _mmu;

  /**
   *  Tau mass squared
   */
  Energy2 _mtau;

  /**
   *  Top mass squared
   */
  Energy2 _mtop;
  //@}
};

}

#endif /* HERWIG_AlphaEM_H */
