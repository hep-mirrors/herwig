// -*- C++ -*-
//
// IIqqxDipole.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_IIqqxDipole_H
#define HERWIG_IIqqxDipole_H
//
// This is the declaration of the IIqqxDipole class.
//

#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief IIqqxDipole implements the D^{q,qbar;b} subtraction dipole.
 *
 */
class IIqqxDipole: public SubtractionDipole {

public:

  /**
   * The default constructor.
   */
  IIqqxDipole();

public:

  /**
   * Return true, if this dipole can possibly handle the indicated
   * emitter.
   */
  virtual bool canHandleEmitter(const cPDVector& partons, int emitter) const {
    return emitter < 2 && abs(partons[emitter]->id()) < 7;
  }

  /**
   * Return true, if this dipole can possibly handle the indicated
   * splitting.
   */
  virtual bool canHandleSplitting(const cPDVector& partons, int emitter, int emission) const {
    return 
      canHandleEmitter(partons,emitter) &&
      partons[emitter]->id() == partons[emission]->id();
  }

  /**
   * Return true, if this dipole can possibly handle the indicated
   * spectator.
   */
  virtual bool canHandleSpectator(const cPDVector& partons, int spectator) const {
    return spectator < 2 && partons[spectator]->coloured();
  }

  /**
   * Return true, if this dipole applies to the selected
   * configuration.
   */
  virtual bool canHandle(const cPDVector& partons,
			 int emitter, int emission, int spectator) const;

  /**
   *  How to sample the z-distribution.
   *  FlatZ = 1
   *  OneOverZ = 2
   *  OneOverOneMinusZ = 3
   *  OneOverZOneMinusZ = 4
   */

  virtual int samplingZ() const {return 2;}

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const;

  /**
   * Return the matrix element averaged over spin correlations.
   */
  virtual double me2Avg(double ccme2) const;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  IIqqxDipole & operator=(const IIqqxDipole &) = delete;

};

}

#endif /* HERWIG_IIqqxDipole_H */
