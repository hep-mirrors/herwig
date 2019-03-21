// -*- C++ -*-
//
// IILightInvertedTildeKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_IILightInvertedTildeKinematics_H
#define HERWIG_IILightInvertedTildeKinematics_H
//
// This is the declaration of the IILightInvertedTildeKinematics class.
//

#include "Herwig/MatrixElement/Matchbox/Phasespace/InvertedTildeKinematics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief IILightInvertedTildeKinematics inverts the final-final tilde
 * kinematics.
 *
 */
class IILightInvertedTildeKinematics: public Herwig::InvertedTildeKinematics {

public:

  /**
   * Perform the mapping of the tilde kinematics for the
   * last selected process and store all dimensionless
   * variables in the subtractionParameters() vector.
   * Return false, if the calculation of the real
   * kinematics was impossible for the selected configuration
   * and true on success.
   */
  virtual bool doMap(const double *);

  /**
   * Return the pt associated to the last generated splitting.
   */
  virtual Energy lastPt() const;

  /**
   * Return the momentum fraction associated to the last splitting.
   */
  virtual double lastZ() const;

  /**
   * Return the upper bound on pt
   */
  virtual Energy ptMax() const;

  /**
   * Given a pt, return the boundaries on z
   */
  virtual pair<double,double> zBounds(Energy pt, Energy hardPt = ZERO) const;

  /**
   * Return true, if this InvertedTildeKinematics object needs to transform
   * all other particles in the process except the emitter, emission and spectator
   */
  virtual bool doesTransform() const { return true; }

  /**
   * If this InvertedTildeKinematics object needs to transform all other particles
   * in the process except the emitter, emission and spectator, return the transformed
   * momentum.
   */
  virtual Lorentz5Momentum transform(const Lorentz5Momentum& p) const {
    return p-(2.*(KplusKtilde*p)/KplusKtilde2)*KplusKtilde+(2.*(Ktilde*p)/K2)*K;
  }

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

  Lorentz5Momentum K;
  Energy2 K2;
  Lorentz5Momentum Ktilde;
  Lorentz5Momentum KplusKtilde;
  Energy2 KplusKtilde2;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  IILightInvertedTildeKinematics & operator=(const IILightInvertedTildeKinematics &) = delete;

};

}

#endif /* HERWIG_IILightInvertedTildeKinematics_H */
