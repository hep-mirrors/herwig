// -*- C++ -*-
//
// IFMassiveTildeKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_IFMassiveTildeKinematics_H
#define HERWIG_IFMassiveTildeKinematics_H
//
// This is the declaration of the IFMassiveTildeKinematics class.
//

#include "Herwig/MatrixElement/Matchbox/Phasespace/TildeKinematics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Martin Stoll
 *
 * \brief IFMassiveTildeKinematics implements the 'tilde' kinematics for
 * a initial-final subtraction dipole.
 *
 */
class IFMassiveTildeKinematics: public TildeKinematics {

public:

  /**
   * Perform the mapping to the tilde kinematics for the
   * last selected process and store all dimensionless
   * variables in the subtractionParameters() vector.
   * Return false, if the calculation of the tilde
   * kinematics was impossible for the selected configuration
   * and true on success.
   */
  virtual bool doMap();

  /**
   * Return the pt associated to the last merged splitting.
   */
  virtual Energy lastPt() const;
  
  /**
   * Return the pt associated to emitter emission and sppectator momentum.
   */
  virtual Energy lastPt(Lorentz5Momentum,Lorentz5Momentum,Lorentz5Momentum) const ;

  /**
   * Given a pt, return the boundaries on z
   */
  virtual pair<double,double> zBounds(Energy pt, Energy hardPt ) const;

  /**
   * Return the momentum fraction associated to the last splitting.
   */
  virtual double lastZ() const;

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
  IFMassiveTildeKinematics & operator=(const IFMassiveTildeKinematics &) = delete;

};

}

#endif /* HERWIG_IFMassiveTildeKinematics_H */
