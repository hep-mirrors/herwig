// -*- C++ -*-
//
// FIqgxDipole.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_FIqgxDipole_H
#define HERWIG_FIqgxDipole_H
//
// This is the declaration of the FIqgxDipole class.
//

#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief FIqgxDipole implements the D_{q,g}^a subtraction dipole.
 *
 */
class FIqgxDipole: public SubtractionDipole {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FIqgxDipole();

  /**
   * The destructor.
   */
  virtual ~FIqgxDipole();
  //@}

public:

  /**
   * Return true, if this dipole applies to the selected
   * configuration.
   */
  virtual bool canHandle(const cPDVector& partons,
			 int emitter, int emission, int spectator) const;

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
  FIqgxDipole & operator=(const FIqgxDipole &);

};

}

#endif /* HERWIG_FIqgxDipole_H */
