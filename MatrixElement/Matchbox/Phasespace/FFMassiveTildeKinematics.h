// -*- C++ -*-
//
// FFMassiveTildeKinematics.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_FFMassiveTildeKinematics_H
#define HERWIG_FFMassiveTildeKinematics_H
//
// This is the declaration of the FFMassiveTildeKinematics class.
//

#include "Herwig++/MatrixElement/Matchbox/Phasespace/TildeKinematics.h"
    
namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief FFMassiveTildeKinematics implements the 'tilde' kinematics for
 * a final-final subtraction dipole.
 *
 */
class FFMassiveTildeKinematics: public TildeKinematics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FFMassiveTildeKinematics();

  /**
   * The destructor.
   */
  virtual ~FFMassiveTildeKinematics();
  //@}

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

public:
  
  /**
   * Triangular / Kallen function
   */
  template <class T>
  inline T rootOfKallen (T a, T b, T c) const {
    return sqrt( a*a + b*b + c*c - 2.*( a*b+a*c+b*c ) ); }

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
  FFMassiveTildeKinematics & operator=(const FFMassiveTildeKinematics &);

};

}

#endif /* HERWIG_FFMassiveTildeKinematics_H */
